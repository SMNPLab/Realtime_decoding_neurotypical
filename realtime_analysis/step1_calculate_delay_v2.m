% step1_calculate_delay_v2.m

%to calculate realtime delay for state-dependent tms

%73-FDI, 74-APB- input channel numbers for muscles in EEG amplifiers
clc
clearvars
delete(gcp('nocreate')); 
parpool('local',5,'IdleTimeout', inf);
user= 'realtime';
addpath(genpath(['C:/Users/' user '/Box/SMNP_Lab/LSL_workingfolder/functions']))
addpath(genpath(['C:/Users/' user '/Box/SMNP_Lab/LSL_workingfolder/mvpa/functions']))
addpath(genpath(['C:/Users/' user '/Box/SMNP_Lab/LSL_workingfolder/mvpa/ttl']));
addpath(['C:/Users/' user '/Box/SMNP_Lab/script_library/toolboxes/fieldtrip-20201214']);
addpath(genpath(['C:/Users/' user '/Box/SMNP_Lab/script_library/toolboxes/MVPA-Light-master']));
addpath(genpath(['C:/Users/' user '/Box/SMNP_Lab/LSL_workingfolder/mvpa/classifiers']));
ft_defaults;
startup_MVPA_Light;

% define stimulation parameters
subjname=input('Enter subject ID: ','s');
state=input('Enter state to target [low, high]: ','s');
trig=input('Enter number of triggers: ');
isi=3;

% load calibration classifier
%cd to classifier location
load([uigetfile]); %get the classifier file and load it, can be obtained from step3_offline_classify code
freqs=[4:0.25:35];
lda_ensemble=offline_classify.lda;
ranks=offline_classify.ranks; % indices of features to keep
csp_params=offline_classify.csp_params;
c1=offline_classify.class1_dval_cutoff;
c2=offline_classify.class2_dval_cutoff;

% instantiate the library
disp('Loading the library...');
lib = lsl_loadlib();

% resolve a stream
disp('Resolving an EEG stream...');
result = {};
while isempty(result)
    result = lsl_resolve_byprop(lib,'name','NeuroneStream');
end

% create a new inlet
disp('Opening an inlet...');
inlet = lsl_inlet(result{1});

% user defined local parameters
num_chan = 64;
sampling_rate = 5000; % native sampling rate of eeg system
processing_rate = 1000; % processing rate of realtime algorithm
resample_rate = round(sampling_rate/processing_rate);
params.pow_freq=freqs;
params.max_triggers=trig; % max number of pulses given

total_triggers= 0; % variable that counts number of pulses given
buffer_length = 0.5; % seconds
num_triggers = 3; % trigger value
params.elec_interest=[1:62]; % electrodes of interest (i.e., eeg electrodes only)

% standard local parameters
address = hex2dec('7FE8'); % startech breakout card parallel port address
eeg_address = hex2dec('2FE8');

% standard global parameters
if strcmpi('high',state)
    params.state = repmat(2,params.max_triggers,1);
    params.eeg_trig_val=repmat(200,params.max_triggers,1);
elseif strcmpi('low',state)
    params.state = repmat(1,params.max_triggers,1);
    params.eeg_trig_val=repmat(100,params.max_triggers,1); 
end

params.fsamp = sampling_rate;
params.isi = isi + rand(params.max_triggers,1); % interstimulus interval + random jitter
params.state_estimation = true; % logical variable to enable/disable phase estimation
stream_eeg = nan(buffer_length.*params.fsamp, num_chan);

if (params.state_estimation)
    params.pulse_num = num_triggers;
    params.trigger_timer = tic;
    
    params.theta = (find(params.pow_freq == 4)):(find(params.pow_freq == 7.75));
    params.alpha1 = (find(params.pow_freq == 8)):(find(params.pow_freq == 10));
    params.alpha2 = (find(params.pow_freq == 10.25)):(find(params.pow_freq == 13));
    params.beta1 = (find(params.pow_freq == 13.25)):(find(params.pow_freq == 20));
    params.beta2 = (find(params.pow_freq == 20.25)):(find(params.pow_freq == 35));
    
    params.state_duration = 0.010; %10 segments after the original 500ms segment-uk check this (2/19/23)
    params.address = address;
    params.eeg_address = eeg_address; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize parallel port %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
config_io;
global cogent; % verify inpoutx64 driver has been successfully been initialized
if( cogent.io.status ~= 0 )
    error('inp/outp installation failed');
end
outp(address,0); % flush parallel port buffer to prepare for stimulation

%%%%%%%%%%%%%%%
% main loop %
%%%%%%%%%%%%%%%

disp('Now receiving chunked data...');
ongoing_state=[];
start_time = tic;
while true
    % get chunk from the inlet
    [new_samples, ts] = inlet.pull_chunk();
    n_samples = length(ts);
    
    % buffer data
    if (n_samples > size(stream_eeg, 1))
        stream_eeg(1:end,:) = new_samples(:, end-size(stream_eeg, 1)+1:end)';
        
    elseif (n_samples > 0)
        stream_eeg(1:(end-n_samples),:) = stream_eeg((n_samples+1):end,:);
        stream_eeg((end-n_samples+1):end,:) = new_samples';
    end
    
    % resampling from native freq to processing freq
    if (sampling_rate ~= processing_rate)
        buffer_eeg = stream_eeg(1:resample_rate:end, :); % sample x chan
    else
        buffer_eeg = stream_eeg;
    end
    
    % select relevant electrodes
    buffer_eeg = buffer_eeg(:,params.elec_interest);
    
    if (~any(isnan(buffer_eeg(:))) && params.state_estimation && not(isempty(ts)))
        
        % this loop must execute faster than 1000 / processing_rate to maintain realtime performance
        if (toc(params.trigger_timer) > params.isi(total_triggers+1) && params.state_estimation)
            
            % rereference data to CAR
            [chunk,dat]=ft_preproc_rereference(buffer_eeg'); % input data must be nchans x ntime
            chunk=chunk'; % convert back to ntime x nchans
            
            % detrend/demean
            chunk=detrend(chunk,'constant');
            chunk=detrend(chunk,'linear'); % samp x chan 
            rr=reshape(chunk,[1,500,62]); % convert to trl(1) x sample x trial - this is working properly 
            csp_chunk=permute(rr,[1,3,2]); % trl x chan x sample
         
            parfor kk=1:5
                [ff, ~]= calc_features_rt(csp_params{kk},csp_chunk, params.pow_freq, processing_rate, params.theta, params.alpha1, params.alpha2, params.beta1, params.beta2);
                [clabel{kk},~]=test_lda_realtime(lda_ensemble{kk},ff(ranks{kk}),c1(kk),c2(kk),1);
            end
            
            if nanmean(cell2mat(clabel)) == 1.5
                state_now=NaN;
            else
                state_now=round(nanmean(cell2mat(clabel)));
            end
            
            ongoing_state=[ongoing_state; state_now];
            
            % if target state is present & maintained for minimum specified time interval
            if length(ongoing_state) > (params.state_duration*processing_rate) 
                if state_now == params.state(total_triggers+1) && sum(ongoing_state(end-(params.state_duration*processing_rate)+1:end)) == params.state(total_triggers+1)*(params.state_duration*1000)
                    outp(params.address,params.pulse_num);
                    outp(params.eeg_address,params.eeg_trig_val(total_triggers+1)); 

                    params.trigger_timer = tic;
                    total_triggers=total_triggers+1;
                    
                    if total_triggers>=params.max_triggers
                        break
                    end
                    
                    % flush the output buffer (parallel port)
                    disp(['trigger ' num2str(total_triggers) ' sent.']);
                    outp(params.address,0); 
                    outp(params.eeg_address,0); 
                    
                    % clear ongoing state to prepare for next trial
                    ongoing_state=[]; 
                    
                end
            end
            
        end
    end
    
end

% all channels of interest 
streamed=buffer_eeg;

% clear unecessary variables before calculating technical delay
clearvars -except streamed params subjname

% compare recorded to streamed data
% cd('C:/Users/realtime/Box/SMNP_Lab/LSL_workingfolder/mvpa/calibration');
% fn=uigetfile;
% load([fn]);
filename=input('Enter filename: ','s');

% load recorded data
cd(uigetdir);
cfg=[];
cfg.channel='Fp1';
cfg.dataset=[filename '.eeg'];
hdr=[filename '.vhdr'];
all_raw=ft_preprocessing(cfg);
recording_time=(numel(all_raw.time{1}) / 5000 ) / 60;

% segment recorded data
cfg=[];
cfg.dataset=[filename '.vhdr'];
cfg.trialdef.eventtype='Stimulus';
cfg.trialdef.eventvalue='A';
cfg.trialdef.prestim=0.500;
cfg.trialdef.poststim=0;
cfg=ft_definetrial(cfg);
trl=cfg.trl;
seg=ft_redefinetrial(cfg,all_raw);

% downsample recorded data to 1000 Hz to match resampling rate of streamed data
saved=downsample(seg.trial{params.max_triggers}(1,:),5);
saved=saved';
streamed=streamed(:,1);

% plot streamed and recorded data
figure;
plot(saved,'-k'); hold on; plot(streamed,'-r');
legend('saved','streamed');

% calculate delay
delay=input('Enter time point 1 and 2, in brackets : '); % sampling rate == 1000 hz, so delay in samples is the same as delay in ms
total_delay = delay(2) - delay(1);
close all;

% print results
figure;
plot(saved,'-k'); hold on; plot(streamed(delay(2)-delay(1)+1:end),'-r');
disp(['delay = ' num2str(total_delay) ' ms']);
disp(['recording time = ' num2str(recording_time) ' minutes']);

user='realtime';
% save technical delay
%cd to data sharing location
t=datenum(date);
save([num2str(subjname) '_technical_delay_' num2str(t) '.mat'],'total_delay');




