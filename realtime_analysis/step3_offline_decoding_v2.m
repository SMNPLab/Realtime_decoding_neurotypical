% step3_offline_decoding_v2.m

% decode high vs low MEP states from EEG/EMG offline


%% step 1: housekeeping
clc
clearvars
user='uuk58-admin';
addpath(['C:/Users/' user '/Box/SMNP_Lab/script_library/toolboxes/fieldtrip-20201214']);
addpath(genpath(['C:/Users/' user '/Box/SMNP_Lab/script_library/toolboxes/MVPA-Light-master']));
addpath(genpath(['C:/Users/' user '/Box/SMNP_Lab/LSL_workingfolder/mvpa/ttl']));
addpath(['C:/Users/' user '/Box/SMNP_Lab/LSL_workingfolder/mvpa/functions']);
%cd to data location
ft_defaults;
startup_MVPA_Light;

%% step 2: preprocessing eeg / emg data

subjname=input('Enter subject ID: ','s'); 
%cd to subject folder
technical_delay=input('Enter technical delay: ');
trial_nums=zeros(7,1);
for ii=[1,2,3,4,5,6]
    num=ii*100;
    eeg_file=['sub-' subjname '_task-trainingblock_' num2str(num) '_eeg.eeg'];
    vhdr_file=['sub-' subjname '_task-trainingblock_' num2str(num) '_eeg.vhdr'];
    vmrk_file=['sub-' subjname '_task-trainingblock_' num2str(num) '_eeg.vmrk'];
  

    % load eeg data
    cfg=[];
    cfg.dataset=eeg_file;
    hdr=vhdr_file;
    cfg.channel={'all'};
    raw=ft_preprocessing(cfg);

    % segment eeg and emg data together into 1000ms before and after the pulse 
    cfg=[];
    cfg.dataset=vhdr_file;
    cfg.trialdef.eventtype='Stimulus';
    cfg.trialdef.eventvalue='A';
    cfg.trialdef.prestim= 1.00; 
    cfg.trialdef.poststim= 1.00;
    cfg=ft_definetrial(cfg);
    trl=cfg.trl;
    segmented{ii}=ft_redefinetrial(cfg,raw);
    idex=ii+1;
    trial_nums(idex)= size(trl,1);

end


cfg=[];
seg=ft_appenddata(cfg,segmented{1},segmented{2},segmented{3},segmented{4},segmented{5},segmented{6});


% select eeg data only
cfg=[];
cfg.latency=[(-technical_delay/1000)-0.500 (-technical_delay/1000)]; % -5xx ms to -xx ms
cfg.channel={'all', '-EMG1', '-EMG2', '-EMG3', '-EMG4', '-EMG5'}; 
eeg=ft_selectdata(cfg,seg);

% rereference eeg data
cfg=[];
cfg.reref='yes';
cfg.refchannel='all'; % common average reference 
cfg.demean='yes';
cfg.detrend='yes';
reref=ft_preprocessing(cfg,eeg);

% reformat, downsample to 1000 Hz
for trl=1:length(reref.trial)
    for chan=1:length(reref.label)
        eeg_mat(:,trl,chan)=reref.trial{trl}(chan,:);
        eeg_down(:,trl,chan)=downsample(eeg_mat(:,trl,chan),5);
    end
end

% select emg data only 
cfg=[];
cfg.latency=[-0.100 0.400]; % -100 ms to +400 ms 
cfg.channel={'EMG1'};
emg=ft_selectdata(cfg,seg);

% reformat emg data, demean/detrend 
for trl=1:size(emg.trial,2)
    emg_mat(:,trl)=emg.trial{trl}(1,:);
    emg_mat(:,trl)=detrend(emg_mat(:,trl),'constant');
    emg_mat(:,trl)=detrend(emg_mat(:,trl),'linear');
end

% subtract baseline emg offset
clear emg
pre=mean(emg_mat(1:375,:),1); % pre-tms background emg levels
for trl=1:size(emg_mat,2)
    emg(:,trl)=emg_mat(:,trl)-pre(trl);
end

% calculating subject-specific bemg noise level - SJH updated, uttara should check  
pre_rms=rms(emg(1:375,:),1); 
rms_sd=std(pre_rms); 
emg_thresh=mean(pre_rms)+(2*rms_sd); 
bad=[];
for trl=1:size(emg,2)
    bemg=rms(emg(1:375,trl)); 
    if bemg > emg_thresh
        bad=[bad; trl];
    end
end

% plot all meps
figure
plot(emg);
win=input('Enter MEP analysis window [x1 x2]: ');
close all;

% calculate mep amplitudes
for trl=1:size(emg,2)
    pp(trl,1)=max(emg(win(1):win(2),trl))-min(emg(win(1):win(2),trl));
end

% remove bad trials from eeg and meps 
keep=[1:numel(pp)];
keep(bad)=[];
mep=pp(keep); %pp w/o bad trials
eeg=eeg_down(:,keep,:); %timeseries x trials x channel

% blockwise demean/detrend/rescale of mep amplitudes (theoretically, bad
% trials should have been removed after detrending, but this is not a huge
% issue - sjh noticed this on 5/23/23 - do not change mid-data collection ) 

plot(mep) % plot good mep amps 
new_pp=[]; 
for nums=1:6
    tr_meps_ini=keep(find(keep<=nums*100));
    tr_meps=tr_meps_ini(find(tr_meps_ini>(nums-1)*100));
    numel(tr_meps)
    
    tmp_pp=pp(tr_meps); 
    tmp_pp=detrend(tmp_pp,'constant');
    tmp_pp=detrend(tmp_pp,'linear');
    tmp_pp=zscore(tmp_pp);
    tmp_pp=rescale(tmp_pp); 
    
    new_pp=[new_pp; tmp_pp]; 
    
end

plot(new_pp) %plot detrended mep amplitudes

%check that emg / eeg trial numbers match
if length(new_pp) ~= size(eeg,2)
    error('emg and eeg trial mismatch!');
end

% combine and save data
preproc.mep.amplitude=new_pp; % blockwise detrended meps, w/o bad trials
preproc.raw.emg=emg; % fs = 5000 hz, all with bad trials
preproc.bemg=pre_rms; %all bemg with bad trials
preproc.raw.all_eeg_trials=eeg_down; %fs=1000 hz with bad trials
preproc.mep.all_meps_no_detrend=pp; %all meps including bad trials, no blockwise detrend
preproc.bad=bad; %bad trial index
preproc.keep=keep; % trial indices to keep

% save all preprocessed data
t=datenum(date);
save([subjname '_preproc_' num2str(t) '.mat'],'preproc','-v7.3');

% clear unecessary variables before performing classification
clearvars -except filename new_pp bad technical_delay subjname keep eeg user

%% step 3: building classifiers

% remove bad trials from mep and downsampled eeg data
mep=new_pp;

% check that meps and features have same number of trials
if length(mep) ~= size(eeg,2)
    error('dimension mismatch');
end

% % dichotimize meps
thresh=median(mep);
for trl=1:length(mep)
    if mep(trl) >= thresh
        label(trl)=2;
    elseif mep(trl) < thresh
        label(trl)=1;
    end
end

% reshape downsampled eeg data for csp calculation
raw_csp=permute(eeg,[2,3,1]); % trial x channel x samples

% build classifier using csp+psd approach
[perf,perf_all,predlabel,runtime,ki,lda,pparam_csp,dval,plabel] = manual_classify_csp_psd(5,raw_csp,label); % perf = k x featurenum x lambdaval

% select best performing classifier without rounding
grid_perf=squeeze(mean(perf,1)); % average across k to get featurenum x lambdaval
f1_nonround=max(max(grid_perf)); 
[fk,lk]=find(grid_perf==f1_nonround); % identify indices for optimal model configuration
lt=length(lk); % take highest lambda value to reduce overfitting
featkeep_nonround=fk(1);
lambdakeep_nonround=lk(lt);

% select best performing classifier with rounding
grid_perf_round=round(squeeze(mean(perf,1)),2); 
f1_round=max(max(grid_perf_round));
[fk_round,lk_round]=find(grid_perf_round==f1_round);
lt_round=length(lk_round);
featkeep_round=fk_round(1);
lambdakeep_round=lk_round(lt_round);

% choose rounded or non rounded classifier
if featkeep_round < featkeep_nonround && (f1_nonround - f1_round) < 0.05 
    featkeep=featkeep_round;
    lambdakeep=lambdakeep_round;
    f1=f1_round;
else
    featkeep=featkeep_nonround;
    lambdakeep=lambdakeep_nonround;
    f1=f1_nonround;
end

% get performance / classifier info without thresholding
grid_perf_all=squeeze(mean(perf_all,1)); % average across k to get featurenum x lambdaval
f1_nonround_all=max(max(grid_perf_all));
[fkall,lkall]=find(grid_perf_all==f1_nonround_all); % identify indices for optimal model configuration
ltall=length(lkall); % take highest lambda value to reduce overfitting
featkeep_nonround_all=fkall(1);
lambdakeep_nonround_all=lkall(ltall);

% select best performing classifier with rounding
grid_perf_round_all=round(squeeze(mean(perf_all,1)),2);
f1_round_all=max(max(grid_perf_round_all));
[fk_round_all,lk_round_all]=find(grid_perf_round_all==f1_round_all);
lt_round_all=length(lk_round_all);
featkeep_round_all=fk_round_all(1);
lambdakeep_round_all=lk_round_all(lt_round_all);

% choose rounded or non rounded classifier
if featkeep_round_all < featkeep_nonround_all && (f1_nonround_all - f1_round_all) < 0.05 
    featkeep_all=featkeep_round_all;
    lambdakeep_all=lambdakeep_round_all;
    f1_all=f1_round_all;
else
    featkeep_all=featkeep_nonround_all;
    lambdakeep_all=lambdakeep_nonround_all;
    f1_all=f1_nonround_all;
end

% select all items needed for realtime classifier
for k=1:5
    lda_ensemble{k}=lda{featkeep,lambdakeep,k}; % feature x lambda x k
    ensemble_ranks{k}=ki{k,featkeep}; % ranks of fold k when there are featkeep features included - i.e., feature indices for each lda per fold
    
    dval_test{k}=dval{featkeep,lambdakeep,k};
    plabel_new{k}=plabel{featkeep,lambdakeep,k};
    dval_low{k}=dval_test{k}(find(plabel_new{k}==1)); % find dvals associated with low states
    dval_high{k}=dval_test{k}(find(plabel_new{k}==2)); % find dvals associated with high states

    dval_low_perc{k}=prctile(dval_low{k},50,"all");
    dval_high_perc{k}=prctile(dval_high{k},50,"all");

    lda_ensemble_all{k}=lda{featkeep_all,lambdakeep_all,k}; % feature x lambda x k
    ensemble_ranks_all{k}=ki{k,featkeep_all}; % ranks of fold k when there are featkeep features included - i.e., feature indices for each lda per fold

end

ss=linspace(1e-10,1); % repeat this from inside function
optimal_p=mv_get_hyperparameter('lda'); % get optimal hyperparams
optimal_p.lambda=ss(lambdakeep);
optimal_p.reg='shrink';
optimal_p.prob=1;  % added to calculate probabilities
optimal_p.form='primal'; % added to calculate probabilities

disp(['f1 = ' num2str(f1)]);

offline_classify.f1=f1;
offline_classify.lda=lda_ensemble;
offline_classify.ranks=ensemble_ranks; % indices of features used for final classifier
offline_classify.hyperparams=optimal_p; % optimal hyperparameters used for final classifier
offline_classify.runtime=runtime;
offline_classify.grid_perf=grid_perf; % k averaged performance across all possible model configurations
offline_classify.technical_delay=technical_delay;
offline_classify.csp_params=pparam_csp;
offline_classify.class1_dval_cutoff=cell2mat(dval_low_perc);
offline_classify.class2_dval_cutoff=cell2mat(dval_high_perc);

offline_classify.f1_all=f1_all;
offline_classify.lda_all=lda_ensemble_all;
offline_classify.ranks_all=ensemble_ranks_all;

t=datenum(date);
save([ subjname '_offline_classify_' num2str(t) '.mat'],'offline_classify','-v7.3');

