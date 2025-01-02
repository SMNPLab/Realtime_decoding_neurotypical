
%% load the defaults 
%load dependencies
clc
clearvars
user='';
%path to fieldtrip
%path to mvpa light
%path to custom written functions
%path to subject data folder
ft_defaults;
startup_MVPA_Light;

%% define the subject related details
subjname={'01','02','03','04','05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19'};
states={'high', 'low', 'random'};
intensities= {'120', '110'};
emgchannels={'EMG1', 'EMG2'}; % the EMG channels of interest (EMG1- FDI, EMG2- APB) 

%% loop across subjects and states to get the mep amplitudes and bemg data

for jjj=1:numel(subjname)
    % cd to subject data folder
    for intensity=1:2 %looping through intensities (120,110) 

        for state= 1:3 %looping through states (high, low, random)

            % define filename based on the condition it belongs to (for weirdly named files) 

            if strcmpi(subjname{jjj}, '18')==1 && strcmpi(states{state}, 'high') && strcmpi(intensities{intensity}, '110')
                filename= (['sub-' subjname{jjj} '_task-testinghigh110RMTfile1_eeg']);
            elseif strcmpi(subjname{jjj}, '14')==1 && strcmpi(states{state}, 'high') && strcmpi(intensities{intensity}, '110')
                filename= (['sub-' subjname{jjj} '_task-testinghigh110RMTfile1_eeg']);
            else 
                filename= (['sub-' subjname{jjj} '_task-testing' states{state} intensities{intensity} 'RMT_eeg']);
            end

            eeg_file=[filename '.eeg'];

            if exist(eeg_file, 'file' )==2  %putting an if condition to confirm if that file exists (if that testing session could be conducted)

                if strcmpi(subjname{jjj}, '18')==1 && strcmpi(states{state}, 'high') && strcmpi(intensities{intensity}, '110')

                    for block=1:3 % for special cases where data was saved in multiple file so need to be appended together here 
                        filename= (['sub-' subjname{jjj} '_task-testing' states{state} intensities{intensity} 'RMTfile' num2str(block) '_eeg']);

                        % define eeg filenames 
                        eeg_file=[filename '.eeg'];
                        vhdr_file=[filename '.vhdr'];
                        vmrk_file=[filename '.vmrk'];

                        % read data from eeg file
                        cfg=[];
                        cfg.dataset=eeg_file;
                        hdr=vhdr_file;
                        raw=ft_preprocessing(cfg); %raw eeg data file

                        %segment data 500 msec before and after the tms pulse for
                        %each trial
                        cfg=[];
                        cfg.dataset=vhdr_file;
                        cfg.trialdef.eventtype='Stimulus';
                        cfg.trialdef.eventvalue='A';
                        cfg.trialdef.prestim=0.500;
                        cfg.trialdef.poststim=0.500;
                        cfg=ft_definetrial(cfg);
                        trl=cfg.trl;
                        seg{block}=ft_redefinetrial(cfg,raw);

                    end

                    cfg=[];
                    seg_eeg=ft_appenddata(cfg,seg{1},seg{2},seg{3});
                    
                elseif strcmpi(subjname{jjj}, '14')==1 && strcmpi(states{state}, 'high') && strcmpi(intensities{intensity}, '110')

                    fname_blocks={'sub-14_task-testing_high110RMTfile1_eeg','sub-14_task-testinghigh110RMTfile2_eeg' }

                    for blocks=1:2
                         % define eeg filenames 
                        eeg_file=[fname_blocks{blocks} '.eeg'];
                        vhdr_file=[fname_blocks{blocks} '.vhdr'];
                        vmrk_file=[fname_blocks{blocks} '.vmrk'];


                        % load eeg data
                        cfg=[];
                        cfg.dataset=eeg_file;
                        hdr=vhdr_file;
                        raw=ft_preprocessing(cfg);

                        % segment eeg data
                        cfg=[];
                        cfg.dataset=vhdr_file;
                        cfg.trialdef.eventtype='Stimulus';
                        cfg.trialdef.eventvalue='A';
                        cfg.trialdef.prestim= 0.500; 
                        cfg.trialdef.poststim= 0.500; 
                        cfg=ft_definetrial(cfg);
                        trl=cfg.trl;
                        seg{blocks}=ft_redefinetrial(cfg,raw);
                    end
                    cfg=[];
                    seg_eeg=ft_appenddata(cfg,seg{1},seg{2});

                else % regular cases where nothing was weird 
                     % define eeg filenames 
                        eeg_file=[filename '.eeg'];
                        vhdr_file=[filename '.vhdr'];
                        vmrk_file=[filename '.vmrk'];

                        % read data from eeg file
                        cfg=[];
                        cfg.dataset=eeg_file;
                        hdr=vhdr_file;
                        raw=ft_preprocessing(cfg); %raw eeg data file

                        %segment data 500 msec before and after the tms pulse for
                        %each trial
                        cfg=[];
                        cfg.dataset=vhdr_file;
                        cfg.trialdef.eventtype='Stimulus';
                        cfg.trialdef.eventvalue='A';
                        cfg.trialdef.prestim=0.500;
                        cfg.trialdef.poststim=0.500;
                        cfg=ft_definetrial(cfg);
                        trl=cfg.trl;
                        seg_eeg=ft_redefinetrial(cfg,raw);
                end

                % select only emg channels and segment data around pulse 
                cfg=[];
                cfg.latency=[-0.100 0.400]; 
                cfg.channel= {'EMG1', 'EMG2'};
                emg_struc=ft_selectdata(cfg,seg_eeg);

                % reformat by detrending constant and linear on a trial and channel by basis
                for chan=1:numel(emg_struc.label) %gives the number of channels
                    for trl=1:length(emg_struc.trial) % gives the number of trials
                        emg(:,chan,trl)=emg_struc.trial{trl}(chan,:); % loads the emg data segment for that trial and emg channel 
                        emg(:,chan,trl)=detrend(emg(:,chan,trl),'constant'); 
                        emg(:,chan,trl)=detrend(emg(:,chan,trl),'linear'); %samples*channels*trials
                    end
                end

                %emg stored as samples*channels*trials

                % subtract the baseline emg offset (background emg) from
                % the emg signal for each channel and trial-- not filtering cuz here we just want to remove the emg offset

                pre=mean(emg(1:375,:,:),1); %% emg signal from 1-75msec from start of emg segment (25msec before the pulse), 1msec has 5 samples (sample frequency = 5kHz)

                for trl=1:size(pre,3)
                    for chan=1:size(pre,2) %pre--1*channel*trial
                        emg_wo_offset(:,chan,trl)=emg(:,chan,trl)-pre(1,chan,trl); %emg-samples*channels*trials, pre- one value per trial, channel and subject
                    end
                end
                %emg_wo_offset--samples*channels*trials

                %plot all meps from both the channels and take input for mep window selection
                figure
                subplot(1,2,1)
                plot(squeeze(emg_wo_offset(:,1,:))); title('FDI'); hold on;
                subplot(1,2,2)
                plot(squeeze(emg_wo_offset(:,2,:))); title('APB')

                %load the windows you have already run from the mep_windows
                %file
                win(1,:)=input('Enter MEP analysis window for fdi [x1 x2]: ');
                win(2,:)=input('Enter MEP analysis window for apb [x1 x2]: ');
                close all;

                % calculate mep amplitudes by calculating maximum and min from the selected mep window
                for trl=1:size(emg_wo_offset,3)
                    for chan=1:size(emg_wo_offset,2)
                        mep_amp(chan,trl)=max(emg_wo_offset(win(chan,1):win(chan,2),chan,trl))-min(emg_wo_offset(win(chan,1):win(chan,2),chan,trl)); %samples*channels*trials
                    end
                end

                %find weird mep trials and reject them
                for chan=1:size(mep_amp,1)
                    figure
                    sgt=sgtitle([ subjname{jjj} ' ' intensities{intensity} ' ' states{state} ' ' emgchannels{chan}])
                    sgt.FontSize=12
                    for trl=1:size(mep_amp,2) 
                        subplot(5,6,trl)
                        plot(squeeze(emg_wo_offset(win(chan,1):win(chan,2),chan,trl))); yline(max(emg_wo_offset(win(chan,1):win(chan,2),chan,trl))); yline(min(emg_wo_offset(win(chan,1):win(chan,2),chan,trl))); box off; title(num2str(trl)); hold on;
                    end
                    weird_mep_bad{chan}=input("enter all trials with weird meps : "); %enter all the trials for that channel with weird meps
                end
                close all;
    %-----------------------------------------------   --------------------------------------------------------------------------- 
                %BEMG calcultions 

                %segment data 100-25msec before the pusle for bemg rms calculations-- need to do this to identify bad trials and reject them
                cfg=[];
                cfg.latency=[-0.100 -0.025];
                cfg.channel= {'EMG1', 'EMG2'}; %fdi and apb
                prestim=ft_selectdata(cfg,seg_eeg);

                % filter, demean and detrend
                cfg=[];
                cfg.demean='yes';
                cfg.detrend='yes';
                cfg.dftfilter='yes';
                cfg.dftfreq=[30 60 90 120 150 180];
                prestim_dd=ft_preprocessing(cfg,prestim); %dd stands for demeaned and detrended

                for chan=1:numel(prestim.label) %gives the number of emg channels
                    for trl=1:numel(prestim.trial) %gives the number of trials
                        bemg_rms(chan,trl)=rms(prestim_dd.trial{trl}(chan,:)); % this just gives one value per channel and trial
                    end
                end

%                 %bemg_rms-- channel*trials
%     
                %find the bad trials threshold based on bemg value and reject bad trials
                for chan=1:numel(prestim.label)
                    bemg_thresh(chan)=mean(bemg_rms(chan,:))+2*std(bemg_rms(chan,:)); % this takes the mean and std across trials for each channel 
                end
                %bemg_thresh- one value per channel

                %find the bad trials using the bemg threshold mep_amp_goodmethod- you want to find these but not store these 
                % not actually using this in statistics, including all
                % trials and bemg as covariate 
                for chan=1:numel(prestim.label)
                    bad=[];
                    for trl=1:numel(prestim.trial)  
                        if bemg_rms(chan,trl)> bemg_thresh(chan)
                            bad=[bad trl];
                        end  
                    end
                    if numel(bad)>0
                        bad_bemg{chan}=bad; %stores the bad trials for all the channels in one structure
                    else 
                        bad_bemg{chan}=NaN;
                    end

                end
% 
%                 %%-----------------------------------------------
                %calculations for emg analysis

                %reject bad trials- based on weird mep waveform, unclear
                %signal to noise
                mep_amp_good=mep_amp;
                for chan=1:size(mep_amp,1)
                    % if ~isnan(bad_bemg)
                    %     mep_amp_good(chan,bad_bemg(chan,:))=NaN;
                    % end
                    if ~isnan(weird_mep_bad{chan})
                        mep_amp_good(chan,weird_mep_bad{chan})=NaN;
                    end
                end

                %create final data structure
                results.mep_amp{jjj,state,intensity}=mep_amp; %mep amplitude for each channel
                results.bemg_rms{jjj,state,intensity}=bemg_rms; %bemg rms value for each channel and trial 
                results.bad_trials_bemg{jjj,state,intensity}=bad_bemg; %bad trials for all channels (different channels might have different bad trials)
                results.emg_wo_offset{jjj,state,intensity}=emg_wo_offset; %emg signals without the bemg offset
                results.bad_trials_weirdmeps{jjj,state,intensity}=weird_mep_bad; %bad trials weird meps for all channels (different channels might have different bad trials)
                results.emg_with_offset{jjj,state,intensity}=emg; %detrended and demeaned emg              
                results.mep_window{jjj,state,intensity}=win; %first row is fdi and second row is apb
                results.state{jjj,state,intensity}=state;
                results.mep_amp_good{jjj,state,intensity}=mep_amp_good; %channels*trials - only rejects weirdly shaped/noisy meps, not ones with high bemg - this is used for main analysis 
                results.intensity{jjj,state,intensity}=intensity;
                results.uniquenum{jjj,state,intensity}=uniquenum{jjj};
                results.subjname{jjj,state,intensity}=subjname{jjj};

            else %if data is not available for subj-intensity-muscle combo 

                results.mep_amp{jjj,state,intensity}=NaN(2,25);
                results.bemg_rms{jjj,state,intensity}=NaN(2,25);
                results.bad_trials_bemg{jjj,state,intensity}=NaN;
                results.emg_wo_offset{jjj,state,intensity}=NaN;
                results.bad_trials_weirdmeps{jjj,state,intensity}=NaN;
                results.mep_window{jjj,state,intensity}=NaN;
                results.mep_amp_good{jjj,state,intensity}=NaN(2,25); %channels*trials%bad trials weird meps for all channels (different channels might have different bad trials)
                results.state{jjj,state,intensity}=state;
                results.intensity{jjj,state,intensity}=intensity;
                results.uniquenum{jjj,state,intensity}=uniquenum{jjj};
                results.subjname{jjj,state,intensity}=subjname{jjj};

            end

            clearvars -except subjname uniquenum user intensity emgchannels results jjj states intensities mep_win

        end

        % calculate normalized mep amplitudes
        %  group all meps from all states within one intensity
        % normalize the mep amplitudes-- do it subject wise and channel wise
        % we are looping over all data for each muscle/state over intensities right now 
        mean_meps_combined=nanmean(horzcat(results.mep_amp_good{jjj,:,intensity}),2);%subject,states,intensity,chan - gives two means, one per channel 
%       % ^ gives us mep normalization factor per subject, intensity and muscle 

        for chan=1:2
            for state=1:3
                meps_norm{chan,state}= results.mep_amp_good{jjj,state,intensity}(chan,:) ./mean_meps_combined(chan);
                mean_meps(chan,state)=nanmean(meps_norm{chan,state});
                results.normalized_meps{jjj,intensity,chan,state}= meps_norm{chan,state}; %subj*intensity*chan*state
            end
        end
% 
        results.normalized_mean_meps{jjj,intensity}=mean_meps; %chan*state
        results.mean_meps_combined{jjj,intensity}=mean_meps_combined; %chan - normalization factor 
      
        %calculate mep mod
        for chan=1:2
            mod_high_low(chan)= mean_meps(chan,1) ./ mean_meps(chan,2);
            mod_high_random(chan) = mean_meps(chan,1) ./ mean_meps(chan,3);
            %mod_random_low(chan) = mean_meps(chan,3) ./ mean_meps(chan,2);
        end 
        
        results.mod_high_low{jjj,intensity}=mod_high_low;
        results.mod_high_random{jjj,intensity}=mod_high_random;
        results.mod_random_low{jjj,intensity}=mod_random_low; % same, incorrect but never used 

        %calculate coefficient of variation
        for chan=1:2
            for state=1:3
                coeff_of_var(chan,state)= nanstd(meps_norm{chan,state}) ./nanmean(meps_norm{chan,state});
            end
        end
        results.coeff_of_var{jjj,intensity}=coeff_of_var; %chan,state

        clearvars -except subjname uniquenum user intensity emgchannels results jjj states intensities mep_win

        %cd to saving location
        save('test_emg_data.mat','results','-v7.3');

    end

end

% 
%variables saved in results
% mep_amp- mep amplitude for each chanel
% bemg_rms- bemg rms value for each channel and trial 
% bad_trials_bemg- bad trials based on high bemg for all channels
% emg_wo_offset- emg signal detrended and demeaned prestim bemg subtracted
% emg_without_offset- detrended and demeaned emg
% bad_trial_weird_meps- bad trials rejected because of weird meps
% mep_amp_good- mep amplitudes only of good trials (removed bad trials from weird meps)
% states- high,low and random strings
% intensities- intensitiy values
% uniquenum- uniquenum array
% subjname- subjname array
% coeff_of_var- coefficient of variation array
% mod_high_low- mep mod for high/low states
% mod_high_random- mep mod for high/random states
% mod_random_low- mep mod for random/low states
% results.normalized_meps- subject*intensity*channel*state cell array
% normalized_mean_meps- mep amplitudes normalized with mean of mep amplitudes for all states and each intensity
% mean_meps_combined- mean mep of all the states combined per intensity and muscle

%% plot the meps for each intensity, condition and muscle 
clc
clearvars
user='';
%cd to data location
load('test_emg_data.mat');
subjname=results.subjname;
intensities={'120','110'};
states={'strong', 'weak', 'random'};
%for muscles 1=fdi, 2=apb

mean_meps=results.normalized_mean_meps; % mep amplitudes normalized with mean of mep amplitudes across states for each intensity {subject, intensity}(muscle,state)

%convert the mean meps from cell array to a matrix-
for jjj=1:19
    for iii=1:2
        for ccc=1:2
            for sss=1:3
                mean_meps_arr(jjj,ccc,sss,iii)=mean_meps{jjj,iii}(ccc,sss); %(subject,channel,state,intensity)
            end
        end
    end
end

%calculate sem across subjects on the matrix for non nan blocks only
for iii=1:2
    for ccc=1:2
        for sss=1:3
            mep_arr_used=mean_meps_arr(:,ccc,sss,iii);
            mep_sem(ccc,sss,iii)= nanstd(mep_arr_used)./sqrt(size(mep_arr_used(~isnan(mep_arr_used)),1));
        end
    end
end

% plot the average meps for state and intensity
col=["#ab590d", "#e39952", '#810290','#DA80E4']; 
muscle={'FDI', 'APB'}
counter=0;
for chan=1:2
    fig=figure;
    for intensity=1:2
        counter=counter+1;
        p(intensity)=subplot(1,2,intensity)  
        box off;
        for jjj=1:19
            plot([1,2,3],[mean_meps{jjj,intensity}(chan,2) mean_meps{jjj,intensity}(chan,1) mean_meps{jjj,intensity}(chan,3)],'color', col(counter) ,'linestyle',':', 'LineWidth',0.2); hold on; 
            scatter([1,2,3],[mean_meps{jjj,intensity}(chan,2) mean_meps{jjj,intensity}(chan,1) mean_meps{jjj,intensity}(chan,3)],'MarkerFaceColor',col(counter),'MarkerEdgeColor',col(counter),'SizeData',50,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); hold on; 
        end
        scatter([1,2,3],[nanmean(mean_meps_arr(:,chan,2,intensity)) nanmean(mean_meps_arr(:,chan,1,intensity)) nanmean(mean_meps_arr(:,chan,3,intensity))],'square','MarkerFaceColor',col(counter),'MarkerEdgeColor',col(counter),'SizeData',250,'LineWidth',1);
        errorbar([1,2,3],[nanmean(mean_meps_arr(:,chan,2,intensity)) nanmean(mean_meps_arr(:,chan,1,intensity)) nanmean(mean_meps_arr(:,chan,3,intensity))],[mep_sem(chan,2,intensity),mep_sem(chan,1,intensity),mep_sem(chan,3,intensity)],'Color',col(counter),'LineWidth',2.5); hold on; 
        xticks([1,2,3])
        xlim([0.75 3.25])
        ax = gca;
        ax.YAxis.FontSize = 17;
        ax.XAxis.FontSize=17;
        box off;
        %title([muscle{chan} ',' intensities{intensity}],'FontSize',19,'FontWeight','bold')

    end
%     han=axes(fig,'visible','off'); 
%     han.YLabel.Visible='on';
%     han.Title.Visible='on';
%     ylabel(han,'Normalized MEP amplitudes','FontSize',19,'FontWeight','bold');
    ylim tight
    linkaxes([p(1) p(2)], 'y');
end

%for the average MEP figure
col=["#ab590d", "#e39952", '#810290','#DA80E4'];
sha=["square", "square", "square","square"];
%mean_meps_arr- %(subject,channel,state,intensity)
figure
%for fdi 120  %(subject,1,state,1)
plot([1, 2, 3],[nanmean(mean_meps_arr(:,1,2,1)) nanmean(mean_meps_arr(:,1,1,1)) nanmean(mean_meps_arr(:,1,3,1))], '-','Color',col{1},'LineWidth',2); hold on;
scatter([2],[nanmean(mean_meps_arr(:,1,1,1))],sha(1),'MarkerFaceColor',col{1},'MarkerEdgeColor',col{1},'SizeData',300); hold on; 
scatter([1],[nanmean(mean_meps_arr(:,1,2,1))],sha(1),'MarkerFaceColor',col{1},'MarkerEdgeColor',col{1},'SizeData',300); hold on; 
scatter([3],[nanmean(mean_meps_arr(:,1,3,1))],sha(1),'MarkerFaceColor',	col{1},'MarkerEdgeColor',col{1},'SizeData',300); hold on;
errorbar([1],[nanmean(mean_meps_arr(:,1,2,1))],[mep_sem(1,2,1)],'Color',col{1},'LineWidth',2,'CapSize',16); hold on; 
errorbar([2],[nanmean(mean_meps_arr(:,1,1,1))],[mep_sem(1,1,1)],'Color',col{1},'LineWidth',2,'CapSize',16); hold on; 
errorbar([3],[nanmean(mean_meps_arr(:,1,3,1))],[mep_sem(1,3,1)],'Color',col{1},'LineWidth',2,'CapSize',16); hold on;

%for fdi 110 %(subject,1,state,2)
plot([1, 2, 3],[nanmean(mean_meps_arr(:,1,2,2)) nanmean(mean_meps_arr(:,1,1,2)) nanmean(mean_meps_arr(:,1,3,2))], '-','Color',col{2},'LineWidth',2)
scatter([2],[nanmean(mean_meps_arr(:,1,1,2))],sha(2),'MarkerFaceColor',col{2},'MarkerEdgeColor',col{2},'MarkerFaceAlpha',1,'SizeData',300,'MarkerEdgeAlpha',1); hold on; 
scatter([1],[nanmean(mean_meps_arr(:,1,2,2))],sha(2),'MarkerFaceColor',col{2},'MarkerEdgeColor',col{2},'MarkerFaceAlpha',1,'SizeData',300,'MarkerEdgeAlpha',1); hold on; 
scatter([3],[nanmean(mean_meps_arr(:,1,3,2))],sha(2),'MarkerFaceColor',	col{2},'MarkerEdgeColor',col{2},'MarkerFaceAlpha',1,'SizeData',300,'MarkerEdgeAlpha',1); hold on;
errorbar([1],[nanmean(mean_meps_arr(:,1,2,2))],[mep_sem(1,2,2)],'Color',col{2},'LineWidth',2,'CapSize',16); hold on; 
errorbar([2],[nanmean(mean_meps_arr(:,1,1,2))],[mep_sem(1,1,2)],'Color',col{2},'LineWidth',2,'CapSize',16); hold on; 
errorbar([3],[nanmean(mean_meps_arr(:,1,3,2))],[mep_sem(1,3,2)],'Color',col{2},'LineWidth',2,'CapSize',16); hold on;

%for apb 120 %(subject,2,state,1)
plot([1 2 3],[nanmean(mean_meps_arr(:,2,2,1)) nanmean(mean_meps_arr(:,2,1,1)) nanmean(mean_meps_arr(:,2,3,1))], '-','Color',col{3},'LineWidth',2)
scatter([2],[nanmean(mean_meps_arr(:,2,1,1))],sha(3),'MarkerFaceColor',col{3},'MarkerEdgeColor',col{3},'SizeData',300); hold on; 
scatter([1],[nanmean(mean_meps_arr(:,2,2,1))],sha(3),'MarkerFaceColor',col{3},'MarkerEdgeColor',col{3},'SizeData',300); hold on; 
scatter([3],[nanmean(mean_meps_arr(:,2,3,1))],sha(3),'MarkerFaceColor',col{3},'MarkerEdgeColor',col{3},'SizeData',300); hold on;
errorbar([1],[nanmean(mean_meps_arr(:,2,2,1))],[mep_sem(2,2,1)],'Color',col{3},'LineWidth',2,'CapSize',16); hold on; 
errorbar([2],[nanmean(mean_meps_arr(:,2,1,1))],[mep_sem(2,1,1)],'Color',col{3},'LineWidth',2,'CapSize',16); hold on; 
errorbar([3],[nanmean(mean_meps_arr(:,2,3,1))],[mep_sem(2,3,1)],'Color',col{3},'LineWidth',2,'CapSize',16); hold on;

%for apb 110  %(subject,2,state,2)
plot([1 2 3],[nanmean(mean_meps_arr(:,2,2,2)) nanmean(mean_meps_arr(:,2,1,2)) nanmean(mean_meps_arr(:,2,3,2))], '-','Color',col{4},'LineWidth',2)
scatter([2],[nanmean(mean_meps_arr(:,2,1,2))],sha(4),'MarkerFaceColor',col{4},'MarkerEdgeColor',col{4},'MarkerFaceAlpha',1,'SizeData',300,'MarkerEdgeAlpha',1); hold on; 
scatter([1],[nanmean(mean_meps_arr(:,2,2,2))],sha(4),'MarkerFaceColor',col{4},'MarkerEdgeColor',col{4},'MarkerFaceAlpha',1,'SizeData',300,'MarkerEdgeAlpha',1); hold on; 
scatter([3],[nanmean(mean_meps_arr(:,2,3,2))],sha(4),'MarkerFaceColor',col{4},'MarkerEdgeColor',col{4},'MarkerFaceAlpha',1,'SizeData',300,'MarkerEdgeAlpha',1); hold on;
errorbar([1],[nanmean(mean_meps_arr(:,2,2,2))],[mep_sem(2,2,2)],'Color',col{4},'LineWidth',2,'CapSize',16); hold on; 
errorbar([2],[nanmean(mean_meps_arr(:,2,1,2))],[mep_sem(2,1,2)],'Color',col{4},'LineWidth',2,'CapSize',16); hold on; 
errorbar([3],[nanmean(mean_meps_arr(:,2,3,2))],[mep_sem(2,3,2)],'Color',col{4},'LineWidth',2,'CapSize',16); hold on;

xticks([1,2,3])
xticklabels({'weak','strong','random'});
xtickangle(0)
xlim([0.75 3.25])
ax = gca;
ax.YAxis.FontSize = 14;
ax.XAxis.FontSize=14;
ax.XAxis.FontWeight='bold';
l1 = plot(nan, nan,sha(1),'markerfacecolor',col{1},'markeredgecolor',col{1},'markersize',150);
l2 = plot(nan, nan,sha(2),'markerfacecolor',col{2},'markeredgecolor',col{2}, 'markersize',150);
l3 = plot(nan, nan,sha(3),'markerfacecolor',col{3},'markeredgecolor',col{3},'markersize',150);
l4 = plot(nan, nan,sha(4),'markerfacecolor',col{4},'markeredgecolor',col{4},'markersize',150);
lgd=legend([l1, l2, l3, l4], {'FDI, 120% RMT', 'FDI, 110% RMT', 'APB, 120% RMT', 'APB, 110% RMT'});
fontsize(lgd,13,'points')
legend boxoff

%% coefficient of variation figure 
clc
clearvars
user='';
%cd to data location
load('test_emg_data.mat');
subjname=results.subjname;
intensities={'120','110'};
states={'strong', 'weak', 'random'};
%muscle 1=fdi, 2= apb

coeff_of_var=results.coeff_of_var; %{subject,intensity}(chan,state)

for iii=1:2
    for ccc=1:2
        for sss=1:3
            for jjj=1:19
                cv(jjj,iii,ccc,sss)= coeff_of_var{jjj,iii}(ccc,sss); %(subject,intensity,channel,state)
            end
            cv_used=cv(:,iii,ccc,sss);
            mean_cv(iii,ccc,sss)= nanmean(cv_used); %mean across subjects for a intensity, state and channel
            sem_cv(iii,ccc,sss)= nanstd(cv(:,iii,ccc,sss))./sqrt(size(cv_used(~isnan(cv_used)),1)); %intensity*channel*state, cv_used and cv(:,iii,ccc,sss) are the same

        end
    end
end
  
col=["#ab590d", "#e39952", '#810290','#DA80E4']; 
muscle={'FDI', 'APB'}
counter=0;
for chan=1:2
    fig=figure;
    for intensity=1:2
        counter=counter+1;
        p(intensity)=subplot(1,2,intensity)  
        box off;
        for jjj=1:19
            plot([1,2,3],[coeff_of_var{jjj,intensity}(chan,2) coeff_of_var{jjj,intensity}(chan,1) coeff_of_var{jjj,intensity}(chan,3)],'color', col(counter) ,'linestyle',':', 'LineWidth',0.2); hold on; 
            scatter([1,2,3],[coeff_of_var{jjj,intensity}(chan,2) coeff_of_var{jjj,intensity}(chan,1) coeff_of_var{jjj,intensity}(chan,3)],'MarkerFaceColor',col(counter),'MarkerEdgeColor',col(counter),'SizeData',50,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); hold on; 
        end
        scatter([1,2,3],[mean_cv(intensity,chan,2)  mean_cv(intensity,chan,1)  mean_cv(intensity,chan,3)],'square','MarkerFaceColor',col(counter),'MarkerEdgeColor',col(counter),'SizeData',250,'LineWidth',1);
        %errorbar([1,2,3],[mean_cv(intensity,chan,2)  mean_cv(intensity,chan,1)  mean_cv(intensity,chan,3)],[sem_cv(intensity,chan,2) sem_cv(intensity,chan,1) sem_cv(intensity,chan,2)],'Color',col(counter),'LineWidth',2.5); hold on; -
        errorbar([1,2,3],[mean_cv(intensity,chan,2)  mean_cv(intensity,chan,1)  mean_cv(intensity,chan,3)],[sem_cv(intensity,chan,2) sem_cv(intensity,chan,1) sem_cv(intensity,chan,3)],'Color',col(counter),'LineWidth',2.5); hold on; 
        xticks([1,2,3])
        xlim([0.75 3.25])
        ax = gca;
        ax.YAxis.FontSize = 16;
%         ax.XAxis.FontSize=17;
        box off;
        %title([muscle{chan} ',' intensities{intensity}],'FontSize',19,'FontWeight','bold')

    end
    ylim tight
    linkaxes([p(1) p(2)], 'y');
end

col={"#ab590d", "#e39952", '#810290','#DA80E4'}; 
sha=["square", "square", "square","square"];
%for the average CV figure-- (intensity,channel,state)
figure
%for fdi 120 -- (1,1,state)
plot([1, 2, 3],[mean_cv(1,1,2) mean_cv(1,1,1) mean_cv(1,1,3)],'-','Color',col{1},'LineWidth',2); hold on;
scatter([2],[mean_cv(1,1,1)],sha(1),'MarkerFaceColor',col{1},'MarkerEdgeColor',col{1},'SizeData',300); hold on; 
scatter([1],[mean_cv(1,1,2)],sha(1),'MarkerFaceColor',col{1},'MarkerEdgeColor',col{1},'SizeData',300); hold on; 
scatter([3],[mean_cv(1,1,3)],sha(1),'MarkerFaceColor',	col{1},'MarkerEdgeColor',col{1},'SizeData',300); hold on;
errorbar([1],[mean_cv(1,1,2)],[sem_cv(1,1,2)],'Color',col{1},'LineWidth',2,'CapSize',16); hold on; 
errorbar([2],[mean_cv(1,1,1)],[sem_cv(1,1,1)],'Color',col{1},'LineWidth',2,'CapSize',16); hold on; 
errorbar([3],[mean_cv(1,1,3)],[sem_cv(1,1,3)],'Color',col{1},'LineWidth',2,'CapSize',16); hold on;

%for fdi 110-- (2,1,state)
plot([1, 2, 3],[mean_cv(2,1,2) mean_cv(2,1,1) mean_cv(2,1,3)], '-','Color',col{2},'Linewidth',2)
scatter([2],[mean_cv(2,1,1)],sha(2),'MarkerFaceColor',col{2},'MarkerEdgeColor',col{2},'MarkerFaceAlpha',1,'SizeData',300,'MarkerEdgeAlpha',1); hold on; 
scatter([1],[mean_cv(2,1,2)],sha(2),'MarkerFaceColor',col{2},'MarkerEdgeColor',col{2},'MarkerFaceAlpha',1,'SizeData',300,'MarkerEdgeAlpha',1); hold on; 
scatter([3],[mean_cv(2,1,3)],sha(2),'MarkerFaceColor',	col{2},'MarkerEdgeColor',col{2},'MarkerFaceAlpha',1,'SizeData',300,'MarkerEdgeAlpha',1); hold on;
errorbar([1],[mean_cv(2,1,2)],[sem_cv(2,1,2)],'Color',col{2},'LineWidth',2,'CapSize',16); hold on; 
errorbar([2],[mean_cv(2,1,1)],[sem_cv(2,1,1)],'Color',col{2},'LineWidth',2,'CapSize',16); hold on; 
errorbar([3],[mean_cv(2,1,3)],[sem_cv(2,1,3)],'Color',col{2},'LineWidth',2,'CapSize',16); hold on;

%for apb 120-- (1,2,state)
plot([1, 2, 3],[mean_cv(1,2,2) mean_cv(1,2,1) mean_cv(1,2,3)], '-','Color',col{3},'Linewidth',2)
scatter([2],[mean_cv(1,2,1)],sha(2),'MarkerFaceColor',col{3},'MarkerEdgeColor',col{3},'MarkerFaceAlpha',1,'SizeData',300,'MarkerEdgeAlpha',1); hold on; 
scatter([1],[mean_cv(1,2,2)],sha(2),'MarkerFaceColor',col{3},'MarkerEdgeColor',col{3},'MarkerFaceAlpha',1,'SizeData',300,'MarkerEdgeAlpha',1); hold on; 
scatter([3],[mean_cv(1,2,3)],sha(2),'MarkerFaceColor',	col{3},'MarkerEdgeColor',col{3},'MarkerFaceAlpha',1,'SizeData',300,'MarkerEdgeAlpha',1); hold on;
errorbar([1],[mean_cv(1,2,2)],[sem_cv(1,2,2)],'Color',col{3},'LineWidth',2,'CapSize',16); hold on; 
errorbar([2],[mean_cv(1,2,1)],[sem_cv(1,2,1)],'Color',col{3},'LineWidth',2,'CapSize',16); hold on; 
errorbar([3],[mean_cv(1,2,3)],[sem_cv(1,2,3)],'Color',col{3},'LineWidth',2,'CapSize',16); hold on;

%for apb 110-- (2,2,state)
plot([1, 2, 3],[mean_cv(2,2,2) mean_cv(2,2,1) mean_cv(2,2,3)], '-','Color',col{4},'Linewidth',2)
scatter([2],[mean_cv(2,2,1)],sha(2),'MarkerFaceColor',col{4},'MarkerEdgeColor',col{4},'MarkerFaceAlpha',1,'SizeData',300,'MarkerEdgeAlpha',1); hold on; 
scatter([1],[mean_cv(2,2,2)],sha(2),'MarkerFaceColor',col{4},'MarkerEdgeColor',col{4},'MarkerFaceAlpha',1,'SizeData',300,'MarkerEdgeAlpha',1); hold on; 
scatter([3],[mean_cv(2,2,3)],sha(2),'MarkerFaceColor',	col{4},'MarkerEdgeColor',col{4},'MarkerFaceAlpha',1,'SizeData',300,'MarkerEdgeAlpha',1); hold on;
errorbar([1],[mean_cv(2,2,2)],[sem_cv(2,2,2)],'Color',col{4},'LineWidth',2,'CapSize',16); hold on; 
errorbar([2],[mean_cv(2,2,1)],[sem_cv(2,2,1)],'Color',col{4},'LineWidth',2,'CapSize',16); hold on; 
errorbar([3],[mean_cv(2,2,3)],[sem_cv(2,2,3)],'Color',col{4},'LineWidth',2,'CapSize',16); hold on;

xticks([1,2,3])
xticklabels({'weak','strong','random'});
xtickangle(0)
xlim([0.75 3.25])
ax = gca;
ax.YAxis.FontSize = 14;
ax.XAxis.FontSize=14;
ax.XAxis.FontWeight='bold';
l1 = plot(nan, nan,sha(1),'markerfacecolor',col{1},'markeredgecolor',col{1},'markersize',150);
l2 = plot(nan, nan,sha(2),'markerfacecolor',col{2},'markeredgecolor',col{2}, 'markersize',150);
l3 = plot(nan, nan,sha(3),'markerfacecolor',col{3},'markeredgecolor',col{3},'markersize',150);
l4 = plot(nan, nan,sha(4),'markerfacecolor',col{4},'markeredgecolor',col{4},'markersize',150);
lgd=legend([l1, l2, l3, l4], {'FDI, 120% RMT', 'FDI, 110% RMT', 'APB, 120% RMT', 'APB, 110% RMT'});
fontsize(lgd,13,'points')
legend boxoff

