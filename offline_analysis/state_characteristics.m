 %% load the defaults

% use existing realtime lda function, but run twice and pull pred high and pred low depending on dvals for those subjects who need
% diff dvalues for states at 120%RMT
clc
clearvars
user='';
%path to fieldtrip
%path to mvpa light
%path to custom written functions
%path to subject data folder
ft_defaults;
startup_MVPA_Light;

%% subject list
subjname={'01','02','03','04','05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19'};


%% get the chunks of eeg data for the resting data sliding window every x milliseconds
% 1min- 60 sec, 5min- 300 sec
% 1 sliding window - 0.01 sec, number of sliding windows- 300/0.01-- 30000,
% example for
%cd to data location
load(['state_change_data_for_fig.mat'])
predicted_states=results.pred_state_dval1; %keeps updating
predicted_states_dval0=results.pred_state_dval0; 

win_size= 50; %ms 


for jjj=1:numel(subjname)
    
    %cd to data location
    disp(subjname{jjj})

    % load classifier
    disp(['Please select the offline classify file for subject ' subjname{jjj}])
    load(uigetfile); % load classifier for that subject
    freqs=[4:0.25:35];
    lda_ensemble=offline_classify.lda; 
    ranks=offline_classify.ranks; % indices of features to keep 
    csp_params=offline_classify.csp_params;
    c1=offline_classify.class1_dval_cutoff;
    c2=offline_classify.class2_dval_cutoff;

    if exist(['sub-' subjname{jjj} '_task-resting_eeg.eeg'], 'file' )
        eeg_file= (['sub-' subjname{jjj} '_task-resting_eeg.eeg']);
        vhdr_file= (['sub-' subjname{jjj} '_task-resting_eeg.vhdr']);
        vmrk_file= (['sub-' subjname{jjj} '_task-resting_eeg.vmrk']);

        % load eeg data
        cfg=[];
        cfg.dataset=eeg_file;
        hdr=vhdr_file;
        cfg.channel={'all','-EMG' ,'-FDI', '-APB', '-ADM', '-ECR', '-EDC'};
        eeg=ft_preprocessing(cfg);

        % rereference eeg data
        cfg=[];
        cfg.reref='yes';
        cfg.refchannel='all';
        reref=ft_preprocessing(cfg,eeg);
 
        %converting eeg trial data from a cell array to a matrix
        eeg_mat=cell2mat(reref.trial); % chan x sample 
        threemins=60*3*5000; % s x min x samples  
        eeg_mat=eeg_mat(:,1:threemins);  %taking only the first 3 min to standardize across people

        % reformat and downsample to 1000 Hz
        for chan=1:size(eeg_mat,1)
            eeg_down(chan,:)=downsample(eeg_mat(chan,:),5); % downsample to 1000hz 
        end

        clear eeg_mat

        for trl=1:(fix(size(eeg_down,2)/(win_size))-9)

            % make the seg
            seg_num=[((trl-1)*win_size)+1 (((trl-1)*win_size)+500)]; % identifies sample numbers 500 ms segments in steps of 50 ms 

            seg_eeg(trl,:,:)= eeg_down(:, seg_num(1):seg_num(2)); %chan*sample -> trial*channels*samples

        end 

        %converting from trial*channels*samples to samples*trials*channels
        seg_eeg= permute(seg_eeg,[3,1,2]);

        % detrend / demean seg 
        seg_detrend=detrend(seg_eeg,'constant'); %demean
        seg_detrend=detrend(seg_detrend,'linear'); %linearly detrend

        results.seg_eeg{jjj}=seg_detrend;

        clear seg_eeg reref eeg_down

        %permuting the data for csp processing
        csp_chunk=permute(seg_detrend,[2,3,1]); %trl*channels*500

        dval=1; %ran with dval=0 for s006 and s019 for high condition

        for k=1:5  

            % get csp features for realtime collected data (i.e., a type of test set)
            csp_params{k}.is_train_set=0;
            [~, segrt, ~] = mv_preprocess_csp(csp_params{k}, csp_chunk,[]);

            % compute psd from csp features
            [~, realtest_features]= calc_ps_rt(segrt); % SJH: confirm that this function doesnt just do one segment at a time 

            rng(42);
            tmp_features=realtest_features(:,ranks{k});

            % run this twice, once for each dval condition, selecting
            % predlabels for the relevant dval condition for relevant
            % subjects 
            [clabel(:,k)]=test_lda_realtime(lda_ensemble{k},tmp_features,c1(k),c2(k),dval);

        end
        clabel(clabel==0)=NaN;

        for trl=1:size(clabel,1)
            if nanmean(clabel(trl,:)) == 1.5
                pred_state(trl)=NaN;
            else
                pred_state(trl)=round(nanmean(clabel(trl,:))); 
            end
        end

        results.predicted_states{jjj}=pred_state;
        results.dval_high{jjj}= 1;
        results.dval_low{jjj}= 1;
        results.sliding_window_length{jjj}= win_size;
        results.subjname{jjj}= subjname{jjj};

      % adding in data previously analyzed 
        pred_state=predicted_states{jjj}; % pulls predictions for subjects where dval =1 for certain states
        pred_state_dval0=predicted_states_dval0{jjj}; % pulls predictions for subjects where dval =0 for certain states

        %pulls high states depending upon dval and calculates their number
        %and proportion
        if strcmpi(subjname{jjj},'04')==1 | strcmpi(subjname{jjj},'17')==1 
           results.instances(jjj,1)= numel(find(pred_state_dval0==2)); %number of high predictions for dval=0
           results.hs_freq(jjj)= numel(find(pred_state_dval0==2))./ numel(pred_state_dval0); %proportion of high predictions for dval=0
           high_states_diff=diff(find(pred_state_dval0==2)); % finds temporal difference (i.e., number of windows) between neighbouring high states for dval=0
        else 
           results.instances(jjj,1)= numel(find(pred_state==2));  %number of high predictions for dval=1
           results.hs_freq(jjj)= numel(find(pred_state==2))./ numel(pred_state); %proportion of high states for dval =1 
           high_states_diff=diff(find(pred_state==2)); % finds temporal difference (i.e., number of windows) between neighbouring high states for dval =1
        end

        %number of low windows
        results.instances(jjj,2)= numel(find(pred_state==1)); % number of low predictions
        results.instances(jjj,3)= numel(pred_state); %number of all predictions
        results.instances(jjj,4)= numel(find(isnan(pred_state))); % number of underconfident predictions
        
        %calculate the proportions of low and underconfident predictions
        results.ls_freq(jjj)= numel(find(pred_state==1))./ numel(pred_state); %proportion of low predictions
        results.us_freq(jjj)= numel(find(isnan(pred_state)))./ numel(pred_state); %proportion of underconfident predictions

        %FOR STATE DURATION %

        %for finding high state duration 
        hs_dur=[];
        dur=[];
        for hs=1:numel(high_states_diff)
            if  high_states_diff(hs) ==1 %if there are consecutive high states, then keep adding 1 to "dur" list.
                dur=[dur 1]
                if hs == numel(high_states_diff) %if it the end of the high_states list, then add the number of high states along with their duration to the list hs_dur
                    hs_dur=[hs_dur (500 + sum(dur)*50) ./1000]; 
                end
            elseif high_states_diff(hs) ~=1 %if the next state in high_states_diff is not a high state, then add the duration of the most recent high states to hs_dur
                 hs_dur=[hs_dur (500 + sum(dur)*50) ./1000];  
                 dur=[] %dur is emptied at the end of high state
            end
        end
        
        %for low_states duration-- the flow of logic is same as for high
        %state duration code
        low_states_diff=diff(find(pred_state==1));
        ls_dur=[];
        dur=[];
        for ls=1:numel(low_states_diff)
            if  low_states_diff(ls) ==1
                dur=[dur 1]
                if ls == numel(low_states_diff)
                    ls_dur=[ls_dur (500 + sum(dur)*50) ./1000] 
                end
            elseif low_states_diff(ls) ~=1 
                 ls_dur=[ls_dur (500 + sum(dur)*50) ./1000]  
                 dur=[]
            end
        end

        %for unconfident states- same flow of logic as unconf_states
        unconf_states_diff=diff(find(isnan(pred_state)));
        us_dur=[];
        dur=[];
        for us=1:numel(unconf_states_diff)
            if  unconf_states_diff(us) ==1
                dur=[dur 1]
                if us == numel(unconf_states_diff)
                    us_dur=[ls_dur (500 + sum(dur)*50) ./1000] 
                end
            elseif unconf_states_diff(us) ~=1 
                 us_dur=[us_dur (500 + sum(dur)*50) ./1000];
                 dur=[]
            end
        end
    
        %time between continuous states
        high_states_diff=high_states_diff-1;
        hs_time= high_states_diff(find(high_states_diff>0)); %finding the number of windows between consecutive same states
        hs_time=hs_time*0.05; % converting the number of windows between consecutive states to milliseconds
        hs_time_mean=nanmean(hs_time);
    
        low_states_diff=low_states_diff-1;
        ls_time= low_states_diff(find(low_states_diff>0)); %finding the number of windows between consecutive same states
        ls_time=ls_time*0.05; % converting the number of windows between consecutive states to milliseconds
        ls_time_mean=nanmean(ls_time);
    
        unconf_states_diff=unconf_states_diff-1;
        us_time= unconf_states_diff(find(unconf_states_diff>0)); %finding the number of windows between consecutive same states
        us_time=us_time*0.05; % converting the number of windows between consecutive states to milliseconds
        us_time_mean=nanmean(us_time);
    
        results.high_duration{jjj}=hs_dur;
        results.high_dur_mean(jjj)=nanmean(hs_dur);
        results.low_duration{jjj}=ls_dur;
        results.low_dur_mean(jjj)=nanmean(ls_dur);
        results.unconfident_duration{jjj}=us_dur;
        results.unconf_dur_mean(jjj)=nanmean(us_dur);
        results.hs_time{jjj}=hs_time;
        results.ls_time{jjj}=ls_time;
        results.us_time{jjj}=us_time;
        results.hs_time_mean{jjj}=hs_time_mean;
        results.ls_time_mean{jjj}=ls_time_mean;
        results.us_time_mean{jjj}=us_time_mean;
        results.subjname=subjname; 
        results.uniquenum=uniquenum;
        results.pred_state_dval1{jjj}=pred_state;
        results.pred_state_dval0{jjj}=pred_state_dval0;
        disp(subjname{jjj})
        
        %cd to data location
        save(['state_change_data_for_fig_v2.mat'],'results','-v7.3');  
    
        clearvars -except subjname uniquenum results win_size dval user predicted_states predicted_states_dval0
    end 
end
    

%% plotting the state duration-
% if no high or low states detected, use 0  
clc
clearvars
user='';
%add path to data location
load('state_change_data_for_fig_v2.mat')
% 
%load the arrays to be plotted
hs_dur=results.high_dur_mean;
ls_dur=results.low_dur_mean; 
us_dur=results.unconf_dur_mean; 

%replace the nans in this array as 0 and then calculate the mean and sem
hs_dur(isnan(hs_dur))=0;
ls_dur(isnan(ls_dur))=0;
us_dur(isnan(us_dur))=0;

sem_dur=zeros(3,1);
sem_dur(1)= std(ls_dur)./sqrt(numel(ls_dur));
sem_dur(2)= std(hs_dur)./sqrt(numel(hs_dur)); 
sem_dur(3)= std(us_dur)./sqrt(numel(us_dur)); 

col=["#ab590d"];
figure;

for jjj=1:numel(hs_dur)
    scatter([2],[hs_dur(jjj)],'MarkerFaceColor',col,'MarkerEdgeColor',col, 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3,'SizeData',70); hold on; 
    scatter([1],[ls_dur(jjj)],'MarkerFaceColor',col,'MarkerEdgeColor',col,'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3,'SizeData',70); hold on; 
    scatter([3],[us_dur(jjj)],'MarkerFaceColor',col,'MarkerEdgeColor',col,'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3,'SizeData',70); hold on; 
    plot([1,2,3],[ls_dur(jjj) hs_dur(jjj) us_dur(jjj)],'color', "#c79c73" ,'linestyle',':', 'LineWidth',0.35); hold on; 
end

scatter([2],[nanmean(hs_dur)],'square','MarkerFaceColor',col,'MarkerEdgeColor',col,'SizeData',600); hold on; %technically should be mean instead of nanmean, but nans removed lines 317-319
scatter([1],[nanmean(ls_dur)],'square','MarkerFaceColor',col,'MarkerEdgeColor',col,'SizeData',600); hold on; 
scatter([3],[nanmean(us_dur)],'square','MarkerFaceColor',col,'MarkerEdgeColor',col,'SizeData',600); hold on; 
errorbar([1,2,3],[nanmean(ls_dur) nanmean(hs_dur) nanmean(us_dur)], sem_dur,'k','LineWidth',2,'CapSize', 19); hold on; 
xticks([1,2,3])
xticklabels({'weak', 'strong', 'unconfident'});
ax = gca;
ax.YAxis.FontSize = 17;
ax.XAxis.FontSize=17;
% ax.XAxis.FontWeight='bold';
xlim([0.75 3.25])
ylim([-0.1 23]);
% breakplot([1,2,3],[nanmean(ls_dur) nanmean(hs_dur) nanmean(us_dur)],4,22,'Line')
% ylabel('average duration of predicted states (in seconds)','FontSize',17,'FontWeight','bold');
breakyaxis([3.45 22]);



%% plot abundance of states
clc
clearvars
user='';
%add path to data location
load('state_change_data_for_fig_v2.mat')
% 
%load the arrays to be plotted
hs_freq=results.hs_freq*100; %mean- 33.63, sem- 5.9
ls_freq=results.ls_freq*100; %mean- 48.4 , sem- 6.5 
us_freq=results.us_freq*100; %mean- 17.8 , sem- 2.6

sem_freq(1)= std(ls_freq)./sqrt(19) %changed by uuk on 7/31/2024 from nanstd to std-- should not work if nans are not converted
sem_freq(2)= std(hs_freq)./sqrt(19)
sem_freq(3)= std(us_freq)./sqrt(19)

figure;
col=["#ab590d"];
for jjj=1:numel(hs_freq)
    scatter([1],[ls_freq(jjj)],'MarkerFaceColor',col,'MarkerEdgeColor',col, 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',0.3,'SizeData',70); hold on; 
    scatter([2],[hs_freq(jjj)],'MarkerFaceColor',col,'MarkerEdgeColor',col,'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3,'SizeData',70); hold on; 
    scatter([3],[us_freq(jjj)],'MarkerFaceColor',col,'MarkerEdgeColor',col,'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3,'SizeData',70); hold on; 
    plot([1,2,3],[ls_freq(jjj) hs_freq(jjj) us_freq(jjj)],'color', "#c79c73",'Linestyle',":",'LineWidth',0.35); hold on; 
end

scatter([1],[nanmean(ls_freq)],'square','MarkerFaceColor',col,'MarkerEdgeColor',col,'SizeData',600); hold on; %tehcnically nanmean should not be used, but this array contains no nans
scatter([2],[nanmean(hs_freq)],'square','MarkerFaceColor',col,'MarkerEdgeColor',col,'SizeData',600); hold on; 
scatter([3],[nanmean(us_freq)],'square','MarkerFaceColor',col,'MarkerEdgeColor',col,'SizeData',600); hold on; 
errorbar([1,2,3],[nanmean(ls_freq) nanmean(hs_freq) nanmean(us_freq)],sem_freq,'k','LineWidth',2,'CapSize', 19); hold on; 
xticks([1,2,3])
xticklabels({'weak', 'strong', 'unconfident'});
%xlabel("predicted CST states",'FontWeight','bold');
ax = gca;
ax.YAxis.FontSize = 17;
ax.XAxis.FontSize=17;
xlim([0.75 3.25])
ylim tight;

%ylabel('% of time spent in a state ','FontSize',17,'FontWeight','bold');

%% plot the average duration between two consecutive states- 
clc
clearvars
user='';
%add path to data location
load('state_change_data_for_fig_v2.mat') % replaced this with v2 after code review check and change on 7/31/2024
col=["#ab590d"];

%load the arrays to be plotted
hs_time=cell2mat(results.hs_time_mean); %v2 mean- 0.93  , sem- 0.16
ls_time=cell2mat(results.ls_time_mean); %v2 mean- 1.00 , sem- 0.45
us_time=cell2mat(results.us_time_mean); %v2 mean- 1.97  , sem- 1.28

%changed by uuk on 7/31/2024 from nanstd to std-- should not work if nay nans are not converted
sem_time(1)= nanstd(ls_time)./sqrt(size(ls_time(~isnan(ls_time)),2)) %sem for time between low states, size(array,2), ls_time is 1x19
sem_time(2)= nanstd(hs_time)./sqrt(size(hs_time(~isnan(hs_time)),2)) %for time between high states
sem_time(3)= nanstd(us_time)./sqrt(size(us_time(~isnan(us_time)),2)) %for time between underconfident states

figure;
for jjj=1:numel(hs_time)
    scatter([1],[ls_time(jjj)],'MarkerFaceColor',col,'MarkerEdgeColor',col, 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3,'SizeData',70); hold on; 
    scatter([2],[hs_time(jjj)],'MarkerFaceColor',col,'MarkerEdgeColor',col,'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3,'SizeData',70); hold on; 
    scatter([3],[us_time(jjj)],'MarkerFaceColor',col,'MarkerEdgeColor',col,'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3,'SizeData',70); hold on; 
end

% 4- us_time is 24.99
% 15- ls_time is 8.55
for jjj=setdiff([1:19],[4,15])
    plot([1,2,3],[ls_time(jjj) hs_time(jjj) us_time(jjj)],'color', "#c79c73" ,'linestyle',':', 'LineWidth',0.35); hold on; 
end

scatter([1],[nanmean(ls_time)],'square','MarkerFaceColor',col,'MarkerEdgeColor',col,'SizeData',600); hold on; 
scatter([2],[nanmean(hs_time)],'square','MarkerFaceColor',col,'MarkerEdgeColor',col,'SizeData',600); hold on; 
scatter([3],[nanmean(us_time)],'square','MarkerFaceColor',col,'MarkerEdgeColor',col,'SizeData',600); hold on; 
errorbar([1,2,3],[nanmean(ls_time) nanmean(hs_time) nanmean(us_time)],sem_time,'k','LineWidth',2,'CapSize', 19); hold on; 
xticks([1,2,3])
xticklabels({'weak', 'strong', 'under-confident'});
ax = gca;
ax.YAxis.FontSize = 17;
ax.XAxis.FontSize=17;
xlim([0.75 3.25])
xlabel("predicted CST states", 'FontSize',17,'FontWeight','bold');
ylabel('Time between consecutive states (in seconds) ','FontSize',17,'FontWeight','bold');
ylim([0,9])
breakyaxis([3.4 8]);
