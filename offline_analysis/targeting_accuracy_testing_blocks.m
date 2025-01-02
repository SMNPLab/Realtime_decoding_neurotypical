%targeting_accuracy for realtime_decoding manuscript 
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
%% 
% get data paths, subject name, filename 
subjname={'01','02','03','04','05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19'};
states={'high', 'low'};
intensities={'120', '110'};

%% the main loop
for jjj=1:length(subjname)
    tdelay=input(['Enter technical delay for subject ' subjname{jjj} '(in ms):']);
    for state=1:length(states)
    
        for intensity=1:length(intensities)
        
            % cd to subject data folder

            if strcmpi(subjname{jjj}, '18')==1 && strcmpi(states{state}, 'high') && strcmpi(intensities{intensity}, '110')
                filename= (['sub-' subjname{jjj} '_task-testinghigh110RMTfile1_eeg']);
            elseif strcmpi(subjname{jjj}, '14')==1 && strcmpi(states{state}, 'high') && strcmpi(intensities{intensity}, '110')
                filename= (['sub-' subjname{jjj} '_task-testinghigh110RMTfile1_eeg']);
            else 
                filename= (['sub-' subjname{jjj} '_task-testing' states{state} intensities{intensity} 'RMT_eeg']);
            end

            eeg_file=[filename '.eeg'];
            
            if exist(eeg_file, 'file') == 2
              
                % load classifier
                
                load(uigetfile); % load classifier 
                freqs=[4:0.25:35];
                lda_ensemble=offline_classify.lda;
                ranks=offline_classify.ranks; % indices of features to keep
                
                csp_params=offline_classify.csp_params;
                c1=offline_classify.class1_dval_cutoff;
                c2=offline_classify.class2_dval_cutoff;

                if strcmpi(subjname{jjj}, '18')==1 && strcmpi(states{state}, 'high') && strcmpi(intensities{intensity}, '110')

                    for block=1:3
                        filename= (['sub-' subjname{jjj} '_task-testing' states{state} intensities{intensity} 'RMTfile' num2str(block) '_eeg']);
                        % define eeg filenames 
                        eeg_file=[filename '.eeg'];
                        vhdr_file=[filename '.vhdr'];
                        vmrk_file=[filename '.vmrk'];
                        
                            
                        % load eeg data
                        cfg=[];
                        cfg.dataset=eeg_file;
                        hdr=vhdr_file;
                        cfg.channel={'all','-EMG','-FDI','-APB','-EDC','-ADM','-ECR'};
                        eeg=ft_preprocessing(cfg);
                        
                        % segment eeg data (500 ms segment that ends X ms before tms pulse, where X = technical delay)
                        tdelay=total_delay;
                        cfg=[];
                        cfg.dataset=vhdr_file;
                        cfg.trialdef.eventtype='Stimulus';
                        cfg.trialdef.eventvalue='A';
                        cfg.trialdef.prestim= (tdelay/1000) + 0.500; 
                        cfg.trialdef.poststim= (-tdelay/1000); 
                        cfg=ft_definetrial(cfg);
                        trl=cfg.trl;
                        seg{block}=ft_redefinetrial(cfg,eeg);
                    end

                    cfg=[];
                    seg_eeg=ft_appenddata(cfg,seg{1},seg{2},seg{3});



                elseif strcmpi(subjname{jjj}, '14')==1 && strcmpi(states{state}, 'high') && strcmpi(intensities{intensity}, '110')

                    fname_blocks={'sub-14_task-testinghigh110RMTfile1_eeg','sub-14_task-testinghigh110RMTfile2_eeg' }

                    for blocks=1:2

                        % define eeg filenames 
                        eeg_file=[fname_blocks{blocks} '.eeg'];
                        vhdr_file=[fname_blocks{blocks} '.vhdr'];
                        vmrk_file=[fname_blocks{blocks} '.vmrk'];
                            
                        % load eeg data
                        cfg=[];
                        cfg.dataset=eeg_file;
                        hdr=vhdr_file;
                        cfg.channel={'all','-EMG','-FDI','-APB','-EDC','-ADM','-ECR'};
                        eeg=ft_preprocessing(cfg);
                        
                        % segment eeg data (500 ms segment that ends X ms before tms pulse, where X = technical delay)
                        tdelay=total_delay;
                        cfg=[];
                        cfg.dataset=vhdr_file;
                        cfg.trialdef.eventtype='Stimulus';
                        cfg.trialdef.eventvalue='A';
                        cfg.trialdef.prestim= (tdelay/1000) + 0.500; 
                        cfg.trialdef.poststim= (-tdelay/1000); 
                        cfg=ft_definetrial(cfg);
                        trl=cfg.trl;
                        seg{blocks}=ft_redefinetrial(cfg,eeg);
                    end
                    seg_eeg=ft_appenddata(cfg,seg{1},seg{2});


                else 

                    % define eeg filenames 
                    eeg_file=[filename '.eeg'];
                    vhdr_file=[filename '.vhdr'];
                    vmrk_file=[filename '.vmrk'];
                    
                        
                    % load eeg data
                    cfg=[];
                    cfg.dataset=eeg_file;
                    hdr=vhdr_file;
                    cfg.channel={'all','-EMG','-FDI','-APB','-EDC','-ADM','-ECR'};
                    eeg=ft_preprocessing(cfg);
                    
                    % segment eeg data (500 ms segment that ends X ms before tms pulse, where X = technical delay)
                    tdelay=total_delay;
                    cfg=[];
                    cfg.dataset=vhdr_file;
                    cfg.trialdef.eventtype='Stimulus';
                    cfg.trialdef.eventvalue='A';
                    cfg.trialdef.prestim= (tdelay/1000) + 0.500; 
                    cfg.trialdef.poststim= (-tdelay/1000); 
                    cfg=ft_definetrial(cfg);
                    trl=cfg.trl;
                    seg_eeg=ft_redefinetrial(cfg,eeg);

                end

                % rereference eeg data
                cfg=[];
                cfg.reref='yes';
                cfg.refchannel='all';
                reref=ft_preprocessing(cfg,seg_eeg);
                
                % reformat and downsample to 1000 Hz 
                for trl=1:length(reref.trial)
                    for chan=1:length(reref.label)
                        eeg_mat(:,trl,chan)=reref.trial{trl}(chan,:);
                        eeg_down(:,trl,chan)=downsample(eeg_mat(:,trl,chan),5); % downsample to 1000hz  samples*trials*channels
                        eeg_down(:,trl,chan)=detrend(eeg_down(:,trl,chan),'constant'); % added by sjh on 5/6
                        eeg_down(:,trl,chan)=detrend(eeg_down(:,trl,chan),'linear'); % added by sjh on 5/6
                    end
                end
                csp_chunk=permute(eeg_down,[2,3,1]); %samples*trials*channels ---> trials*channels*samples
    
                %add conditions with dval criteria that are 0
                if strcmpi(subjname{jjj},'04')==1 && strcmpi(states{state},'high')==1
                    dval=0;
                elseif strcmpi(subjname{jjj},'15')==1 && strcmpi(states{state},'low')==1 && strcmpi(intensities{intensity},'110')==1
                    dval=0;
                elseif strcmpi(subjname{jjj},'17')==1 && strcmpi(states{state},'high')==1
                    dval=0;
                else
                    dval=1;
                end 

    
                for k=1:5
                    % get csp features for realtime collected data (i.e., a type of test set)
                    csp_params{k}.is_train_set=0;
                    [~, segrt, ~] = mv_preprocess_csp(csp_params{k}, csp_chunk,[]);
                    
                    % compute psd from csp features
                    [~, realtest_features]= calc_ps_rt(segrt);
                    
                    rng(42);
                    tmp_features=realtest_features(:,ranks{k});
                    [clabel(1:size(csp_chunk,1),k)]=test_lda_realtime(lda_ensemble{k},tmp_features,c1(k),c2(k),dval);
                end
                clabel(clabel==0)=NaN; 
                
                if strcmpi(states{state},'high')
                    targeted_states=repmat(2,numel(clabel),1);
                else
                    targeted_states=repmat(1,numel(clabel),1);
                end
                
                for k=1:size(clabel,1)
                    if nanmean(clabel(k,:)) == 1.5
                       pred_state(k)=NaN;
                    else
                       pred_state(k)=round(nanmean(clabel(k,:))); 
                    end
                end
                accurate=0;
                inaccurate=0;
                
                for t=1:size(clabel,1)
                    if pred_state(t)==targeted_states(t)   
                        accurate=accurate+1;
                    elseif pred_state(t)~=targeted_states(t)
                        inaccurate=inaccurate+1;
                    end
                end
                
                correct=sum(accurate)/numel(targeted_states(1:size(clabel,1))) * 100;
                accuracy.accuracy(jjj,state,intensity)=correct;
                accuracy.dval(jjj,state,intensity)=dval;
                accuracy.state(jjj,state,intensity)=state;
                accuracy.intensity(jjj,state,intensity)=intensity;
                clearvars -except subjname states accuracy intensities uniquenum user state jjj

            else
                accuracy.accuracy(jjj,state,intensity)=NaN;
                accuracy.dval(jjj,state,intensity)=NaN;
                accuracy.state(jjj,state,intensity)=state;
                accuracy.intensity(jjj,state,intensity)=intensity;
                clearvars -except subjname states accuracy intensities uniquenum user state jjj
  
            end
        end
       
    end
   
end
        
%path to folder to save data
save(['targeting_accuracy_updated.mat'],'accuracy','-v7.3'); 


%% calculate sem for accuracy values
clc
clearvars
user='';
%cd to data location
acc=load('targeting_accuracy_updated.mat').accuracy.accuracy; % loading accuracy from previously saved array 

 
for state=1:2
    for intensity=1:2
        acc_used=acc(:,state,intensity);
        acc_sem(state,intensity)= nanstd(acc_used)./sqrt(size(acc_used(~isnan(acc_used)),1)); %SJH: need to recalculate to account for actual number of data points rather than max number of data points 
        acc_mean(state,intensity)= nanmean(acc_used);
        m(state,intensity)=state;
    end
end

accuracy_sem.mean=acc_mean; % state*intensity (state=1 is high and state=2 is low), intensity=1 is high and intensity=2 is low
accuracy_sem.sem=acc_sem;
%path to folder to save data
save(['accuracy_sem_updated.mat'],'accuracy_sem','-v7.3'); 


%% load accuracies and make a plot out of it 
clc
clearvars
user='';
%cd to data location
acc=load('targeting_accuracy_updated.mat').accuracy.accuracy; %subject*state*intensity, states={'high', 'low'}; intensities={'120', '110'};
acc_sem=load('accuracy_sem_updated.mat').accuracy_sem.sem; %state*intensity (state 1=high, 2=low, intensity 1=120, 2=110)

figure;
for subj=1:size(acc,1) %subjects
    plot([1,2], [acc(subj,2,1) acc(subj,1,1)],'color', '#ab590d','linestyle',':', 'LineWidth',0.25); hold on; %plots low,120 and high 120
    plot([2.5 3.5], [acc(subj,2,2) acc(subj,1,2)],'color', '#e39952' ,'linestyle',':', 'LineWidth',0.25); hold on; %plots low,110 and high 110
    scatter([1,2], [acc(subj,2,1) acc(subj,1,1)],'markerfacecolor', '#ab590d', 'markeredgecolor','#ab590d', 'SizeData',80,'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3); hold on;
    scatter([2.5 3.5], [acc(subj,2,2) acc(subj,1,2)],'markerfacecolor', '#e39952', 'markeredgecolor','#e39952','SizeData',80,'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3); hold on;
end
scatter([1 2],[nanmean(acc(:,2,1)) nanmean(acc(:,1,1))],'square','MarkerFaceColor','#ab590d','MarkerEdgeColor','#ab590d','SizeData',400,'LineWidth',1) %plotting the subject wise mean for low120 and high 120
scatter([2.5 3.5],[nanmean(acc(:,2,2)) nanmean(acc(:,1,2))],'square','MarkerFaceColor','#e39952','MarkerEdgeColor','#e39952','SizeData',400,'LineWidth',1)  %plotting the subject wise mean for low110 and high 110
errorbar([1,2],[nanmean(acc(:,2,1))  nanmean(acc(:,1,1))],[acc_sem(2,1) acc_sem(1,1)],'Color','k','LineWidth',1.5,'CapSize',15); hold on; %for 120% RMT
errorbar([2.5,3.5],[nanmean(acc(:,2,2))  nanmean(acc(:,1,2)) ],[acc_sem(2,2) acc_sem(1,2)],'color','k','LineWidth',1.5,'CapSize',15); hold on; %for 110% RMT
xticks([1,2,2.5,3.5])
xticklabels({'weak','strong','weak','strong'});
ax = gca;
ax.YAxis.FontSize = 14;
ax.XAxis.FontSize=14;
xlim([0.75 3.75])
ylim([0 105]);
xlabel('CST states','FontSize',17,'FontWeight','bold')
ylabel('state-targeting accuracy (%)','FontSize',17,'FontWeight','bold');
L1 = plot(nan, nan,'square','markerfacecolor','#ab590d','markeredgecolor','#ab590d','markersize',150);
L2 = plot(nan, nan,'square','markerfacecolor','#e39952','markeredgecolor','#e39952','markersize',150);
lgd=legend([L1, L2], {'120% RMT', '110% RMT'});
fontsize(lgd,13,'points')
legend boxoff
