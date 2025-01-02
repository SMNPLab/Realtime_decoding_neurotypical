
clc; 
clearvars;
user='';
%path to fieldtrip
%path to mvpa light
%path to custom written functions
%path to subject data folder
ft_defaults;
startup_MVPA_Light;

%load the subjname and the uniquenum arrays
subjname={'01','02','03','04','05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19'};

states={'high','low','random'};
intensities={'120', '110'};

%% obtaining the difference between power spectral feastures during classifier predicted strong and weak states 
for jjj=1:length(subjname)
    % cd to subject folder
    tdelay=input(['Enter technical for subject ' subjname{jjj} ' (in ms) : ']);
    for intensity=1:2
        for state=1:3

            if strcmpi(subjname{jjj}, '18')==1 && strcmpi(states{state}, 'high') && strcmpi(intensities{intensity}, '110')
                filename= (['sub-' subjname{jjj} '_task-testinghigh110RMTfile1_eeg']);
            elseif strcmpi(subjname{jjj}, '14')==1 && strcmpi(states{state}, 'high') && strcmpi(intensities{intensity}, '110')
                filename= (['sub-' subjname{jjj} '_task-testinghigh110RMTfile1_eeg']);
            else 
                filename= (['sub-' subjname{jjj} '_task-testing' states{state} intensities{intensity} 'RMT_eeg']);
            end

        
            eeg_file=[filename '.eeg'];

            if exist(eeg_file, 'file') == 2 % if file exists 
             
                
                if strcmpi(subjname{jjj}, '18')==1 && strcmpi(states{state}, 'high') && strcmpi(intensities{intensity}, '110') % special case for Subj 20 

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

                elseif strcmpi(subjname{jjj}, '14')==1 && strcmpi(states{state}, 'high') && strcmpi(intensities{intensity}, '110') % special case 
                    
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

                else % all other cases 

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
                        eeg_down(:,trl,chan)=detrend(eeg_down(:,trl,chan),'constant'); 
                        eeg_down(:,trl,chan)=detrend(eeg_down(:,trl,chan),'linear'); 
                    end
                end
                eeg_struc{jjj,intensity,state}= eeg_down;
            else 
                eeg_struc{jjj,intensity,state}= NaN(500,25,62); % if no file available 
            end
        end
    end

    %combine statewise data from both intensities
    high_eeg= cat(2,eeg_struc{jjj,1,1},eeg_struc{jjj,2,1});
    low_eeg= cat(2,eeg_struc{jjj,1,2},eeg_struc{jjj,2,2});
    random_eeg= cat(2,eeg_struc{jjj,1,3},eeg_struc{jjj,2,3});

    %compute the power spectral values and perform log here-- for high, low and random
    %states
    for chan=1:size(high_eeg,3)
        for trl=1:size(high_eeg,2)
            high_pow{jjj}(:,trl,chan)=log(pwelch(high_eeg(:,trl,chan),500,[],freqs,1000)); % freq x trl x chan x time
        end
    end
    
    %for low    
    for chan=1:size(low_eeg,3)
        for trl=1:size(low_eeg,2)
            low_pow{jjj}(:,trl,chan)=log(pwelch(low_eeg(:,trl,chan),500,[],freqs,1000)); % freq x trl x chan x time
        end
    end

    % for random
    for chan=1:size(random_eeg,3)
        for trl=1:size(random_eeg,2)
            random_pow{jjj}(:,trl,chan)=log(pwelch(random_eeg(:,trl,chan),500,[],freqs,1000)); % freq x trl x chan x time
        end
    end

    %contrast high/low univariate comparisons - these dont end up getting
    %used 
    hilo_pow{jjj}=nanmean(high_pow{jjj},2)./nanmean(low_pow{jjj},2);
    hirand_pow{jjj}=nanmean(high_pow{jjj},2)./nanmean(random_pow{jjj},2);

    results.eeg_struc=eeg_struc; 
    results.high_pow=high_pow;
    results.low_pow=low_pow; 
    results.random_pow=random_pow; 
    results.hilo_pow=hilo_pow; % ratios
    results.hirand_pow=hirand_pow; % ratios 
    %cd to data location
    save(['univariate.mat'],'results','-v7.3');

    clearvars -except eeg_struc user freqs subjname uniquenum states intensities hilo_pow hirand_pow results high_pow low_pow random_pow

end

%% combine the univariate data into separate frequency bands to be used for topoplots
clc
clearvars
user='';
%cd to data location
load(['univariate.mat']);
high_pow=results.high_pow; %{jjj}(freq*trials*channels)
low_pow=results.low_pow;

%compute average high and low power for every person- this is same as
%before and picked from the univariate analysis put in the initially submitted
%manuscript
for jjj=1:size(high_pow,2)
    uni(jjj,:,:)=squeeze(mean(high_pow{jjj},2))-squeeze(mean(low_pow{jjj},2)); %obtain the univariate contrasts for each frequency and channel
end

%obtain the frequency band specific univariate data
freqs=[4:0.25:35];
theta=[find(freqs==4):find(freqs==7.75)];
alpha1=[find(freqs==8):find(freqs==10)];
alpha2=[find(freqs==10.25):find(freqs==13)]; 
beta1=[find(freqs==13.25):find(freqs==20)];
beta2=[find(freqs==20.25):find(freqs==35)];


univariate.theta_pow=uni(:,theta,:); %jjj*freq*chan, here obtaining mean over the frequency band
univariate.alpha1_pow=uni(:,alpha1,:);
univariate.alpha2_pow=uni(:,alpha2,:);
univariate.beta1_pow=uni(:,beta1,:);
univariate.beta2_pow=uni(:,beta2,:);

%mean across the frequency band
univariate.theta_pow_mean=squeeze(mean(uni(:,theta,:),2)); %jjj*freq*chan, here obtaining mean over the frequency band
univariate.alpha1_pow_mean=squeeze(mean(uni(:,alpha1,:),2));
univariate.alpha2_pow_mean=squeeze(mean(uni(:,alpha2,:),2));
univariate.beta1_pow_mean=squeeze(mean(uni(:,beta1,:),2));
univariate.beta2_pow_mean=squeeze(mean(uni(:,beta2,:),2));
univariate.univariate_all_freqs=uni;
univariate.channel_labels=results.chan_labels;

%cd to saving location
save(['univariate_topo.mat'],'univariate','-v7.3');


%% main plotting starts
%load univariate data for topoplots
clc
clearvars
user='';
%cd to data location
load(['univariate_topo.mat']);

%load chan_labels
chan_labels=univariate.channel_labels;

data_to_plot(:,:,1) = univariate.theta_pow_mean; 
data_to_plot(:,:,2) = univariate.alpha1_pow_mean;
data_to_plot(:,:,3) = univariate.alpha2_pow_mean;
data_to_plot(:,:,4) = univariate.beta1_pow_mean;
data_to_plot(:,:,5) = univariate.beta2_pow_mean; 

titles = {'Theta (4-8 Hz)', 'Alpha1 (8-10 Hz)', 'Alpha2 (10-13 Hz)', 'Beta1 (13-20 Hz)', 'Beta2 (20-35 Hz)'};


for jjj=1:19
    f = figure;
    f.Position = [0 0 2400 450];
    f.WindowState = 'maximized';   
    
    for col = 1:5
        subplot_col(col)=subplot(1,5,col);
        data = [];
        data.label = chan_labels; 
        data.dimord = 'chan'; 
        data.time = 0; 
        data.avg=data_to_plot(jjj,:,col);

        cfg = [];
        cfg.layout = 'elec1005.lay';  
        cfg.marker = 'on'; 
        cfg.gridscale = 300; 
        cfg.zlim = [min(data_to_plot(jjj,:,col)) max(data_to_plot(jjj,:,col))]; 
        cfg.interactive = 'yes';
        cfg.comment = 'no';
        
        % Plot
        format bank
        ft_topoplotER(cfg, data); c(col)=colorbar('FontSize', 12); c(col).Location='southoutside';
        c(col).Label.String = 'power(ln(\muV^2))';  c(col).Label.FontWeight='bold';
        c(col).Label.FontSize=15; 
        title(titles{col}, 'FontSize', 18, 'FontWeight','bold')
        c(col).TicksMode="manual";
        c(col).Ticks = round(linspace(min(data_to_plot(jjj,:,col)),max(data_to_plot(jjj,:,col)),4),2);
        c(col).TickLabels =num2cell(round(linspace(min(data_to_plot(jjj,:,col)),max(data_to_plot(jjj,:,col)),4),2)); 
        c(col).LimitsMode="manual";
        c(col).Limits=[round(min(data_to_plot(jjj,:,col)),2) round(max(data_to_plot(jjj,:,col)),2)];
        c(col).Position(1)=c(col).Position(1)-0.035;
        c(col).Position(2)=0.32;
        c(col).Position(3)=0.12;
        c(col).Position(4)=0.03;
        clearvars -except subplot_col pos1 new_pos1 f user chan_labels data_to_plot jjj titles c univariate
    end

    for col=1:5
        pos1(col,:)= get(subplot_col(col), 'Position'); % gives the position of current sub-plot
        pos1(col,[2,3,4])=[0.1100 0.1237 0.8150];
        new_pos1(col,:) = pos1(col,:) +[-0.03 0 0.015 0.03];
        set(subplot_col(col), 'Position',new_pos1(col,:));

    end
    annotation('textbox', [0.4, 0.75, 0.2, 0.1], 'String', ['Participant ' num2str(jjj)], 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontWeight','bold','FontSize', 18.5);
    
    %cd to data location
    saveas(f,[num2str(jjj) 'topoplot.fig']);
    exportgraphics(f, [num2str(jjj) 'topoplot.tif'], 'Resolution', 600);

end


%% for the group average plots
data_to_plot_avg(:,1) = mean(univariate.theta_pow_mean,1);
data_to_plot_avg(:,2) = mean(univariate.alpha1_pow_mean,1);
data_to_plot_avg(:,3) = mean(univariate.alpha2_pow_mean,1);
data_to_plot_avg(:,4) = mean(univariate.beta1_pow_mean,1);
data_to_plot_avg(:,5) = mean(univariate.beta2_pow_mean,1); 

f = figure;
f.Position = [0 0 2400 450];
f.WindowState = 'maximized';   

for col = 1:5
    subplot_col(col)=subplot(1,5,col);
    data = [];
    data.label = chan_labels; 
    data.dimord = 'chan'; 
    data.time = 0; 
    data.avg=data_to_plot_avg(:,col);

    cfg = [];
    cfg.layout = 'elec1005.lay';  
    cfg.marker = 'on'; 
    cfg.gridscale = 300; 
    cfg.zlim = [min(data_to_plot_avg(:,col)) max(data_to_plot_avg(:,col))]; 
    cfg.interactive = 'yes';
    cfg.comment = 'no';
    
    % Plot
    format bank
    ft_topoplotER(cfg, data); c(col)=colorbar('FontSize', 14); c(col).Location='southoutside';
    c(col).Label.String = 'power(ln(\muV^2))';  c(col).Label.FontWeight='bold';
    c(col).Label.FontSize=15; 
    title(titles{col}, 'FontSize', 18, 'FontWeight','bold')
    c(col).TicksMode="manual";
    c(col).Ticks = round(linspace(min(data_to_plot_avg(:,col)),max(data_to_plot_avg(:,col)),4),2);
    c(col).TickLabels =num2cell(round(linspace(min(data_to_plot_avg(:,col)),max(data_to_plot_avg(:,col)),4),2)); 
    c(col).LimitsMode="manual";
    c(col).Limits=[round(min(data_to_plot_avg(:,col)),2) round(max(data_to_plot_avg(:,col)),2)];
    c(col).Position(1)=c(col).Position(1)-0.035;
    c(col).Position(2)=0.32;
    c(col).Position(3)=0.12;
    c(col).Position(4)=0.03;
    clearvars -except subplot_col pos1 new_pos1 f user chan_labels data_to_plot_avg jjj titles c
end

for col=1:5
    pos1(col,:)= get(subplot_col(col), 'Position'); % gives the position of current sub-plot
    pos1(col,[2,3,4])=[0.1100 0.1237 0.8150];
    new_pos1(col,:) = pos1(col,:) +[-0.03 -0.025 0.015 0.03];
    set(subplot_col(col), 'Position',new_pos1(col,:));

end

%cd to data location
saveas(f,'group_avg_topoplot.fig');
exportgraphics(f, ('group_topoplot.tif'), 'Resolution', 600);

 