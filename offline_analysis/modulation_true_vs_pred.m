%code to check for targeting accruacy wrt to true labels
clc
clearvars
user='';
%path to data
%path to codes
subjname={'01','02','03','04','05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19'};

%% load the mep amplitudes from the training set
%load all the mep amplitudes for states for fdi at 120% rmt
%path to data location
load('test_emg_data.mat');
for jjj=1:19
    mep_test_high{jjj}=results.normalized_meps{jjj,1,1,1};
    mep_test_low{jjj}=results.normalized_meps{jjj,1,1,2};
    pred_mod(jjj)= nanmean(mep_test_high{jjj}./nanmean(mep_test_low{jjj}));
end    

%for true mod calculation
for jjj=1:numel(subjname)
all_pred_meps_withnan=cat(2,mep_test_high{jjj}, mep_test_low{jjj});
all_pred_meps_no_nan=all_pred_meps_withnan(~isnan(all_pred_meps_withnan));
all_pred_meps{jjj}=all_pred_meps_no_nan;
    high_pred=[];
    low_pred=[];
    %get true labels
    for trl=1:numel(all_pred_meps{jjj})
        if all_pred_meps{jjj}(trl)>=median(all_pred_meps{jjj})
            high_pred=[high_pred all_pred_meps{jjj}(trl)];
        elseif all_pred_meps{jjj}(trl)< median(all_pred_meps{jjj})
            low_pred=[low_pred all_pred_meps{jjj}(trl)];
        end
    end

    high_meps_true{jjj}=high_pred;
    low_meps_true{jjj}=low_pred;

    true_mod(jjj)=nanmean(high_meps_true{jjj}) ./nanmean(low_meps_true{jjj});

    clearvars -except true_mod high_meps_true low_meps_true all_pred_meps mep_test_high mep_test_low pred_mod user subjname uniquenum 
end
results.true_mod=true_mod;
results.pred_mod=pred_mod;
results.high_meps_true=high_meps_true;
results.low_meps_true=low_meps_true;
results.mep_high_pred=mep_test_high;
results.mep_low_pred=mep_test_low;

sem_pred_mod=nanstd(pred_mod(1:18))./sqrt(numel(pred_mod(~isnan(pred_mod(1:18))))) %did not include subject 19 since one MEPs for one of the states was not available
sem_true_mod=nanstd(true_mod(1:18))./sqrt(numel(true_mod(~isnan(true_mod(1:18)))))
%% for plotting the graph

col="#ab590d";

figure;
for jjj=[1 3:18]
    scatter([1],[pred_mod(jjj)],'MarkerFaceColor',col,'MarkerEdgeColor',col, 'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2,'SizeData',80); hold on; 
    scatter([2],[true_mod(jjj)],'MarkerFaceColor',col,'MarkerEdgeColor',col,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2,'SizeData',80); hold on; 
    plot([1,2],[pred_mod(jjj) true_mod(jjj)],'color', '#e39952' ,'linestyle',':', 'LineWidth',0.25); hold on; 
end
scatter([1],[pred_mod(2)],'MarkerFaceColor',col,'MarkerEdgeColor',col, 'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2,'SizeData',80); hold on; 
scatter([2],[true_mod(2)],'MarkerFaceColor',col,'MarkerEdgeColor',col,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2,'SizeData',80); hold on; 



scatter([1],[nanmean(pred_mod(1:18))],'square','MarkerFaceColor',col,'MarkerEdgeColor',col,'SizeData',400); hold on; 
scatter([2],[nanmean(true_mod(1:18))],'square','MarkerFaceColor',col,'MarkerEdgeColor',col,'SizeData',400); hold on; 
errorbar([1,2],[nanmean(pred_mod(1:18)) nanmean(true_mod(1:18))],[sem_pred_mod sem_true_mod],'k','LineWidth',2,'Capsize', 17); hold on; %updated the sem values here
xticks([1,2])
xticklabels({'predicted', 'true'});
ax = gca;
ax.YAxis.FontSize = 17;
ax.XAxis.FontSize=17;
ax.XAxis.FontWeight="bold";
xlim([0.75 2.25])
ylim tight;
breakyaxis([3 5.5]);
ylabel('(strong/weak) CST output ','FontSize',17,'FontWeight','bold');

