function [perf,perf_all,predlabel,runtime,ki,lda,pparam_csp,dval_all,plabel,test_indices] = manual_classify_csp_psd(numk,raw,label) 

% the artist formerly known as manual_optimize_classify_rt_uk 

% checked by sjh 10/13/2022 prior to uttara beginning full realtime mvpa experiment 

% classifies data for mvpa_tmseeg data by hand
% performs feature selection and hyperparameter tuning simultaneously for each fold
% also includes csp

% INPUT ARGUMENTS:
%   - numk = number of folds
%   - raw = raw data for csp calculation; trl x chan x sample 
%   - label = vector of true class labels

% OUTPUT ARGUMENTS: confirm with uttara these are correct 
%   - perf = matrix of auc values with dval cutoff (k x feature x lambda)
%   - perf_all = matrix of auc values without dval cutoff (k x feature x lambda) 
%   - predlabel = cell array of all predicted labels (feature x lambda)
%   - plabel = cell array of predicted labels per fold (feature x lambda x k)
%   - runtime = scalar value indicating runtime for entire function
%   - ki = cell array of feature ranks (k x feature)
%   - lda =  cell array of trained classifiers (feature x lambda x k)
%   - pparam_csp = csp filters (one per fold)
%   - prob = probability values (k x feature x lambda)
%   - dval_all = all dvals, regardless of whether they meet threshold (k x feature x lambda) 

% partition data into numk folds using sprinkly approach
tic;
rng(42);
c=cvpartition(label,'KFold',numk,'Stratify',true);
for k=1:numk
    folds{k}=test(c,k);
end

% define lambda search space - this is the regularization rate
ss=linspace(1e-10,1); 

% preallocate predlabel cell array so it can be filled during cross validation
predlabel=cell(310,numel(ss)); % feature x lambda
for a=1:size(predlabel,1)
    for b=1:size(predlabel,2)
        predlabel{a,b}=zeros(length(label),1);
    end
end

% set up csp parameters
pparam=[];
pparam.n=31; % for a total of nchan components 
pparam.calculate_spatial_pattern=1;
pparam.feature_dimension=2; % channels
pparam.target_dimension=3; % time
pparam.lambda=10^-10;
pparam.calculate_variance=0;
pparam.calculate_log=0;

% train and test lda - all optimization done per fold
for k=1:numk
    
    disp(['... fold ' num2str(k) ' ...']);
    
    % select train and test indices
    train_idx=find(folds{k}==0);
    test_idx=find(folds{k}==1);
    test_indices{k}=test_idx;
    
    % select train and test labels
    train_klabel=label(train_idx);
    test_klabel=label(test_idx);
    
    % select train and test raw data
    train_raw=raw(train_idx,:,:);
    test_raw=raw(test_idx,:,:);
    
    % apply csp to training data and get weights/spatial filter
    rng(42); 
    pparam.is_train_set=1;
    [pparam_csp{k}, segtrain, ~] = mv_preprocess_csp(pparam, train_raw, train_klabel); %%% sjh: replaced last output argument with wiggly, we dont need it 
    
    %calculate power from csp_features
    [~,train_kfeatures]=calc_ps_rt(segtrain);
    
    %test features
    rng(42); 
    pparam_csp{k}.is_train_set=0;
    [~, segtest, ~] = mv_preprocess_csp(pparam_csp{k}, test_raw, []);
    [~,test_kfeatures] = calc_ps_rt(segtest);

    % rank training features
    rng(42);
    ranks(k,:)=fscchi2(squeeze(train_kfeatures),train_klabel);
    
    % simultaneously optimize feature # (based on fold-specific feature ranking) and lambda value
    for ff=1:size(train_kfeatures,2)
        
        % select topmost features (feature num = ff)
        ki{k,ff}=ranks(k,1:ff);
        rel_train_features{ff,k}=train_kfeatures(:,ki{k,ff});
        rel_test_features{ff,k}=test_kfeatures(:,ki{k,ff});
        
        % search all predefined lambda values
        for ii=1:length(ss)
            
            pp{ii}=mv_get_hyperparameter('lda'); 
            pp{ii}.lambda=ss(ii); % regularization coefficient 
            pp{ii}.reg='shrink';
            pp{ii}.prob=1;  % added to calculate probabilities
            pp{ii}.form='primal'; % added to calculate probabilities
            
            % train lda using topmost features (feature num = ff)
            rng(42);
            warning('');
            lda{ff,ii,k}=train_lda(pp{ii},rel_train_features{ff,k},train_klabel); % train and test
            
            % test lda using topmost features (feature num = ff)
            rng(42);
            [predlabel{ff,ii}(test_idx),dval,prob{ff,ii,k}]=test_lda(lda{ff,ii,k},rel_test_features{ff,k});
            
            % select high confidence predictions (top 50 percentile of dvals per class)
            lodval=find(dval<0);
            hidval=find(dval>=0); 
            dind=[find(dval<=prctile(dval(lodval),50)); find(dval>=prctile(dval(hidval),50))]; 
            %dind= [find(dval<=median(dval(lodval))); find(dval>=median(dval(hidval))) ];
            dind=sort(dind); % sjh - why sorting here? 
            
            % exit code if singularity error occurs
            if ~isempty(lastwarn)
                error('singularity error.');
            end
            
            % performance
            perf(k,ff,ii)=mv_calculate_performance('f1','dval',dval(dind),test_klabel(dind)'); % performance for high confidence predictions 
            perf_all(k,ff,ii)=mv_calculate_performance('f1','dval',dval,test_klabel'); % performance for all predictions 
            dval_all{ff,ii,k}=dval;
            plabel{ff,ii,k}=predlabel{ff,ii}(test_idx);
        end
        
    end
    runtime=toc;
    clearvars -except plabel dval_all perf test_indices predlabel runtime ki lda pparam_csp  folds pparam ss label numk raw perf_all
     
end

end
