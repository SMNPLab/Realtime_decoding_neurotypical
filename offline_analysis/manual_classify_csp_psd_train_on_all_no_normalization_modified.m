function [perf,perf_all,runtime,ki,lda,pparam_csp,dval_all,plabel] = manual_classify_csp_psd_train_on_all_no_normalization_modified(numk,raw,label,k) 

% classifies data for mvpa_tmseeg data 
% performs feature selection and hyperparameter tuning simultaneously for each fold
% also includes csp

% INPUT ARGUMENTS:
%   - numk = number of folds
%   - raw = raw data for csp calculation; trl x chan x sample 
%   - label = vector of true class labels
%   - k = subject to use as test set 

% OUTPUT ARGUMENTS: 
%   - perf = matrix of auc values with dval cutoff (k x feature x lambda)
%   - perf_all = matrix of auc values without dval cutoff (k x feature x lambda) 
%   - plabel = cell array of predicted labels per fold (feature x lambda x k)
%   - runtime = scalar value indicating runtime for entire function
%   - ki = cell array of feature ranks (k x feature)
%   - lda =  cell array of trained classifiers (feature x lambda x k)
%   - pparam_csp = csp filters (one per fold)
%   - dval_all = all dvals, regardless of whether they meet threshold (k x feature x lambda) 

% partition data into numk folds using sprinkly approach
tic;
rng(42);

% define lambda search space - this is the regularization rate
ss=linspace(1e-10,1); 

%preallocate predlabel cell array so it can be filled during cross validation

% set up csp parameters
pparam=[];
pparam.n=31; % for a total of nchan components 
pparam.calculate_spatial_pattern=1;
pparam.feature_dimension=2; % channels
pparam.target_dimension=3; % time
pparam.lambda=10^-10;
pparam.calculate_variance=0;
pparam.calculate_log=0;

predlabel=cell(310,numel(ss)); % feature x lambda- cell array for predicted states

disp(['... fold ' num2str(k) ' ...']);

% set training indices, excluding test indices 
train_array=[1:numk];
train_array(k)=[]; 

% select train and test raw data    
train_raw=cell2mat(raw(1,train_array)'); %all but one subject
test_raw=cell2mat(raw(1,k)); %the one subject used for test

train_klabel= cell2mat(label(1,train_array));
test_klabel= cell2mat(label(1,k));

% apply csp to training data and get weights/spatial filter
rng(42); 
pparam.is_train_set=1;
[pparam_csp, segtrain, ~] = mv_preprocess_csp(pparam, train_raw, train_klabel); 

%calculate power from csp_features
[~,train_kfeatures]=calc_ps_rt(segtrain);

% rank training features
rng(42);
ranks=fscchi2(squeeze(train_kfeatures),train_klabel);

%test features
rng(42); 
pparam_csp.is_train_set=0;
[~, segtest, ~] = mv_preprocess_csp(pparam_csp, test_raw, []);
[~,test_kfeatures] = calc_ps_rt(segtest);

% simultaneously optimize feature # (based on fold-specific feature ranking) and lambda value
for ff=1:size(train_kfeatures,2)

    % select topmost features (feature num = ff)
    ki{ff}=ranks(1:ff);
    rel_train_features{ff}=train_kfeatures(:,ki{ff});
    rel_test_features{ff}=test_kfeatures(:,ki{ff});

% search all predefined lambda values
    for ii=1:100 
        
        pp{ii}=mv_get_hyperparameter('lda'); 
        pp{ii}.lambda=ss(ii); % regularization coefficient 
        pp{ii}.reg='shrink';
        pp{ii}.prob=1;  % added to calculate probabilities
        pp{ii}.form='primal'; % added to calculate probabilities
        
        % train lda using topmost features (feature num = ff)
        rng(42);
        warning('');
        lda{ff,ii}=train_lda(pp{ii},rel_train_features{ff},train_klabel); % train and test
        
        % test lda using topmost features (feature num = ff)
        rng(42);
        [predlabel{ff,ii},dval,prob{ff,ii}]=test_lda(lda{ff,ii},rel_test_features{ff});
        
        % select high confidence predictions (top 50 percentile of dvals per class)
        lodval=find(dval<0);
        hidval=find(dval>=0); 
        dind=[find(dval<=prctile(dval(lodval),50)); find(dval>=prctile(dval(hidval),50))]; 
        dind=sort(dind); 
        
        % exit code if singularity error occurs
        if ~isempty(lastwarn)
            error('singularity error.');
        end
        
        % performance
        perf(ff,ii)=mv_calculate_performance('f1','dval',dval(dind),test_klabel(dind)'); % performance for high confidence predictions 
        perf_all(ff,ii)=mv_calculate_performance('f1','dval',dval,test_klabel'); % performance for all predictions 
        dval_all{ff,ii}=dval;
        plabel{ff,ii}=predlabel{ff,ii};
    
        runtime=toc;
    end
    
end


end
