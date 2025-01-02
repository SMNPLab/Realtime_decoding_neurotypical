function [clabel,dval] = test_lda_realtime(cf,X,c1,c2,dvalcriteria)
% Applies an LDA classifier to test data and produces class labels,
% decision values, and posterior probabilities.
% updated by SJH to use pre-defined d-value cutoffs for class decisions

% checked by sjh on 10/13/2022 prior to uttara beginning full realtime mvpa experiment
%
% Usage:
% [clabel,dval] = test_lda(cf,X)
%
%Parameters:
% cf             - struct describing the classifier obtained from training
%                  data. Must contain at least the fields w and b,
%                  see train_lda
% X              - [number of samples x number of features] matrix of
%                  test samples
% c1             - dvalue threshold for class 1 predictions (should be a positive number)
% c2             - dvalue threshold for class 2 predictions (should be a negative number)
% dvalcriteria   - mention if you want the dval criteria to be applied (for
% confident predictions) 0= no dval criteria wanted, 1= dval criteria
% wanted
%
%Output:
% clabel        - predicted class labels (1's and 2's, 0's reflect no prediciton)

dval = X*cf.w + cf.b;
if dvalcriteria==1
    clabel = double(dval >= c1) + 2*double(dval <= c2);
elseif dvalcriteria==0
    clabel = double(dval >0) + 2*double(dval <0);
end


end
