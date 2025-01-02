function [pxx, features]= calc_ps_rt(seg)

% written by uuk to calculate power spectral features from csp data for step3 / offline decoding 
% this is called in manual_classify_csp_psd 

%INPUT
% data= trials*channels*timepoints

%OUTPUT
%pxx= power spectrum values (freq*trial*channel)
%feature matrix = ps feature matrix (trial * (channel* number of freq bands))

% checked by sjh on 10/13/2022 prior to uttara beginning full realtime mvpa experiment 

seg= permute(seg,[3,1,2]); % original = trl, chan, sample --> new = sample, trl, chan 
freqs=[4:0.25:35];
for chan=1:size(seg,3)
    for trl=1:size(seg,2)
        pxx(:,trl,chan)=pwelch(seg(:,trl,chan),500,[],freqs,1000); % freq x trl x chan x time
    end
end
theta=[find(freqs==4):find(freqs==7.75)];
alpha1=[find(freqs==8):find(freqs==10)];
alpha2=[find(freqs==10.25):find(freqs==13)]; 
beta1=[find(freqs==13.25):find(freqs==20)];
beta2=[find(freqs==20.25):find(freqs==35)];

% compute power timeseries - averaging over frequency (dimension 1)
theta_pow=squeeze(mean(pxx(theta,:,:),1)); % trial x channel x time
alpha1_pow=squeeze(mean(pxx(alpha1,:,:),1));
alpha2_pow=squeeze(mean(pxx(alpha2,:,:),1));
beta1_pow=squeeze(mean(pxx(beta1,:,:),1));
beta2_pow=squeeze(mean(pxx(beta2,:,:),1));

features=cat(2,theta_pow,alpha1_pow,alpha2_pow,beta1_pow,beta2_pow); 
    
end    