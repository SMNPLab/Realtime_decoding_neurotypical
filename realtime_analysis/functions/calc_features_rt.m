function [ff, pxx]= calc_features_rt(csp_params,csp_chunk, pow_freq, processing_rate, theta, alpha1, alpha2, beta1, beta2)

% written by uuk to calculate power spectra for a single trial in realtime code 
% used in step4 / realtime decoding 

% checked by sjh on 10/13/2022 prior to uttara starting full realtime mvpa exp

[~,segrt,~]=mv_preprocess_csp(csp_params,csp_chunk,[]); % segrt should match dim of csp_chunk
seg=squeeze(permute(segrt,[3,1,2])); % trial x channel x time --> time x trial x channel 

pxx = pwelch(seg, size(seg,1),[],pow_freq, processing_rate); % freq x chan 

% compute features
theta_pow=squeeze(mean(pxx(theta,:,1))); % trial x channel x time
alpha1_pow=squeeze(mean(pxx(alpha1,:,1)));
alpha2_pow=squeeze(mean(pxx(alpha2,:,1)));
beta1_pow=squeeze(mean(pxx(beta1,:,1)));
beta2_pow=squeeze(mean(pxx(beta2,:,1)));
ff=[theta_pow alpha1_pow alpha2_pow beta1_pow beta2_pow];

end