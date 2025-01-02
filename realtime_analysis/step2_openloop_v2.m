% step2_openloop_v2.m 
% delivers single-pulse open-loop TMS at a desired ISI with random jitter
% queries user to indicate the number of pulses and the desired ISI

% housekeeping 
clc
clearvars
addpath(genpath('C:/Users/realtime/Box/SMNP_Lab/LSL_workingfolder/functions'));
addpath(genpath('C:/Users/realtime/Box/SMNP_Lab/LSL_workingfolder/mvpa/ttl'));

% set up stimulation parameters 
max_triggers=input('Enter the number of pulses you want to deliver: '); % desired number of pulses 
isi=3; 

trigger_val=3; % value to write to parallel port
total_triggers=0; % counts number of pulses iteratively 
isi=isi + rand(max_triggers,1); % inter stimulus interval + random variation 
address = hex2dec('7FE8'); % parallel port address 

% initialize parallel port 
config_io;
global cogent; 
if( cogent.io.status ~= 0 )
    error('inp/outp installation failed');
end
outp(address,0); % flush parallel port 

% main loop
for trigger=1:max_triggers
    
    outp(address, trigger_val); % send pulse 
    total_triggers=total_triggers+1;                     
    disp(['trigger ' num2str(total_triggers) ' sent.']);
    
    pause(isi(trigger)); 
    outp(address,0); % flush buffer 

end

disp('stimulation complete.') 

