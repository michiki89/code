%% Granger Causality in Time Domain

clearvars
close all
clc


for ii = 1:18
%% Get RRi and Resp Data

% and the time vector + Fs
% data should be in columns

% define Path and File Name
if ii <10
    fileName = ['PB0',mat2str(ii)];
else
    if ii == 18
        fileName = 'PB19';
    else
        fileName = ['PB',mat2str(ii)];
    end
end

pathName = '/Volumes/beaujolais/Benutzer/mk572/Documents/Data/Paced Respiration Study/Data/';

whichInterval = '0.2';    % In Paced Respiration Study you can choose: 0.1 - 0.6 Hz or spont breathing
whichLead = 1;
whichInterp = 1;          % whichInterp = 0;    --> Interp at Respiration sampling
                          % whichInterp = 1;    --> Interp at 4Hz
% Load data from .mat file
[time,respSig,~,time_RRi,RRiSig,~,Fs] = chooseIntLead(whichInterval,whichLead,fileName,pathName);

% Build column vectors
time_Resp = time';
time_RRi = time_RRi';
RRiSig = (RRiSig' - mean(RRiSig))/std(RRiSig);      % no mean and unit variance signal
respSig = (respSig - mean(respSig))/std(respSig);   % no mean and unit variance signal

% Interpolate Data
switch whichInterp
    case 0      % interpolation to Fs
        RRiSig_int = interp1(time_RRi,RRiSig,time_Resp,'pchip');
        respSig_int = respSig;      % signals with _int will be used for analysis
        time = time_Resp;
        Fs_int = Fs;
    case 1      % interpolation to 4Hz
        time_start = min([time_Resp(1),time_RRi(1)]);
        time_stop = max([time_Resp(end),time_RRi(end)]);
        Fs_int = 4;
        time = [time_start:1/Fs_int:time_stop]';   

        RRiSig_int = interp1(time_RRi,RRiSig,time,'pchip'); % signals with _int will be used for analysis
        respSig_int = interp1(time_Resp,respSig,time,'pchip');
end
        
[F(ii,1),c_v(ii,1),Fmax(ii,1)] = granger_cause_withNorm(RRiSig_int,respSig_int,0.1,20,0.0001,10);
[F(ii,2),c_v(ii,2),Fmax(ii,2)] = granger_cause_withNorm(respSig_int,RRiSig_int,0.1,20,0.0001,10);
end
%% Estimate AR-Modell of RRi

N_RRi   = 1;
N_Resp = 0;

RRi_AR = armax_LS(respSig_int,RRiSig_int,time,Fs_int,N_Resp,N_RRi);

%% Estimate AR-Modell of RRi

N_Resp   = 1;
N_RRi = 0;

Resp_AR = armax_LS(RRiSig_int,respSig_int,time,Fs_int,N_RRi,N_Resp);

%% Estimate ARMAX of RRi with Resp as Input
%% Estimate AR-Modell of RRi

N_RRi   = 1;
N_Resp = 1;

RRi_ARX = armax_LS(respSig_int,RRiSig_int,time,Fs_int,N_Resp,N_RRi);


%% Estimate ARMAX of Resp with RRi as Input
N_Resp   = 1;
N_RRi = 1;

Resp_ARX = armax_LS(RRiSig_int,respSig_int,time,Fs_int,N_RRi,N_Resp);

%% Calculate Variance of Residuum
res_RRi_AR = RRiSig_int - RRi_AR;
res_RRi_ARX = RRiSig_int - RRi_ARX;
res_Resp_AR = respSig_int - Resp_AR;
res_Resp_ARX = respSig_int - Resp_ARX;

var_RRi_AR = var(res_RRi_AR);
var_RRi_ARX = var(res_RRi_ARX);
var_Resp_AR = var(res_Resp_AR);
var_Resp_ARX = var(res_Resp_ARX);

F_RRi = log(var_RRi_AR/var_RRi_ARX);
F_Resp = log(var_Resp_AR/var_Resp_ARX);



%% PLot Interpolation of functions
figure(1)
ax1 = subplot(2,1,1);
plot(time_Resp,respSig)
hold on
plot(time,respSig_int)

ax2 = subplot(2,1,2);
plot(time_RRi,RRiSig)
hold on
plot(time,RRiSig_int)

figure(2)
subplot(2,1,1)
plot(time,RRiSig_int)
hold on
plot(time,RRi_AR);
plot(time,RRi_ARX);

subplot(2,1,2)
plot(time,respSig_int)
hold on
plot(time,Resp_AR);
plot(time,Resp_ARX);


