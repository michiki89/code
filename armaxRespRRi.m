% --------------------------------------------------------- %
% ARMA RRi Resp
% --------------------------------------------------------- %
clearvars
close all
clc

%% Load Data
whichInterval = 'spont';
fileName = 'PB02';
pathName = '/Volumes/beaujolais/Benutzer/mk572/Documents/Data/Paced Respiration Study/Data/';

%% LOAD DATA

% LOAD ORIGINAL DATA ------------------------------------------- %
whichLead = 1;
[time,respSig,signal,t_RRi,RRi,~,Fs] = chooseIntLead(whichInterval,whichLead,fileName,pathName);
respSig = (respSig - mean(respSig))/std(respSig);
% FPT = Process_ECG(signal',Fs);
% FPT = FPTtoStruct(FPT);

%% Interpolate RRi to get equidistant sampling points
t_RRi_int = time;%t_RRi(1):1/4:t_RRi(end);
RRi_int = interp1(t_RRi,RRi,time,'pchip');
RRi_int = (RRi_int - mean(RRi_int))/std(RRi_int);

%% LS Estimation of ARMAX Parameters
output = RRi_int';
input = respSig;
N_ex = 1;
N_y = 0;
[RRi_resp,RRi_res,sysTF,coeff_armax] = armax_LS(input,output,time,Fs,N_ex,N_y);


%% RRi-ARMAX
% sigRRi = iddata(RRi_int',[],1/Fs);
% 
% N_RRi = 0;
% N_n = 100;
% % 
% % armaRRi = armax(sigRRi,[N_RRi,N_n]);
% [armaRRi,errRRi,sysRRi,coeff_RRi] = armax_LS(input,output,time,Fs,N_n,N_RRi);
% 
% 
% %% RESP-ARMA
% % sigResp = iddata(respSig,[],1/Fs);
% % 
% N_Resp = 0;
% N_n = 100;
% % 
% % 
% % 
% % armaResp = armax(sigResp,[N_Resp,N_n]);
% [armaResp,errResp,sysResp,coeff_Resp] = armax_LS(output,input,time,Fs,N_n,N_Resp);

sigRRi = iddata(RRi_int',respSig,1/Fs);

N_RRi = 1;
N_Resp = 0;
N_n = 10;
% 
armaRRi = armax(sigRRi,[N_RRi,N_Resp,N_n,1]);
% [armaRRi,errRRi,sysRRi,coeff_RRi] = armax_LS(input,output,time,Fs,N_n,N_RRi);

%% Plotting
subplot(2,1,1)
plot(time,respSig,'-b','DisplayName','Respiration')
hold on
plot(time,RRi_int,'-r','DisplayName','original')
plot(time,RRi_resp,'-g','DisplayName','Armax')
plot(time,RRi_res,'--k','DisplayName','Residual')
legend show
box on
grid on

subplot(2,1,2)
N_fft = 2^nextpow2(length(time));
plot(0:Fs/N_fft:Fs - Fs/N_fft,abs(fft(RRi_int,N_fft)),'-r','DisplayName','original')
hold on
plot(0:Fs/N_fft:Fs - Fs/N_fft,abs(fft(RRi_resp,N_fft)),'-g','DisplayName','Armax')
plot(0:Fs/N_fft:Fs - Fs/N_fft,abs(fft(RRi_res,N_fft)),'--k','DisplayName','Residual')
legend show
box on
grid on
xlim([0,1])