% Max. Peak Power Coupling Measure


clearvars -except xcorr_nf xcorr_f
close all
clc


for ii= 1:18

    
%% Load Data
whichDataset = 0;
switch whichDataset
    case 0
        whichInterval = 'spont';%mat2str(round(0.1*jj,1));
        if ii<10
            fileName = ['PB0',mat2str(ii)];
        else
            if ii == 18
                fileName = ['PB19'];
            else
                fileName = ['PB',mat2str(ii)];
            end
        end
        pathName = '/Volumes/beaujolais/Benutzer/mk572/Documents/Data/Paced Respiration Study/Data/';
        
        
        %% LOAD DATA
        
        % LOAD ORIGINAL DATA ------------------------------------------- %
        whichLead = 1;
        [time,respSig,signal,t_RRi,RRi,~,Fs] = chooseIntLead(whichInterval,whichLead,fileName,pathName);

        if isempty(time)
            xcorr_f(ii,jj) = NaN;
            xcorr_nf(ii,jj) = NaN;
            powerMeas(ii,jj) = NaN;
            continue
        end
    case 1
        fileName = 'f1y01';
        pathName = '/Volumes/beaujolais/Benutzer/mk572/Documents/Data/Fantasia/';
        load([pathName,fileName]);
        

end

respSig = (respSig - mean(respSig))/std(respSig);

%% Interpolate RRi to get equidistant sampling points
whichInterp = 1;

switch whichInterp
    case 0
        t_RRi_int = time;
        RRi_int = interp1(t_RRi,RRi,time,'pchip');
    case 1
        Fs = 4;
        t_start = min([time(1),t_RRi(1)]);
        t_end   = max([time(end),t_RRi(end)]);
        t_RRi_int = t_start:1/Fs:t_end;
        RRi_int = interp1(t_RRi,RRi,t_RRi_int,'pchip');
        respSig = interp1(time,respSig,t_RRi_int,'pchip');
        time = t_RRi_int;
end




RRi_int = (RRi_int - mean(RRi_int))/std(RRi_int);
N_fft = 8*2^nextpow2(length(time));

welchInterval = length(time);
welchOverlap = 0;
welchZP = 1;
spec_RRi = fft(RRi_int,N_fft*welchZP);



Pxx_RRi = pwelch(RRi_int,welchInterval,welchOverlap,welchZP*N_fft,Fs);

%% Compute cross_spec
Pxx_Resp = pwelch(respSig,welchInterval,welchOverlap,welchZP*N_fft,Fs);
[CPSD,f] = cpsd(respSig,RRi_int',welchInterval,welchOverlap,welchZP*N_fft,Fs,'onesided');


%% Find Maximal Peak and its frequency in certain frequency band
lowerBound = 0.02;
upperBound = 1;  % Bandwidth of area to search for peak in cross spec
thresFact = 1/2;

coupleMeas = 0;
i=0;

indUpper = find(f <=upperBound,1,'last');
indLower = find(f >= lowerBound ,1,'first');

[peaks,indPeaks] = findpeaks(abs(CPSD(indLower:indUpper)));
[peaks,rearr] = sort(peaks,1,'descend');
indPeaks = indPeaks(rearr) + indLower-1;

while coupleMeas <= 0
    i=i+1;

    f_max(i) = f(indPeaks(i));
    
    threshold = thresFact*peaks(i);
    
    CPSD_red = abs(CPSD);
    CPSD_red(CPSD_red < threshold) = 0;
    ind = find(CPSD_red(1:indPeaks(i)) == 0,1,'last');
    if isempty(ind)
        ind = 0;
    end
    indStart(i) = ind +1;
    indEnd(i)   = find(CPSD_red(indPeaks(i):end) == 0,1,'first');
    indEnd(i)   = indPeaks(i) + indEnd(i) -1;
    
    coupleMeas = coupleMeas + sum(abs(CPSD(indStart:indEnd)))/sum(abs(CPSD(indLower:indUpper)))*100;
    peakWidth(i) = f(indEnd(i)) - f(indStart(i));
end    



%% Filter
for jj = 4
whichDomain = 'time';
whichFilter = 'adaptiveArmax';
switch whichDomain
    case 'freq'
        switch whichFilter
            case 'gauss'
                %% Gauss Spectral Filter
                [filt_sig,filt_spec,spec_fil] = RRiGaussFilter(RRi_int,time,Fs,N_fft,f_max(i),peakWidth(i));
                
                
            case 'lsfilter'
                
                [filt_sig,filt_spec,spec_fil] = RRiLSFitFilter(RRi_int,time,Fs,N_fft,f_max(1),indPeaks(1));
                
               
        end
                
    case 'time'
        switch whichFilter
            case 'notch'
                [filt_sig,filt_spec,spec_fil] = RRiNotchFilter(RRi_int,time,Fs,N_fft,f_max(i),peakWidth(i),10);

            case 'lsfilter'
                [filt_sig,filt_spec,spec_fil] = lsGauss_timeDomain(RRi_int,respSig,N_fft,welchInterval,0,4,f_max(1));
                
                  case 'adaptiveArmax'
                [filt_sig,RRiSig_resp,resVar,sqError] = adaptiveFilter_TD(RRi_int,respSig,10);
                filt_spec = fft(filt_sig,N_fft);
                spec_fil = NaN(N_fft/2,1);
           
        end
end

%[filtered_signal,gauss] = filtGauss(CPSD,Fs,4*(f(indEnd(1))-f(indStart(1))),f_max(1));


% time_new = linspace(time(1),time(2),length(filt_sig));
% filt_sig = interp1(time_new,abs(filt_sig)',time);

[CPSD_new,f_new] = cpsd(respSig,filt_sig,welchInterval,welchOverlap,welchZP*N_fft,Fs);
Pxx_RRi_new = pwelch(filt_sig,welchInterval,welchOverlap,welchZP*N_fft,Fs);
% Pxx_RRi_new = Pxx_RRi_new(1:end/2,1) + flip(Pxx_RRi_new(end/2+1:end,1));


%% Cross Correlation of respiration with RRi and filtered RRi
xcorr_nf(ii,jj) = max(xcorr(respSig,RRi_int,'coeff'));
xcorr_f(ii,jj)  = max(xcorr(respSig,filt_sig,'coeff'));

powerMeas(ii,jj) = coupleMeas;
end
[F(ii,2),c_v(ii,2),Fmax(ii,2)] = granger_cause_withNorm(filt_sig,respSig,0.1,5,0.0001,2);
[F(ii,1),c_v(ii,1),Fmax(ii,1)] = granger_cause_withNorm(respSig,filt_sig,0.1,5,0.0001,2);
[F_raw(ii,1),c_v(ii,1),Fmax(ii,1)] = granger_cause_withNorm(respSig,RRi_int,0.1,5,0.0001,2);
[F_raw(ii,2),c_v(ii,1),Fmax(ii,2)] = granger_cause_withNorm(RRi_int,respSig,0.1,5,0.0001,2);
end
%% Plot maxPeak Measure
figure(1)
plot(f,abs(CPSD),'DisplayName','Cross Spectrum')
hold on
for ii=1:i
plot(f_max(ii),peaks(ii),'or','DisplayName','maximal Frequency')
% plot([f(indStart),f(indStart)],[0,round(maxPeak/10)*10],'--r')
% plot([f(indEnd),f(indEnd)],[0,round(maxPeak/10)*10],'--r')
area(f(indStart(ii):indEnd(ii)),abs(CPSD(indStart(ii):indEnd(ii))),'FaceColor',[1,0,0],'DisplayName',[mat2str(round(100*coupleMeas)/100) '% of Power' ])
end

xlim([0:1])
grid on
box on
legend show


%% Plot Welch RRi with and without Resp
figure(2)
ax1 = subplot(2,1,1);
plot(time,respSig)
hold on
plot(time,RRi_int)
plot(time,real(filt_sig))

grid on
box on
xlabel('time in s')
ylabel('RRi signals')

ax2 = subplot(2,1,2);
plot(f,abs(Pxx_RRi))
hold on
plot(f_new,abs(Pxx_RRi_new))
xlim([0,1])
grid on
box on
xlabel('frequency in Hz')
ylabel('Welch Pxx')

%% Plot FFT RRi with and without Resp and Filter
figure(3)
plot(f,abs(spec_RRi(1:end/2+1))/max(abs(spec_RRi(1:end/2+1))))
hold on
plot(f(1:end-1),abs(spec_fil))
plot(f(1:end-1),abs(filt_spec(1:end/2))/max(abs(spec_RRi)))
xlim([0,1])
grid on 
box on
xlabel('frequency in Hz')
ylabel('FFt Spectrum')

%% Plot CPSD before and after Filtering
figure(4)
plot(f,abs(CPSD))
hold on
plot(f_new,abs(CPSD_new))
xlim([0,1])

ylabel('Pxy Welch')


