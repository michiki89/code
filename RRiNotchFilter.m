function [filt_sig,filt_spec,spec_fil] = RRiNotchFilter(signal,time,Fs,N_fft,f_max,peakWidth,filterOrder)


notchFilter = fdesign.notch('N,F0,BW',filterOrder,2*f_max(1)/Fs,2*peakWidth/Fs);
        filterNotch = design(notchFilter);
%         [b,a] = sos2tf(filterNotch.sosMatrix);
        spec_fil = freqz(filterNotch,N_fft/2);
        
        input = signal;%[flip(signal),signal,flip(signal)];
        filt_sig = filtfilt(filterNotch.sosMatrix,filterNotch.ScaleValues,input)';
%         filt_sig = filt_sig(length(time)+1:2*length(time));
        filt_spec = fft(filt_sig,N_fft);