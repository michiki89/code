function [filt_sig,filt_spec,spec_fil] = RRiGaussFilter(RRi_int,time,Fs,N_fft,f_m,peakWidth)
               

f_filter = -Fs/2:Fs/N_fft:Fs/2-Fs/N_fft;
                spec_RRi = fft(RRi_int,N_fft);
%                 f_m = f_max(1);
                FWHM_f = peakWidth/2;
                sigma_f = FWHM_f/2*sqrt(log(2)*2);              
                
                
                %% Gauss Spectral Domain
                spec_fil = 1 - (exp( -(f_filter-f_m).^2/(2*(sigma_f.^2))) + exp( -(f_filter+f_m).^2/(2*(sigma_f.^2))));
                
                
                spectral_filter = [spec_fil(end/2+1:end),spec_fil(1:end/2)];

                filt_spec = spectral_filter .* spec_RRi;

                filt_sig = ifft(filt_spec);
                filt_sig = filt_sig(1:length(time));
                spec_fil = spec_fil(1,end/2+1:end);
%                 filt_spec = filt_spec(1,1:end/2);