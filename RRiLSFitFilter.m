function [filt_sig,filt_spec,spec_fil] = RRiLSFitFilter(RRi_int,time,Fs,N_fft,f_m,maxInd);
                f_filter = [0:Fs/N_fft:Fs-Fs/N_fft]';
                Func1 = fft(RRi_int',N_fft);
%                 f_m = f_max;
%                 maxInd = indPeaks;
                N_samp = length(time);
                %% Aquire points for LS-Gauss Fit
                % search min next to fmax and use points in between maximum and minimus for
                % gauss fit
                [~,minInd] = findpeaks(-abs(Func1(1:maxInd,1)));
                startInd = minInd(end);
                
                [~,minInd] = findpeaks(-abs(Func1(maxInd+1:end,1)));
                endInd = minInd(1)+maxInd;
                
                
                %% Create linear Model for LS Estimation
                % y = B*exp(-A*f) --> take natural logarithm on both sides and substitute
                % ln(y) = y_tilde and ln(B) = C --> y_tilde(f) = [-(f-fmax)^2 1]* [A;C]
                % with [-(f-fmax)^2 1] = P and theta = [A;C]
                
                y = abs(Func1(startInd:endInd,1));
                y_tilde = log(y);
                
                P = [-(f_filter(startInd:endInd,1)-f_m).^2,ones(length(f_filter(startInd:endInd,1)),1)];
                
                theta = inv(P'*P)*P'*y_tilde;
                
                A = theta(1);
                B = exp(theta(2));
                
                gaussFilter = B*exp(-A*(f_filter(1:end/2,1)-f_m).^2);
                gaussFilter = 1 -gaussFilter/max(gaussFilter);
                gaussFilter = [gaussFilter;1;flip(gaussFilter(2:end,1))];
                
                %% Filter signal spectrally
                Filt_Func1 = Func1.*gaussFilter;
                
                filt_func1 = ifft(Filt_Func1);
                filt_func1 = filt_func1(1:N_samp,1);
                
                filt_sig = filt_func1;
                spec_fil = gaussFilter(1:end/2,1);
                filt_spec = Filt_Func1;