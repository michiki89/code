function [RRiSig_filt,RRiSig_filt_spec,gauss_filter_spec] = lsGauss_timeDomain(RRiSig,respSig,N_fft,N_window,N_overlap,Fs,fmax)
%% Function to calculate a time domain filter, designed in the spectral domain with a LS Approximation
% The function estimates the CPSD between RRi and Respiration signal, which
% have to have same the sampling rate. The CPSD is estimated and the
% gaussian spectral filter is fitted to the peak around fmax. after fitting
% the filter, the time domain representation of the filter is calculated
% and the RRi is filtered in time domain with the forward backward filtering
% (zero phase filtering)
%
% ----- Inputs -----
% RRiSig:       column vector of RRi time series sampling points (e.g. 4Hz interpolated
%               data or respiration sampling time)
% respSig:      column vector of RRi time series sampling points (e.g. 4Hz interpolated
%               data or respiration sampling time)
% N_fft:        zeropadding / sampling points of CPSD
% N_window:     window length in samples for Welch CPSD estimation
% N_overlap:    window overlap in samples for Welch CPSD estimation
% Fs:           sampling rate of respiration and RRi
% fmax:         frequency of maximum in cross PSD

%----- Outputs -----
% RRiSig_filt: filtered RRi time series

%% If RRi and Resp not column vector --> transposed
%Make sure RRi is a column vector
[a,b] = size(RRiSig);
if (b>a)
    %x is a row vector -- fix this
    RRiSig = RRiSig';
end

%Make sure resp is a column vector
[a,b] = size(respSig);
if (b>a)
    %x is a row vector -- fix this
    respSig = respSig';
end


%% Sampling Parameter Definition
N_samp = length(RRiSig);    % number of samples of RRi and Respiration time series
dt = 1/Fs;
df = Fs/N_fft;              % spectral resolution
f = [0:df:Fs-df]';          % spectral sampling points


%% CPSD estimation
xspec = cpsd(RRiSig,respSig,N_window,N_overlap,N_fft,'twosided');   % cross spectrum


%% LS sampling point extraction
[~,Nmax] = min(abs(f-fmax));        % find sample, where fmax in cross spectrum

% Find minima next to main lobe (assumption: zeropadding and welch averaging is large enough
% that PSD is smoothed good enough for robust minima detection)
[~,help] = findpeaks(-abs(xspec(1:Nmax,1)));   % Find minimum to the left of the main peak
Nmin_left = help(end);
[~,help] = findpeaks(-abs(xspec(Nmax+1:end/2,1)));   % Find minimum to the right of the main peak
Nmin_right = Nmax + help(1);

% Check if Nmin_right larger than Nmin_left
if Nmin_right < Nmin_right
    error('LS sampling point extraction failed');
end


%% LS estimation
% the parameters of a gaussian bell are estimated:
% y(n) = B*exp(-A*n^2) --> y_tilde(n) = ln(y(n)) = ln(B) - A*n^2 = C - A*n^2

% setting up the linear model for the LS estimation
y_tilde = log(abs(xspec(Nmin_left:Nmin_right,1)));
freq    = (f(Nmin_left:Nmin_right,1)-fmax).^2;
X = [-freq,ones(length(freq),1)];

% Actual LS estimation
lambda = inv(X'*X)*X'*y_tilde;

% Parameter extraction
B = exp(lambda(2,1));
A = lambda(1,1);

% get spectral filter curve
gaussSpec = B*exp(-A*(f(1:end/2,1)-fmax).^2);
gaussSpec = [gaussSpec;0;flip(gaussSpec(2:end,1))];     % flip for negative frequencies, otherwise filtered signal complex


%% Calculate sigma in time domain to be able to filter in time domain
sigma_f = 1/sqrt(2*A);              % sigma in frequency domain: gauss = exp(-0.5*(f-f0)^2/(sigma^2))
N_sigma_f = ceil(sigma_f/df);       % calculate number of bins of sigma_f


% sigma_f * sigma_t = 1/(2*pi) [Kiencke11]
sigma_t = 1/(2*pi*sigma_f);         % sigma in time domain
N_sigma_t = ceil(sigma_t/dt);       % calculate number of bins of sigma_t


%% construct gaussian Filter in time Domain
filter_length = 8*N_sigma_t + 1;   % length of filter in time domain, odd number of samples to make it symmetric around dirac
gaussSig = fspecial('gaussian',[filter_length,1],N_sigma_t);    % create gaussian with calculated sigma
% gaussSig = (gaussSig-min(gaussSig))/(max(gaussSig)-min(gaussSig));
t = [dt:dt:(filter_length-1)/2*dt]';        % construct time vector
df = Fs/filter_length;
f_filter = [-Fs/2 + df/2:df:Fs/2 - df/2]';

cosSig = [flip(cos(2*pi*fmax*t));1;cos(2*pi*fmax*t)];
diracSig = [zeros((filter_length-1)/2,1);1;zeros((filter_length-1)/2,1)];
gaussFilt = diracSig - 2*gaussSig.*cosSig;          % frequency shifting by fmax



RRiSig_filt = conv(RRiSig,gaussFilt,'same');
% RRiSig_filt = conv(flip(RRiSig_filt),gaussFilt,'same');

RRiSig_filt_spec = fft(RRiSig_filt,N_fft);

%% Plotting
dirac_spec =fftshift(fft(diracSig));
gauss_spec = fftshift(fft(gaussSig));
cos_spec  = fftshift(fft(cosSig));
gauss_filter_spec = fft(gaussFilt,N_fft);
gauss_filter_spec = gauss_filter_spec(1:end/2,1);

% figure(1)
% subplot(4,1,1)
% plot(f_filter,abs(gauss_spec))
% subplot(4,1,2)
% plot(f_filter,abs(dirac_spec))
% subplot(4,1,3)
% plot(f_filter,abs(cos_spec))
% subplot(4,1,4)
% plot(f_filter,abs(gauss_filter_spec))

