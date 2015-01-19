function [filter,b,a] = lsGauss_timeDomain(RRiSig,respSig,N_fft,N_window,N_overlap,Fs,fmax)
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


%% Sampling Parameter Definition
N_samp = length(RRiSig);    % number of samples of RRi and Respiration time series
df = Fs/N_fft;              % spectral resolution
f = [0:df:Fs-df]';          % spectral sampling points

%% CPSD estimation
xspec = cpsd(RRiSig,respSig,N_window,N_overlap,N_fft,'twosided');   % cross spectrum

%% LS sampling point extraction
[~,Nmax] = min(abs(f-fmax));        % find sample, where fmax in cross spectrum

% Find minima next to main lobe (assumption: zeropadding and welch averaging is large enough
% that PSD is smoothed good enough for robust minima detection)
[~,Nmin_left] = findpeaks(-abs(xspec(1:Nmax,1)),1,'last');   % Find minimum to the left of the main peak
[~,Nmin_right] = findpeaks(-abs(xspec(N_max+1:end/2,1)),1,'first');   % FInd minimum to the right of the main peak

% Check if Nmin_right larger than Nmin_left
if Nmin_right < Nmin_right
    error('LS sampling point extraction failed');
end

%% LS estimation
% the parameters of a gaussian bell are estimated:
% y(n) = B*exp(-A*n) --> y_tilde(n) = ln(y(n)) = ln(B) - A*n = C - A*n

% setting up the linear model for the LS estimation
y_tilde = ln(abs(xspec(Nmin_left:Nmin_right,1)));
freq    = f(Nmin_left:Nmin_right,1);
X = [-freq,ones(length(freq),1)];

% Actual LS estimation
lambda = inv(X'*X)*X'*y_tilde;

% Parameter extraction
B = exp(lambda(2,1));
A = lambda(1,1);

% get spectral filter curve
gaussFilt = B*exp(-A*f(1:end/2,1));
gaussFilt = [gaussFilt;0;flip(gaussFilt(1:end-1,1))];










