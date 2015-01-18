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
% [~,Nmin_left] = findpeaks(



