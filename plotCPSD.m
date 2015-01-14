function [figHandle,h,f] = plotCPSD(y,Fs,N_fft,N_window,N_overlap)
%% Function to plot all cross spectra of a multivariate stochastic process

%---- Inputs ----%
% y:    Function to be plotted (#observations X #variables)
% Fs: sampling frequency
% Zeropadding/# of frequency samples: if not chosen N_fft = 2^nextpow2(N_window)
% N_window: length of window for welch PSD estimation
% N_overlap: # of samples which overlap between windows

% default use: [figHandle,h,f] = plotCPSD(y,1,512,length(t)/10,0);


%% Define frequency parameters
df = 1/Fs;      % frequency resolution
N_samp = length(y(:,1));       % # of samples/observations
N_variables = length(y(1,:));  % # of variables


%% Calculate and Plot CPSD
k = 0;
figHandle = figure(1);
for ii = 1:N_variables
    for jj = 1:N_variables
        k = k +1;
        
        [h(:,k),f] = cpsd(y(:,ii),y(:,jj),N_window,N_overlap,N_fft,Fs,'twosided');
        
        subplot(N_variables,N_variables,k);
        plot(f,abs(h(:,k)))
        grid on
        box on
        xlabel('frequency in Hz')
        ylabel('|PSD|^2')
        title(['S_{',mat2str(ii),mat2str(jj),'}'])
        
    end
end


    