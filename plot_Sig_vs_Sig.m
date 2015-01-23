%% Function to plot time series and phase of time series against each other

function [signalPlot] = plot_Sig_vs_Sig(signal1,signal2,time)
%% Make sure RRi & Resp are the same length
if (length(signal1) ~= length(signal2))
    error('Time series must have the same length');
end

signal1Ang = angle(hilbert(signal1));
signal2Ang = angle(hilbert(signal2));


signalPlot = figure(1);


ax1 = subplot(3,2,[1,2]);
plot(time,signal1)
xlabel('time [s]')
ylabel('Amplitude of signal 1')
grid on
box on

ax2 = subplot(3,2,[3,4]);
plot(time,signal2)
xlabel('time [s]')
ylabel('Amplitude of signal 2')
grid on
box on

ax3 = subplot(3,2,5);
plot(signal1,signal2)
xlabel('Amplitude of signal 1')
ylabel('Amplitude of signal 2')
grid on
box on

ax4 = subplot(3,2,6);
plot(signal1Ang,signal2Ang,'or')
xlabel('Phase of signal 1 [rad]')
ylabel('Phase of signal 2 [rad]')
grid on
box on

linkaxes([ax1,ax2],'x')
