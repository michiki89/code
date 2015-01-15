function [y,t,A1,A2,B] = getDing06(N_samp,Fs,y0)
%% Function to build outputs of Arma model in
% Faes10 : Testing Frequency-Domain Causality in Multivariate Time Series

% The Arma Model has an order of two and has three outputs: 
% [y1(n);y2(n);y3(n)] = A1*[y1(n-1);y2(n-1);y3(n-1);] + A2*[y1(n-2);y2(n-2);y3(n-2)] + B*[w1(n);w2(n);w3(n)];

% ----- Inputs -----
% yx0:      starting values yx0 = [yx(-2);yx(-1)];
% paramteres of Arma Model: roh,f1,f2
% N_samp: # of samples to simulate time series for

% ----- Outputs -----
% y:      time function of ARMA outputs (#Observations X #Outputs);
% ARMA coefficients: A1,A2,B


%% Calculate ARMA Matrices
A1 = [   0.9,    0;
        0.16,    0.8];
    
A2 = [  -0.5,      0;
        -0.2,   -0.5];
    
    
B = [1 0;0,1];

%% Define Noise

w1 = 1*randn(N_samp,1);
w2 = 1*randn(N_samp,1);


%% Calculate First two Steps with starting values
y(:,1) = A1*y0(:,2) + A2*y0(:,1) + B*[w1(1,1);w2(1,1)];
y(:,2) = A1*y(:,1) + A2*y0(:,2) + B*[w1(2,1);w2(2,1)];

for ii = 3:N_samp
    y(:,ii) = A1*y(:,ii-1) + A2*y(:,ii-2) + B*[w1(ii,1);w2(ii,1)];
end
    
y = y';
t = [0:N_samp-1]'/Fs;
