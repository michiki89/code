function [y,t,A,B] = getBaccala(N_samp,Fs,y0)
%% Function to generate the ARMA-Model from
% Baccala00: Partial Directed Coherence: a new concept in neural science
% determination

%---- Input ----
% N_samp,Fs: sampling points and sampling frequency
% y0:        start values of output functions

%---- Output ----
% y: output function (1 X #output functions)
% t:        start values of output functions
% A,B:  ARMA Matrices

%% Define Noise
w1 = 1*randn(N_samp,1);
w2 = 1*randn(N_samp,1);
w3 = 1*randn(N_samp,1);

%% Define Matrices
A  = [  0.5,    0.3,    0.4;
        -0.5,   0.3,    1;
        0,      -0.3,   -0.2];
    
B = [1,0,0;
     0,1,0;
     0,0,1];
 

%% Create Output y

y(:,1) = A*y0(1,:)' + B*[w1(1,1);w2(1,1);w3(1,1)];

for ii = 2:N_samp
    
    y(:,ii) = A*y(:,ii-1) + B*[w1(ii,1);w2(ii,1);w3(ii,1)];
    
end

y = y';
t = [0:N_samp-1]'/Fs;