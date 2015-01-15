function [Y_est,Residual,sys,coeff] = armax_LS(input,output,time,Fs,N_ex,N_y)
%% Function to calculate the ARMA coefficients, using Least Squares

% --- Inputs ---
% input: exogeneous input of the ARMA model (column vector)
% output: autoregressive part (column vector)
% time, Fs: sampling points and sampling frequency
% N_ex,N_y: order of exogeneous part (N_ex) and AR part (N_y)



N_max = max([N_ex,N_y]);
Y = zeros(length(output)-N_max,1);
Y = Y + output(N_max + 1:end);
t_Y = time(N_max + 1:end);

for ii = N_max + 1: length(time)
    if N_y == 0
            X(ii-N_max,:) = [flip(input(ii-N_ex:ii-1,1))];   
    elseif N_ex == 0
            X(ii-N_max,:) = [flip(output(ii-N_y:ii-1,1))];
    else
            X(ii-N_max,:) = [flip(output(ii-N_y:ii-1,1));flip(input(ii-N_ex:ii-1,1))];
    end
end

coeff = inv(X'*X)*X'*Y;
Y_est = X*coeff;
Residual = Y - Y_est;

Y_est = interp1(t_Y,Y_est,time,'pchip');
Residual = interp1(t_Y,Residual,time,'pchip');


    if N_y == 0
            sys = tf(flip(coeff)',1,1/Fs);   
    elseif N_ex == 0
            sys = tf(1,flip(coeff)',1/Fs);
    else
            sys = tf(flip(coeff(N_y + 1:end))',flip(coeff(1:N_y))',1/Fs);
    end


