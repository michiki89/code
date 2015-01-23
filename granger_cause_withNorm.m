function [F,c_v,Fmax] = granger_cause_withNorm(x,y,alpha,max_lag,min_var,normShift)

%Make sure x & y are the same length
if (length(x) ~= length(y))
    error('x and y must be the same length');
end

%Make sure x is a column vector
[a,b] = size(x);
if (b>a)
    %x is a row vector -- fix this
    x = x';
end

%Make sure y is a column vector
[a,b] = size(y);
if (b>a)
    %y is a row vector -- fix this
    y = y';
end


%% Shift input vector to find max coupling
% signal_power = sum(x.^2);
% noise_power  = signal_power/(10^(SNR/10));
x_shifted = [x(end-normShift+1:end,1);x(1:end-normShift)]+ 0.01*randn(length(x),1);

[Fmax,~] = granger_cause(x_shifted,x,alpha,max_lag,min_var);

[F,c_v] = granger_cause(x,y,alpha,max_lag,min_var);

F = F/Fmax;




    
    
    
    

