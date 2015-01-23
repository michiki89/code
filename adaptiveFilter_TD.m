function [RRiSig_filt,RRiSig_resp,resVar,sqError] = adaptiveFilter_TD(RRiSig,respSig,max_lag)
%% Adaptive Filter in Time Domain



%% Make sure RRi & Resp are the same length
if (length(RRiSig) ~= length(respSig))
    error('HR and Respiration timeseries must be the same length');
end

%% Make sure RRi is a column vector
[a,b] = size(RRiSig);
if (b>a)
    %x is a row vector -- fix this
    RRiSig = RRiSig';
end

%% Make sure Resp is a column vector
[a,b] = size(respSig);
if (b>a)
    %y is a row vector -- fix this
    respSig = respSig';
end

%% Definitions
N_samp = length(RRiSig);
x = respSig;
y = RRiSig;

%% Regression using respiration data as input and RRi as data to be estimated
max_lag_AR = 0;
max_lag_ex = max_lag;
max_lag = max([max_lag_AR,max_lag_ex]);


X = [zeros(N_samp-max_lag,max_lag_AR + max_lag_ex)];
indStart = max_lag+1;

for ii = indStart:N_samp
    y_in = flip(y(ii-max_lag_AR:ii-1,1));
    x_in = flip(x(ii-max_lag_ex:ii-1,1));
    for jj=1:max_lag_AR
        X(ii-max_lag,jj) = y_in(jj,1);
    end
    
    for jj=1:max_lag_ex
        X(ii-max_lag,jj+max_lag_AR) = x_in(jj,1);
    end
    
end


[coeff,~,res] = regress(y(indStart:end,1),X);
y_tilde = zeros(N_samp-max_lag,1);

for ii = indStart:N_samp
    y_in = flip(y(ii-max_lag_AR:ii-1,1));
    x_in = flip(x(ii-max_lag_ex:ii-1,1));
    for jj=1:max_lag_AR
        y_tilde(ii-max_lag,1) = y_tilde(ii-max_lag,1) + coeff(jj,1)*y_in(jj,1);
    end
    
    for jj=1:max_lag_ex
        y_tilde(ii-max_lag,1) = y_tilde(ii-max_lag,1) + coeff(jj+max_lag_AR,1)*x_in(jj,1);
    end
end


%% Create Outputs
RRiSig_resp = [mean(y_tilde)*ones(max_lag,1);y_tilde];%interp1([max_lag+1:1:N_samp]',y_tilde,[1:N_samp],'pchip','extrap');
RRiSig_filt = [mean(res)*ones(max_lag,1);res];%interp1([max_lag+1:1:N_samp]',res,[1:N_samp],'pchip','extrap');
resVar = var(RRiSig_filt);
sqError = mean((y(max_lag+1:end,1)-y_tilde).^2);

