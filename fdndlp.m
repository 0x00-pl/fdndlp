function y = fdndlp(x, cfg)
%
% *****************************************************
%
% This program is the implementation of Variance-Normalizied Delayed Linear
% Prediction in time-frequency domain, which is aimed at speech
% dereverberation, known as weighted prediction error(WPE) method.
%
% Reference:
% [1] Nakatani T, Yoshioka T, Kinoshita K, et al. Speech Dereverberation 
% Based on Variance-Normalized Delayed Linear Prediction[J]. 
% IEEE Transactions on Audio Speech & Language Processing, 2010, 18(7):1717-1731.
%
% 
% Main parameters:
%
% mic_num                  the number of channels
% K                        the number of subbands
% F                        over-sampling rate
% N                        decimation factor
% D1                       subband preditction delay
% Lc1                      subband prediction order 
% eps                      lower bound of rho^2(Normalizaton factor)
%
%
% ***************************************************
% Created by Teng Xiang at Oct-14-2017 
% ***************************************************


% ***************************************************
% Load Parameters
% ***************************************************

num_mic = cfg.num_mic;
num_out = cfg.num_out;

K = cfg.K;                         % the number of subbands
N = cfg.N;                         % decimation factor

D1 = cfg.D1;                       % subband preditction delay                     
Lc1 = cfg.Lc1;                     % subband prediction order 

eps = cfg.eps;                     % lower bound of rho(Normalizaton factor)

len = length(x);


% ***************************************************
% Subband Decomposition 
% ***************************************************

xk = cell(num_mic, 1);

for m = 1 : num_mic
    tmp = stftanalysis(x(:,m), K, N);
    xk{m} = tmp;
end

LEN = size(xk{1}, 1);
dk = zeros(LEN, K, num_out);



% ***************************************************
% Subband NDLP
% ***************************************************
% 因为时域信号是实数，所以频域只需要考虑一半的频率
for k = 1 : K/2 + 1
%     fprintf('iRe=%d\n',k);
    xk_tmp = zeros(LEN+Lc1, num_mic);
    for m = 1 : num_mic 
        xk_tmp(Lc1+1:end,m) = xk{m}(:,k);
    end
    xk_tmp = xk_tmp.';
    x_buf = xk_tmp(1:num_out,Lc1+1:end).';
    X_BUF = zeros(num_mic * Lc1, LEN);
    for ii = 1 : LEN-D1
        xn_D = xk_tmp(:,ii+Lc1:-1:ii+1);
        X_BUF(:,ii+D1) = xn_D(:); 
    end 
    rho2 = max(mean(abs(x_buf(:,1:num_out)).^2, 2), eps);
    c_Err = cfg.iterations;

    while (c_Err > 1e-2)
        Lambda = diag(1 ./ rho2);
        Phi = X_BUF*Lambda*X_BUF';
        p = X_BUF*conj(x_buf./rho2(:,ones(1,num_out)));
        c = pinv(Phi)*p;
        dk(:,k,:) = (x_buf.' - c'*X_BUF).';
        rho2 = max(mean(squeeze(abs(dk(:,k,:)).^2),2), eps);
        c_Err = c_Err - 1;
    end
end
dk(:,K/2+2:end,:) = conj(dk(:,K/2:-1:2,:));


% ***************************************************
% Synthesis
% ***************************************************

y = [];
for m = 1 : num_out
    y_tmp = stftsynthesis(dk(:,:,m), K, N);
    y = [y, y_tmp];
end
y = y(1 : len, :);

