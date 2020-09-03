function [x,errbar] = bcs_vb(Phi,y,hyper,options)
% This function is to solve the Bayesian Compressive Sensing via
% Variational Bayesian EM approach
%
% Input:
%   Phi     - Sensing Matrix, with size M x N;
%   y       - Measurements, with size M;
%
% Output:
%   x       - Solutions, with size N;
%   hyper - Some parameters possibily used in this programme.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 时间：2010年9月
% 作者：余磊
% 单位：武汉大学电子信息学院
% 邮箱：yuleiwhu@gmail.com
%
% 声明：本代码学习使用，商用请联系作者。若出现任何问题和建议，欢迎来信讨论！
% Copyright (C) LEI YU 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% 参数设置
if nargin < 2
    error('大哥，你至少给我两个参数撒！！！');
end
% Done, problem dimensions
[M,N] = size(Phi);

if nargin < 3
    hyper.a = 1e-6;
    hyper.b = 1e-6;
    hyper.c = ones(N,1)*1e-6;
    hyper.d = ones(N,1)*1e-6;
    hyper.alpha0 = 1;
    hyper.alpha =  ones(N,1)*inf;%ones(N,1)*0.1;
end

if nargin < 4
    options.tol = 1e-4;
    options.iter = 1000;
end

if ~isfield(hyper,'a'); hyper.a = 1e-6; end
if ~isfield(hyper,'b'); hyper.b = 1e-6; end
if ~isfield(hyper,'c'); hyper.c = ones(N,1)*1e-6; end
if ~isfield(hyper,'d'); hyper.d = ones(N,1)*1e-6; end
if ~isfield(hyper,'alpha0'); hyper.alpha0 = 1; end
if ~isfield(hyper,'alpha'); hyper.alpha = ones(N,1)*inf; end


%% 模型参数初始化
a = hyper.a; b = hyper.b;
c = hyper.c; d = hyper.d;
alpha0 = hyper.alpha0;
alpha = hyper.alpha*inf;

%% 迭代过程中要用到的变量
PTP = Phi'*Phi;
Pty = Phi'*y;
Yy = y'*y;
mu = Phi'*y;%ones(size())
Sigma = zeros(N,N);
res = inf*ones(options.iter+1,1);
%% VB迭代过程
for i = 1:options.iter;
    
    %% alpha
    cnew = c+1/2; dnew = d+ (mu.^2+diag(Sigma))/2;
    alpha = cnew./dnew;
       
    %% alpha0
    bb = Sigma + mu*mu';
    anew = a+M/2; bnew = b+ ( Yy - 2*Pty'*mu + ones(1,N)*(bb.*PTP)*ones(N,1))/2;
    
    alpha0 = anew/bnew;
    %% mu
    Sigma = inv(diag(alpha) + alpha0*PTP);
    mu = alpha0*Sigma*Pty;
    %% 迭代终止条件

    res(i+1) = max(abs(mu));
    if abs(res(i)-res(i+1))/res(i) < options.tol
        break;
    end
%      x_r_image=reshape(abs(mu),MM,NN);
%  if mod(i,5)==0
%      figure
%      imagesc((abs(x_r_image)))
%      xlabel('D cell')
%      ylabel('Range cell')
%      axis tight;
%   end
    %% display
%         pause(0.1);
%         figure(1000)
%         plot(mu);
end

errbar = sqrt(diag(Sigma));
x = mu;
end

