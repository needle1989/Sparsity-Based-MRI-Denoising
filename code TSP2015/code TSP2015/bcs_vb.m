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
% ʱ�䣺2010��9��
% ���ߣ�����
% ��λ���人��ѧ������ϢѧԺ
% ���䣺yuleiwhu@gmail.com
%
% ������������ѧϰʹ�ã���������ϵ���ߡ��������κ�����ͽ��飬��ӭ�������ۣ�
% Copyright (C) LEI YU 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% ��������
if nargin < 2
    error('��磬�����ٸ�������������������');
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


%% ģ�Ͳ�����ʼ��
a = hyper.a; b = hyper.b;
c = hyper.c; d = hyper.d;
alpha0 = hyper.alpha0;
alpha = hyper.alpha*inf;

%% ����������Ҫ�õ��ı���
PTP = Phi'*Phi;
Pty = Phi'*y;
Yy = y'*y;
mu = Phi'*y;%ones(size())
Sigma = zeros(N,N);
res = inf*ones(options.iter+1,1);
%% VB��������
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
    %% ������ֹ����

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

