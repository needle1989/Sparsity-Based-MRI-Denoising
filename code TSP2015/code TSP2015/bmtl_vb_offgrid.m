function [X beta]=bmtl_vb_offgrid(Phi,B,Y,resolution)

%  Bayesian Multi-task Leaning based on DP
%  Implemented by Variational Bayes methods.
%
% Input:
%   Phi     - Sensing Matrix, with size Nd x N;
%   Y       - Measurements data matrix, with size Nd x M;
%
% Output:
%   X       - image data matirx, with size N x M;
%   hyper - Some parameters possibily used in this programme.
% nd number of ULA; M measurements; N number of grids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Done, problem dimensions
[Mnd,N] = size(Phi); % Dimension of measurement matrix
[nd,M] = size(Y);   % Dimension of measurement signal
    hyper.truncateL = M;
    hyper.a = (1/std(Y(:,1)))^2;%11e-6
    hyper.b = 1e-2;%1;
%     hyper.c = ones(N,hyper.truncateL)*1e-6;
%     hyper.d = ones(N,hyper.truncateL)*1e-6;
    hyper.c =1e-6 ;%1
    hyper.d = 1e-6;
    hyper.e = 1e-6;
    hyper.f = 1e-6;
    hyper.tao1 = ones(hyper.truncateL,1);
    hyper.tao2 = ones(hyper.truncateL,1);
    hyper.alpha0 = hyper.a / hyper.b;%1 / 1e-1;hyper.a / hyper.b 
%     [hyper.alpha,temp] = kmeans(hyper.alpha.',hyper.truncateL);
    hyper.z = zeros(hyper.truncateL,M);
    hyper.lambda = 1;
    hyper.pi = zeros(hyper.truncateL,1);
    hyper.mu = zeros(N,M);%zeros(N,M);%Phi' * Y;
    hyper.alpha=zeros(N,1);
    for ii = 1:M
        hyper.sigma(:,:,ii) = ones(N,N);
        hyper.alpha=  hyper.alpha+(abs(Phi((ii-1)*nd+1:ii*nd,:)' * Y(:,ii))); 
    end
    hyper.alpha=M./hyper.alpha;

    options.tol = 1e-5;% 0.005 may be suitable
    options.iter = 300;
%% Initialization
a = hyper.a; b = hyper.b;
c = hyper.c; d = hyper.d;
alpha0 = hyper.alpha0;
alpha = hyper.alpha;
mu = hyper.mu;
beta = zeros(N,1);
Sigma = hyper.sigma;
r=resolution*pi/180;
N_DOA=4;
%% Main Iterations

for mmm=1:options.iter;
   mu_old=mu; 
fprintf(1,'This is the %dth iteration.\n',mmm);  % iteration number 

%% Update mu
for ii = 1 : M;
    C = 1 / alpha0 * eye(nd) + (Phi((ii-1)*nd+1:ii*nd,:)+B((ii-1)*nd+1:ii*nd,:)*diag(beta)) * inv(diag(alpha)) * (Phi((ii-1)*nd+1:ii*nd,:)+B((ii-1)*nd+1:ii*nd,:)*diag(beta))';
    Cinv = inv(C);
    Sigma(:,:,ii) = inv(diag(  alpha )) - inv( diag( alpha ) ) * (Phi((ii-1)*nd+1:ii*nd,:)+B((ii-1)*nd+1:ii*nd,:)*diag(beta))' * Cinv * (Phi((ii-1)*nd+1:ii*nd,:)+B((ii-1)*nd+1:ii*nd,:)*diag(beta)) * inv(diag( alpha ));
    
    mu(:,ii) = alpha0 * Sigma(:,:,ii) * (Phi((ii-1)*nd+1:ii*nd,:)+B((ii-1)*nd+1:ii*nd,:)*diag(beta))' * Y( : , ii );
end   
% if mod(mmm,60) == 0
%     beta*180/pi
%     figure; imagesc(abs(mu))
% end

%% Update alpha
cnew = c+M;
dnew = M*mean(abs( mu.^2),2);
temp=zeros(N,1);
for jj = 1 : M
    temp = temp+diag(Sigma(: , : , jj));
end
dnew = dnew + d+temp;
alpha = abs(cnew ./ (dnew + eps));
p_alpha=abs(1 ./ alpha);
%% Upate \alpha_0
anew = a + M * nd;
bnew = b;
for ii = 1 : M
    bnew = bnew +abs(norm( Y(:,ii) - (Phi((ii-1)*nd+1:ii*nd,:)+B((ii-1)*nd+1:ii*nd,:)*diag(beta)) * mu(:,ii), 'fro')^2)+ real(trace((Phi((ii-1)*nd+1:ii*nd,:)+B((ii-1)*nd+1:ii*nd,:)*diag(beta))  * Sigma(:,:,ii) * (Phi((ii-1)*nd+1:ii*nd,:)+B((ii-1)*nd+1:ii*nd,:)*diag(beta))'));
end

alpha0 = abs(anew / bnew)/1.3;%1

%% update beta
    [temp, idx] = sort(p_alpha, 'descend');
    idx = sort(idx(1:N_DOA),'ascend');
    diff_idx = diff(idx);
    if any(diff_idx)==1
        pos_diff=find(diff_idx==1);
        for nnn=1:length(pos_diff)
        if p_alpha(idx(pos_diff(nnn)))>=p_alpha(idx(pos_diff(nnn)+1))
           idx(pos_diff(nnn)+1)= idx(pos_diff(nnn));
        else
           idx(pos_diff(nnn))= idx(pos_diff(nnn)+1);
        end
        end
    end
    idx = unique(idx);
    diff_idx = diff(idx);
     if any(diff_idx)==1
        pos_diff=find(diff_idx==1);
        for nnn=1:length(pos_diff)
        if p_alpha(idx(pos_diff(nnn)))>=p_alpha(idx(pos_diff(nnn)+1))
           idx(pos_diff(nnn)+1)= idx(pos_diff(nnn));
        else
           idx(pos_diff(nnn))= idx(pos_diff(nnn)+1);
        end
        end
     end
    idx = unique(idx);
     diff_idx = diff(idx);
     if any(diff_idx)==1
        pos_diff=find(diff_idx==1);
        for nnn=1:length(pos_diff)
        if p_alpha(idx(pos_diff(nnn)))>=p_alpha(idx(pos_diff(nnn)+1))
           idx(pos_diff(nnn)+1)= idx(pos_diff(nnn));
        else
           idx(pos_diff(nnn))= idx(pos_diff(nnn)+1);
        end
        end
     end
    idx = unique(idx);
    N_DOA_used=length(idx);
    temp = beta;
    beta = zeros(N,1);
    beta(idx) = temp(idx);
    temp_p_beta = zeros(N_DOA_used,N_DOA_used);
    temp_v_beta = zeros(N_DOA_used,1);
    for ii= 1 : M
        BHB=B((ii-1)*nd+1:ii*nd,idx)'*B((ii-1)*nd+1:ii*nd,idx);
        temp_p_beta = temp_p_beta + conj(BHB).*(mu(idx,ii)*mu(idx,ii)');
        temp_v_beta = temp_v_beta + conj(mu(idx,ii)).*(B((ii-1)*nd+1:ii*nd,idx)'*(Y( : , ii)-Phi((ii-1)*nd+1:ii*nd,:)*mu(:,ii)));
    end
%  temp1 = real(temp_p_beta)\ real( temp_v_beta);

    if any(abs(diag(temp_p_beta)) <0.000001)
            for n = 1:N_DOA_used               
                if abs (temp_p_beta(n,n)) < 0.000001
                   beta(idx(n)) = 0;
                else
                temp_beta = beta(idx);
                temp_beta(n) = 0;
                beta(idx(n)) = (temp_v_beta(n) - temp_p_beta(n,:) * temp_beta) / temp_p_beta(n,n);
                end
            end
    else
       beta = zeros(N,1);
       beta(idx) = real(temp_p_beta)\ real( temp_v_beta);
    end
    for n=1:N_DOA_used
                if beta(idx(n)) > r/2
                    beta(idx(n)) = r/2;
                end
                if beta(idx(n)) < -r/2
                    beta(idx(n)) = -r/2;
                end
    end
    beta= real(beta);

%% stop criterion

if max(max(abs(mu-mu_old))) < options.tol
    break;
end

end
X = mu;

