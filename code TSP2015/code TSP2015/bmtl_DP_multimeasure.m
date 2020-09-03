function [X alpha_inv alpha0_inv]=bmtl_DP_multimeasure(Phi,Y,hyper)

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
% nd number of ULA; M subbands; N number of grids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
if nargin < 2
    error('error!');
end

% Done, problem dimensions
[Mnd,N] = size(Phi); % Dimension of measurement matrix
 N_measure = length(Y);
[nd,M] = size(Y{1});   % Dimension of measurement signal

if nargin < 3
    hyper.truncateL = M;
    hyper.a = 1e-4;
    hyper.b = 1e-4;
%     hyper.c = ones(N,hyper.truncateL)*1e-6;
%     hyper.d = ones(N,hyper.truncateL)*1e-6;
    hyper.c = 1e-4;
    hyper.d = 1e-4;
    hyper.e = 1e-4;
    hyper.f = 1e-4;
    hyper.tao1 = ones(hyper.truncateL,1);
    hyper.tao2 = ones(hyper.truncateL,1);
    hyper.alpha0 = 1;% (1/std(Y(:)))^2*1e3/ 1e-1;hyper.a / hyper.b 
%     [hyper.alpha,temp] = kmeans(hyper.alpha.',hyper.truncateL);
    hyper.z = zeros(hyper.truncateL,M);
    hyper.lambda = 1;
    hyper.pi = zeros(hyper.truncateL,1);
    
    for ii = 1:M
        hyper.sigma(:,:,ii) = ones(N,N);
        hyper.alpha(:,ii) =  (1 ./ abs(Phi((ii-1)*nd+1:ii*nd,:)' * Y{1}(:,ii)));%ones(N,1)*inf
        for nnn=1:N_measure
            hyper.mu{nnn}(:,ii)=Phi((ii-1)*nd+1:ii*nd,:)' * Y{nnn}(:,ii);
        end
    end
    hyper.kappa = abs(ones(hyper.truncateL,M));
    hyper.kappa = hyper.kappa ./ repmat(sum(hyper.kappa,1),hyper.truncateL,1);
%     hyper.kappa = zeros(hyper.truncateL,M);
%     hyper.kappa(1,:) = ones(1,M);
end

if nargin < 4
    options.tol = 0.005;% 0.005 may be suitable
    options.iter = 500;
end

%% Initialization
K = hyper.truncateL;
a = hyper.a; b = hyper.b;
c = hyper.c; d = hyper.d;
cnew = c * ones( N , K ); dnew = d * ones( N , K );
% cnew = ones(hyper.truncateL,1); dnew = ones(N,hyper.truncateL);
e = hyper.e; f = hyper.f;
z = hyper.z; Pi = hyper.pi;
tao1_new = hyper.tao1; tao2_new =  hyper.tao2;
lambda = hyper.lambda;
alpha0 = hyper.alpha0;
alpha = hyper.alpha;
mu = hyper.mu;
% figure
% imagesc(abs(mu));
Sigma = hyper.sigma;
kappa = hyper.kappa;
%% Main Iterations

for mmm=1:options.iter;
    mu_old=mu; 
%     fprintf(1,'This is the %dth iteration.\n',mmm);  % iteration number 

%% Update mu

for ii = 1 : M;
    C = 1 / alpha0 * eye(nd) + Phi((ii-1)*nd+1:ii*nd,:) * inv(diag( sum (  alpha * diag( kappa(:,ii) ), 2))) * Phi((ii-1)*nd+1:ii*nd,:)';
    Cinv = inv(C);
    Sigma(:,:,ii) = inv(diag( sum (  alpha * diag( kappa(:,ii) ), 2))) - inv( diag( sum (  alpha * diag( kappa(:,ii) ), 2)) ) * Phi((ii-1)*nd+1:ii*nd,:)' * Cinv * Phi((ii-1)*nd+1:ii*nd,:) * inv(diag( sum (  alpha * diag( kappa(:,ii) ), 2)));
    for nnn=1:N_measure
    mu{nnn}(:,ii) = alpha0 * Sigma(:,:,ii) * Phi((ii-1)*nd+1:ii*nd,:)' * Y{nnn}( : , ii );
    end
end  
%% Update alpha
cnew = abs( repmat ( ( c + sum(kappa,2)*N_measure ).' , N , 1 ));
temp = zeros(N,M);
for jj = 1 : M
    for nnn=1:N_measure
    temp(: , jj) = temp(: , jj)+abs( mu{nnn}(: , jj).^2) + diag(Sigma(: , : , jj));
    end
end
for jj = 1 : K
    dnew(: , jj) = sum(temp * diag(kappa( jj , : )) , 2);
end
dnew = dnew + d;
alpha = abs(cnew ./ (dnew + eps));

% for ii = 1 : M
% z_r = iszeros(mu(:,ii),1/alpha0/3);
% alpha (z_r) = 1;
% end

% figure; imagesc(abs(temp));
% alpha
%     if mod(mmm,50) == 0
% %     figure; imagesc(abs(mu{1}));
%     figure; imagesc(abs(1 ./ alpha));
%     end
% if mmm > 200
%     for ii=1:M
%     pos=find(1>alpha(:,ii) & alpha(:,ii)>0.01);
%     alpha(pos,ii)=1;
%     end
% end
%% Update z

for jj = 1 : M;
    for ii = 1 : K;
        if ii == 1;
        temp1(ii,jj) = 0;
        temp2(ii,jj) = psi(tao1_new(1) + eps) - psi(tao1_new(1) + tao2_new(1) + eps);
        elseif ii == K
        temp1(ii,jj) = sum ( psi( tao2_new(1 : ii-1) + eps ) - psi( tao1_new(1 : ii-1) + tao2_new(1 : ii-1) + eps) );
        temp2(ii,jj) = 0;
        else
        temp1(ii,jj) = sum ( psi( tao2_new(1 : ii-1) + eps ) - psi( tao1_new(1 : ii-1) + tao2_new(1 : ii-1) ) + eps );
        temp2(ii,jj) = psi(tao1_new(ii)+ eps ) - psi(tao1_new(ii) + tao2_new(ii) + eps );            
        end
        temp3(ii,jj) = N_measure*real( sum(log((cnew(:,ii)+eps) ) - log((dnew(:,ii)) ) ));
%         psi((cnew(:,ii)) )
%         log((cnew(:,ii)) 
        temp4(ii,jj)=0;
        for nnn=1:N_measure
        temp4(ii,jj) = temp4(ii,jj)-real( trace( Sigma(:,:,jj) * diag( alpha(:,ii) ) ) + mu{nnn}(:,jj)' * diag(alpha(:,ii)) * mu{nnn}(:,jj) );
        end
        temp5(ii,jj) = temp1(ii,jj) + temp2(ii,jj) + temp3(ii,jj) + temp4(ii,jj);

    end
end
%         psi((cnew) )
%         log((cnew) )
temp6 = exp( temp5 - repmat( max( temp5,[],1 ) , K , 1));
kappa = temp6 ./ (repmat( sum( temp6 , 1) , K , 1 ) );
% kappa
% figure; imagesc(abs(kappa))

%% Update pi
tao1 = 1;
tao2 = lambda;
for ii = 1 : K-1;
    tao1_new(ii) = sum(kappa(ii,:));
%     tao2_new(ii) = (K - ii) * sum(kappa(ii,:));
    tao2_new(ii) = sum(sum(kappa(ii+1:end,:)));
end
tao1_new = tao1_new + tao1;
tao2_new = tao2_new + tao2;
tao1_new(ii + 1) = 1;
tao2_new(ii + 1) = 1;

%% Update \lambda
enew = e + K - 1;
fnew = f - sum( psi(tao2_new(1:K-1) + eps) - psi(tao1_new(1:K-1) + tao2_new(1:K-1) + eps));
lambda = enew / (fnew + eps)*1.1;
%% Upate \alpha_0
anew = a + M * nd * N_measure;
bnew = b;    
for nnn=1:N_measure
for ii = 1 : M
    bnew = bnew +abs(norm( Y{nnn}(:,ii) - Phi((ii-1)*nd+1:ii*nd,:) * mu{nnn}(:,ii), 'fro')^2)+ trace(Phi((ii-1)*nd+1:ii*nd,:)  * Sigma(:,:,ii) * Phi((ii-1)*nd+1:ii*nd,:)');
end
end
alpha0 = abs(anew / bnew)/1.2;
%% stop criterion
for nnn=1:N_measure
    cri(nnn)=max(max(abs(mu{nnn}-mu_old{nnn})));
end

if max(cri) < options.tol
    for ii=1:M
    pos=find(1>alpha(:,ii) & alpha(:,ii)>0.001);
    alpha(pos,ii)=1;
    end
    break;
end

end
for ii=1:M
    pos=find(1>alpha(:,ii) & alpha(:,ii)>0.001);
    alpha(pos,ii)=1;
end
for ii = 1 : M;
    C = 1 / alpha0 * eye(nd) + Phi((ii-1)*nd+1:ii*nd,:) * inv(diag( sum (  alpha * diag( kappa(:,ii) ), 2))) * Phi((ii-1)*nd+1:ii*nd,:)';
    Cinv = inv(C);
    Sigma(:,:,ii) = inv(diag( sum (  alpha * diag( kappa(:,ii) ), 2))) - inv( diag( sum (  alpha * diag( kappa(:,ii) ), 2)) ) * Phi((ii-1)*nd+1:ii*nd,:)' * Cinv * Phi((ii-1)*nd+1:ii*nd,:) * inv(diag( sum (  alpha * diag( kappa(:,ii) ), 2)));
    for nnn=1:N_measure
    mu{nnn}(:,ii) = alpha0 * Sigma(:,:,ii) * Phi((ii-1)*nd+1:ii*nd,:)' * Y{nnn}( : , ii );
    end
end  
X = mu;
pos_temp=find(sum(alpha)~= N);
alpha_inv=1./abs(alpha(:,pos_temp));
alpha0_inv=1/abs(alpha0);
end
% function value = iszeros(theta,noise_level)
% theta = abs(theta);
% value = theta<=noise_level;
% end

