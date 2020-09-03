function [X Beta]=bmtl_DP_offgrid(A,B,Y,hyper)

%  Bayesian Multi-task Leaning based on DP
%  Implemented by Variational Bayes methods.
%
% Input:
%   A     - on grid Sensing Matrix, with size nd x N;
%   B       - derivation of A, with size nd x N;
%   Y       - Measurements data matrix, with size nd x M;
%
% Output:
%   X       - image data matirx, with size N x M;
%   hyper   - Some parameters possibily used in this programme.
% nd number of ULA; M measurements; N number of grids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
if nargin < 3
    error('error!');
end

% Done, problem dimensions
[Mnd,N] = size(A); % Dimension of measurement matrix
[nd,M] = size(Y);   % Dimension of measurement signal

if nargin < 4
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
%     hyper.pi = zeros(hyper.truncateL,1);
    hyper.mu_x = zeros(N,M);%zeros(N,M);%A' * Y;
    hyper.mu_beta = zeros(N,hyper.truncateL);
    for ii = 1:M
        hyper.sigma_x(:,:,ii) = eye(N,N);
        hyper.sigma_beta (:,:,ii) = zeros(N,N);
        hyper.alpha(:,ii) =  (1 ./ abs(A((ii-1)*nd+1:ii*nd,:)' * Y(:,ii)));%ones(N,1)*inf
        hyper.mu_x(:,ii)=A((ii-1)*nd+1:ii*nd,:)' * Y(:,ii);
    end
   
    hyper.kappa = abs(ones(hyper.truncateL,M));
    hyper.kappa = hyper.kappa ./ repmat(sum(hyper.kappa,1),hyper.truncateL,1);
%     hyper.kappa = zeros(hyper.truncateL,M);
%     hyper.kappa(1,:) = ones(1,M);
end

if nargin < 4
    options.tol = 5e-3;% 0.005 may be suitable
    options.iter = 100;
end

%% Initialization
K = hyper.truncateL;
a = hyper.a; b = hyper.b;
c = hyper.c; d = hyper.d;
cnew = c * ones( N , K ); dnew = d * ones( N , K );
% cnew = ones(hyper.truncateL,1); dnew = ones(N,hyper.truncateL);
e = hyper.e; f = hyper.f;
z = hyper.z; 
% Pi = hyper.pi;
tao1_new = hyper.tao1; tao2_new =  hyper.tao2;
anew = a;
bnew = b;
lambda = hyper.lambda;
alpha0 = hyper.alpha0;
alpha = hyper.alpha;
mu_x = hyper.mu_x;
mu_beta=hyper.mu_beta;
Sigma_beta = hyper.sigma_beta;
Sigma_x = hyper.sigma_x;
kappa = hyper.kappa;
%% Main Iterations
for mmm=1:options.iter;
    mu_old=mu_x; 
    fprintf(1,'This is the %dth iteration.\n',mmm);  % iteration number 
%% Update mu
for ii = 1 : M;
    temp_sigma_x = zeros(N,N);
    temp_mu_x = zeros(N,1);
    BHB=B((ii-1)*nd+1:ii*nd,:)'*B((ii-1)*nd+1:ii*nd,:);
    for jj= 1 : K
        Phi_k = A((ii-1)*nd+1:ii*nd,:)+B((ii-1)*nd+1:ii*nd,:)*diag(mu_beta(:,jj));
        PhikHPhik = Phi_k'*Phi_k;
        temp_sigma_x = temp_sigma_x + kappa(jj,ii)*(diag(alpha(:,jj))+alpha0*(PhikHPhik+BHB.*Sigma_beta(:,:,jj)));
        temp_mu_x = temp_mu_x + kappa(jj,ii)* Phi_k'* Y( : , ii );
    end
    Sigma_x(:,:,ii) = inv(temp_sigma_x);
    mu_x(:,ii) = alpha0 * Sigma_x(:,:,ii) *  temp_mu_x;
end  
%% Update beta
for jj = 1 : K;
%     temp_sigma_beta = zeros(N,N);
%     temp_mu_beta = zeros(N,1);
%     for ii= 1 : M
%         BHB=B((ii-1)*nd+1:ii*nd,:)'*B((ii-1)*nd+1:ii*nd,:);
%         BHA=B((ii-1)*nd+1:ii*nd,:)'*A((ii-1)*nd+1:ii*nd,:);
%         temp_sigma_beta = temp_sigma_beta + kappa(jj,ii)*(BHB.*(mu_x(:,ii)*mu_x(:,ii)'+Sigma_x(:,:,ii)));
%         temp_mu_beta = temp_mu_beta + kappa(jj,ii)* (diag(mu_x(:,ii))*B((ii-1)*nd+1:ii*nd,:)'*(Y( : , ii )-A((ii-1)*nd+1:ii*nd,:)*mu_x(:,ii))-diag(BHA*Sigma_x(:,:,ii)));
%     end
%     Sigma_beta(:,:,jj) = inv(alpha0*temp_sigma_beta+diag(alpha(:,jj)));
%     mu_beta(:,jj) = alpha0 * real(Sigma_beta(:,:,jj) )* real( temp_mu_beta);
    mu_beta(:,jj) = zeros(N,1); 
    Sigma_beta(:,:,jj)= zeros(N,N);
end 

%% Update alpha
cnew = abs( repmat ( ( c + sum(kappa,2)).' , N , 1 ));
for jj = 1 : M
    temp(: , jj) = abs( mu_x(: , jj).^2) + real(diag(Sigma_x(: , : , jj)));
end
for jj = 1 : K
    dnew(: , jj) = sum(temp * diag(kappa( jj , : )) , 2);%+abs( mu_beta(: , jj).^2) + real(diag(Sigma_beta(: , : , jj)));
end
dnew = dnew + d;
alpha = abs(cnew ./ (dnew + eps));
%% plot
if mod(mmm,20)==0
  figure; imagesc(abs(mu_x));
  figure; imagesc(abs(1 ./ alpha));
end
%% Update z
ita = zeros(K,M);
for ii = 1 : M;
    AHA = A((ii-1)*nd+1:ii*nd,:)'*A((ii-1)*nd+1:ii*nd,:);
    BHB = B((ii-1)*nd+1:ii*nd,:)'*B((ii-1)*nd+1:ii*nd,:);
    AHB = A((ii-1)*nd+1:ii*nd,:)'*B((ii-1)*nd+1:ii*nd,:);
    ita_1 = norm(Y( : , ii )-A((ii-1)*nd+1:ii*nd,:)*mu_x(:,ii),'fro')^2+trace(AHA*Sigma_x(:,:,ii));
    ita_2 = Sigma_x(:,:,ii);%mu_x(:,ii)*mu_x(:,ii)'+
    for jj = 1 : K;
%         ita(jj,ii) = ita_1+trace((BHB.*(mu_beta(:,jj)*mu_beta(:,jj)'+Sigma_beta(:,:,jj)))*(ita_2))...
%             + 2*real(Y( : , ii )'*B((ii-1)*nd+1:ii*nd,:)*diag(mu_beta(:,jj))*mu_x(:,ii)-trace(AHB*diag(mu_beta(:,jj))*(ita_2)));
        ita(jj,ii) = norm(Y( : , ii )-(A((ii-1)*nd+1:ii*nd,:)+B((ii-1)*nd+1:ii*nd,:)*diag(mu_beta(:,jj)))*mu_x(:,ii),'fro')^2 ...
        + trace((A((ii-1)*nd+1:ii*nd,:)+B((ii-1)*nd+1:ii*nd,:)*diag(mu_beta(:,jj)))'*(A((ii-1)*nd+1:ii*nd,:)+B((ii-1)*nd+1:ii*nd,:)*diag(mu_beta(:,jj)))*ita_2);
        if jj == 1;
        temp1(jj,ii) = 0;
        temp2(jj,ii) = psi(tao1_new(1) + eps) - psi(tao1_new(1) + tao2_new(1) + eps);
        elseif jj == K
        temp1(jj,ii) = sum ( psi( tao2_new(1 : jj-1) + eps ) - psi( tao1_new(1 : jj-1) + tao2_new(1 : jj-1) + eps) );
        temp2(jj,ii) = 0;
        else
        temp1(jj,ii) = sum ( psi( tao2_new(1 : jj-1) + eps ) - psi( tao1_new(1 : jj-1) + tao2_new(1 : jj-1) ) + eps );
        temp2(jj,ii) = psi(tao1_new(jj)+ eps ) - psi(tao1_new(jj) + tao2_new(jj) + eps );            
        end
        temp3(jj,ii) = real(nd*(log(anew)-log(bnew)))-alpha0*real(ita(jj,ii));
        temp3(jj,ii) = 0;

        temp4(jj,ii) = real( sum( log((cnew(:,jj)+eps) ) - log((dnew(:,jj)) ) ))-real( trace( Sigma_x(:,:,ii) * diag( alpha(:,jj) ) ) + mu_x(:,ii)' * diag(alpha(:,jj)) * mu_x(:,ii) );
        temp5(jj,ii) = temp1(jj,ii) + temp2(jj,ii) + temp3(jj,ii) + temp4(jj,ii);
    end
end
%         psi((cnew) )
%         log((cnew) )
temp6 = exp( temp5 - repmat( max( temp5,[],1 ) , K , 1));
kappa = temp6 ./ (repmat( sum( temp6 , 1) , K , 1 ) );
% kappa
% figure; imagesc(abs(kappa))

%% Update tao
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
lambda = enew / (fnew + eps)*1.2;
%% Upate \alpha_0
anew = a + M * nd;
bnew = b +(sum(sum(kappa.*real(ita),1)));
alpha0 = abs(anew / bnew)/1.2;
%% stop criterion
% max(max(abs(mu_x-mu_old))) 
% if max(max(abs(mu_x-mu_old))) < options.tol
%     for ii=1:M
%     pos=find(1>alpha(:,ii) & alpha(:,ii)>0.01);
%     alpha(pos,ii)=1;
%     end
%     break;
% end

end
%  for ii=1:M
%     pos=find(1>alpha(:,ii) & alpha(:,ii)>0.01);
%     alpha(pos,ii)=1;
%  end
% for ii = 1 : M;
%     temp_sigma_x = zeros(N,N);
%     BHB=B((ii-1)*nd+1:ii*nd,:)'*B((ii-1)*nd+1:ii*nd,:);
%     for jj= 1 : K
%         Phi_k = A((ii-1)*nd+1:ii*nd,:)+B((ii-1)*nd+1:ii*nd,:)*diag(mu_beta(:,jj));
%         PhikHPhik = Phi_k'*Phi_k;
%         temp_sigma_x = temp_sigma_x + kappa(jj,ii)*(diag(alpha(:,jj))+alpha0*(PhikHPhik+BHB.*Sigma_beta(:,:,jj)));
%         temp_mu_x = temp_mu_x + kappa(jj,ii)* Phi_k'* Y( : , ii );
%     end
%     Sigma_x(:,:,ii) = temp_sigma_x;
%     mu_x(:,ii) = alpha0 * inv(Sigma_x(:,:,ii)) *  temp_mu_x;
% end  
X = mu_x;
Beta = mu_beta;
end
% function value = iszeros(theta,noise_level)
% theta = abs(theta);
% value = theta<=noise_level;
% end

