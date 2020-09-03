function [Xhat Theta]=DL_MTR(Y,M,grid,fs,T,d);
%% ----------initialization-----------
f=-fs/2:fs/T:fs/2-fs/T;
c=3*10^8;
d=c/fs/2;
N=length(grid);
A = zeros(T,M,N);
B = zeros(T,M,N);
iter=0;
tol=0.001;
maxiter=1000;
for m = 1:M
    for n = 1:N
        for t=1:T
        A(t,m,n) = exp(2*sqrt(-1) * pi * f(t)*d*(m-(M+1)/2) * cos(grid(n)/180*pi)/c);
        B(t,m,n) = -2*sqrt(-1) * pi * f(t)*d*(m-(M+1)/2) * sin(grid(n)/180*pi) * temp;
        end
    end
end
y=reshape(Y,N*T,1);
Xhat=zeros(N*T,1);
a = 1e-6;
b = 1e-6;
c = ones(N*T,1)*1e-6;
d = ones(N*T,1)*1e-6;
mu_old=zeros(N*T,1);
alpha0 = std(y)^2*1e-1;
for t=1:T
alpha((t-1)*N+1:t*N) = ones(N,1)*inf;
Sigma((t-1)*N+1:t*N,(t-1)*N+1:t*N) = zeros(N,N);
mu_old((t-1)*N+1:t*N) = A(t,:,:)'*y((t-1)*N+1:t*N);
end
theta_old=grid;
while ~converged
iter=iter+1;
% %% ------------------pruning variables as they go to zero-----------------
%     if (min(gamma) < PRUNE_GAMMA)
%         index = find(gamma > PRUNE_GAMMA);
%         usedNum = length(index);
%         keep_list = keep_list(index);
%         blkStartLoc = blkStartLoc(index);
%         blkLenList = blkLenList(index);
%         
%         % prune gamma and associated components in Sigma0 
%         gamma = gamma(index);
%         temp = Sigma0;
%         Sigma0 = [];
%         for k = 1 : usedNum
%             Sigma0{k} = temp{index(k)};
%         end
%         
%         % construct new Phi
%         temp = [];
%         for k = 1 : usedNum
%             temp = [temp, Phi0(:,blkStartLoc(k):blkStartLoc(k)+blkLenList(k)-1)];
%         end
%         Phi = temp;
%         %clear temp;
%     end

%% -------learning the sparse coefficient for each channel t=1:T--------
mu_old = mu;
for t=1:T   
PTP((t-1)*N+1:t*N,(t-1)*N+1:t*N) = A(t,:,:)'*A(t,:,:);
Pty((t-1)*N+1:t*N) = A(t,:,:)'*y((t-1)*N+1:t*N);

%% alpha
cnew = c+1/2; dnew = d+ (mu((t-1)*N+1:t*N).^2+diag(Sigma((t-1)*N+1:t*N,(t-1)*N+1:t*N)))/2;
alpha((t-1)*N+1:t*N) = cnew./dnew;
       
%% mu
Sigma((t-1)*N+1:t*N,(t-1)*N+1:t*N)  = inv(diag(alpha(((t-1)*N+1:t*N))) + alpha0*PTP((t-1)*N+1:t*N,(t-1)*N+1:t*N));
mu ((t-1)*N+1:t*N)= alpha0*Sigma((t-1)*N+1:t*N,(t-1)*N+1:t*N)*Pty((t-1)*N+1:t*N);
Xhat((t-1)*N+1:t*N)=mu((t-1)*N+1:t*N);
end
%% alpha0

bb = Sigma + mu*mu';
anew = a+M*T/2; bnew = b+ ( y'*y - 2*Pty'*mu + ones(1,N*T)*(bb.*PTP)*ones(N*T,1))/2;    
alpha0 = anew/bnew;
%% -------learning the directions Theta--------
D_theta_old=;
theta=theta_old-D_theta_old;
%% -------stopping---------------------
    if norm(alpha - alpha_last)/norm(alpha_last) < tol || iter >= maxiter
        converged = true;
    end
end

end