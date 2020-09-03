clear all;
close all;
clc;
%% %%%%%%%%%----------------parameters for Radar------------------%%%%%%%%%%
% %% ------------------------array signal--------------------------------- %%
% C=3*10^8;
% M=10; %number of ULA 
% K=2;
% theta = [68.78 108.42]';
% resolution = 2;
% grid = (0:resolution:180)';
% N = length(grid);
% i=sqrt(-1);
% %% ----------------------wideband signal generation--------------------- %%
% T=256;%signal length
% C_N=256;%number of channels
% fs=4096*10^6;
% B=[485 485]*10^6;
% r=B./(T/fs);
% t=(0:T-1)/fs;
% N_t=length(t);
% f0=[544 1024]*10^6;
% amp = [1 1]';
% S=zeros(K,T);
% for k=1:K
% S(k,:) = amp(k)*exp(2*i* pi *(f0(k)*t+0.5*r(k)*t.^2));
% end
% % S(1,:) = exp(2*i* pi *(100*t))+exp(2*i* pi *(200*t))+exp(2*i* pi *(300*t))+exp(2*i* pi *(400*t));%+exp(2*i* pi *(500*t))
% % S(2,:) = exp(2*i* pi *(145*t))+exp(2*i* pi *(2*145*t))+exp(2*i* pi *(3*145*t))+exp(2*i* pi *(4*145*t));%+exp(2*i* pi *(5*145*t))
% 
% % observed signal
% Y0 = (fft(S,C_N,2));
% f=[fs/C_N:fs/C_N:fs/2];
% C_N=C_N/2;
% Y00 = Y0(:,1:C_N);
% figure
% plot(f/10^6,abs(Y00(1,:)),'b*-')
% hold on
% plot(f/10^6,abs(Y00(2,:)),'ro-')
% xlabel('frequency channels (MHz)')
% ylabel('magnitude')
% axis tight;
% figure
% plot(abs(Y00(1,:)),'b*-')
% hold on
% plot(abs(Y00(2,:)),'ro-')
% lambda=C./f;
% f_ref=640*10^6;
% d=0.5*C/f_ref;
% x_true=zeros(N,C_N);
% x_true(29,34:64)=1;
% x_true(55,64:98)=1;
% figure
% imagesc(f/10^6,grid,x_true)
% xlabel('frequency channels (MHz)')
% ylabel('direction')
%% %%%%%%%%----------------parameters for Sonar-------------%%%%%%%%%%%%%%%
%% ------------------------array signal--------------------------------- %%
C=1500;
M=10; %number of ULA 
K=2;
theta = [58.5 108.5]';
%theta = [68.78 108.42]';
resolution = 1;
grid = (0:resolution:180)';
N = length(grid);
i=sqrt(-1);
%% ----------------------wideband signal generation--------------------- %%
T=256;%signal length
C_N=256;%number of channels
fs=4096;
C_f=fs/C_N;
B=[1024 1024];
r=B./(T/fs);
t=(0:T-1)/fs;
N_t=length(t);
f0=[896 896];
%f0=[1024-256+128 1024-256+128];
%f0=[1024-256+128 1408];
amp = [1 1]';
S=zeros(K,T);
for k=1:K
S(k,:) = amp(k)*exp(2*i* pi *(f0(k)*t+0.5*r(k)*t.^2));
end
% S(1,:) = exp(2*i* pi *(100*t))+exp(2*i* pi *(200*t))+exp(2*i* pi *(300*t))+exp(2*i* pi *(400*t));%+exp(2*i* pi *(500*t))
% S(2,:) = exp(2*i* pi *(145*t))+exp(2*i* pi *(2*145*t))+exp(2*i* pi *(3*145*t))+exp(2*i* pi *(4*145*t));%+exp(2*i* pi *(5*145*t))
Y0 = (fft(S,C_N,2));
f=[C_f:C_f:fs/2];
C_N=C_N/2;
Y00 = Y0(:,1:C_N);
% figure
% plot(f,abs(Y00(1,:)),'b*-')
% hold on
% plot(f,abs(Y00(2,:)),'ro-')
% xlabel('frequency channels (Hz)')
% ylabel('magnitude')
% axis tight;

lambda=C./f;
f_ref=1000;
d=0.5*C/f_ref;
x_true=zeros(N,C_N);
S11 = f0(1)/C_f;
S12 = (f0(1)+B(1))/C_f;
S21 = f0(2)/C_f;
S22 = (f0(2)+B(2))/C_f;
pos_true= [fix(theta(1)/resolution)+1 fix(theta(2)/resolution)+1]';
x_true(pos_true(1),S11:S12)=1;
x_true(pos_true(2),S21:S22)=1;
% figure
% imagesc(f,grid,x_true)
% xlabel('frequency channels (Hz)')
% ylabel('direction')

%% --------uniform linear array (ULA), with the origin at the middle-----%%
A_true = zeros(C_N*M,K);
for t=1:C_N
    for m = 1:M
       for k = 1:K
       A_true((t-1)*M+m,k) = exp(-2*i * pi *d*(m-(M+1)/2) * cos(theta(k)/180*pi)/lambda(t));
        end
    end
end
%% --------------------observation generation----------------------------%%
Y=zeros(M,C_N);
for t=1:C_N
    Y(:,t) = A_true((t-1)*M+1:t*M,:)*Y00(:,t);
end
%% -----------------------Dictionary construct---------------------------%%
A = zeros(C_N*M,N);
B = zeros(C_N*M,N);
% for t=1:C_N
%     for m = 1:M
%        for n = 1:N
%        A((t-1)*M+m,n) = exp(-2*i * pi *d*(m-(M+1)/2) * cos(grid(n)/180*pi)/lambda(t));
%        B((t-1)*M+m,n) = 2*i * pi *d*(m-(M+1)/2) * sin(grid(n)/180*pi)/lambda(t)* A((t-1)*M+m,n);
%         end
%     end
% end
pos=[-(M-1)/2:(M-1)/2].';
for t=1:C_N
       A((t-1)*M+1:t*M,:) = exp(-2*i * pi *d*(pos) * cos(grid.'/180*pi)/lambda(t));
       A((t-1)*M+1:t*M,:) = A((t-1)*M+1:t*M,:)* diag(1 ./ diag(A((t-1)*M+1:t*M,:)'*A((t-1)*M+1:t*M,:)));
       B((t-1)*M+1:t*M,:) = 2*i * pi *d*(pos) * sin(grid.'/180*pi)/lambda(t).* A((t-1)*M+1:t*M,:);
end
%%  ----------------add noise-------------------------
snr = [0];
N_exp = 5;

Err_DP=zeros(K,length(snr));

Err_VB=zeros(K,length(snr));
for nnn=1:length(snr)
    nnn
    Err_1=zeros(K,1);
    Err_2=zeros(K,1);
    for mmm=1:N_exp
        mmm
Y_n=awgn(Y,snr(nnn),'measured');
Sigma_a = sqrt(norm(Y_n - Y,'fro')^2/M/T);
%% ------------------------conventioal beamforming-----------------------%%
Xbeam_hat=zeros(N,C_N);
for t=S11:S22
    Xbeam_hat(:,t)=A((t-1)*M+1:(t)*M,:)'*Y_n(:,t);
end
figure
imagesc(f,grid,abs(Xbeam_hat))
xlabel('frequency channels (Hz)')
% ylabel('direction')
%% ----------------------------Lasso+averaging---------------------------%%

%% ------------------------------fused Lasso-----------------------------%%
% Xbp_hat_h = zeros(N,C_N);
% C_N_1=length(S11:S22);
% D_oper = zeros( (C_N_1 - 1) , N * C_N_1 );
% for t = 1 : C_N_1-1 
%     D_oper( t , : ) = [ zeros(1 , (t-1)*N) , ones(1,N) ,  -1*ones(1,N) , zeros(1 , N*C_N_1-(2+(t-1))*N)];
% end
% lam1=1;
% lam2=1;
% tic;
% cvx_begin
%         variable  Xbp(N*C_N_1) complex;
%         minimize(lam1*norm(Xbp,1)+lam2*norm(D_oper*Xbp,1));  
%         temp=0;
%         for t=1:C_N_1
%         temp=temp+norm(Y_n(:,t+S11-1)-A((S11-1+t-1)*M+1:(t+S11-1)*M,:)*Xbp((t-1)*N+1:t*N),2);
%         end
%         subject to
%         temp<=350*(Sigma_a);
% cvx_end 
% toc
% Xbp_hat=reshape(Xbp,N,C_N_1);
% Xbp_hat_h (:,S11:S22)=Xbp_hat;
% figure
% imagesc(f,grid,abs(Xbp_hat_h));
% xlabel('frequency channels (Hz)')
% ylabel('direction')
%% --------------------SBL without joint sparsity------------------------%%
% Xbcs_s_hat=zeros(N,C_N);
% for t=S11:S22
% Xbcs_s_hat(:,t)=bcs_vb( A((t-1)*M+1:t*M,:),Y_n(:,t));    
% end
% figure
% imagesc(f,grid,abs(Xbcs_s_hat))
% xlabel('frequency channels (Hz)')
% ylabel('direction')
%% --------------------SBL with joint sparsity---------------------------%% 
Xbcs_g_hat_h=zeros(N,C_N); 
% Xbcs_g_hat=bmtl_vb(A((S11-1)*M+1:S22*M,:),Y_n(:,S11:S22))
[Xbcs_g_hat beta]=bmtl_vb_offgrid(A((S11-1)*M+1:S22*M,:),B((S11-1)*M+1:S22*M,:),Y_n(:,S11:S22),1);
Xbcs_g_hat_h(:,S11:S22)=Xbcs_g_hat;
figure
imagesc(f,grid,abs(Xbcs_g_hat_h))
xlabel('frequency channels (Hz)')
ylabel('direction')
% 
%% ------------------------MTL-CS-DP-(on-grid)---------------------------%% 
% Xdp_hat_h=zeros(N,C_N);
% Xdp_hat=bmtl_DP(A((S11-1)*M+1:S22*M,:),Y_n(:,S11:S22));
% Xdp_hat_h(:,S11:S22)=Xdp_hat;
% figure
% imagesc(f,grid,abs(Xdp_hat_h))
% xlabel('frequency channels (Hz)')
% ylabel('direction')

%% -----------------------MTL-CS-DP-(off-grid)---------------------------%% 
Xdp_offgrid_hat_h=zeros(N,C_N);
[Xdp_offgrid_hat X_Beta]=bmtl_DP_offgrid_2(A((S11-1)*M+1:S22*M,:),B((S11-1)*M+1:S22*M,:),Y_n(:,S11:S22),resolution);
Xdp_offgrid_hat_h(:,S11:S22)=Xdp_offgrid_hat;
% figure
% imagesc(f,grid,abs(Xdp_offgrid_hat_h))
% xlabel('frequency channels (Hz)')
% ylabel('direction')

pos_a=find(sum(abs(X_Beta),2)>0.001);
temp_Beta = X_Beta;
beta=zeros(N,1);
for n=1:length(pos_a)
    temp_pos=find(abs(X_Beta(pos_a(n),:))>0.001);
    beta(pos_a(n))=mean(X_Beta(pos_a(n),temp_pos));
end
% xp_rec = grid + beta;
% figure(100),plot(theta, 10*log10(amp.^2), 'bo', xp_rec, 10*log10(abs(Xdp_offgrid_hat_h(:,S11+16)).^2/max(abs(Xdp_offgrid_hat_h(:,S11+16)).^2)), 'r*-');
% hold on
% plot(xp_rec, 10*log10(abs(Xdp_offgrid_hat_h(:,S22-16)).^2/max(abs(Xdp_offgrid_hat_h(:,S22-16)).^2)), 'rx-');
% axis([0,180,min([10*log10(amp.^2); 10*log10(abs(Xdp_offgrid_hat_h(:,S11+16)).^2/max(abs(Xdp_offgrid_hat_h(:,S11+16)).^2))]),max([10*log10(amp.^2); 10*log10(abs(Xdp_offgrid_hat_h(:,S11+16)).^2/max(abs(Xdp_offgrid_hat_h(:,S11+16)).^2))])+3]);
% xlabel('DOA (degrees)', 'fontsize',12); ylabel('Power (dB)','fontsize',12);
% legend('True DOAs','target-1','target-2');
% axis tight;

est_theta=grid(pos_a)+beta(pos_a);

Err_1(1)  = Err_1(1)+min(abs(est_theta-theta(1)));
Err_1(2)  = Err_1(2)+min(abs(est_theta-theta(2)));
% figure(1001),plot(theta, 10*log10(amp.^2), 'bo', grid, 10*log10(abs(Xdp_hat_h(:,S11+16)).^2/max(abs(Xdp_hat_h(:,S11+16)).^2)), 'r*-');
% hold on
% plot(grid, 10*log10(abs(Xdp_hat_h(:,S22-16)).^2/max(abs(Xdp_hat_h(:,S22-16)).^2)), 'rx-');
% axis([0,180,min([10*log10(amp.^2); 10*log10(abs(Xdp_hat_h(:,S22-16)).^2/max(abs(Xdp_hat_h(:,S22-16)).^2));10*log10(abs(Xdp_hat_h(:,S11+16)).^2/max(abs(Xdp_hat_h(:,S11+16)).^2))]),max([10*log10(amp.^2); 10*log10(abs(Xdp_hat_h(:,S22-16)).^2/max(abs(Xdp_hat_h(:,S22-16)).^2));10*log10(abs(Xdp_hat_h(:,S11+16)).^2/max(abs(Xdp_hat_h(:,S11+16)).^2))])+3]);
% xlabel('DOA (degrees)', 'fontsize',12); ylabel('Power (dB)','fontsize',12);
% legend('True DOAs','target-1','target-2');
% axis tight;
% figure(1002),plot(theta, 10*log10(amp.^2), 'bo', grid, 10*log10(abs(Xbcs_g_hat_h(:,S11+16)).^2/max(abs(Xbcs_g_hat_h(:,S11+16)).^2)), 'r*-');
% hold on
% plot(grid, 10*log10(abs(Xbcs_g_hat_h(:,S22-16)).^2/max(abs(Xbcs_g_hat_h(:,S22-16)).^2)), 'rx-');
% axis([0,180,min([10*log10(amp.^2); 10*log10(abs(Xbcs_g_hat_h(:,S22-16)).^2/max(abs(Xbcs_g_hat_h(:,S22-16)).^2));10*log10(abs(Xdp_hat_h(:,S11+16)).^2/max(abs(Xdp_hat_h(:,S11+16)).^2))]),max([10*log10(amp.^2); 10*log10(abs(Xbcs_g_hat_h(:,S22-16)).^2/max(abs(Xbcs_g_hat_h(:,S22-16)).^2));10*log10(abs(Xbcs_g_hat_h(:,S11+16)).^2/max(abs(Xbcs_g_hat_h(:,S11+16)).^2))])+3]);
% xlabel('DOA (degrees)', 'fontsize',12); ylabel('Power (dB)','fontsize',12);
% legend('True DOAs','target-1','target-2');
% axis tight;
% % %% --------------CSSM-----------------%%%%%%%%%%%%%%%%%%%%%%%%
% %  Xff = [];
% %         for k = S11:S22
% %             Xff = [Xff Y_n(:,k)];
% %         end
% %         R_f = (Xff-mean(Xff,2)*ones(1,(S22-S11+1)*1))*(Xff-mean(Xff,2)*ones(1,(S22-S11+1)*1)'/((S22-S11+1)*1);
% %         [Uaa Vaa] = eigs(R_f,L);
% %         
% %         % Initial estimates for CSSM
% %         
% %         for tta = 1:length(grid)
% %             Pmvdr(tta) = 10*log10(1./abs(diag(sv(:,tta)'*inv(R_f)*sv(:,tta))).^1);
% %         end
% %         
% %         bb = findpeaks1(Pmvdr);
% %         azi = [];
% %         for kb = 1:length(bb)
% %             azi = [azi Tta(bb(kb)) Tta(bb(kb))-0.125*beamwidth Tta(bb(kb))+0.125*beamwidth];
% %         end
%       azi = theta-[0.1 0.1]';
%       Rn  = zeros(M); Ry = zeros(M);
%       fc_n= (S11+S22)/2;
%         for kkk = S11:S22
%             Rx = (Y_n(:,kkk)-mean(Y_n(:,kkk)))*(Y_n(:,kkk)-mean(Y_n(:,kkk)))';
%             % For cssm 
%             ksv=zeros(M,length(azi));
%             ksv0=zeros(M,length(azi));
%                 for m = 1:M
%                     for n = 1:length(azi)
%                        ksv(m,n)=exp(-2*i * pi *d*(m-(M+1)/2) * cos(azi(n)/180*pi)/lambda(kkk));
%                        ksv0(m,n) = exp(-2*i * pi *d*(m-(M+1)/2) * cos(azi(n)/180*pi)/lambda(fc_n));
%                     end
%                 end
% %             ksv = exp(-2*i * pi *d*(m-(M+1)/2) * cos(azi/180*pi)/lambda(kkk));%exp(-i*2*pi*f(kkk)/c*AG'*zia);
% %             ksv0 = exp(-2*i * pi *d*(m-(M+1)/2) * cos(azi/180*pi)/lambda(fc_n));%exp(-i*2*pi*f(fc_n)/c*AG'*zia);
%             prod = ksv*ksv0';
%             [Uf D Vf] = svd(prod);
%             TT = Vf*Uf';
%             Ryf{kkk} = TT*Rx*TT';
%             [aa bb] = eigs(Rx,M);
%             db = diag(bb);
%             sigv(kkk) = mean(db(K+1:end));
%             Ry = Ry + Ryf{kkk}/(98-34+1);%;
%             Rn = Rn + sigv(kkk)*TT*TT';  
%             Fa{kkk} = aa(:,1:K);
%             W{kkk} = aa(:,K+1:end);   
%         end 
%         F1 = Fa{S22-S11+1};
%         % For cssm
%         [Ucssm Vcssm] = eigs(inv(Rn)*Ry,M);
%         Uncssm = Ucssm(:,K+1:end);    
%         
%         grid_test=(0:0.1:180)';
%         AA=zeros(M,length(grid_test));
%            
%           for n = 1:length(grid_test)
%              for m = 1:M
%              AA(m,n) = exp(-2*i * pi *d*(m-(M+1)/2) * cos(grid_test(n)/180*pi)/lambda(fc_n));
%              end
%           end
%    
%         P_cssm=zeros(1,length(grid_test));
%         P_tops=zeros(1,length(grid_test));
%         for tta = 1:length(grid_test);
%             P_cssm(tta) = (1./abs(AA(:,tta)'*Uncssm*Uncssm'*AA(:,tta)));
%             % For TOPS
%             Du = []; Dudash = [];
%             for kk = S11:S22
%             a_tops=zeros(M,1);
%             b_tops=zeros(M,1);
%             for m = 1:M
%             a_tops(m) = exp(-2*i * pi*d*(m-(M+1)/2) * cos(grid_test(tta)/180*pi)/lambda(kk));
%             end   
%             proj = eye(M) - a_tops*a_tops';
%             for m = 1:M
%             b_tops(m) = exp(-2*i *(f(kk)-f(fc_n))/C* pi*d*(m-(M+1)/2) * cos(grid_test(tta)/180*pi));
%             end 
%             Phi = diag(b_tops);
%             Ui = Phi*F1;
%             Du = [Du Ui'*W{kk}];
%             Uidash = proj*Ui;
%             Dudash = [Dudash Uidash'*W{kk}];
%             end
%             [temp1 temp2] = eigs(Dudash*Dudash',K);
%             P_tops(tta) = (1/temp2(K,K)); 
%         end
%         figure('name','P_CSSM')
%         plot(grid_test,10*log10(abs(P_cssm)/max(abs(P_cssm))));
%         figure('name','P_TOPS')
%         plot(grid_test,10*log10(abs(P_tops)/max(abs(P_tops))));
%% ---------------------------------------------------------%%
Xvb_offgrid_hat_h=zeros(N,C_N);
[Xvb_offgrid_hat X_Beta_vb]=bmtl_VB_offgrid_2(A((S11-1)*M+1:S22*M,:),B((S11-1)*M+1:S22*M,:),Y_n(:,S11:S22),resolution);
Xvb_offgrid_hat_h(:,S11:S22)=Xvb_offgrid_hat;
% figure
% imagesc(f,grid,abs(Xvb_offgrid_hat_h))
% xlabel('frequency channels (Hz)')
% ylabel('direction')

pos_a_vb=find(sum(abs(X_Beta_vb),2)>0.001);
beta_vb=zeros(N,1);
for n=1:length(pos_a_vb)
    temp_pos=find(abs(X_Beta_vb(pos_a_vb(n),:))>0.001);
    beta_vb(pos_a_vb(n))=mean(X_Beta_vb(pos_a_vb(n),temp_pos));
end
% xp_rec = grid + beta_vb;
% figure(1000),plot(theta, 10*log10(amp.^2), 'bo', xp_rec, 10*log10(abs(Xvb_offgrid_hat_h(:,S11+16)).^2/max(abs(Xvb_offgrid_hat_h(:,S11+16)).^2)), 'r*-');
% hold on
% plot(xp_rec, 10*log10(abs(Xvb_offgrid_hat_h(:,S22-16)).^2/max(abs(Xvb_offgrid_hat_h(:,S22-16)).^2)), 'rx-');
% axis([0,180,min([10*log10(amp.^2); 10*log10(abs(Xvb_offgrid_hat_h(:,S11+16)).^2/max(abs(Xvb_offgrid_hat_h(:,S11+16)).^2))]),max([10*log10(amp.^2); 10*log10(abs(Xvb_offgrid_hat_h(:,S11+16)).^2/max(abs(Xvb_offgrid_hat_h(:,S11+16)).^2))])+3]);
% xlabel('DOA (degrees)', 'fontsize',12); ylabel('Power (dB)','fontsize',12);
% legend('True DOAs','target-1','target-2');
% axis tight;
est_theta_vb=grid(pos_a_vb)+beta_vb(pos_a_vb);
Err_2(1)  = Err_2(1)+min(abs(est_theta_vb-theta(1)));
Err_2(2)  = Err_2(2)+min(abs(est_theta_vb-theta(2)));
    end
    Err_DP(:,nnn) = Err_1/N_exp;
    Err_VB(:,nnn) = Err_2/N_exp;
end
figure
plot(snr,Err_DP(1,:),'*-');hold on;plot(snr,Err_DP(2,:),'o-');hold on;plot(snr,Err_VB(1,:),'--');hold on;
plot(snr,Err_VB(2,:),'-.')