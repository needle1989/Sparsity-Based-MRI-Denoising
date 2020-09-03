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
resolution = 2;
grid = (0:resolution:180)';
N = length(grid);
i=sqrt(-1);
%% ----------------------wideband signal generation--------------------- %%
N_measurements=[2 4 6 10 15 20 30];
T=256;%signal length
C_N=128;%number of channels
fs=4096;
C_f=fs/T;
f=[C_f:C_f:fs/2];
lambda=C./f;
f_ref=500;
d=0.5*C/f_ref;
% d_err = 0.2*randn(M,1);
% d_err(1) = 0;
d_err = zeros(M,1);

%% -----------------------Dictionary construct---------------------------%%
A = zeros(C_N*M,N);
B = zeros(C_N*M,N);
for tn=1:C_N
    for m = 1:M
       for n = 1:N
       A((tn-1)*M+m,n) = exp(-2*i * pi *(d*(m-(M+1)/2)) * cos(grid(n)/180*pi)/lambda(tn));
       B((tn-1)*M+m,n) = 2*i * pi *(d*(m-(M+1)/2)) * sin(grid(n)/180*pi)/lambda(tn)* A((tn-1)*M+m,n);
        end
    end
end
%%  
Y=[];Y_n=[];
for nnnn=1:length(N_measurements)
t=(nnnn-1)*T+(0:T-1)/fs;
N_t=length(t);
Order=10;
f0=[151 197];
amp = [1 1]';
S=zeros(K,T);
for k=1:K
    for n=1:Order
       S(k,:) = S(k,:)+amp(k)*exp(2*i* pi *(n*f0(k)*t));
    end
end
Y0 = (fft(S,T,2));
Y00 = Y0(:,1:C_N);
% figure
% plot(f,abs(Y00(1,:)),'b*-')
% hold on
% plot(f,abs(Y00(2,:)),'ro-')
% hold on
% plot(f,abs(Y00(3,:)),'k.-')
% hold on
% plot(f,abs(Y00(4,:)),'y--')
% xlabel('frequency channels (Hz)')
% ylabel('magnitude')
% axis tight;
pos_channel_1=find(abs(Y00(1,:))>80);
pos_channel_2=find(abs(Y00(2,:))>80);
% pos_channel_3=find(abs(Y00(3,:))>80);
% pos_channel_4=find(abs(Y00(4,:))>80);
pos_channel=unique([pos_channel_1 pos_channel_2]);%pos_channel_3 pos_channel_4
theta = [58.77 58.77+N_measurements(nnnn)]';
%% --------uniform linear array (ULA), with the origin at the middle-----%%
A_true = zeros(C_N*M,K);
for tn=1:C_N
    for m = 1:M
       for k = 1:K
       A_true((tn-1)*M+m,k) = exp(-2*i * pi *(d*(m-(M+1)/2)+d_err(m))* cos(theta(k)/180*pi)/lambda(tn));
        end
    end
end
x_true=zeros(N,C_N);
pos_true= [fix(theta(1)/resolution)+1 fix(theta(2)/resolution)+1]'; %fix(theta(3)/resolution)+1 fix(theta(4)/resolution)+1
x_true(pos_true(1),pos_channel_1)=1;
x_true(pos_true(2),pos_channel_2)=1;
% x_true(pos_true(3),pos_channel_3)=1;
% x_true(pos_true(4),pos_channel_4)=1;
% figure
% imagesc(f,grid,x_true)
% xlabel('frequency channels (Hz)')
% ylabel('direction')

%% --------------------observation generation----------------------------%%
for tn = 1:C_N
    Y{nnnn}(:,tn) = A_true((tn-1)*M+1:tn*M,:)*Y00(:,tn);
end
end
%%----------------add noise-------------------------
%snr = [-10 -15 -18];
snr = -44.7;

Err_DP_I=zeros(K,length(N_measurements));
Err_DP_S=zeros(K,length(N_measurements));
Err_RVM_S=zeros(K,length(N_measurements));
Err_CSSM=zeros(K,length(N_measurements));
%N_exp = 10;
N_exp = 10;
Sigma_a=sqrt(10^(-snr/10));
for nnn=1:length(N_measurements)
    nnn
    theta = [58.77 58.77+N_measurements(nnn)]';
    Y_n=[];
    for mmm=1:N_exp
        mmm
         Noise=Sigma_a*randn(M,C_N)/sqrt(2)+sqrt(-1)*Sigma_a*randn(M,C_N)/sqrt(2);
         Y_n = Y{nnn}+ Noise; 
%% ------------------------conventioal beamforming-----------------------%%
Xbeam_hat=zeros(N,C_N);
for tn=1:C_N
    Xbeam_hat(:,tn)=Xbeam_hat(:,tn)+abs(A((tn-1)*M+1:(tn)*M,:)'*Y_n(:,tn));
end
figure
imagesc(f,grid,abs(Xbeam_hat))
xlabel('frequency channels (Hz)')
ylabel('direction')
temp=sort(max(abs(Xbeam_hat(:,10:end))));
temp_a=temp(89);
pos_channel = find(max(abs(Xbeam_hat(:,10:end)))>temp_a)+9;
Xbeam_hat_channel=zeros(size(Xbeam_hat));
Xbeam_hat_channel(:,pos_channel)=Xbeam_hat(:,pos_channel);
% figure
% imagesc(f,grid,abs(Xbeam_hat_channel))
%%
A_temp=zeros(M*length(pos_channel),N);
% Y_temp=zeros(M,length(pos_channel));
Y_temp=zeros(M,length(pos_channel));
for nn=1:length(pos_channel)
    A_temp((nn-1)*M+1:nn*M,:)=A((pos_channel(nn)-1)*M+1:pos_channel(nn)*M,:);
    B_temp((nn-1)*M+1:nn*M,:)=B((pos_channel(nn)-1)*M+1:pos_channel(nn)*M,:);
    Y_temp(:,nn)=Y_n(:,pos_channel(nn));
end
%% ------------------------MTL-CS-DP-(on-grid)---------------------------%%
Xdp_hat_h=zeros(N,C_N);
[Xdp_hat alpha_m alpha0_DP]=bmtl_DP(A_temp,Y_temp);
Xdp_hat_h(:,pos_channel)=Xdp_hat;
alpha_m=alpha_m/max(max(alpha_m));

figure
imagesc(f,grid,abs(Xdp_hat_h))
xlabel('frequency channels (Hz)')
ylabel('direction')
[ind alpha_DP]=kmeans(alpha_m.',2,'EmptyAction','drop');
theta_ex = [58 60];
est_theta_DP_S_1=DPSeperatedSearch(Xdp_hat,alpha_DP,theta_ex,alpha0_DP,pos_channel,Y_temp,A_temp,theta(1));
theta_ex = theta_ex+N_measurements(nnn);
est_theta_DP_S_2=DPSeperatedSearch(Xdp_hat,alpha_DP,theta_ex,alpha0_DP,pos_channel,Y_temp,A_temp,theta(2));
Err_DP_S(1,nnn) = Err_DP_S(1,nnn) + (abs(est_theta_DP_S_1-theta(1)));
Err_DP_S(2,nnn) = Err_DP_S(2,nnn) + (abs(est_theta_DP_S_2-theta(2)));
% 
% %% ------------------------RVM zhangmeng Liu---------------------
% X_RVM_hat_h=zeros(N,C_N);
X_RVM_hat_h=zeros(N,C_N);
[X_RVM_hat alpha alpha0] = mt_CS(A_temp,Y_temp);
X_RVM_hat_h(:,pos_channel)=X_RVM_hat;
figure
imagesc(f,grid,abs(X_RVM_hat_h))
xlabel('frequency channels (Hz)')
ylabel('direction')
%% testing
theta_ex = [56 58 60 62];
est_theta_RVM_1=RVMSeperatedSearch(theta_ex,alpha,alpha0,pos_channel,Y_temp,A_temp,theta(1));
theta_ex=theta_ex+N_measurements(nnn);
est_theta_RVM_2=RVMSeperatedSearch(theta_ex,alpha,alpha0,pos_channel,Y_temp,A_temp,theta(2));
Err_RVM_S(1,nnn) = Err_RVM_S(1,nnn) + (abs(est_theta_RVM_1-theta(1)));
Err_RVM_S(2,nnn) = Err_RVM_S(2,nnn) + (abs(est_theta_RVM_2-theta(2)));

%% -----------------------MTL-CS-DP-(off-grid)---------------------------%% 
Xdp_offgrid_hat_h=zeros(N,C_N);
[Xdp_offgrid_hat X_Beta]=bmtl_DP_offgrid_2(A_temp,B_temp,Y_temp,resolution);
Xdp_offgrid_hat_h(:,pos_channel)=Xdp_offgrid_hat;
figure
imagesc(f,grid,abs(Xdp_offgrid_hat_h))
xlabel('frequency channels (Hz)')
ylabel('direction')


pos_a=find(sum(abs(X_Beta),2)>0.001);
temp_Beta = X_Beta;
beta=zeros(N,1);
for n=1:length(pos_a)
    temp_pos=find(abs(X_Beta(pos_a(n),:))>0.001);
    beta(pos_a(n))=mean(X_Beta(pos_a(n),temp_pos));
end
% xp_rec = grid + beta;
% Xdp_offgrid_hat_ref=max(abs(Xdp_offgrid_hat_h(:,pos_channel_1(2))).^2,abs(Xdp_offgrid_hat_h(:,pos_channel_2(1))).^2);
% figure(100),plot(theta, 10*log10(amp.^2), 'bo', xp_rec, 10*log10(abs(Xdp_offgrid_hat_h(:,pos_channel_1(1))).^2/Xdp_offgrid_hat_ref), 'r*-');
% hold on
% plot(xp_rec, 10*log10(abs(Xdp_offgrid_hat_h(:,pos_channel_2(1))).^2/Xdp_offgrid_hat_ref), 'rx-');
% % axis([0,180,min([10*log10(amp.^2); 10*log10(abs(Xdp_offgrid_hat_h(:,pos_channel_2(1))).^2/Xdp_offgrid_hat_ref)]),...
% %     max([10*log10(amp.^2); 10*log10(abs(Xdp_offgrid_hat_h(:,pos_channel_2(1))).^2/Xdp_offgrid_hat_ref)])+3]);
% xlabel('DOA (degrees)', 'fontsize',12); ylabel('Power (dB)','fontsize',12);
% legend('True DOAs','target-1','target-2');
% axis tight;

est_theta=grid(pos_a)+beta(pos_a);

if min(abs(est_theta-theta(1)))>2
   Err_DP_I(1,nnn) = Err_DP_I(1,nnn)+2;
else
  Err_DP_I(1,nnn) = Err_DP_I(1,nnn)+min(abs(est_theta-theta(1)));  
end
if min(abs(est_theta-theta(2)))>2
   Err_DP_I(2,nnn) = Err_DP_I(2,nnn)+2;
else
  Err_DP_I(2,nnn) = Err_DP_I(2,nnn)+min(abs(est_theta-theta(2)));  
end
%% ------------------------CSSM-------------------------------------------
% XCSSM_hat=CSSM_multimeasure(A_temp,Y_temp);
% 
% Err_CSSM(1,nnn) = Err_CSSM(1,nnn)+min(abs(XCSSM_hat-theta(1)));
% Err_CSSM(2,nnn) = Err_CSSM(2,nnn)+min(abs(XCSSM_hat-theta(2)));

    end
Err_DP_I(1,nnn) = Err_DP_I(1,nnn)/N_exp;
Err_DP_I(2,nnn) = Err_DP_I(2,nnn)/N_exp;
Err_DP_S(1,nnn) = Err_DP_S(1,nnn)/N_exp;
Err_DP_S(2,nnn) = Err_DP_S(2,nnn)/N_exp;
Err_RVM_S(1,nnn) = Err_RVM_S(1,nnn)/N_exp;
Err_RVM_S(2,nnn) = Err_RVM_S(2,nnn)/N_exp;
% save Err_DP_I.mat Err_DP_I;
% save Err_DP_S.mat Err_DP_S;
% save Err_RVM_S.mat Err_RVM_S;
% Err_CSSM(1,nnn) = Err_CSSM(1,nnn)/N_exp;
% Err_CSSM(2,nnn) = Err_CSSM(2,nnn)/N_exp;
end
%save Err_Seperation_25dB.mat Err_DP_I Err_DP_S Err_RVM_S N_measurements;
figure
plot(N_measurements,Err_DP_I(1,:),'r-^',N_measurements,Err_DP_S(1,:),'r-.*',N_measurements,Err_RVM_S(1,:),'b--o');
hold on
plot(N_measurements,Err_DP_I(2,:),'r-^',N_measurements,Err_DP_S(2,:),'r-.*',N_measurements,Err_RVM_S(2,:),'b-o');
figure
plot(N_measurements,sum(Err_RVM_S)/K,'b--*',N_measurements,sum(Err_DP_S)/K,'r--o',N_measurements,sum(Err_DP_I)/K,'r-.^');
axis([4,30,0,1.8]);
xlabel('Source seperation(Deg.)');ylabel('RMSE(Deg.)')
% % figure(1001),plot(theta, 10*log10(amp.^2), 'bo', grid, 10*log10(abs(Xdp_hat_h(:,S11+16)).^2/max(abs(Xdp_hat_h(:,S11+16)).^2)), 'r*-');
% % hold on
% % plot(grid, 10*log10(abs(Xdp_hat_h(:,S22-16)).^2/max(abs(Xdp_hat_h(:,S22-16)).^2)), 'rx-');
% % axis([0,180,min([10*log10(amp.^2); 10*log10(abs(Xdp_hat_h(:,S22-16)).^2/max(abs(Xdp_hat_h(:,S22-16)).^2));10*log10(abs(Xdp_hat_h(:,S11+16)).^2/max(abs(Xdp_hat_h(:,S11+16)).^2))]),max([10*log10(amp.^2); 10*log10(abs(Xdp_hat_h(:,S22-16)).^2/max(abs(Xdp_hat_h(:,S22-16)).^2));10*log10(abs(Xdp_hat_h(:,S11+16)).^2/max(abs(Xdp_hat_h(:,S11+16)).^2))])+3]);
% % xlabel('DOA (degrees)', 'fontsize',12); ylabel('Power (dB)','fontsize',12);
% % legend('True DOAs','target-1','target-2');
% % axis tight;
% % figure(1002),plot(theta, 10*log10(amp.^2), 'bo', grid, 10*log10(abs(Xbcs_g_hat_h(:,S11+16)).^2/max(abs(Xbcs_g_hat_h(:,S11+16)).^2)), 'r*-');
% % hold on
% % plot(grid, 10*log10(abs(Xbcs_g_hat_h(:,S22-16)).^2/max(abs(Xbcs_g_hat_h(:,S22-16)).^2)), 'rx-');
% % axis([0,180,min([10*log10(amp.^2); 10*log10(abs(Xbcs_g_hat_h(:,S22-16)).^2/max(abs(Xbcs_g_hat_h(:,S22-16)).^2));10*log10(abs(Xdp_hat_h(:,S11+16)).^2/max(abs(Xdp_hat_h(:,S11+16)).^2))]),max([10*log10(amp.^2); 10*log10(abs(Xbcs_g_hat_h(:,S22-16)).^2/max(abs(Xbcs_g_hat_h(:,S22-16)).^2));10*log10(abs(Xbcs_g_hat_h(:,S11+16)).^2/max(abs(Xbcs_g_hat_h(:,S11+16)).^2))])+3]);
% % xlabel('DOA (degrees)', 'fontsize',12); ylabel('Power (dB)','fontsize',12);
% % legend('True DOAs','target-1','target-2');
% % axis tight;
% %% ---------------------------------------------------------%%
% Xvb_offgrid_hat_h=zeros(N,C_N);
% [Xvb_offgrid_hat X_Beta_vb]=bmtl_VB_offgrid_2(A((S11-1)*M+1:S22*M,:),B((S11-1)*M+1:S22*M,:),Y_n(:,S11:S22),resolution);
% Xvb_offgrid_hat_h(:,S11:S22)=Xvb_offgrid_hat;
% % figure
% % imagesc(f,grid,abs(Xvb_offgrid_hat_h))
% % xlabel('frequency channels (Hz)')
% % ylabel('direction')
% 
% pos_a_vb=find(sum(abs(X_Beta_vb),2)>0.001);
% beta_vb=zeros(N,1);
% for n=1:length(pos_a_vb)
%     temp_pos=find(abs(X_Beta_vb(pos_a_vb(n),:))>0.001);
%     beta_vb(pos_a_vb(n))=mean(X_Beta_vb(pos_a_vb(n),temp_pos));
% end
% % xp_rec = grid + beta_vb;
% % figure(1000),plot(theta, 10*log10(amp.^2), 'bo', xp_rec, 10*log10(abs(Xvb_offgrid_hat_h(:,S11+16)).^2/max(abs(Xvb_offgrid_hat_h(:,S11+16)).^2)), 'r*-');
% % hold on
% % plot(xp_rec, 10*log10(abs(Xvb_offgrid_hat_h(:,S22-16)).^2/max(abs(Xvb_offgrid_hat_h(:,S22-16)).^2)), 'rx-');
% % axis([0,180,min([10*log10(amp.^2); 10*log10(abs(Xvb_offgrid_hat_h(:,S11+16)).^2/max(abs(Xvb_offgrid_hat_h(:,S11+16)).^2))]),max([10*log10(amp.^2); 10*log10(abs(Xvb_offgrid_hat_h(:,S11+16)).^2/max(abs(Xvb_offgrid_hat_h(:,S11+16)).^2))])+3]);
% % xlabel('DOA (degrees)', 'fontsize',12); ylabel('Power (dB)','fontsize',12);
% % legend('True DOAs','target-1','target-2');
% % axis tight;
% est_theta_vb=grid(pos_a_vb)+beta_vb(pos_a_vb);
% Err_2(1)  = Err_2(1)+min(abs(est_theta_vb-theta(1)));
% Err_2(2)  = Err_2(2)+min(abs(est_theta_vb-theta(2)));