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
theta = [58.77 108.42]';
%theta = [68.78 108.42]';
resolution = 2;
grid = (0:resolution:180)';
N = length(grid);
i=sqrt(-1);
%% ----------------------wideband signal generation--------------------- %%
T=256;%signal length
C_N=256;%number of channels
fs=4096;
C_f=fs/C_N;
t=(0:T-1)/fs;
N_t=length(t);
Order=10;
f0=[151 197];%180
amp = [1 1]';
S=zeros(K,T);
for k=1:K
    for n=1:Order
       S(k,:) = S(k,:)+amp(k)*exp(2*i* pi *(n*f0(k)*t));
    end
end
Y0 = (fft(S,C_N,2));
f=[C_f:C_f:fs/2];
C_N=C_N/2;
Y00 = Y0(:,1:C_N);
figure
plot(f,abs(Y00(1,:)),'b*-')
hold on
plot(f,abs(Y00(2,:)),'ro-')
% hold on
% plot(f,abs(Y00(3,:)),'k.-')
% hold on
% plot(f,abs(Y00(4,:)),'y--')
xlabel('frequency channels (Hz)')
ylabel('magnitude')
axis tight;
pos_channel_1=find(abs(Y00(1,:))>80);
pos_channel_2=find(abs(Y00(2,:))>80);
% pos_channel_3=find(abs(Y00(3,:))>80);
% pos_channel_4=find(abs(Y00(4,:))>80);
pos_channel=unique([pos_channel_1 pos_channel_2]);%pos_channel_3 pos_channel_4

lambda=C./f;
f_ref=500;
d=0.5*C/f_ref;
x_true=zeros(N,C_N);
pos_true= [fix(theta(1)/resolution)+1 fix(theta(2)/resolution)+1]'; %fix(theta(3)/resolution)+1 fix(theta(4)/resolution)+1
x_true(pos_true(1),pos_channel_1)=abs(Y00(1,pos_channel_1));
x_true(pos_true(2),pos_channel_2)=abs(Y00(2,pos_channel_2));
% x_true(pos_true(3),pos_channel_3)=1;
% x_true(pos_true(4),pos_channel_4)=1;
figure
imagesc(f,grid,x_true)
xlabel('frequency channels (Hz)')
ylabel('direction')
% d_err = 0.2*randn(M,1);
% d_err(1) = 0;
d_err = zeros(M,1);
%% --------uniform linear array (ULA), with the origin at the middle-----%%
A_true = zeros(C_N*M,K);
for tn=1:C_N
    for m = 1:M
       for k = 1:K
       A_true((tn-1)*M+m,k) = exp(-2*i * pi *(d*(m-(M+1)/2)+d_err(m))* cos(theta(k)/180*pi)/lambda(tn));
        end
    end
end
%% --------------------observation generation----------------------------%%
Y=zeros(M,C_N);
for tn = 1:C_N
    Y(:,tn) = A_true((tn-1)*M+1:tn*M,:)*Y00(:,tn);
end
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
%%  ----------------add noise-------------------------
%snr = [-10 -15 -18];
snr = -35;
Err_1=zeros(length(snr),1);
Err_2=zeros(length(snr),1);
%N_exp = 10;
N_exp = 1;
for nnn=1:length(snr)
    nnn
    for mmm=1:N_exp
Sigma_a=sqrt(10^(-snr(nnn)/10));
Noise=Sigma_a*randn(size(Y))/sqrt(2)+sqrt(-1)*Sigma_a*randn(size(Y))/sqrt(2);
Y_n = Y + Noise;
SNR=0;
for nnnn=1:10 
    SNR=SNR+10*log10(norm(Y(nnnn,:))^2/norm(Noise(nnnn,:))^2);
end
SNR=SNR/20
%% ------------------------conventioal beamforming-----------------------%%
Xbeam_hat=zeros(N,C_N);
for tn=1:C_N
    Xbeam_hat(:,tn)=A((tn-1)*M+1:(tn)*M,:)'*Y_n(:,tn);
end
figure
imagesc(f,grid,abs(Xbeam_hat))
xlabel('frequency channels (Hz)')
ylabel('direction')
temp=sort(mean(abs(Xbeam_hat(:,10:end))));
temp_a=temp(88);
pos_channel = find(mean(abs(Xbeam_hat(:,10:end)))>temp_a)+9;
Xbeam_hat_channel=zeros(size(Xbeam_hat));
Xbeam_hat_channel(:,pos_channel)=Xbeam_hat(:,pos_channel);
figure
imagesc(f,grid,abs(Xbeam_hat_channel))
xlabel('frequency channels (Hz)')
ylabel('direction')
%%
A_temp=zeros(M*length(pos_channel),N);
Y_temp=zeros(M,length(pos_channel));
for nn=1:length(pos_channel)
    A_temp((nn-1)*M+1:nn*M,:)=A((pos_channel(nn)-1)*M+1:pos_channel(nn)*M,:);
    B_temp((nn-1)*M+1:nn*M,:)=B((pos_channel(nn)-1)*M+1:pos_channel(nn)*M,:);
    Y_temp(:,nn)=Y_n(:,pos_channel(nn));
end
%% ----------------------------BP+averaging---------------------------%%
% Xbp_hat_h = zeros(N,C_N);
% lam1=0.8;
% Xbp_hat=zeros(N,length(pos_channel));
% for nn=1:length(pos_channel)
% cvx_begin
%         variable  Xbp(N) complex;
%         minimize(norm(Xbp,1)+lam1*norm(Y_temp(:,nn)-A_temp((nn-1)*M+1:(nn)*M,:)*Xbp));  
% cvx_end 
% Xbp_hat(:,nn)=Xbp;
% end
% Xbp_hat_h (:,pos_channel)=Xbp_hat;
% figure
% imagesc(f,grid,abs(Xbp_hat_h));
% xlabel('frequency channels (Hz)')
% ylabel('direction')
% %% ----------------------------joint BPDN---------------------------%%
% Xbp_hat_h_joint = zeros(N,C_N);
% C_N_1=length(pos_channel);
% lam2=0.1;
% 
% cvx_begin
%      variable xgbp(N*C_N_1) complex;
%      expression temp(N); 
%      expression x_temp(C_N_1);  
%      for t=1:N     
%      for nn=1:C_N_1
%          x_temp(nn)=xgbp(t+(nn-1)*N);
%      end
%          temp(t)=norm(x_temp,2);  
%      end
%        temp_e=0;
%         for t=1:C_N_1
%         temp_e=temp_e+norm(Y_temp(:,t)-A_temp((t-1)*M+1:(t)*M,:)*xgbp((t-1)*N+1:t*N),2);
%         end
%      minimize(sum(temp,1)+lam2*temp_e);    
% cvx_end
% 
% Xbp_hat_j=reshape(xgbp,N,C_N_1);
% Xbp_hat_h_joint (:,pos_channel)=Xbp_hat_j;
% figure
% imagesc(f,grid,abs(Xbp_hat_h_joint));
% xlabel('frequency channels (Hz)')
% ylabel('direction')
%% ------------------------------fused Lasso-----------------------------%%
% Xbp_hat_h = zeros(N,C_N);
% C_N_1=length(pos_channel);
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
%         for t=1:length(pos_channel)
%         temp=temp+norm(Y_n(:,pos_channel(t))-A((pos_channel(t)-1)*M+1:(pos_channel(t))*M,:)*Xbp((t-1)*N+1:t*N),2);
%         end
%         subject to
%         temp<=350*(Sigma_a);
% cvx_end 
% toc
% Xbp_hat=reshape(Xbp,N,C_N_1);
% Xbp_hat_h (:,pos_channel)=Xbp_hat;
% figure
% imagesc(f,grid,abs(Xbp_hat_h));
% xlabel('frequency channels (Hz)')
% ylabel('direction')
%% --------------------SBL without joint sparsity------------------------%%
% Xbcs_s_hat=zeros(N,C_N);
% for t=1:length(pos_channel)
% Xbcs_s_hat(:,pos_channel(t))=bcs_vb( A((pos_channel(t)-1)*M+1:pos_channel(t)*M,:),Y_n(:,pos_channel(t)));    
% end
% figure
% imagesc(f,grid,abs(Xbcs_s_hat))
% xlabel('frequency channels (Hz)')
% ylabel('direction')
%% --------------------SBL with joint sparsity---------------------------%% 
% Xbcs_g_hat_h=zeros(N,C_N); 
% Xbcs_g_hat=bmtl_vb(A_temp,Y_temp);
% Xbcs_g_hat_h(:,pos_channel)=Xbcs_g_hat;
% figure
% imagesc(f,grid,abs(Xbcs_g_hat_h))
% xlabel('frequency channels (Hz)')
% ylabel('direction')
% 
%% ------------------------MTL-CS-DP-(on-grid)---------------------------%%
Xdp_hat_h=zeros(N,C_N);
[Xdp_hat alpha_m alpha0_DP]=bmtl_DP(A_temp,Y_temp);
Xdp_hat_h(:,pos_channel)=Xdp_hat;
alpha_m=alpha_m/max(max(alpha_m));
pos_alpha_m=find(sum(alpha_m)>0.5);
figure
imagesc(f,grid,abs(Xdp_hat_h))
xlabel('frequency channels (Hz)')
ylabel('direction')
[ind alpha_DP]=kmeans(alpha_m(:,pos_alpha_m).',2,'EmptyAction','drop');


theta_ex=[58 60];
%theta_ex=[108 110];
pos_theta_ex = fix(theta_ex/resolution)+1;

pos_DP=find(alpha_DP(:,pos_theta_ex)>0.01);
alpha_DP_1 = alpha_DP(pos_DP(1),:);
if pos_DP(1)==2;
alpha_DP_2 = alpha_DP(1,:);
else
alpha_DP_2 = alpha_DP(2,:);  
end
theta_test = [min(theta_ex):0.01:max(theta_ex)];
pos_test=[1:min(pos_theta_ex)-1 max(pos_theta_ex)+1:N];
temp_func = zeros(length(theta_test),1);
temp_func_inv = zeros(length(theta_test),1);
[temp bbb_temp]=max(Xdp_hat);
%f_test_temp_1=  pos_channel_1;%[10 20 29 39 48 58 67 77 86 95]
%f_test_temp_1=[13 26 38 51 63 76 88 100 112 125];
%f_test_temp_1=[13 26 38 51 63 76 88 100 112 125];
% f_test_temp_1 = pos_channel;
% pos_test_temp = zeros(length(f_test_temp_1),1);
% for nn=1:length(f_test_temp_1)
%     if isempty(find(pos_channel==f_test_temp_1(nn)))
%        pos_test_temp(nn) = 0;
%     else
%     pos_test_temp(nn)=find(pos_channel==f_test_temp_1(nn));
%     end
% end
% pos_test_temp_temp=find(pos_test_temp);
% pos_test_temp = pos_test_temp(pos_test_temp_temp);
pos_test_temp_1 = find(bbb_temp==pos_theta_ex(1));
pos_test_temp_2 = find(bbb_temp==pos_theta_ex(2));
if length(pos_test_temp_1)>=length(pos_test_temp_2)
pos_test_temp = pos_test_temp_1;
else
pos_test_temp = pos_test_temp_2; 
end
f_temp = f(pos_channel);
a_steer=zeros(M,1);
d_steer=zeros(M,1);
d_M=d*([1:M]-(M+1)/2)';
for kk = 1:length(theta_test);
 for nn=1:length(pos_test_temp)
%  nn=11;
 f_t = f_temp (pos_test_temp(nn));
 [aa_temp bb_pos]=max(Xdp_hat(:,pos_test_temp(nn)));
 if abs(bb_pos-pos_theta_ex(1))<=1;
     alpha_DP_temp=alpha_DP_1;
 else
     alpha_DP_temp=alpha_DP_2; 
 end
 lambda_temp = C./f_t;
 a_steer = exp(-2*i * pi *(d_M)* cos(theta_test(kk)/180*pi)/lambda_temp);
 d_steer = 2*i * pi *(d_M) * sin(theta_test(kk)/180*pi)/lambda_temp.*a_steer;
 A_temp_test = A_temp((pos_test_temp(nn)-1)*M+1:pos_test_temp(nn)*M,pos_test);
 R_y_f = Y_temp(:,pos_test_temp(nn))*Y_temp(:,pos_test_temp(nn))';
 Sigma_test = A_temp_test*diag(alpha_DP_temp(pos_test))* A_temp_test' + Sigma_a^2*eye(M);%alpha0_DP/T*sqrt(2)
 temp_func(kk) = temp_func(kk) + ...
     (a_steer'*inv(Sigma_test)*(a_steer*a_steer'*inv(Sigma_test)*R_y_f-R_y_f*inv(Sigma_test)*a_steer*a_steer')*inv(Sigma_test)*d_steer)/(a_steer'*inv(Sigma_test)*R_y_f*inv(Sigma_test)*a_steer);
 %temp_func(kk) = (a_steer'*inv(Sigma_test)*(a_steer*a_steer'*inv(Sigma_test)*R_y_f-R_y_f*inv(Sigma_test)*a_steer*a_steer')*inv(Sigma_test)*d_steer)/(a_steer'*inv(Sigma_test)*R_y_f*inv(Sigma_test)*a_steer);
 end
 temp_func_inv(kk)=1/(real(temp_func(kk)));
end
% figure
% plot(theta_test,real(temp_func));
figure
plot(theta_test,10*log(abs(temp_func_inv)/max((abs(temp_func_inv)))));
est_theta_pos_1=find(abs(temp_func_inv)==max(abs(temp_func_inv)));
est_theta_DP = theta_test(est_theta_pos_1)

theta_ex=[108 110];
pos_theta_ex = fix(theta_ex/resolution)+1;

pos_DP=find(alpha_DP(:,pos_theta_ex)>0.01);
alpha_DP_1 = alpha_DP(pos_DP(1),:);
if pos_DP(1)==2;
alpha_DP_2 = alpha_DP(1,:);
else
alpha_DP_2 = alpha_DP(2,:);  
end
theta_test = [min(theta_ex):0.01:max(theta_ex)];
pos_test=[1:min(pos_theta_ex)-1 max(pos_theta_ex)+1:N];
temp_func = zeros(length(theta_test),1);
temp_func_inv = zeros(length(theta_test),1);
[temp bbb_temp]=max(Xdp_hat);
%f_test_temp_1=  pos_channel_1;%[10 20 29 39 48 58 67 77 86 95]
%f_test_temp_1=[13 26 38 51 63 76 88 100 112 125];
%f_test_temp_1=[13 26 38 51 63 76 88 100 112 125];
% f_test_temp_1 = pos_channel;
% pos_test_temp = zeros(length(f_test_temp_1),1);
% for nn=1:length(f_test_temp_1)
%     if isempty(find(pos_channel==f_test_temp_1(nn)))
%        pos_test_temp(nn) = 0;
%     else
%     pos_test_temp(nn)=find(pos_channel==f_test_temp_1(nn));
%     end
% end
% pos_test_temp_temp=find(pos_test_temp);
% pos_test_temp = pos_test_temp(pos_test_temp_temp);
pos_test_temp_1 = find(bbb_temp==pos_theta_ex(1));
pos_test_temp_2 = find(bbb_temp==pos_theta_ex(2));
if length(pos_test_temp_1)>=length(pos_test_temp_2)
pos_test_temp = pos_test_temp_1;
else
pos_test_temp = pos_test_temp_2; 
end
f_temp = f(pos_channel);
a_steer=zeros(M,1);
d_steer=zeros(M,1);
d_M=d*([1:M]-(M+1)/2)';
for kk = 1:length(theta_test);
 for nn=1:length(pos_test_temp)
%  nn=11;
 f_t = f_temp (pos_test_temp(nn));
 [aa_temp bb_pos]=max(Xdp_hat(:,pos_test_temp(nn)));
 if abs(bb_pos-pos_theta_ex(1))<=1;
     alpha_DP_temp=alpha_DP_1;
 else
     alpha_DP_temp=alpha_DP_2; 
 end
 lambda_temp = C./f_t;
 a_steer = exp(-2*i * pi *(d_M)* cos(theta_test(kk)/180*pi)/lambda_temp);
 d_steer = 2*i * pi *(d_M) * sin(theta_test(kk)/180*pi)/lambda_temp.*a_steer;
 A_temp_test = A_temp((pos_test_temp(nn)-1)*M+1:pos_test_temp(nn)*M,pos_test);
 R_y_f = Y_temp(:,pos_test_temp(nn))*Y_temp(:,pos_test_temp(nn))';
 Sigma_test = A_temp_test*diag(alpha_DP_temp(pos_test))* A_temp_test' + Sigma_a^2*eye(M);%alpha0_DP/T*sqrt(2)
 temp_func(kk) = temp_func(kk) + ...
     (a_steer'*inv(Sigma_test)*(a_steer*a_steer'*inv(Sigma_test)*R_y_f-R_y_f*inv(Sigma_test)*a_steer*a_steer')*inv(Sigma_test)*d_steer)/(a_steer'*inv(Sigma_test)*R_y_f*inv(Sigma_test)*a_steer);
 %temp_func(kk) = (a_steer'*inv(Sigma_test)*(a_steer*a_steer'*inv(Sigma_test)*R_y_f-R_y_f*inv(Sigma_test)*a_steer*a_steer')*inv(Sigma_test)*d_steer)/(a_steer'*inv(Sigma_test)*R_y_f*inv(Sigma_test)*a_steer);
 end
 temp_func_inv(kk)=1/(real(temp_func(kk)));
end
% figure
% plot(theta_test,real(temp_func));
figure
plot(theta_test,10*log(abs(temp_func_inv)/max((abs(temp_func_inv)))));
est_theta_pos_1=find(abs(temp_func_inv)==max(abs(temp_func_inv)));
est_theta_DP = theta_test(est_theta_pos_1)


% xp_rec = grid;
% Xdp_hat_ref=max(abs(Xdp_hat_h(:,pos_channel_1(2))).^2,abs(Xdp_hat_h(:,pos_channel_2(1))).^2);
% figure(101),plot(theta, 10*log10(amp.^2), 'bo', xp_rec, 10*log10(abs(Xdp_hat_h(:,pos_channel_1(2))).^2/Xdp_hat_ref), 'r*-');
% hold on
% plot(xp_rec, 10*log10(abs(Xdp_hat_h(:,pos_channel_2(1))).^2/Xdp_hat_ref), 'rx-');
% %axis([0,180,min([10*log10(amp.^2); 10*log10(abs(Xdp_hat_h(:,pos_channel_2(1))).^2/max(abs(Xdp_hat_h(:,pos_channel_2(1))).^2))]),max([10*log10(amp.^2); 10*log10(abs(Xdp_hat_h(:,pos_channel_2(1))).^2/max(abs(Xdp_hat_h(:,pos_channel_2(1))).^2))])+3]);
% xlabel('DOA (degrees)', 'fontsize',12); ylabel('Power (dB)','fontsize',12);
% legend('True DOAs','target-1','target-2');
% axis tight;
%% ------------------------RVM zhangmeng Liu---------------------
X_RVM_hat_h=zeros(N,C_N);
[X_RVM_hat alpha alpha0] = mt_CS(A_temp,Y_temp);
X_RVM_hat_h(:,pos_channel)=X_RVM_hat;
figure
imagesc(f,grid,abs(X_RVM_hat_h))
xlabel('frequency channels (Hz)')
ylabel('direction')
%% testing
theta_ex=[56 58 60];
%theta_ex=[106 108 110];
pos_theta_ex = fix(theta_ex/resolution)+1;
theta_test = [min(theta_ex):0.01:max(theta_ex)];
pos_test=[1:min(pos_theta_ex)-1 max(pos_theta_ex)+1:N];
temp_func = zeros(length(theta_test),1);
temp_func_inv = zeros(length(theta_test),1);
%f_test_temp_1=  [10 20 29 39 48 58 67 77 86 95];
%f_test_temp_1=[13 26 38 51 63 76 88 100 112 125];
%f_test_temp_1=[13 26 38 51 63 76 88 100 112 125];
f_test_temp_1 = pos_channel;
pos_test_temp = zeros(length(f_test_temp_1),1);
for nn=1:length(f_test_temp_1)
    if isempty(find(pos_channel==f_test_temp_1(nn)))
       pos_test_temp(nn) = 0;
    else
    pos_test_temp(nn)=find(pos_channel==f_test_temp_1(nn));
    end
end
pos_test_temp_temp=find(pos_test_temp);
pos_test_temp = pos_test_temp(pos_test_temp_temp);
f_temp = f(pos_channel);
a_steer=zeros(M,1);
d_steer=zeros(M,1);
d_M=d*([1:M]-(M+1)/2)';
for kk = 1:length(theta_test);
 for nn=1:length(pos_test_temp)
%  nn=11;
 f_t = f_temp (pos_test_temp(nn));
 lambda_temp = C./f_t;
 a_steer = exp(-2*i * pi *(d_M)* cos(theta_test(kk)/180*pi)/lambda_temp);
 d_steer = 2*i * pi *(d_M) * sin(theta_test(kk)/180*pi)/lambda_temp.*a_steer;
 A_temp_test = A_temp((pos_test_temp(nn)-1)*M+1:pos_test_temp(nn)*M,pos_test);
 R_y_f = Y_temp(:,pos_test_temp(nn))*Y_temp(:,pos_test_temp(nn))';
 Sigma_test = A_temp_test*diag(alpha(pos_test))* A_temp_test' + alpha0/T*sqrt(2)*eye(M);
 temp_func(kk) = temp_func(kk) + ...
     (a_steer'*inv(Sigma_test)*(a_steer*a_steer'*inv(Sigma_test)*R_y_f-R_y_f*inv(Sigma_test)*a_steer*a_steer')*inv(Sigma_test)*d_steer)/(a_steer'*inv(Sigma_test)*R_y_f*inv(Sigma_test)*a_steer);
 %temp_func(kk) = (a_steer'*inv(Sigma_test)*(a_steer*a_steer'*inv(Sigma_test)*R_y_f-R_y_f*inv(Sigma_test)*a_steer*a_steer')*inv(Sigma_test)*d_steer)/(a_steer'*inv(Sigma_test)*R_y_f*inv(Sigma_test)*a_steer);
 end
 temp_func_inv(kk)=1/(real(temp_func(kk)));
end
% figure
% plot(theta_test,real(temp_func));
figure
plot(theta_test,10*log(abs(temp_func_inv)/max((abs(temp_func_inv)))));
est_theta_pos_1=find(abs(temp_func_inv)==max(abs(temp_func_inv)));
est_theta_1 = theta_test(est_theta_pos_1)

theta_ex=[106 108 110];
pos_theta_ex = fix(theta_ex/resolution)+1;
theta_test = [min(theta_ex):0.01:max(theta_ex)];
pos_test=[1:min(pos_theta_ex)-1 max(pos_theta_ex)+1:N];
temp_func = zeros(length(theta_test),1);
temp_func_inv = zeros(length(theta_test),1);
%f_test_temp_1= pos_channel_2;
%f_test_temp_1=[13 26 38 51 63 76 88 100 112 125];
%f_test_temp_1=[13 26 38 51 63 76 88 100 112 125];
f_test_temp_1 = pos_channel;
pos_test_temp = zeros(length(f_test_temp_1),1);
for nn=1:length(f_test_temp_1)
    if isempty(find(pos_channel==f_test_temp_1(nn)))
       pos_test_temp(nn) = 0;
    else
    pos_test_temp(nn)=find(pos_channel==f_test_temp_1(nn));
    end
end
pos_test_temp_temp=find(pos_test_temp);
pos_test_temp = pos_test_temp(pos_test_temp_temp);
f_temp = f(pos_channel);
a_steer=zeros(M,1);
d_steer=zeros(M,1);
d_M=d*([1:M]-(M+1)/2)';
for kk = 1:length(theta_test);
 for nn=1:length(pos_test_temp)
%  nn=11;
 f_t = f_temp (pos_test_temp(nn));
 lambda_temp = C./f_t;
 a_steer = exp(-2*i * pi *(d_M)* cos(theta_test(kk)/180*pi)/lambda_temp);
 d_steer = 2*i * pi *(d_M) * sin(theta_test(kk)/180*pi)/lambda_temp.*a_steer;
 A_temp_test = A_temp((pos_test_temp(nn)-1)*M+1:pos_test_temp(nn)*M,pos_test);
 R_y_f = Y_temp(:,pos_test_temp(nn))*Y_temp(:,pos_test_temp(nn))';
 Sigma_test = A_temp_test*diag(alpha(pos_test))* A_temp_test' + alpha0/T*sqrt(2)*eye(M);
 temp_func(kk) = temp_func(kk) + ...
     (a_steer'*inv(Sigma_test)*(a_steer*a_steer'*inv(Sigma_test)*R_y_f-R_y_f*inv(Sigma_test)*a_steer*a_steer')*inv(Sigma_test)*d_steer)/(a_steer'*inv(Sigma_test)*R_y_f*inv(Sigma_test)*a_steer);
 %temp_func(kk) = (a_steer'*inv(Sigma_test)*(a_steer*a_steer'*inv(Sigma_test)*R_y_f-R_y_f*inv(Sigma_test)*a_steer*a_steer')*inv(Sigma_test)*d_steer)/(a_steer'*inv(Sigma_test)*R_y_f*inv(Sigma_test)*a_steer);
 end
 temp_func_inv(kk)=1/(real(temp_func(kk)));
end
% figure
% plot(theta_test,real(temp_func));
figure
plot(theta_test,10*log(abs(temp_func_inv)/max((abs(temp_func_inv)))));
est_theta_pos_1=find(abs(temp_func_inv)==max(abs(temp_func_inv)));
est_theta_2 = theta_test(est_theta_pos_1)
Error_noise_power = abs(alpha0/T*sqrt(2)-Sigma_a^2)/Sigma_a^2
% pos_theta_1=find(theta_test==theta(1));
% pos_theta_2=find(theta_test==theta(2));
% hold on;
% stem(theta_test(pos_theta_1),0,'r-*')
% hold on;
% stem(theta_test(pos_theta_2),0,'r-*')
% 
%% -----------------------MTL-CS-DP-(off-grid)---------------------------%% 
% Xdp_offgrid_hat_h=zeros(N,C_N);
% [Xdp_offgrid_hat X_Beta]=bmtl_DP_offgrid_2(A_temp,B_temp,Y_temp,resolution);
% Xdp_offgrid_hat_h(:,pos_channel)=Xdp_offgrid_hat;
% 
% figure
% imagesc(f,grid,abs(Xdp_offgrid_hat_h))
% xlabel('frequency channels (Hz)')
% ylabel('direction')
% 
% pos_a=find(sum(abs(X_Beta),2)>0.001);
% temp_Beta = X_Beta;
% beta=zeros(N,1);
% for n=1:length(pos_a)
%     temp_pos=find(abs(X_Beta(pos_a(n),:))>0.001);
%     beta(pos_a(n))=mean(X_Beta(pos_a(n),temp_pos));
% end
% % xp_rec = grid + beta;
% % Xdp_offgrid_hat_ref=max(abs(Xdp_offgrid_hat_h(:,pos_channel_1(2))).^2,abs(Xdp_offgrid_hat_h(:,pos_channel_2(1))).^2);
% % figure(100),plot(theta, 10*log10(amp.^2), 'bo', xp_rec, 10*log10(abs(Xdp_offgrid_hat_h(:,pos_channel_1(1))).^2/Xdp_offgrid_hat_ref), 'r*-');
% % hold on
% % plot(xp_rec, 10*log10(abs(Xdp_offgrid_hat_h(:,pos_channel_2(1))).^2/Xdp_offgrid_hat_ref), 'rx-');
% % % axis([0,180,min([10*log10(amp.^2); 10*log10(abs(Xdp_offgrid_hat_h(:,pos_channel_2(1))).^2/Xdp_offgrid_hat_ref)]),...
% % %     max([10*log10(amp.^2); 10*log10(abs(Xdp_offgrid_hat_h(:,pos_channel_2(1))).^2/Xdp_offgrid_hat_ref)])+3]);
% % xlabel('DOA (degrees)', 'fontsize',12); ylabel('Power (dB)','fontsize',12);
% % legend('True DOAs','target-1','target-2');
% % axis tight;
% 
% est_theta=grid(pos_a)+beta(pos_a);
% 
% Err_1(nnn) = Err_1(nnn)+min(abs(est_theta-theta(1)));
% Err_2(nnn) = Err_2(nnn)+min(abs(est_theta-theta(2)));
%% ------------------------------CSSM------------------------
azi = theta+2*rand(2,1);
Rn  = zeros(M); Ry = zeros(M);
fc_n= 16*pos_channel(fix(length(pos_channel)/2));
lambda_temp_1=C./f_temp;
for kkk = 1:length(pos_channel)
Rx = (Y_temp(:,kkk)-mean(Y_temp(:,kkk)))*(Y_temp(:,kkk)-mean(Y_temp(:,kkk)))';
% For cssm
ksv=zeros(M,length(azi));
ksv0=zeros(M,length(azi));
for m = 1:M
for n = 1:length(azi)
ksv(m,n)=exp(-2*i * pi *d*(m-(M+1)/2) * cos(azi(n)/180*pi)/lambda_temp_1(kkk));
ksv0(m,n) = exp(-2*i * pi *d*(m-(M+1)/2) * cos(azi(n)/180*pi)/lambda_temp_1(fix(length(pos_channel)/2)));
end
end
prod = ksv*ksv0';
[Uf D Vf] = svd(prod);
TT = Vf*Uf';
Ryf{kkk} = TT*Rx*TT';
[aa bb] = eigs(Rx,M);
db = diag(bb);
sigv(kkk) = mean(db(K+1:end));
Ry = Ry + Ryf{kkk}/(length(pos_channel));%;
Rn = Rn + sigv(kkk)*TT*TT';
Fa{kkk} = aa(:,1:K);
W{kkk} = aa(:,K+1:end);
end
F1 = Fa;
[Ucssm Vcssm] = eigs(inv(Rn)*Ry,M);
Uncssm = Ucssm(:,K+1:end);
grid_test=(0:0.01:180)';
AA=zeros(M,length(grid_test));
for n = 1:length(grid_test)
for m = 1:M
AA(m,n) = exp(-2*i * pi *d*(m-(M+1)/2) * cos(grid_test(n)/180*pi)/lambda_temp_1(fix(length(pos_channel)/2)));
end
end
P_cssm=zeros(1,length(grid_test));
P_tops=zeros(1,length(grid_test));
for tta = 1:length(grid_test);
P_cssm(tta) = (1./abs(AA(:,tta)'*Uncssm*Uncssm'*AA(:,tta)));
end
figure('name','P_CSSM')
plot(grid_test,10*log10(abs(P_cssm)/max(abs(P_cssm))),'r-');
hold on;plot(theta, 10*log10([1 1].^2), 'bo');
xlabel('DOA (Deg.)'); ylabel('Gain (dB)');
    end
Err_1(nnn) = Err_1(nnn)/N_exp;
Err_2(nnn) = Err_2(nnn)/N_exp;
end
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