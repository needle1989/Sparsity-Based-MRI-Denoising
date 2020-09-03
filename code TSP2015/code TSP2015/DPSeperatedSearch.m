function est_theta_DP_S_1=DPSeperatedSearch(Xdp_hat,alpha_DP,theta_ex,alpha0_DP,pos_channel,Y_temp,A_temp,theta);
i=sqrt(-1);
T=256;
M=10;
C=1500;
d=1.5;
f=[1:T/2]*16;
resolution=2;
pos_theta_ex = fix(theta_ex/resolution)+1;
N=length(0:resolution:180);
[temp pos_DP]=find(max(sum(alpha_DP(:,pos_theta_ex))));
alpha_DP_1 = alpha_DP(pos_DP(1),:);
theta_test = [min(theta_ex):0.01:max(theta_ex)];
pos_test=[1:min(pos_theta_ex)-1 max(pos_theta_ex)+1:N];
temp_func = zeros(length(theta_test),1);
temp_func_inv = zeros(length(theta_test),1);
[temp bbb_temp]=max(Xdp_hat);
pos_test_temp_1 = find(bbb_temp==pos_theta_ex(1));
pos_test_temp_2 = find(bbb_temp==pos_theta_ex(2));
% 
% if length(pos_test_temp_1)>=length(pos_test_temp_2)
% pos_test_temp = pos_test_temp_1;
% else
% pos_test_temp = pos_test_temp_2; 
% end
pos_test_temp=[pos_test_temp_1 pos_test_temp_2];
f_temp = f(pos_channel);
a_steer=zeros(M,1);
d_steer=zeros(M,1);
d_M=d*([1:M]-(M+1)/2)';
for kk = 1:length(theta_test);
 for nn=1:length(pos_test_temp)
%  nn=11;
 f_t = f_temp (pos_test_temp(nn));
 alpha_DP_temp=alpha_DP_1;
 lambda_temp = C./f_t;
 a_steer = exp(-2*i * pi *(d_M)* cos(theta_test(kk)/180*pi)/lambda_temp);
 d_steer = 2*i * pi *(d_M) * sin(theta_test(kk)/180*pi)/lambda_temp.*a_steer;
 A_temp_test = A_temp((pos_test_temp(nn)-1)*M+1:pos_test_temp(nn)*M,pos_test);
 R_y_f = Y_temp(:,pos_test_temp(nn))*Y_temp(:,pos_test_temp(nn))';
 Sigma_test = A_temp_test*diag(alpha_DP_temp(pos_test))* A_temp_test' + alpha0_DP/T*sqrt(2)*eye(M);%Sigma_a^2
 temp_func(kk) = temp_func(kk) + ...
     (a_steer'*inv(Sigma_test)*(a_steer*a_steer'*inv(Sigma_test)*R_y_f-R_y_f*inv(Sigma_test)*a_steer*a_steer')*inv(Sigma_test)*d_steer)/(a_steer'*inv(Sigma_test)*R_y_f*inv(Sigma_test)*a_steer);
 %temp_func(kk) = (a_steer'*inv(Sigma_test)*(a_steer*a_steer'*inv(Sigma_test)*R_y_f-R_y_f*inv(Sigma_test)*a_steer*a_steer')*inv(Sigma_test)*d_steer)/(a_steer'*inv(Sigma_test)*R_y_f*inv(Sigma_test)*a_steer);
 end
 temp_func_inv(kk)=1/(real(temp_func(kk)));
end
% figure
% plot(theta_test,real(temp_func));
est_theta_pos_1=find(abs(temp_func_inv)==max(abs(temp_func_inv)));
est_theta_DP_S_1 = theta_test(est_theta_pos_1);
% pos_temp=find(abs(theta-est_theta_DP_S_1)<2);
% figure
% plot(theta_test,10*log(abs(temp_func_inv)/max((abs(temp_func_inv)))));
% hold on;plot(theta, 10*log10([1].^2), 'bo');
% xlabel('DOA (Deg.)'); ylabel('Gain (dB)');

