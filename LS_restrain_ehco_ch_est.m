%ֻ�Խ����źŽ��й��ƣ�Ȼ�󷵻��������������Ƶ����źż�ȥ��Ȼ����ʹ��ZF���������Զ���źţ�����������������������
%����治ͬ�㻹��û�й� FIR filter!!
%r = Xh + exp(i*w0)FYg + n;���Ǳ������н�h ,g��������1
%X,Y ����ȥ��cp,ֻ�������ݲ��ֵĴ�������Ϣ����
%X = [X1 X2],�ֱ�Ϊ��������1�Ϸ��͵��źţ��ͽ�������2�Ϸ��͵��ź�
%Y = [Y1 Y2],�ֱ�Ϊ��������1�Ϸ��͵��źţ��ͽ�������2�Ϸ��͵��ź�
%FӦ��Ϊÿ���ز��Ͼ��е���λ���ʽΪ
%��ʼ��λ��Ҳ��Ϊ��w0 = 0;
%w = 0; coffients =i*pi*w;               
% F = diag([1,exp(coffients*1),exp(coffients*2),exp(coffients*3),exp(coffients*4),exp(coffients*5),exp(coffients*6),exp(coffients*7),exp(coffients*8),exp(coffients*(N_sym - 1))]); %�����ź���ÿһ�����ŵ���λ��
%FӦ��Ҳ��һ�����ƣ������Ȱ�����������£�����ȡ����0����ʵ����w = ��f_carrier_near -f_carrier_far��*T/R_rx(TΪ���������R_rx:�������½������ߵĹ�������)
%������Ҫ�ڶ������ٶ�h_near���й��ƣ������ڻز����ƣ����Եõ���׼ȷ��Զ�˽ڵ��źŵĹ���[h(��׼ȷ)��exp(i*w0)g(g��Զ�˵�Ƶ��)]
%function [Y_restraint_ehco_est,X_restraint_ehco_est,R_new_temp,MSE_far_after_ch_est,h_mse_all_ch] = LS_restrain_ehco_ch_est(transmit_signal,transmit_signal_near,recv_signal,PrefixRatio,N_subc,N_sym,N_Rx_ant,F)
function [X_restraint_ehco_est,MSE_near_after_ch_est,h_mse_all_near_ch] = LS_restrain_ehco_ch_est(transmit_signal,transmit_signal_near,recv_signal,PrefixRatio,N_subc,N_sym,N_Rx_ant,F)
%cp���ȣ�������Զ���źŵ����ݷ��ͳ�����ͬ
cp_length = PrefixRatio * N_subc;
data_length = length(transmit_signal);

% cp_far = reshape(transmit_signal,data_length/N_sym,N_sym,N_Rx_ant);   
% cp_far = cp_far(1:cp_length,:,:); %257��ȥ��cp ���ݿ�ʼ��λ��
            
cp_near = reshape(transmit_signal_near,data_length/N_sym,N_sym,N_Rx_ant);
cp_near =cp_near(1:cp_length,:,:);
%%Զ�˽ڵ���ź�
% Y_far = reshape(transmit_signal,data_length/N_sym,N_sym,N_Rx_ant);   
% Y_far = Y_far((cp_length + 1):end,:,:); %257��ȥ��cp ���ݿ�ʼ��λ��
% Y1 = Y_far(:,:,(N_Rx_ant - 1));
% Y2 = Y_far(:,:,N_Rx_ant);

%%���˽ڵ���ź�(��ȥ��Ƶ������)
X_near = reshape(transmit_signal_near,data_length/N_sym,N_sym,N_Rx_ant);  %�����transmit_signal_near ���ڽ���ʱ���Ѿ�����100��
%X_near = X_near((cp_length + 1):end,:,:);
X1 = X_near(:,:,(N_Rx_ant - 1));  %���˽ڵ�����1�ϵ��ź�
X2 = X_near(:,:,N_Rx_ant);  %���˽ڵ�����2�ϵ��ź�
X_est = [X1 X2];

%%���������ϵ��ź�
R_all = reshape(recv_signal,data_length/N_sym,N_sym,N_Rx_ant);   
%R_all = R_all((cp_length + 1):end,:,:); %257��ȥ��cp ���ݿ�ʼ��λ��
R1 = R_all(:,:,(N_Rx_ant - 1));          %���սڵ�����1�ϵ��ź�
R2 = R_all(:,:,N_Rx_ant);          %���սڵ�����2�ϵ��ź�

%%Զ����������ߡ�������������ߵ��ŵ�����
%������������ߵĵ�һ�׶ε��ŵ�����
 h11 = inv(X1'*X1)*X1'*R1;
 h21 = inv(X2'*X2)*X2'*R1;
 h12 = inv(X1'*X1)*X1'*R2;
 h22 = inv(X2'*X2)*X2'*R2;
 h_near_est_1stage =[h11 h12;h21 h22];


%�����źŵĹ���
r_est_near = X_est * h_near_est_1stage;
R_near_est = zeros(N_subc+cp_length,N_sym,N_Rx_ant); 
R_near_est(:,:,(N_Rx_ant - 1)) = r_est_near(:,1:N_sym);
R_near_est(:,:,N_Rx_ant) = r_est_near(:,(N_sym +1):end);

% %��һ�׶εĻز����ƣ�R_far_ehco ֻ����������������Ҫ��Զ���źţ�
% R_restrain_ehco = R_all -100 * R_near_est;
%  
% %�����˽����źŵĽ����ź�
% R_after_restrain_1 = R_restrain_ehco(:,:,(N_Rx_ant - 1)); 
% R_after_restrain_2 = R_restrain_ehco(:,:,N_Rx_ant);
% 
% %(F*[Y1.' Y2.']).' + n = [R_after_restrain_1  R_after_restrain_2]
% R_after_restrain_1_inv_F = inv(F)*R_after_restrain_1.';
% R_after_restrain_1_inv_F =  R_after_restrain_1_inv_F.';
% R_after_restrain_2_inv_F = inv(F)*R_after_restrain_2.';
% R_after_restrain_2_inv_F =  R_after_restrain_2_inv_F.';
% 
% %%Զ���ŵ��Ĺ���
% h_far11 = inv(Y1'*Y1)*Y1'* R_after_restrain_1_inv_F ;
% h_far21 = inv(Y2'*Y2)*Y2'* R_after_restrain_1_inv_F ;
% h_far12 = inv(Y1'*Y1)*Y1'* R_after_restrain_2_inv_F ;
% h_far22 = inv(Y2'*Y2)*Y2'* R_after_restrain_2_inv_F ;
% h_far_est_1stage = [h_far11 h_far12;h_far21 h_far22];
% 
% %Զ���źŵĹ���
% R_after_restrain_inv_F = [R_after_restrain_1_inv_F  R_after_restrain_2_inv_F];
% 
% Y_est_far = R_after_restrain_inv_F * inv(h_far_est_1stage);
% Y_est_far_1 =  Y_est_far(:,1:N_sym);
% Y_est_far_2 =  Y_est_far(:,(N_sym + 1):end);
% 
% %������Ҫ�ڶ������ٶ�h_near���й��ƣ������ڽ�������Ĳ��裬���Եõ���׼ȷ��Զ�˽ڵ��źŵĹ���
% %���� A = [r_near_est (F*Y_est_far_1.').' (F*Y_est_far_2.').'];
% %��ΪF�ǶԽ�Ϊ1�ľ������Ծ�û�����,
% %X_est = [X1 X2];
% 
% A = [X_est Y_est_far_1 Y_est_far_2];

%%%%��mse 
% %Զ���źŵĹ��ƺ��mse
% mean_Y1 = mean(Y1);
% mean_Y2 = mean(Y2);
% mean_Y_est_far_1 = mean(Y_est_far_1);
% mean_Y_est_far_2 = mean(Y_est_far_2);
% e_1 = mean_Y1 - mean_Y_est_far_1;
% e_2 = mean_Y2 - mean_Y_est_far_2;
% mse_1 = mse(abs(e_1));
% mse_2 = mse(abs(e_2));%����Լ����[1.015e-05;1.868e-05;];
% MSE_far_after_ch_est = [mse_1;mse_2];
%�����źŵĹ��ƺ��mse
mean_X1 = mean(X1);
mean_X2 = mean(X2);
mean_X_est_near_1 = mean(X1);
mean_X_est_near_2 = mean(X2);
e_1 = mean_X1 - mean_X_est_near_1;
e_2 = mean_X2 - mean_X_est_near_2;
mse_1 = mse(abs(e_1));
mse_2 = mse(abs(e_2));%����Լ����[1.015e-05;1.868e-05;];
MSE_near_after_ch_est = [mse_1;mse_2];


% %�����Ƶ���Ҫ��Զ���ź�������װ
% Y_new = cat(3,Y_est_far_1,Y_est_far_2);  %��ά������άΪ�ڼ�������
% Y_new_new = cat(1,cp_far,Y_new);         %��Ϊ���Ƴ�����Զ���źţ�����ֻ�ý�cp_far��Y_new,����������
% Y_restraint_ehco_est = reshape(Y_new_new,1,data_length,N_Rx_ant);%����Ϊ֮ǰ����cp���ź�

%�����Ƶ����Ƶ��Ľ����ź�������װ
X_new = cat(3,X1,X2);
% X_new_new = cat(1,cp_near,X_new);
X_restraint_ehco_est = reshape(X_new,1,data_length,N_Rx_ant);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R_new = cat(3,R_after_restrain_1,R_after_restrain_2);  %��ά������άΪ�ڼ�������
% R_new_new = cat(1,cp_far,R_new);         %��Ϊ���Ƴ�����Զ���źţ�����ֻ�ý�cp_far��Y_new,����������
% R_new_temp = reshape(R_new_new,1,data_length,N_Rx_ant);%����Ϊ֮ǰ����cp���ź�


%�ŵ���mse��
h = zeros((N_sym * N_Rx_ant),(N_sym * N_Rx_ant));
for i =1:(N_sym * N_Rx_ant)
h(i,i) = 1;
end

% h_erro_all_far_ch = mean(h) - mean(h_far_est_1stage);
% %h_erro_all_far_ch = diag(h_erro_all_far_ch);
% %h_erro_all_far_ch = diag(h_erro_all_far_ch);%����ȡ��ȡ�Խǣ����Ľ������һ����
% h_erro_all_far_ch = abs(h_erro_all_far_ch);
% h_mse_all_far_ch = mse(h_erro_all_far_ch); %��̫����
h_erro_all_near_ch = mean(h) - mean(h_near_est_1stage);
h_erro_all_near_ch = abs(h_erro_all_near_ch);
h_mse_all_near_ch= mse(h_erro_all_near_ch); %