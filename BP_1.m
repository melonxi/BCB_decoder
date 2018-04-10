%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic                     
amp=1;                          %����BPSK  amplitude(���)Ϊ1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LDPC����
t = 0.5;                          %�����ϵ��
iteration=5;
load H.mat;
[rows,cols]=size(H);
rate=(cols-rows)/cols;          %����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
someSNR=0.5;               %���������
BER=zeros(1,length(someSNR));
frame = 1000;                  %����֡��
% count=1;
for S_num=1:length(someSNR) 
    total_num=0;
    hard_total_num=0;
    SNR=someSNR(S_num);
    EbNo=10.^(SNR/10);          %��������ɷֱ���ʾת��Ϊ���ʱȱ�ʾ
    sigma=1/sqrt(2*rate*EbNo);
    No=2*sigma^2;
    for i=1:frame
        i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        x=round(rand(1,cols-rows)); 
        c=mod(P*x',2); %  c' = P*s'
        u1=[c' x];      %%  u=[c | s]
        u=reorder_bits(u1,rearranged_cols);
        check1=sum(mod(u*H',2)) %�����Ƿ�У��ɹ�
        %x = (sign(randn(1,size(G,1)))+1)/2; % random bits
        %y = mod(x*G,2);                     % coding 
        s = 2*u-1;                          % BPSK modulation
        len_symbol=length(s);
        Sigma1 = related_matrix(cols,t);
        related_n = sqrtm(Sigma1);
        noise = sigma*randn(1,cols);
        related_noise = noise*related_n; 
        y = s + related_noise;   % AWGN transmission
        scale(1:length(u))=1;
        [s_hat1b,s_hat1s] = decodeLogDomain(y,H,sigma,iteration);
        noise_hat = y - s_hat1b;
        noise_hat_matrix(i,:) = noise_hat;
        y_matrix(i,:) = y;
        u_matrix(i,:) = u;
        errors=length(find(u~=s_hat1b)); 
        total_num=total_num+errors;
        total_num%������������
        
    end  
  BER(S_num)=total_num/(frame*cols);
end
toc
save noise_hat_t0.5_s0.5.txt -ascii noise_hat_matrix;
save y_t0.5_s0.5.txt -ascii y_matrix;
save u_t0.5_s0.5.txt -ascii u_matrix;
BER