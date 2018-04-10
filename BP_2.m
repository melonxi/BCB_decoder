%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;
y_matrix= textread('y_0.5_1.txt');
u_matrix = textread('u_0.5_1.txt');
noise_wave_matrix= textread('noise_wave_0.5_1.txt');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic                     
amp=1;                          %����BPSK  amplitude(���)Ϊ1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LDPC����
t=0.5;
iteration=25;
load H.mat;
[rows,cols]=size(H);
rate=(cols-rows)/cols;          %����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
someSNR=1;               %���������
BER=zeros(1,length(someSNR));
frame = 1000;
total_num=0;
hard_total_num=0;%����֡��
% count=1;
for S_num=1:length(someSNR) 
    
    SNR=someSNR(S_num);
    EbNo=10.^(SNR/10);          %��������ɷֱ���ʾת��Ϊ���ʱȱ�ʾ
    sigma=1/sqrt(2*rate*EbNo);
    No=2*sigma^2;
    for i=1:frame
        i
        y = y_matrix(i,:);
        u = u_matrix(i,:);
        check1=sum(mod(u*H',2))
        noise_wave = noise_wave_matrix(i,:);
        y_hat = y - noise_wave;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [s_hat2b,s_hat2s] = decodeLogDomain(y_hat,H,sigma,iteration);
        
        
        errors=length(find(u~=s_hat2b));
        errors
        total_num=total_num+errors;
        total_num %������������
       
    end
    BER(S_num)=total_num/(frame*cols);
end
toc
BER