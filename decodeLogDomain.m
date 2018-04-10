%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vhat,z] = decodeLogDomain(rx_waveform,H,sigma,iteration)
%  此算法是对数域的和积译码算法
%  rx        : 接收信号矢量
%  H         : LDPC的校验矩阵
%  N0        : 噪声方差
%  iteration : 最大迭代次数
%  vhat      : 译码后的码字
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                
[M N] = size(H);                       
Lci = -(2*rx_waveform)/(sigma^2);       %Qji0/Qji1  似然比
Lrji = zeros(M, N);                    %生成一个cows行cols列的零矩阵
Pibetaij = zeros(M, N);                
L = repmat(Lci,M,1);                   %将Lci复制cows行
Lqij = L.*H;                           %表示校验节点与变量节点所含的信息量
[r, c] = find(H);                                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 1:iteration                    
   alphaij = sign(Lqij);       %符号位        
   betaij = abs(Lqij);         %绝对值       

   for l = 1:length(r)                 
      Pibetaij(r(l), c(l)) = log((exp(betaij(r(l), c(l))) + 1)/(exp(betaij(r(l), c(l))) - 1));    %Uj0                     
   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%校验节点消息处理
   for i = 1:M                 %逐行处理        
      c1 = find(H(i, :));   %找到每一行1的位置
      for k = 1:length(c1)              
         sumOfPibetaij = 0;            
         prodOfalphaij = 1;
         
         sumOfPibetaij = sum(Pibetaij(i, c1)) - Pibetaij(i, c1(k));   %计算第i行除去第k列的信息量   
         
         if sumOfPibetaij < 1e-20
            sumOfPibetaij = 1e-10;
         end       
         
         PiSumOfPibetaij = log((exp(sumOfPibetaij) + 1)/(exp(sumOfPibetaij) - 1));
         prodOfalphaij = prod(alphaij(i, c1))*alphaij(i, c1(k));    %求符号位
         Lrji(i, c1(k)) = prodOfalphaij*PiSumOfPibetaij;
      end % for k
   end % for i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%变量节点消息处理
   for j = 1:N                         
      r1 = find(H(:, j));              
      for k = 1:length(r1)        
         % Update L(qij) by summation of L(rij)\r1(k)
         Lqij(r1(k), j) = Lci(j) + sum(Lrji(r1, j)) - Lrji(r1(k), j);
      end % for k
      LQi = Lci(j) + sum(Lrji(r1, j)); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % 判决
      if LQi < 0
         vhat(j) = 1;
      else
         vhat(j) = 0;
      end
   end % for j
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%加校验语句，如果译码结果满足校验方程，则停止迭代
z = vhat;
z(z==0)=-1;
if (mod(H*vhat',2)==0) 
    break;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
end % for n
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%