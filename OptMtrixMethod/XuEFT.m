function [PhiNew] = XuEFT(Phi,Psi,m,muG,Alpha)
%%          参数说明
%Phi:测量矩阵
%Psi:稀疏矩阵
%m:测量值的长度
%
%%          初始化
%
Iter = 100;         %最大迭代次数
D_hat = Phi_k*Psi;
%列标准化
D_hat = D_hat./repmat(sqrt(sum(D_hat.^2,1)),size(D_hat,1),1);
%构造Gram矩阵
G_k = D_hat'*D_hat;

end