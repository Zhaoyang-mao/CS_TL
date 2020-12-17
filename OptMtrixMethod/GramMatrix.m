%基于特征值分解的测量矩阵优化矩阵的计算
%author:mao
%时间：2019年11月29日
function [A_matrix,Phi_matrix] = GramMatrix(phi,Y)
%参数说明
%phi：测量矩阵为m*n  这里直接选用128*256
%Y：为稀疏矩阵为n*n  这里为图像的DCT变换为256*256

D = phi*Y;              %
%%          对D进行列单位化
[m,n] =size(D);
theat = ((n/m)^2)*(m-1)-n;
D_I = I_vec(D);
%%          Gram矩阵
% G = D_I'*D_I;               %得到Gram矩阵G
%%          循环开始
i=0;
while(1)
    G = D_I'*D_I;%得到Gram矩阵G
    [V,H] = eig(G);     %对矩阵进行特征分解得到，特征值H，和特征向量V
    aa = im2bw(H,1e-14);
    H_nm = aa*(n/m);      %将对角元素中非0的项赋值为n/m,H尖
    L = [sqrt(n/m)*eye(m),zeros(m,n-m)];%L的大小为m*n 这里为128*256
    L = rot90(L,2);                 %对其旋转180度
    D_jian = L*V';
    D_jian_I = I_vec(D_jian);       %对D尖进行单位化
    NewGram = D_jian_I'*D_jian_I;
    SumNewGram = CalSumDiga(NewGram);
    if abs(SumNewGram-theat) < 0.1 || i==50 %迭代终止条件
        A_matrix = D_I*inv(Y);
        Phi_matrix = D_I;
        break
    end
    i= i+1;
    D_I = D_jian_I;
%     disp(['第',num2str(i),'次迭代','误差为：',num2str(abs(SumNewGram-theat))]);
%     pause(0.1);
end
end
%%              列单位化矩阵函数
function Ivec = I_vec(A)
% [m1,n1] = size(A);        %测量矩阵的大小
% vec = ones(1,m1);        %生成1*m全为1的矩阵
% A_d = 1./(sqrt(vec*(A.^2))); %每列的平方和开根号
% A_dd = repmat(A_d,m1,1);      %生成一个m行都是D_d的矩阵
% Ivec = A.*A_dd;              %到此矩阵D的列单位化完毕
Ivec = A*diag(1./sqrt(sum(A.*A)));
end

%%              计算矩阵的非对角元素之和
function SumDiga = CalSumDiga(A)
[m2,n2] = size(A);
NI = [~eye(m2),ones(m2,n2-m2)];
Vect = A.*NI;
SumDiga = sum(Vect(:));
end