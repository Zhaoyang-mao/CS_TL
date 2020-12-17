function [ theta,erro_rn ] = OMP( y,A,K )
N=max(size(A));
M=min(size(A));
theta=zeros(1,N);   %  待重构的向量    
Base_t=[];              %  记录基向量的矩阵
r_n=y;                  %  残差值
for times=1:K;                                    %  迭代次数(有噪声的情况下,该迭代次数为K)
    for col=1:N;                                  %  恢复矩阵的所有列向量
        product(col)=abs(A(:,col)'*r_n);          %  恢复矩阵的列向量和残差的投影系数(内积值) 
    end
    [val,pos]=max(product);                       %  最大投影系数对应的位置，val值，pos位置
    Base_t=[Base_t,A(:,pos)];                       %  矩阵扩充，记录最大投影的基向量
    A(:,pos)=zeros(M,1);                          %  选中的列置零（实质上应该去掉，为了简单我把它置零）
    aug_y=(Base_t'*Base_t)^(-1)*Base_t'*y;   %  最小二乘,使残差最小
    r_n=y-Base_t*aug_y;                            %  残差
    erro_rn(times)=norm(r_n,2);      %迭代误差
    pos_array(times)=pos;                         %  纪录最大投影系数的位置
    if erro_rn(times)<1e-6 %
            break; %跳出for循环
    end
end
theta(pos_array)=aug_y;                           %  重构的向量
end