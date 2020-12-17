%
%参考文献：一种基于迁移学习的感知矩阵优化方法
clc;close all
clear 
%%                  读取并处理数据
addpath('./../Src/','./../RecoverAlgorithm/','./../OptMtrixMethod/');
img = imread('../Data/peppers256.png');
img = imresize(img,[256,256]);      %将图像设置为256*256大小
img = double(img);                  %将数据转换为double类型
[m,n] = size(img);                  %图像的尺寸

%%                  参数设置
Param.Style = 4;            %4种实验
Param.Psnr = [];            %存放信噪比
Param.Time = [];            %存放运行时间
Param.Rate = 0.1:0.1:0.3;  %采样率
Param.Psnr = zeros(Param.Style,length(Param.Rate));
Param.Time = zeros(Param.Style,length(Param.Rate));
for i = 1:1:Param.Style
    for j = 1:1:length(Param.Rate)
%%                  构造稀疏基
        tic                     %计时开始
        Psi = DWT1(n);          %构造稀疏基，小波基
        Psi=Psi*diag(1./sqrt(diag(Psi'*Psi)));
%%                  测量矩阵构造
        mm = floor(n*Param.Rate(j));
        Phi = randn(mm,n);              %构造随机高斯矩阵
        for ii = 1:1:mm
            Phi(ii,:) = Phi(ii,:)/norm(Phi(ii,:));  %归一化处理
        end
        switch i
            case 1
                %随机高斯 不做变换
                disp(['GaussMatrix...']);
            case 2
                %特征分解法
                [~,Phi] = GramMatrix(Phi,Psi);
                disp(['EIGMatrix...']);
            case 3
                %elad方法
                EladIter = 100; dd1=0.8;dd2 = 0.95;     %参数设置
                [Phi] = DesignProjection(Psi,mm,EladIter,dd1,dd2,Phi);
                disp(['EladMatrix...']);
            case 4
                %本文方法
                PsiIter = randn(size(Psi));OurK = 50;Ourlambd = 0.35;
                Ourbet = 0.01;
                [Phi] = TOTL_EIG(img,Phi,Psi,PsiIter,OurK,Ourlambd,Ourbet);
                disp(['OurMatrix...']);
        end


        disp(['Sampling Rate = ',num2str(Param.Rate(j)),'...']);
%%                  重构算法
        %测量值
        s = Psi * img;
        y = Phi *img*Psi';
        %观测矩阵
        A = Phi * Psi';
        %重构
        for jj = 1:1:n
            Yhat(:,jj) = OMP(y(:,jj),A,50);
        end
        Yhat = Psi'*sparse(Yhat)*Psi;
        Yhat = full(Yhat);
        ErrorYhat = sum(sum(abs(Yhat-img).^2));
        PsnrYhat = 10*log10(255*255/(ErrorYhat/m/n));
        Param.Psnr(i,j) = PsnrYhat;
        Param.Time(i,j) = toc;
    end
end
% imshow(Yhat,[]);
for iii = 1:1:length(Param.Psnr)
   hold on
   plot(Param.Rate,Param.Psnr,'LineWidth',3);
   hold off
   legend('GaussMatrix','EIGMatrix','EladMatrix','OurMatrix');
   xlabel('Sampling Rate(%)');
   ylabel('Psnr(dB)')
end
