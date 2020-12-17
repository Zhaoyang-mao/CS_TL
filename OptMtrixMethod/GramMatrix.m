%��������ֵ�ֽ�Ĳ��������Ż�����ļ���
%author:mao
%ʱ�䣺2019��11��29��
function [A_matrix,Phi_matrix] = GramMatrix(phi,Y)
%����˵��
%phi����������Ϊm*n  ����ֱ��ѡ��128*256
%Y��Ϊϡ�����Ϊn*n  ����Ϊͼ���DCT�任Ϊ256*256

D = phi*Y;              %
%%          ��D�����е�λ��
[m,n] =size(D);
theat = ((n/m)^2)*(m-1)-n;
D_I = I_vec(D);
%%          Gram����
% G = D_I'*D_I;               %�õ�Gram����G
%%          ѭ����ʼ
i=0;
while(1)
    G = D_I'*D_I;%�õ�Gram����G
    [V,H] = eig(G);     %�Ծ�����������ֽ�õ�������ֵH������������V
    aa = im2bw(H,1e-14);
    H_nm = aa*(n/m);      %���Խ�Ԫ���з�0���ֵΪn/m,H��
    L = [sqrt(n/m)*eye(m),zeros(m,n-m)];%L�Ĵ�СΪm*n ����Ϊ128*256
    L = rot90(L,2);                 %������ת180��
    D_jian = L*V';
    D_jian_I = I_vec(D_jian);       %��D����е�λ��
    NewGram = D_jian_I'*D_jian_I;
    SumNewGram = CalSumDiga(NewGram);
    if abs(SumNewGram-theat) < 0.1 || i==50 %������ֹ����
        A_matrix = D_I*inv(Y);
        Phi_matrix = D_I;
        break
    end
    i= i+1;
    D_I = D_jian_I;
%     disp(['��',num2str(i),'�ε���','���Ϊ��',num2str(abs(SumNewGram-theat))]);
%     pause(0.1);
end
end
%%              �е�λ��������
function Ivec = I_vec(A)
% [m1,n1] = size(A);        %��������Ĵ�С
% vec = ones(1,m1);        %����1*mȫΪ1�ľ���
% A_d = 1./(sqrt(vec*(A.^2))); %ÿ�е�ƽ���Ϳ�����
% A_dd = repmat(A_d,m1,1);      %����һ��m�ж���D_d�ľ���
% Ivec = A.*A_dd;              %���˾���D���е�λ�����
Ivec = A*diag(1./sqrt(sum(A.*A)));
end

%%              �������ķǶԽ�Ԫ��֮��
function SumDiga = CalSumDiga(A)
[m2,n2] = size(A);
NI = [~eye(m2),ones(m2,n2-m2)];
Vect = A.*NI;
SumDiga = sum(Vect(:));
end