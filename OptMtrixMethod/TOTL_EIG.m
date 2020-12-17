function [PhiNew] = TOTL_EIG(x,Phi,Psi0,PsiIter,k,lambd,bet)
    PsiNew = TOTL(x,Psi0,PsiIter,k,lambd,bet);
    [~,PhiNew] = GramMatrix(Phi,PsiNew);
end

function psi = TOTL(x,psi0,psi1,k,lambda,bet)
%     bet = 0.1;
    for ii = 1:k
        W = psi1*x;
        s_ii = sign(W)*(max(abs(W)-lambda/2,0));
        [U,S,V] = svd(s_ii*(x')+bet*psi0);
        psi1 = U*(V');
    end
    psi = psi1;
end