function [Heq,W] = SVD(H,K,Nt)
    [u,d,v] = svd(H);
    V = v(:,1:K);
    U = u(:,1:K)';
    W = U;
    Heq = H*V;
    for i = 1:K
        W(:,i) = W(:,i)/norm(W(:,i));
    end
end
