function [W] = ZF(H,K)
    W = H'*inv(H*H');

    for i = 1:K
        W(:,i) = W(:,i)/norm(W(:,i));
    end
end
