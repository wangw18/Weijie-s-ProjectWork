function [W] = MRT(H,K)
    W = H';

    for i = 1:K
        W(:,i) = W(:,i)/norm(W(:,i));
    end
end
