function [Heq,W] = PZF(H,K,Nt,lamada)
    F=zeros(Nt,Nt);
    for k1=1:Nt
        F(k1,:)=sqrt(1/Nt)*exp(1j*phase(H(k1,:)'));
    end

    Heq=H*F;
    W = Heq'*inv(Heq*Heq');

    for i = 1:K
        W(:,i) = W(:,i)/norm(W(:,i))*lamada;
    end
end
