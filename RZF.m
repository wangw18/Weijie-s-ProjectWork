function [W] = RZF(H,K,P)
% J=min(K,Nt);
% w=ones(J,1);
% W=zeros(Nt,K);
% 
%     for j1 = 1:J
%         W=W+w(j1)*(H.'*det(H)*inv(H))^j1*H.';
%         %W=W+0.01*w(j1)*(H*H')^j1*H*eye(K);
%     end
    alpha = K*P;
    W = H'*inv(H*H'+alpha*eye(K));
    
    for i = 1:K
        W(:,i) =W(:,i)/norm(W(:,i));
    end
end
