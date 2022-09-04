function [count,sum_rate] = receiver(H,W,Nt,K,Rc,P,n,s)
    yd=zeros(K,Rc);
    yi=zeros(K,Rc);
    y=zeros(K,Rc);
    yr=zeros(K,Rc);
    yi_d=zeros(K,Rc);
    sinr_n=zeros(K,Rc);
    sinr_d=zeros(K,Rc);
    sinr=zeros(K,Rc);
    Q_map=sqrt(0.5)*[1+1j; 1-1j; -1+1j; -1-1j];
    
    for k1=1:K
        yd(k1,:)=sqrt(P/Nt)*H(k1,:)*W(:,k1)*s(k1,:);

        for k2=1:K
            if k1~=k2
                yi(k1,:)=yi(k1,:)+H(k1,:)*W(:,k2)*s(k2,:);
            end
        end 
        yi(k1,:)=sqrt(P/Nt)*yi(k1,:);

        y(k1,:)=yd(k1,:)+yi(k1,:)+sqrt(1/Nt)*n(k1,:);
        [~,pos]=min(abs(y(k1,:)-Q_map));
        for n1=1:Rc
            yr(k1,n1)=Q_map(pos(n1));
        end

        for k3=1:K
            if k1~=k3
                yi_d(k1,:)=yi_d(k1,:)+H(k1,:)*W(:,k3)*W(:,k3)'*H(k1,:)';
            end
        end 
        yi_d(k1,:)=P/Nt*yi_d(k1,:);

        sinr_n(k1,:)=P/Nt*H(k1,:)*W(:,k1)*W(:,k1)'*H(k1,:)';
        sinr_d(k1,:)=yi_d(k1,:)+1;
        sinr(k1,:)=sum(sinr_n(k1,:))/sum(sinr_d(k1,:));
    end
    
    count=sum(sum(abs(yr-s)>0.001));
    sum_rate=real(sum(sum(log2(1+sinr))));
end

