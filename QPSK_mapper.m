function [s] = QPSK_mapper(K,Tc)
    s=sqrt(0.5)*complex(sign(randn(K,Tc)),sign(randn(K,Tc)));
end