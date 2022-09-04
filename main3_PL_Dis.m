clear;
close all;
Nt = 64; % number of trans antennas
K = 64;  % number of users
Rs = 1; % symbol rate
sigma2 = 1.8; % noise density
f0 = 1;% 1GHz
f = f0/6; % multiples relative to 6GHz
Gt = 15;
Gr = 15;

SNRdB=50;

d_n=0.01:0.005:0.08;
sum_rate1=zeros(1,length(d_n));
sum_rate2=zeros(1,length(d_n));
sum_rate3=zeros(1,length(d_n));
sum_rate4=zeros(1,length(d_n));
sum_rate5=zeros(1,length(d_n));
sum_rate6=zeros(1,length(d_n));

BER1=zeros(1,length(d_n));
BER2=zeros(1,length(d_n));
BER3=zeros(1,length(d_n));
BER4=zeros(1,length(d_n));
BER5=zeros(1,length(d_n));
BER6=zeros(1,length(d_n));

for i_d=1:length(d_n)
    for i_SNR=1:length(SNRdB)
        pathloss=148.1 + 37.6*log10(d_n(i_d)) + 20*log10(f);
        var=10^((-pathloss + Gt + Gr)/10);
        SNR=10^(SNRdB(1,i_SNR)/10);
        P=sigma2*SNR;

        error_cnt1=0;
        error_cnt2=0;
        error_cnt3=0;
        error_cnt4=0;
        error_cnt5=0;
        error_cnt6=0;
        data_rate1=0;
        data_rate2=0;
        data_rate3=0;
        data_rate4=0;
        data_rate5=0;
        data_rate6=0;

        loop_num=1000;

        for i_loop=1:loop_num
            s=QPSK_mapper(K,Rs);
            H_real=randn(K,Nt);
            H_imag=randn(K,Nt);
            H=(sqrt(var/2))*complex(H_real,H_imag);
            %H=(sqrt(1/2))*complex(H_real,H_imag);
            n_real=sqrt(sigma2/2)*randn(K,Rs);
            n_imag=sqrt(sigma2/2)*randn(K,Rs);
            n=complex(n_real,n_imag);


            % ZF
            [W1]=ZF(H,K);
            [count_temp,sum_rate_temp] = receiver(H,W1,Nt,K,Rs,P,n,s);
            error_cnt1=error_cnt1+count_temp;
            data_rate1=data_rate1+sum_rate_temp;

            % MRT
            [W2]=MRT(H,K);
            [count_temp,sum_rate_temp] = receiver(H,W2,Nt,K,Rs,P,n,s);
            error_cnt2=error_cnt2+count_temp;
            data_rate2=data_rate2+sum_rate_temp;

            % RZF
            [W3]=RZF(H,K,sigma2/P);
            [count_temp,sum_rate_temp] = receiver(H,W3,Nt,K,Rs,P,n,s);
            error_cnt3=error_cnt3+count_temp;
            data_rate3=data_rate3+sum_rate_temp;

            % PZF
            [H4,W4]=PZF(H,K,Nt,4);
            [count_temp,sum_rate_temp] = receiver(H4,W4,Nt,K,Rs,P,n,s);
            error_cnt4=error_cnt4+count_temp;
            data_rate4=data_rate4+sum_rate_temp;

            % SVD
            [H5,W5]=SVD(H,K,Nt);
            [count_temp,sum_rate_temp] = receiver(H5,W5,Nt,K,Rs,P,n,s);
            error_cnt5=error_cnt5+count_temp;
            data_rate5=data_rate5+sum_rate_temp;
        end

        BER1(1,i_d)=error_cnt1/(loop_num*K*Rs);
        BER2(1,i_d)=error_cnt2/(loop_num*K*Rs);
        BER3(1,i_d)=error_cnt3/(loop_num*K*Rs);
        BER4(1,i_d)=error_cnt4/(loop_num*K*Rs);
        BER5(1,i_d)=error_cnt5/(loop_num*K*Rs);

        sum_rate1(1,i_d)=data_rate1/loop_num/Rs;
        sum_rate2(1,i_d)=data_rate2/loop_num/Rs;
        sum_rate3(1,i_d)=data_rate3/loop_num/Rs;
        sum_rate4(1,i_d)=data_rate4/loop_num/Rs;
        sum_rate5(1,i_d)=data_rate5/loop_num/Rs;

    end
end

d_n = d_n * 1000;
% figure
% semilogy(d_n,BER1(1,:),'r-x');
% hold on;
% semilogy(d_n,BER2(1,:),'b-<');
% hold on;
% semilogy(d_n,BER3(1,:),'c-p');
% hold on;
% semilogy(d_n,BER4(1,:),'k-s');
% hold on;
% semilogy(d_n,BER5(1,:),'g-o');
% hold on;
% legend('ZF','MRT','RZF','PZF','SVD','FontName','Times New Roman');
% xlabel('Distance (m)','FontName','Times New Roman');
% ylabel('SER','FontName','Times New Roman');
% hold on;
% title(['SER vs PL Dis(BS, UE) with K = ', num2str(K), ', Nt = ',num2str(Nt), ', f = ', num2str(f0), ' GHz, SNR = 50 dB'],'FontName','Times New Roman');

figure
plot(d_n,sum_rate1(1,:),'r-x');
hold on;
plot(d_n,sum_rate2(1,:),'b-<');
hold on;
plot(d_n,sum_rate3(1,:),'c-p');
hold on;
plot(d_n,sum_rate4(1,:),'k-s');
hold on;
plot(d_n,sum_rate5(1,:),'g-o');
hold on;
%legend('ZF','MF','SVD','BD','MMSE','SLNR');
legend('ZF','MRT','RZF','PZF','SVD','FontName','Times New Roman');
xlabel('Distance (m)','FontName','Times New Roman');
ylabel('Sum Rate (bit/s/Hz)','FontName','Times New Roman');
hold on;
title(['Sum Rate vs PL Dis(BS, UE) with K = ', num2str(K), ', Nt = ',num2str(Nt), ', f = ', num2str(f0), ' GHz, SNR = 50 dB'],'FontName','Times New Roman');