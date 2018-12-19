%---------------------------Exact mercury/wf-------------------------------
%tol=1e-5;   %Tolerance

c=0;   %Original lower bound for internal bisection
e=0;   %Original lower bound for internal bisection

d=1500;  %Original upper bound for internal bisection
f=1500;  %Original upper bound for internal bisection

Gamma1_dB=-20:25;
Gamma2_dB=Gamma1_dB+10;
Gamma3_dB=Gamma2_dB+10;  %xx-fold gain

Gamma1=10.^(Gamma1_dB/10);
Gamma2=10.^(Gamma2_dB/10);
Gamma3=10.^(Gamma3_dB/10);

p1=zeros(1,length(Gamma1));
p2=zeros(1,length(Gamma1));
p3=zeros(1,length(Gamma1));

p11=zeros(1,length(Gamma1));
p21=zeros(1,length(Gamma1));
p31=zeros(1,length(Gamma1));

C_1=zeros(1,length(Gamma1));
C_2=zeros(1,length(Gamma1));
C_3=zeros(1,length(Gamma1));
C_4=zeros(1,length(Gamma1));
for n=1:length(Gamma1)
    gamma1=Gamma1(n);  %64-QAM
    gamma2=Gamma2(n);
    gamma3=Gamma3(n);
    
    a=gamma1*MMSE_64_QAM_23_new(300*gamma1/111);  %Original lower bound for internal bisection
    b=gamma3*MMSE_64_QAM_23_new(300*gamma1/111);  %Original upper bound for internal bisection
    %max1=-1+ceil((log(b-a)-log(tol))/log(2));   %Number of iterations
    
    %d=Bisection_4_PAM(0,100,1e-5,a,gamma2);  %Original upper bound
    %f=Bisection_4_PAM(0,100,1e-5,a,gamma2);  %Original upper bound
    
    for k=1:20
        eta=(a+b)/2;  %bisection
        
        %rou_1a=Bisection_64_QAM_23_new(c,d,1e-5,a,gamma1);
        %rou_2a=Bisection_64_QAM_23_new(e,f,1e-5,a,gamma2);
        %rou_3a=Bisection_64_QAM_23_new(e,f,1e-5,a,gamma3);
        
        %rou_1b=Bisection_64_QAM_23_new(c,d,1e-5,b,gamma1);
        %rou_2b=Bisection_64_QAM_23_new(e,f,1e-5,b,gamma2);
        %rou_3b=Bisection_64_QAM_23_new(e,f,1e-5,b,gamma3);
        
        rou_1e=Bisection_64_QAM_23_new(c,d,1e-5,eta,gamma1);
        rou_2e=Bisection_64_QAM_23_new(e,f,1e-5,eta,gamma2);
        rou_3e=Bisection_64_QAM_23_new(e,f,1e-5,eta,gamma3);
        
        %f_a=(1/(3*gamma1))*rou_1a+(1/(3*gamma2))*rou_2a+(1/(3*gamma3))*rou_3a-1;
        %f_b=(1/(3*gamma1))*rou_1b+(1/(3*gamma2))*rou_2b+(1/(3*gamma3))*rou_3b-1;
        f_e=(1/(3*gamma1))*rou_1e+(1/(3*gamma2))*rou_2e+(1/(3*gamma3))*rou_3e-1;
        
        if f_e==0   %Find the root
            p1(n)=rou_1e/gamma1;
            p2(n)=rou_2e/gamma2;
            p3(n)=rou_3e/gamma3;
            break
        elseif f_e<0
            b=eta;
        else
            a=eta;
        end
        if abs(f_e)<1e-3
            p1(n)=rou_1e/gamma1;
            p2(n)=rou_2e/gamma2;
            p3(n)=rou_3e/gamma3;
            break
        end
    end
end

%----------------------Constallation Constrained WF------------------------
p=10.^(Gamma1_dB/10);  %16-QAM
Pt=3;  %Total power

g1=1;
g2=10;
g3=100;  %Channel gain

M=64;  %Constellation order
tol=1e-5;   %Tolerance

p_1=zeros(1,length(p));
p_2=zeros(1,length(p));
p_3=zeros(1,length(p));
%-----------------------------Constellation WF-----------------------------
for n=1:length(p)
    a=0;   %Original lower bound
    b=g3*(1-1/M);   %Original upper bound
    max1=-1+ceil((log(b-a)-log(tol))/log(2));   %Number of iterations
    
    for k=1:max1+1
        lamda=(a+b)/2;   %Bisection
        %p_1a=constellation1(p(n),g1,M,a);   %Power allocation at a
        %p_2a=constellation1(p(n),g2,M,a);
        %p_3a=constellation1(p(n),g3,M,a);
        
        %p_1b=constellation1(p(n),g1,M,b);   %Power allocation at b
        %p_2b=constellation1(p(n),g2,M,b);
        %p_3b=constellation1(p(n),g3,M,b);
        
        p_1l=constellation1(p(n),g1,M,lamda);   %Power allocation at lamda
        p_2l=constellation1(p(n),g2,M,lamda);
        p_3l=constellation1(p(n),g3,M,lamda);
        
        %fa=p_1a+p_2a+p_3a-Pt;
        %fb=p_1b+p_2b+p_3b-Pt;
        fl=p_1l+p_2l+p_3l-Pt;
        
        if fl==0   %Find the root
            p_1(n)=p_1l;
            p_2(n)=p_2l;
            p_3(n)=p_3l;
            break
        elseif fl<0
            b=lamda;
        else
            a=lamda;
        end
        if b-a<tol
            p_1(n)=p_1l;
            p_2(n)=p_2l;
            p_3(n)=p_3l;
            break
        end
    end
end

plot(Gamma1_dB,p1)
hold on
grid on
plot(Gamma1_dB,p_1,'--')

plot(Gamma1_dB,p2)
plot(Gamma1_dB,p_2,'--')

plot(Gamma1_dB,p3)
plot(Gamma1_dB,p_3,'--')

xlabel('P/dB')
ylabel('Power Allocation')
ylim([0 3])
legend('Exact mercury/waterfilling','AOPA')