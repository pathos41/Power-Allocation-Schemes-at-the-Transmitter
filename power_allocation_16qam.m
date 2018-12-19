%---------------------------Exact mercury/wf-------------------------------
%tol=1e-5;   %Tolerance

c=0;   %Original lower bound for internal bisection
e=0;   %Original lower bound for internal bisection

d=1000;  %Original upper bound for internal bisection
f=1000;  %Original upper bound for internal bisection

Gamma1_dB=-20:25;
Gamma2_dB=Gamma1_dB+20;

Gamma1=10.^(Gamma1_dB/10);
Gamma2=10.^(Gamma2_dB/10);

p1=zeros(1,length(Gamma1));
p2=zeros(1,length(Gamma1));

p11=zeros(1,length(Gamma1));
p21=zeros(1,length(Gamma1));

for n=1:length(Gamma1)
    gamma1=Gamma1(n);  %16-QAM
    gamma2=Gamma2(n);
    
    a=gamma1*MMSE_16_QAM_23_new(200*gamma1/101);  %Original lower bound for external bisection
    b=gamma2*MMSE_16_QAM_23_new(200*gamma1/101);  %Original upper bound for external bisection
    
    %max1=-1+ceil((log(b-a)-log(tol))/log(2));   %Number of iterations
    
    %d=Bisection_4_PAM(0,100,1e-5,a,gamma2);  %Original upper bound
    %f=Bisection_4_PAM(0,100,1e-5,a,gamma2);  %Original upper bound
    %num=max(max1,1);
    for k=1:20
        eta=(a+b)/2;  %bisection
        
        %rou_1a=Bisection_16_QAM_23_new(c,d,1e-5,a,gamma1);
        %rou_2a=Bisection_16_QAM_23_new(e,f,1e-5,a,gamma2);
        
        %rou_1b=Bisection_16_QAM_23_new(c,d,1e-5,b,gamma1);
        %rou_2b=Bisection_16_QAM_23_new(e,f,1e-5,b,gamma2);
        
        rou_1e=Bisection_16_QAM_23_new(c,d,1e-5,eta,gamma1);
        rou_2e=Bisection_16_QAM_23_new(e,f,1e-5,eta,gamma2);
        
        %f_a=(1/(2*gamma1))*rou_1a+(1/(2*gamma2))*rou_2a-1;
        %f_b=(1/(2*gamma1))*rou_1b+(1/(2*gamma2))*rou_2b-1;
        f_e=(1/(2*gamma1))*rou_1e+(1/(2*gamma2))*rou_2e-1;
        
        if f_e==0   %Find the optimal point
            p1(n)=rou_1e/gamma1;
            p2(n)=rou_2e/gamma2;
            break
        elseif f_e<0
            b=eta;
        else
            a=eta;
        end
        if abs(f_e)<1e-3
            p1(n)=rou_1e/gamma1;
            p2(n)=rou_2e/gamma2;
            break
        end
    end
end
%----------------------Constallation Constrained WF------------------------
p=10.^(Gamma1_dB/10);  %16-QAM
Pt=2;  %Total power

g1=1;
g2=100;  %Channel gain

M=16;  %Constellation order
tol=1e-7;   %Tolerance

p_1=zeros(1,length(p));
p_2=zeros(1,length(p));

%-----------------------------Constellation WF-----------------------------
for n=1:length(p)
    a=0;   %Original lower bound
    b=g2*(1-1/M);   %Original upper bound
    max1=-1+ceil((log(b-a)-log(tol))/log(2));   %Number of iterations
    
    for k=1:1+max1
        lamda=(a+b)/2;   %Bisection
        %p_1a=constellation1(p(n),g1,M,a);   %Power allocation at a
        %p_2a=constellation1(p(n),g2,M,a);
        
        %p_1b=constellation1(p(n),g1,M,b);   %Power allocation at b
        %p_2b=constellation1(p(n),g2,M,b);
        
        p_1l=constellation1(p(n),g1,M,lamda);   %Power allocation at lamda
        p_2l=constellation1(p(n),g2,M,lamda);
        
        %fa=p_1a+p_2a-Pt;
        %fb=p_1b+p_2b-Pt;
        fl=p_1l+p_2l-Pt;
        
        if fl==0   %Find the optimal point
            p_1(n)=p_1l;
            p_2(n)=p_2l;
            break
        elseif fl<0
                b=lamda;
            else
                a=lamda;
        end
        if b-a<tol
            p_1(n)=p_1l;
            p_2(n)=p_2l;
            break
        end
    end
end
%--------------------------------Regular WF--------------------------------
tol=1e-5;   %Tolerance
for n=1:length(p)
    c=0;   %Original lower bound
    d=g2;   %Original upper bound
    max2=-1+ceil((log(d-c)-log(tol))/log(2));  %Number of iterations
    for k=1:max2+1
        lamda_new=(c+d)/2;   %Bisection
        %p_1c=regular_wf1(p(n),g1,c);   %Power allocation at c
        %p_2c=regular_wf1(p(n),g2,c);
        
        %p_1d=regular_wf1(p(n),g1,d);   %Power allocation at d
        %p_2d=regular_wf1(p(n),g2,d);
        
        p_1l_new=regular_wf1(p(n),g1,lamda_new);   %Power allocation at lamda
        p_2l_new=regular_wf1(p(n),g2,lamda_new);
        
        %fc=p_1c+p_2c-Pt;   %Total power is 2
        %fd=p_1d+p_2d-Pt;
        fl_new=p_1l_new+p_2l_new-Pt;
        
        if fl_new==0   %Find the optimal point
            p11(n)=p_1l_new;
            p21(n)=p_2l_new;
            break
        elseif fl_new<0
                d=lamda_new;
            else
                c=lamda_new;
        end
        if d-c<tol
            p11(n)=p_1l_new;
            p21(n)=p_2l_new;
            break
        end
    end
end

plot(Gamma1_dB,p1)
hold on
grid on
plot(Gamma1_dB,p_1,'--')
plot(Gamma1_dB,p11,'-.')

plot(Gamma1_dB,p2)
plot(Gamma1_dB,p_2,'--')
plot(Gamma1_dB,p21,'-.')

xlabel('P/dB')
ylabel('Power Allocation')
ylim([0 2])
legend('Exact mercury/waterfilling','AOPA','WF')