function example3_Sigma_vrs_MSE()

    % Example 3: Sigma vrs MSE

    % Reference:
    %   Paper   = Optimal modeling of nonlinear systems: method of variable injections. (Submitted paper - 2023)
    %   Authors = Soto-Quiros, Pablo and Torokhti, Anatoli

    clc; clear; close all

    er0=[]; er1=[]; er2=[];
    vectSigma=[0.01 0.05 0.1 0.25:0.25:15];
    for sigma=vectSigma
       % Compute Covariance matrices
        m1=10; n=m1;        
        phase=exp(1i*rand(m1,1));
        Ac=diag(phase);
        m=2*m1;  r=floor(m/2); 
        b1=0.4; 
        b=b1*(m+n-1)/(m-1); 
        Exwxw=cov_new(m+n,b);
        A=[real(Ac) -imag(Ac); imag(Ac) real(Ac)];
        Exx=Exwxw(1:m,1:m);
        k=n;    
        Exy=Exx*A'; Eyx=Exy';
        Exw=Exwxw(1:m,m+1:m+k); Ewx=Exw';
        Eww=Exwxw(m+1:m+k,m+1:m+k);
        Ewy=Ewx*A'; Eyw=Ewy';
        Eyy=A*Exx*A'+sigma^2*eye(m);     
        Exp=[Exy Exw]; Epx=Exp';
        Epp=[Eyy Eyw; Ewy Eww];


        % Compute error GBT2
        T0=Exp*pinv(Epp)*Epx;
        vect0=svd(T0);
        error0=trace(Exx)-sum(vect0(1:r));   
        er0=[er0 error0];

        % Compute error T0(v0)
        T1=Exy*pinv(Eyy)*Eyx;
        vect1=svd(T1);
        error1=trace(Exx)-sum(vect1(1:r));
        er1=[er1 error1];

        % Compute error of T1(v0,v1)
        n=1000;
        phase=exp(1i*rand(m1,1));
        Ac=diag(phase);
        m=2*m1;  r1=floor(m/4); r2=floor(m/4);   
        b=b1*(m+n-1)/(m-1); 
        Exwxw=cov_new(m+n,b);
        A=[real(Ac) -imag(Ac); imag(Ac) real(Ac)];
        Exx=Exwxw(1:m,1:m);
        k=n;    
        Exy=Exx*A';
        Exw=Exwxw(1:m,m+1:m+k); Ewx=Exw';
        Eww=Exwxw(m+1:m+k,m+1:m+k);
        Ewy=Ewx*A'; Eyw=Ewy';
        Eyy=A*Exx*A'+sigma^2*eye(m);
        Exv=Exw-Exy*pinv(Eyy)*Eyw; Evx=Exv';
        Evv=Eww-Ewy*pinv(Eyy)*Eyw;       
        T2=Exv*pinv(Evv)*Evx;    
        vect2=svd(T2);
        error2=trace(Exx)-sum(vect1(1:r1))-sum(vect2(1:r2));
        er2=[er2 error2];
    end
    
    vectGBT2=er0;
    vectTv0=er1;
    vectTv1=er2;
    
    figure
    hold on
    plot(vectSigma,vectTv0,'-rx','LineWidth',2,'MarkerSize',8) 
    plot(vectSigma,vectGBT2,'-g+','LineWidth',2,'MarkerSize',8)
    plot(vectSigma,vectTv1,'-b.','LineWidth',2,'MarkerSize',15)           
    xlabel('Sigma', 'FontSize',14)
    ylabel('MSE', 'FontSize',14)
    grid on
    set(gca,'fontsize',14)
    legend('T_0(v_0)','GBT2','T_1(v_0,v_1)')
    
end

function C=cov_new(K,bi)
    C=zeros(K);
    for i=1:K
        for j=1:K
            if abs(i-j)<=(K-1)/bi
                C(i,j)=1-(abs(i-j)*bi)/(K-1);
            end
        end
    end    
end

