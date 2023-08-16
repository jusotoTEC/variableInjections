function example3_Dimension_vrs_MSE()

    % Example 3: Dimension vrs MSE

    % Reference:
    %   Paper   = Optimal modeling of nonlinear systems: method of variable injections. (Submitted paper - 2023)
    %   Authors = Soto-Quiros, Pablo and Torokhti, Anatoli

    clc; clear; close all

    % Compute Covariance matrices
    m1=10; n=10;
    phase=exp(1i*rand(m1,1));
    Ac=diag(phase);
    m=2*m1;  r1=floor(m/2); r2=floor(m/2); r=r1+r2;
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
    sigma=1;
    Eyy=A*Exx*A'+sigma^2*eye(m);     
    Exp=[Exy Exw]; Epx=Exp';
    Epp=[Eyy Eyw; Ewy Eww];
    
    
    % Compute error GBT2
    T0=Exp*pinv(Epp)*Epx;
    vect0=svd(T0);
    error0=trace(Exx)-sum(vect0(1:r));    
    
    
    % Compute error T0(v0)
    T1=Exy*pinv(Eyy)*Eyx;
    vect1=svd(T1);
    error1=trace(Exx)-sum(vect1(1:r));
    
    
    % Compute error of T1(v0,v1)
    er2=[];
    vectDimension=20:100;
    for q=vectDimension
        m1=10; 
        phase=exp(1i*rand(m1,1));
        Ac=diag(phase);
        m=2*m1;  r1=floor(m/2); r2=floor(m/2); 
        b1=0.4; 
        b=b1*(m+q-1)/(m-1);    
        Exwxw=cov_new(m+q,b);
        A=[real(Ac) -imag(Ac); imag(Ac) real(Ac)];
        Exx=Exwxw(1:m,1:m);
        k=q;    
        Exy=Exx*A';
        Exw=Exwxw(1:m,m+1:m+k); Ewx=Exw';
        Eww=Exwxw(m+1:m+k,m+1:m+k);
        Ewy=Ewx*A'; Eyw=Ewy';
        sigma=1;
        Eyy=A*Exx*A'+sigma^2*eye(m);
        Exv=Exw-Exy*pinv(Eyy)*Eyw; Evx=Exv';
        Evv=Eww-Ewy*pinv(Eyy)*Eyw;       

        T2=Exv*pinv(Evv)*Evx;    
        vect2=svd(T2);
        error2=trace(Exx)-sum(vect1(1:r1))-sum(vect2(1:r2));
        er2=[er2 error2];
    end
    
    vectGBT2=error0*ones(length(vectDimension),1);
    vectTv0=error1*ones(length(vectDimension),1);
    vectTv1=er2;
    
    figure
    hold on
    plot(vectDimension,vectTv0,'-rx','LineWidth',2,'MarkerSize',8)       
    plot(vectDimension,vectGBT2,'-g+','LineWidth',2,'MarkerSize',8)
    plot(vectDimension,vectTv1,'-b.','LineWidth',2,'MarkerSize',15)     
    xlabel('Dimension (q)', 'FontSize',14)
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

    %C=C+C';
    
end

