function example2()

    % Example 2

    % Reference:
    %   Paper   = Optimal modeling of nonlinear systems: method of variable injections. (Submitted paper - 2023)
    %   Authors = Soto-Quiros, Pablo and Torokhti, Anatoli


    clc; clear; close all

    %%%%%%% Case 1: q1=25 and q2=500 %%%%%%%  
    fprintf('Case 1: q1=25 and q2=500\n')    
    m=100; 
    % Example T0(v0)
    q0=100; r0=50;
    s=6750;    
    X=rand(m,s);
    Aux1=rand(q0,s); N=Aux1-X; 
    Y=X+N;
    V0=Y;
    Z0=V0;  
    %Covariance matrices    
    Exz0=(1/s)*(X*Z0'); Ez0z0=(1/s)*(Z0*Z0');
    T0=Exz0*pinv(Ez0z0)*(Exz0)';          
    [U,~,~]=svd(T0);
    Ur=U(:,1:r0);
    G0=Ur;
    H0=Ur'*Exz0*pinv(Ez0z0);   
    %Compute error
    Exx=(1/s)*(X*X');
    Exw=Exz0;
    Eww=Ez0z0;
    S=G0*H0;   
    error0_c1=trace(Exx-Exw*S'-S*Exw'+S*Eww*S');
    fprintf(['1) MSE for T(v0) = ', num2str(error0_c1),'\n'])

    % Example T1(v0,v1)
    q0=100; r0=25;
    q1=25; r1=25;
    s=1500;
    X=rand(m,s);
    Aux1=rand(q0,s); N=Aux1-X; 
    Y=X+N;
    V0=Y; V1=rand(q1,s);
    Z0=V0;          
    Ev1z0=(1/s)*(V1*Z0'); Ez0z0=(1/s)*(Z0*Z0'); Z1=V1-Ev1z0*pinv(Ez0z0)*Z0;
    %Covariance matrices    
    Exz0=(1/s)*(X*Z0'); Ez0z0=(1/s)*(Z0*Z0');
    T0=Exz0*pinv(Ez0z0)*(Exz0)';                 
    [U0,~,~]=svd(T0);
    U0r=U0(:,1:r0);
    G0=U0r;
    H0=U0r'*Exz0*pinv(Ez0z0);   

    Exz1=(1/s)*(X*Z1'); Ez1z1=(1/s)*(Z1*Z1');
    T1=Exz1*pinv(Ez1z1)*(Exz1)';                 
    [U1,~,~]=svd(T1);
    U1r=U1(:,1:r1);
    G1=U1r;
    H1=U1r'*Exz1*pinv(Ez1z1);   

    %Compute error
    Exx=(1/s)*(X*X');
    Exw=[Exz0 Exz1];
    Eww=blkdiag(Ez0z0,Ez1z1);
    S=[G0*H0 G1*H1];   
    error1_c1=trace(Exx-Exw*S'-S*Exw'+S*Eww*S');
    fprintf(['2) MSE for T(v0,v1) = ', num2str(error1_c1),'\n'])

    % Example T2(v0,v1,v2)
    q0=100; r0=17;
    q1=25; r1=17;
    q2=500; r2=16;
    s=1300;
    X=rand(m,s);
    Aux1=rand(q0,s); N=Aux1-X; 
    Y=X+N;
    V0=Y; V1=rand(q1,s); V2=rand(q2,s);
    Z0=V0;          
    Ev1z0=(1/s)*(V1*Z0'); Ez0z0=(1/s)*(Z0*Z0'); Z1=V1-Ev1z0*pinv(Ez0z0)*Z0;
    Ev2z0=(1/s)*(V2*Z0'); Ev2z1=(1/s)*(V2*Z1'); Ez1z1=(1/s)*(Z1*Z1'); Z2=V2-Ev2z0*pinv(Ez0z0)*Z0-Ev2z1*pinv(Ez1z1)*Z1;
    %Covariance matrices    
    Exz0=(1/s)*(X*Z0'); Ez0z0=(1/s)*(Z0*Z0');
    T0=Exz0*pinv(Ez0z0)*(Exz0)';                 
    [U0,~,~]=svd(T0);
    U0r=U0(:,1:r0);
    G0=U0r;
    H0=U0r'*Exz0*pinv(Ez0z0);   

    Exz1=(1/s)*(X*Z1'); Ez1z1=(1/s)*(Z1*Z1');
    T1=Exz1*pinv(Ez1z1)*(Exz1)';                 
    [U1,~,~]=svd(T1);
    U1r=U1(:,1:r1);
    G1=U1r;
    H1=U1r'*Exz1*pinv(Ez1z1);   

    Exz2=(1/s)*(X*Z2'); Ez2z2=(1/s)*(Z2*Z2');
    T2=Exz2*pinv(Ez2z2)*(Exz2)';                 
    [U2,~,~]=svd(T2);
    U2r=U2(:,1:r2);
    G2=U2r;
    H2=U2r'*Exz2*pinv(Ez2z2);           

    %Compute error
    Exx=(1/s)*(X*X');
    Exw=[Exz0 Exz1 Exz2];
    Eww=blkdiag(blkdiag(Ez0z0,Ez1z1),Ez2z2);
    S=[G0*H0 G1*H1 G2*H2];   
    error2_c1=trace(Exx-Exw*S'-S*Exw'+S*Eww*S');
    fprintf(['3) MSE for T(v0,v1,v2) = ', num2str(error2_c1),'\n'])    

    
    %%%%%%% Case 2: q1=200 and q2=500 %%%%%%%  
    fprintf('\nCase 2: q1=200 and q2=500\n')  
    
    % Example T0(v0)
    q0=100; r0=50;
    s=6750;    
    X=rand(m,s);
    Aux1=rand(q0,s); N=Aux1-X; 
    Y=X+N;
    V0=Y;
    Z0=V0;  
    %Covariance matrices    
    Exz0=(1/s)*(X*Z0'); Ez0z0=(1/s)*(Z0*Z0');
    T0=Exz0*pinv(Ez0z0)*(Exz0)';          
    [U,~,~]=svd(T0);
    Ur=U(:,1:r0);
    G0=Ur;
    H0=Ur'*Exz0*pinv(Ez0z0);   
    %Compute error
    Exx=(1/s)*(X*X');
    Exw=Exz0;
    Eww=Ez0z0;
    S=G0*H0;   
    error0_c2=trace(Exx-Exw*S'-S*Exw'+S*Eww*S');
    fprintf(['1) MSE for T(v0) = ', num2str(error0_c2),'\n'])

    % Example T1(v0,v1)
    q0=100; r0=25;
    q1=200; r1=25;
    s=1800;
    X=rand(m,s);
    Aux1=rand(q0,s); N=Aux1-X; 
    Y=X+N;
    V0=Y; V1=rand(q1,s);
    Z0=V0;          
    Ev1z0=(1/s)*(V1*Z0'); Ez0z0=(1/s)*(Z0*Z0'); Z1=V1-Ev1z0*pinv(Ez0z0)*Z0;
    %Covariance matrices    
    Exz0=(1/s)*(X*Z0'); Ez0z0=(1/s)*(Z0*Z0');
    T0=Exz0*pinv(Ez0z0)*(Exz0)';                 
    [U0,~,~]=svd(T0);
    U0r=U0(:,1:r0);
    G0=U0r;
    H0=U0r'*Exz0*pinv(Ez0z0);   

    Exz1=(1/s)*(X*Z1'); Ez1z1=(1/s)*(Z1*Z1');
    T1=Exz1*pinv(Ez1z1)*(Exz1)';                 
    [U1,~,~]=svd(T1);
    U1r=U1(:,1:r1);
    G1=U1r;
    H1=U1r'*Exz1*pinv(Ez1z1);   

    %Compute error
    Exx=(1/s)*(X*X');
    Exw=[Exz0 Exz1];
    Eww=blkdiag(Ez0z0,Ez1z1);
    S=[G0*H0 G1*H1];   
    error1_c2=trace(Exx-Exw*S'-S*Exw'+S*Eww*S');
    fprintf(['2) MSE for T(v0,v1) = ', num2str(error1_c2),'\n'])
    
    % Example T2(v0,v1,v2)
    q0=100; r0=17;
    q1=200; r1=17;
    q2=500; r2=16;
    s=1050;
    X=rand(m,s);
    Aux1=rand(q0,s); N=Aux1-X; 
    Y=X+N;
    V0=Y; V1=rand(q1,s); V2=rand(q2,s);
    Z0=V0;          
    Ev1z0=(1/s)*(V1*Z0'); Ez0z0=(1/s)*(Z0*Z0'); Z1=V1-Ev1z0*pinv(Ez0z0)*Z0;
    Ev2z0=(1/s)*(V2*Z0'); Ev2z1=(1/s)*(V2*Z1'); Ez1z1=(1/s)*(Z1*Z1'); Z2=V2-Ev2z0*pinv(Ez0z0)*Z0-Ev2z1*pinv(Ez1z1)*Z1;
    %Covariance matrices    
    Exz0=(1/s)*(X*Z0'); Ez0z0=(1/s)*(Z0*Z0');
    T0=Exz0*pinv(Ez0z0)*(Exz0)';                 
    [U0,~,~]=svd(T0);
    U0r=U0(:,1:r0);
    G0=U0r;
    H0=U0r'*Exz0*pinv(Ez0z0);   

    Exz1=(1/s)*(X*Z1'); Ez1z1=(1/s)*(Z1*Z1');
    T1=Exz1*pinv(Ez1z1)*(Exz1)';                 
    [U1,~,~]=svd(T1);
    U1r=U1(:,1:r1);
    G1=U1r;
    H1=U1r'*Exz1*pinv(Ez1z1);   

    Exz2=(1/s)*(X*Z2'); Ez2z2=(1/s)*(Z2*Z2');
    T2=Exz2*pinv(Ez2z2)*(Exz2)';                 
    [U2,~,~]=svd(T2);
    U2r=U2(:,1:r2);
    G2=U2r;
    H2=U2r'*Exz2*pinv(Ez2z2);           

    %Compute error
    Exx=(1/s)*(X*X');
    Exw=[Exz0 Exz1 Exz2];
    Eww=blkdiag(blkdiag(Ez0z0,Ez1z1),Ez2z2);
    S=[G0*H0 G1*H1 G2*H2];   
    error2_c2=trace(Exx-Exw*S'-S*Exw'+S*Eww*S');
    fprintf(['3) MSE for T(v0,v1,v2) = ', num2str(error2_c2),'\n'])
    
    x=[1 2 3];
    y_c1=[error0_c1 error1_c1 error2_c1];
    y_c2=[error0_c2 error1_c2 error2_c2];
    
    hold on    
    plot(x,y_c2,'-b.','LineWidth',2,'MarkerSize',15)
    plot(x,y_c1,'-rx','LineWidth',2,'MarkerSize',8)
    xticks([1 2 3])
    xticklabels({'T_0(v_0)','T_1(v_0,v_1)','T_2(v_0,v_1,v_2)'})
    xlim([0.75,3.25])
    grid on
    ylabel('MSE')
    legend('q1=25, q2=500', 'q1=200, q2=500')
    
end

