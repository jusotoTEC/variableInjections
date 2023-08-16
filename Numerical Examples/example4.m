function example4()

    % Example 4
    
    % Reference:
    %   Paper   = Optimal modeling of nonlinear systems: method of variable injections. (Submitted paper - 2023)
    %   Authors = Soto-Quiros, Pablo and Torokhti, Anatoli

    clc; clear; close all

    m=50; 

    % Error T0(v0)
    q0=50; 
    s=10000; alpha=1.65;
    X=alpha*rand(m,s);
    V0=rand(q0,s);
    Z0=V0;          
    Exx=(1/s)*(X*X');
    Exz0=(1/s)*(X*Z0'); Ez0z0=(1/s)*(Z0*Z0');
    A0=Exz0*sqrtm(pinv(Ez0z0));                  
    error0=trace(Exx)-norm(A0,'fro')^2;
    fprintf(['MSE for T(v0) = ', num2str(error0),'\n'])

    % Example T1(v0,v1)
    q0=50; q1=50; 
    s=200; alpha=2.1;
    X=alpha*rand(m,s);
    V0=rand(q0,s); V1=rand(q1,s);
    Z0=V0;          
    Ev1z0=(1/s)*(V1*Z0'); Ez0z0=(1/s)*(Z0*Z0'); Z1=V1-Ev1z0*pinv(Ez0z0)*Z0;
    %Covariance matrices    
    Exz0=(1/s)*(X*Z0'); Ez0z0=(1/s)*(Z0*Z0');
    A0=Exz0*sqrtm(pinv(Ez0z0));                  
    Exz1=(1/s)*(X*Z1'); Ez1z1=(1/s)*(Z1*Z1');
    A1=Exz1*sqrtm(pinv(Ez1z1));

    %Compute error
    Exx=(1/s)*(X*X');

    error1=trace(Exx)-norm(A0,'fro')^2-norm(A1,'fro')^2;
    fprintf(['MSE for T(v0,v1) = ', num2str(error1),'\n'])

    % Example T2(v0,v1,v2)
    q0=50; q1=50; q2=50; 
    s=200; alpha=2.75;
    X=alpha*rand(m,s);
    V0=rand(q0,s); V1=rand(q1,s); V2=rand(q2,s);
    Z0=V0;          
    Ev1z0=(1/s)*(V1*Z0'); Ez0z0=(1/s)*(Z0*Z0'); Z1=V1-Ev1z0*pinv(Ez0z0)*Z0;
    Ev2z0=(1/s)*(V2*Z0'); Ev2z1=(1/s)*(V2*Z1'); Ez1z1=(1/s)*(Z1*Z1'); Z2=V2-Ev2z0*pinv(Ez0z0)*Z0-Ev2z1*pinv(Ez1z1)*Z1;
    %Covariance matrices    
    Exz0=(1/s)*(X*Z0'); Ez0z0=(1/s)*(Z0*Z0');
    A0=Exz0*sqrtm(pinv(Ez0z0));         
    Exz1=(1/s)*(X*Z1'); Ez1z1=(1/s)*(Z1*Z1');
    A1=Exz1*sqrtm(pinv(Ez1z1));                 
    Exz2=(1/s)*(X*Z2'); Ez2z2=(1/s)*(Z2*Z2');
    A2=Exz2*sqrtm(pinv(Ez2z2));
    %Compute error
    Exx=(1/s)*(X*X');
    error2=trace(Exx)-norm(A0,'fro')^2-norm(A1,'fro')^2-norm(A2,'fro')^2;
    fprintf(['MSE for T(v0,v1,v2) = ', num2str(error2),'\n'])
    
    x=[1 2 3];
    y=[error0 error1 error2];
    
    p=plot(x,y,'-bx','LineWidth',2,'MarkerSize',20);
    grid on
    xticks([1 2 3])
    xticklabels({'T_0(v_0)','T_1(v_0,v_1)','T_2(v_0,v_1,v_2)'})
    xlim([0.75,3.25])
    ylim([7,12])
    ylabel('MSE')
    p.MarkerEdgeColor = [1 0 0];
      
end

