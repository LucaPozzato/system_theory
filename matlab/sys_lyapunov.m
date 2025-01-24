function f=proj_sys(t,X)
    % Constants
    L = 11;
    C = 11;
    a = 0.04;
    b = 0.18;
    R0 = -0.64; % Nominal resistance value
    omega = 0.02; % Initial guess for omega
    
    % Time varying resistance
    R_t = R0 * (1 + 0.25 * sin(omega * t));

    x=X(1); y=X(2);
    
    Y= [X(3), X(5);
        X(4), X(6)];
    
    f=zeros(4,1);
    
    % System equations
    f(1) = -(R_t * x) / L + y / L;
    f(2) = -x / C + (a * y) / C - (b * y^3) / C;
    
    %Linearized system
    
     Jac=[-R_t / L, 1 / L;
         -1 / C, a / C - (3 * b * y^2) / C];
      
    %Variational equation   
    f(3:6)=Jac*Y;
end

%Output data must be a column vector

[T,Res]=lyapunov(2,@proj_sys,@ode15s,0,0.5,50000,[1, 1],10000);
plot(T,Res);
title('Dynamics of Lyapunov exponents');
xlabel('Time'); ylabel('Lyapunov exponents');