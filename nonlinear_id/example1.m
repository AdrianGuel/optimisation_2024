% nonlinear parameter estimation algorithm using 
% fmincon from MATLAB

close all;
mu=1;
[t,yr] = ode45(@(t,y) f(t,y,mu),[0 20],[2; 0]);
lb = [0,0,0];
ub = [2,5,5];
A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = [];
x0 = (lb + ub)/2;
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point');
xsol = fmincon(@(x) J(t,yr,x),x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

[t,ye] = ode45(@(t,y) f(t,y,xsol(1)),t,[xsol(2); xsol(3)]);
plot(t,yr(:,1),'k-',t,ye(:,1),'b--')
hold on
plot(t,yr(:,2),'r-',t,ye(:,2),'g--')
title('Solution of van der Pol Equation (\mu = 1) with ODE45');
xlabel('Time t');
ylabel('Solution y');
legend('true','estimated')

function dydt = f(t,y,mu)
    %the van der Pol ODEs
    dydt = [y(2); mu*(1-y(1)^2)*y(2)-y(1)];
end

function cost=J(t,yr,theta)
    [~,ye] = ode45(@(t,y) f(t,y,theta(1)),t,[theta(2); theta(3)]);
    cost=norm(yr-ye)^2;
end
