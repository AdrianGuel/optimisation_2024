% PDI control design for complex systems
% fmincon from MATLAB

close all;
clearvars;
%% control initial parameters
T=0.1; tf=30;
t=0:T:tf;
kp=0.001;kd=0.005;ki=0;
q0=kp+kd/T+ki*T;
q1=-kp-2*kd/T;
q2=kd/T;

%% optimiser
lb = [0,-7e1,0];
ub = [6e2,0,6e2];
A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = [];
x0 = [q0,q1,q2];
yr=1;
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point');
xsol = fmincon(@(x) J(t,yr,x),x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
[yc,u]=closed_loop(t,yr,xsol);
yyaxis left
yline(yr,'k-')
hold on
plot(t,yc(1,:),'b--')
hold on
plot(t,yc(2,:),'g--')
xlabel('Time t');
ylabel('Solution y');
legend('goal','x','dx')
yyaxis right
plot(t,u)
xlabel('Time t');
ylabel('Solution u');
title('PID controlled system');


function dydt = f(y,u)
    dydt = [y(2);-0.01*y(1)-y(2)+u];
end

function cost=J(t,yr,gain)
yc=zeros(2,length(t));
u=zeros(1,length(t));
uk=0; uk1=0;
ek1=0; ek2=0;
T=diff(t(1:2));
    for k=2:length(t)
        %% Euler
        %yc(:,k)=yc(:,k-1)+T*f(yc(:,k-1),uk);
        %% 4th order RK
        k1=f(yc(:,k-1),uk);
        k2=f(yc(:,k-1)+T*k1/2,uk);
        k3=f(yc(:,k-1)+T*k2/2,uk);
        k4=f(yc(:,k-1)+T*k3,uk);
        yc(:,k)=yc(:,k-1)+(T/6)*(k1+2*k2+2*k3+k4);
        ek=yr-yc(1,k); %error
	    uk=uk1+gain(1)*ek+gain(2)*ek1+gain(3)*ek2; % Recursive PID
        u(k)=uk;
	    uk1=uk; 
	    ek2=ek1;
	    ek1=ek;        
    end
    cost=norm(yr-yc(1,:))^2+0.08*norm(u)^2;
end

function [yc,u]=closed_loop(t,yr,gain)
yc=zeros(2,length(t));
u=zeros(1,length(t));
uk=0; uk1=0;
ek1=0; ek2=0;
T=diff(t(1:2));
    for k=2:length(t)
        yc(:,k)=yc(:,k-1)+T*f(yc(:,k-1),uk);
        ek=yr-yc(1,k); %error
	    uk=uk1+gain(1)*ek+gain(2)*ek1+gain(3)*ek2; % Recursive PID
        u(k)=uk;
	    uk1=uk; 
	    ek2=ek1;
	    ek1=ek;        
    end
end
