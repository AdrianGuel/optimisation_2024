%% constrained optimisation example
% bar{L}= \frac{1}{2}(u_1^2+u_2^2)-2(u_1+u_2)+\lambda(u_1+u_2-1)
% x = u1, u = u2
% params and functions
close all;
clearvars;

gamma= 1;
dcdx= 1;
dcdu= 1;
c_minus_y= @(x,u) x+u-gamma;
dLdx= @(x,u) x-2;
dLdu= @(x,u) u-2;
lambda=@(x,u) -dLdx(x,u)/dcdx;
dbarLdu=@(x,u) dLdu(x,u)+lambda(x,u)*dcdu;

% numerical method params
num_iter=30;
epsilon1=1e-3;
epsilon2=1e-3;

%% numerical method
figure
f=@(x,u) 0.5*x.^2+0.5*u.^2-2*x-2*u;%goal;
fsurf(f,[-5 5],'FaceAlpha',0.2,'ShowContours','on')
view(0,90)
hold on;
x=zeros(1,num_iter);
u=zeros(1,num_iter);
% Guess x and u
x(1)=4; %u_1
u(1)=5; %u_2
Luu=2;
for k=1:num_iter-1
    if norm(c_minus_y(x(k),u(k)))<epsilon1
        lam=lambda(x(k),u(k));
        if norm(dbarLdu(x(k),u(k)))<epsilon2
            disp("terminate search");
            break;
        else
            u(k+1)=u(k)-dbarLdu(x(k),u(k))/Luu;
        end
    else
        x(k+1)=x(k)-c_minus_y(x(k),u(k))/dcdx;
    end
    plot([x(k+1),x(k)],[u(k+1),u(k)],'k-o','LineWidth',3)
    hold on
    pause(0.5)
end
