close all; clearvars;
L=@(u1,u2) u1.^4+u1.*u2+(1+u2).^2;
g=@(u) [4*u(1)^3+u(2);u(1)+2*(1+u(2))];
G=@(u) [[12*u(1)^2,1];[1,2]];
num_iter=7;
uk=zeros(2,num_iter);
gk=zeros(2,num_iter);

figure
fsurf(L,[-2 2],'FaceAlpha',0.5,'ShowContours','on')
hold on

uk(:,1)=[.0;-1];
for k=1:num_iter-1
    gk(:,k)=g(uk(:,k));
    Gk=G(uk(:,k));
    deltak=-inv(Gk)*gk(:,k);
    uk(:,k+1)=uk(:,k)+deltak;
    plot3([uk(1,k) uk(1,k+1)],[uk(2,k) uk(2,k+1)],...
        [L(uk(1,k),uk(2,k)) L(uk(1,k+1),uk(2,k+1))],'k-o','LineWidth',2)
end

xlabel('$u_1$',Interpreter='latex')
ylabel('$u_2$',Interpreter='latex')
zlabel('$L$',Interpreter='latex')
