close all;
%% Example preliminaries
u=-3:0.1:3;
L=@(u) (u-1).^2;
dL=@(u) 2*(u-1);
figure
plot(u,L(u))
hold on
plot(u,dL(u))
hold on
plot(1,L(1),'p')
legend('L(u)','dL/du','minimum')

%% Example 1
f =@(x1,x2) x1.^2+x2.^2;
df =@(x1,x2) 2*x1+2*x2;
figure
fsurf(f,[-30 30],'FaceAlpha',0.5,'ShowContours','on')
hold on
fsurf(df,[-30 30],'FaceAlpha',0.05,'ShowContours','on')

f =@(x1,x2) -x1.^2-x2.^2;
df =@(x1,x2) -2*x1-2*x2;
figure
fsurf(f,[-30 30],'FaceAlpha',0.5,'ShowContours','on')
hold on
fsurf(df,[-30 30],'FaceAlpha',0.05,'ShowContours','on')

f =@(x1,x2) x1.^2-x2.^2;
df =@(x1,x2) 2*x1-2*x2;
figure
fsurf(f,[-30 30],'FaceAlpha',0.5,'ShowContours','on')
hold on
fsurf(df,[-30 30],'FaceAlpha',0.05,'ShowContours','on')

f =@(x1,x2) -x1.^2+x2.^2;
df =@(x1,x2) -2*x1+2*x2;
figure
fsurf(f,[-30 30],'FaceAlpha',0.5,'ShowContours','on')
hold on
fsurf(df,[-30 30],'FaceAlpha',0.05,'ShowContours','on')
