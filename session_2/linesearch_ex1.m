close all;
clearvars;
f=@(x) 0.5+2*(x-3)^2;
df=@(x) 4*(x-3);

%% Bracketing Algoritm
alpha0=0; alpha1=1; rho=0.25; sigma=0.5; 
num_iter=10;
alpha=zeros(1,num_iter);
a=zeros(1,num_iter);
b=zeros(1,num_iter);
alpha(1)=alpha0;
alpha(2)=alpha1;
fmin=0;
mu= (fmin-f(0))/(rho*df(0)); tau1=9;
for iter=2:num_iter
    val=f(alpha(iter)); %evaluate f(alpha_i)
    if val <= fmin
        fmin=val;
        disp('terminate line search');
    end
    if (val>f(0)+alpha(iter)*rho*df(0)) || ...
            (f(alpha(iter))>= f(alpha(iter-1)))
        a(iter)=alpha(iter-1);
        b(iter)=alpha(iter);
        disp('terminate bracket');
        break;
    end
    val2=df(alpha(iter)); %evaluate f(alpha_i)
    if abs(val2)<=-sigma*df(0)
        fmin=val;
        disp('terminate line search');
    end
    if val2>=0
        a(iter)=alpha(iter-1);
        b(iter)=alpha(iter);
        disp('terminate bracket');
        break;        
    end
    if mu <= (2*alpha(iter)-alpha(iter-1))
        alpha(iter+1)=mu;
    else
        a_inter=2*alpha(iter)-alpha(iter-1);
        b_inter=min(mu,alpha(iter)+tau1*(alpha(iter)-alpha(iter-1)));
        dfz=@(z) (b_inter-a_inter)*df(z);
        zmin=-dfz(a_inter)/(2*(f(b_inter)-f(a_inter)-dfz(a_inter)));
        alpha(iter+1)=a_inter+zmin*(b_inter-a_inter);
    end
end

%% sectioning algorithm
tau2=0.1; tau3=0.5;
for iter2=iter:num_iter
    %choose alpha_j
    a_inter=a(iter2)+tau2*(b(iter2)-a(iter2));
    b_inter=b(iter2)-tau3*(b(iter2)-a(iter2));
    dfz=@(z) (b_inter-a_inter)*df(z);
    zmin=-dfz(a_inter)/(2*(f(b_inter)-f(a_inter)-dfz(a_inter)));
    alpha(iter2)=a_inter+zmin*(b_inter-a_inter);
    %evaluate f(alpha_j)
    val1=f(alpha(iter2));
    if (val1>f(0)+rho*alpha(iter2)*df(0)) || ...
            (f(alpha(iter2))>=f(a(iter2)))
        a(iter2+1)=a(iter2);
        b(iter2+1)=alpha(iter2);
    else
        val2=df(alpha(iter2));
        if abs(val2)<=-sigma*df(0)
            disp('terminate line search')
            break;
        end
        a(iter2+1)=alpha(iter2);
        if (b(iter2)-a(iter2))*df(alpha(iter2))>=0
            b(iter2+1)=a(iter2);
        else
            b(iter2+1)=b(iter2);
        end
    end
end
