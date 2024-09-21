%% RC proportional control example 2024
function J=cost(gains,t,type)
    R=15e3;
    C=100e-9;
    dt=700e-6;
    alpha=dt/(R*C+dt);
    y=zeros(1,length(t));
    u=zeros(1,length(t));
    y2=zeros(1,length(t));
    error=zeros(1,length(t));
    setpoint=[0 1 2 3 0]; %Kp=0.5;
    m=1;
    for k=2:length(t)
        if mod(k,500)==0
            m=m+1;
        end
        error(k)=(setpoint(m)-y2(k-1));
        if k>3
            u(k)=u(k-1)+gains(1)*error(k)+gains(2)*error(k-1)+gains(3)*error(k-2);
        else
            u(k)=u(k-1)+gains(1)*error(k)+gains(2)*error(k-1);
        end
        if u(k)>3.3
            u(k)=3.3;
        end
        if u(k)<0
            u(k)=0;
        end
        %% DAC
        %u(k)=u(k)*4096/3.3;
        %% process
        y(k)=alpha*u(k)+(1-alpha)*y(k-1);
        y2(k)=alpha*y(k)+(1-alpha)*y2(k-1);
    
        if y2(k)>3.3
            y2(k)=3.3;
        end
        if y2(k)<0
            y2(k)=0;
        end
    end
    if type == "QR"
        J=norm(error)^2+norm(gradient(y2))^2;
    elseif type == "QR2"
        J=norm(error)^2+norm(gradient(u))^2;
    elseif type == "QR23"
        J=norm(error)^2+norm(gradient(u))^2+norm(gradient(y2))^2;
    else
        J=norm(error)^2;
    end
end
