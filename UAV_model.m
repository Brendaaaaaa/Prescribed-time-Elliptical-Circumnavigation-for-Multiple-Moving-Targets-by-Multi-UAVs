function [v,psi,x,y] = UAV_model(dt,u_x,u_y,omega,psi,x,y)
v=(u_x^2+u_y^2)^0.5;
psi=omega*dt+psi;
dx=v*cos(psi);
dy=v*sin(psi);
x=x+dx*dt;
y=y+dy*dt;

% 速度约束
v_max=44;v_min=10;
if v>v_max
    v=v_max;
end
if v<v_min
    v=v_min;
end

% 偏航转化
if psi>pi
    psi=psi-2*pi;
end
if psi<-pi
    psi=psi+2*pi;
end

if psi>pi/2
    psi=pi/2;
end
if psi<-pi/2
    psi=-pi/2;
end

end