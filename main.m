%%
clc;
clear all;
%%
% params
global mo s beta_in gamma_in beta_it gamma_it omega_d ep1
mo=3;s=5;
beta_in=10;gamma_in=10;
beta_it=2;gamma_it=2;
omega_d=0.5*pi;
ep1=0.001;
i=1;dt=0.001;
t(i)=0;
t_p=4;
t_end=7;
% initial targets states
xt_1(i)=10;yt_1(i)=3;vt_1(i)=2;psit_1(i)=0;
xt_2(i)=-5;yt_2(i)=-1;vt_2(i)=2;psit_2(i)=0;
xt_3(i)=-6;yt_3(i)=5;vt_3(i)=2;psit_3(i)=0;
xt=[xt_1(i),xt_2(i),xt_3(i)];yt=[yt_1(i),yt_2(i),yt_3(i)];
[a(i),b(i),x0(i),y0(i),theta(i)]=EllipticSelect(xt,yt);
% initial UAV states
x_1(i)=-10;y_1(i)=-10;v_1(i)=20;psi_1(i)=0;
x_2(i)=-20;y_2(i)=6;v_2(i)=20;psi_2(i)=0.5;
x_3(i)=-15;y_3(i)=-10;v_3(i)=20;psi_3(i)=0.1;
x_4(i)=-5;y_4(i)=15;v_4(i)=20;psi_4(i)=-0.5;
x_5(i)=-12;y_5(i)=8;v_5(i)=20;psi_5(i)=-0.1;
x_6(i)=-15;y_6(i)=-15;v_6(i)=20;psi_6(i)=-0.3;

while t(i)<=t_end
%----------  TBG  ----------%
    kappa(i)=TBG(t(i),t_p);
%----------  Elliptic  ----------%
    xt=[xt_1(i),xt_2(i),xt_3(i)];yt=[yt_1(i),yt_2(i),yt_3(i)];
    [a(i),b(i),x0(i),y0(i),theta(i)]=EllipticSelect(xt,yt);
%----------  Targets  ----------%
    vt_1(i)=5;psit_1(i)=0.5*sin(t(i));
    vt_2(i)=5;psit_2(i)=0.5*cos(t(i));
    vt_3(i)=5;psit_3(i)=0;
    
    alpha_1(i)=atan2(y_1(i)-y0(i),x_1(i)-x0(i));
    alpha_2(i)=atan2(y_2(i)-y0(i),x_2(i)-x0(i));
    alpha_3(i)=atan2(y_3(i)-y0(i),x_3(i)-x0(i));
    alpha_4(i)=atan2(y_4(i)-y0(i),x_4(i)-x0(i));
    alpha_5(i)=atan2(y_5(i)-y0(i),x_5(i)-x0(i));
    alpha_6(i)=atan2(y_6(i)-y0(i),x_6(i)-x0(i));
%----------  UAV 1  ----------%
    r_1(i)=((x_1(i)-x0(i))^2+(y_1(i)-y0(i))^2)^0.5;
%     alpha_1(i)=psi_1(i);
    R_1(i)=(a(i)*b(i))/((b(i)^2*cos(alpha_1(i))^2+a(i)^2*sin(alpha_1(i))^2)^0.5);
    dR_1(i)=(-a(i)*b(i)*(a(i)^2-b(i)^2)*sin(alpha_1(i))*cos(alpha_1(i))*omega_d)/(((b(i)^2*cos(alpha_1(i))^2+a(i)^2*sin(alpha_1(i))^2)^0.5)^3);
    e_1n(i)=r_1(i)-R_1(i);
    [u_1n(i)]=CtrlNorm(kappa(i),e_1n(i),dR_1(i));
    e_1t(i)=alpha_1(i)-alpha_2(i)+pi/3;
    if e_1t(i)>pi
        e_1t(i)=e_1t(i)-2*pi;
    elseif e_1t(i)<-pi
        e_1t(i)=e_1t(i)+2*pi;
    else
        e_1t(i)=e_1t(i);
    end
    [omega_1(i),u_1t(i)]=CtrlTan(kappa(i),e_1t(i),r_1(i));
    u_1x(i)=u_1n(i)*cos(alpha_1(i))-u_1t(i)*sin(alpha_1(i));
    u_1y(i)=u_1n(i)*sin(alpha_1(i))+u_1t(i)*cos(alpha_1(i));
%----------  UAV 2  ----------%
    r_2(i)=((x_2(i)-x0(i))^2+(y_2(i)-y0(i))^2)^0.5;
    R_2(i)=(a(i)*b(i))/((b(i)^2*cos(alpha_2(i))^2+a(i)^2*sin(alpha_2(i))^2)^0.5);
    dR_2(i)=(-a(i)*b(i)*(a(i)^2-b(i)^2)*sin(alpha_2(i))*cos(alpha_2(i))*omega_d)/(((b(i)^2*cos(alpha_2(i))^2+a(i)^2*sin(alpha_2(i))^2)^0.5)^3);
    e_2n(i)=r_2(i)-R_2(i);
    [u_2n(i)]=CtrlNorm(kappa(i),e_2n(i),dR_2(i));
    e_2t(i)=alpha_2(i)-alpha_1(i)-pi/3;
    if e_2t(i)>pi
        e_2t(i)=e_2t(i)-2*pi;
    elseif e_2t(i)<-pi
        e_2t(i)=e_2t(i)+2*pi;
    else
        e_2t(i)=e_2t(i);
    end
    [omega_2(i),u_2t(i)]=CtrlTan(kappa(i),e_2t(i),r_2(i));
    u_2x(i)=u_2n(i)*cos(alpha_2(i))-u_2t(i)*sin(alpha_2(i));
    u_2y(i)=u_2n(i)*sin(alpha_2(i))+u_2t(i)*cos(alpha_2(i));
%----------  UAV 3  ----------%
    r_3(i)=((x_3(i)-x0(i))^2+(y_3(i)-y0(i))^2)^0.5;
    R_3(i)=(a(i)*b(i))/((b(i)^2*cos(alpha_3(i))^2+a(i)^2*sin(alpha_3(i))^2)^0.5);
    dR_3(i)=(-a(i)*b(i)*(a(i)^2-b(i)^2)*sin(alpha_3(i))*cos(alpha_3(i))*omega_d)/(((b(i)^2*cos(alpha_3(i))^2+a(i)^2*sin(alpha_3(i))^2)^0.5)^3);
    e_3n(i)=r_3(i)-R_3(i);
    [u_3n(i)]=CtrlNorm(kappa(i),e_3n(i),dR_3(i));
    e_3t(i)=alpha_3(i)-alpha_2(i)-pi/3;
    if e_3t(i)>pi
        e_3t(i)=e_3t(i)-2*pi;
    elseif e_3t(i)<-pi
        e_3t(i)=e_3t(i)+2*pi;
    else
        e_3t(i)=e_3t(i);
    end
    [omega_3(i),u_3t(i)]=CtrlTan(kappa(i),e_3t(i),r_3(i));
    u_3x(i)=u_3n(i)*cos(alpha_3(i))-u_3t(i)*sin(alpha_3(i));
    u_3y(i)=u_3n(i)*sin(alpha_3(i))+u_3t(i)*cos(alpha_3(i));
%----------  UAV 4  ----------%
    r_4(i)=((x_4(i)-x0(i))^2+(y_4(i)-y0(i))^2)^0.5;
    R_4(i)=(a(i)*b(i))/((b(i)^2*cos(alpha_4(i))^2+a(i)^2*sin(alpha_4(i))^2)^0.5);
    dR_4(i)=(-a(i)*b(i)*(a(i)^2-b(i)^2)*sin(alpha_4(i))*cos(alpha_4(i))*omega_d)/(((b(i)^2*cos(alpha_4(i))^2+a(i)^2*sin(alpha_4(i))^2)^0.5)^3);
    e_4n(i)=r_4(i)-R_4(i);
    [u_4n(i)]=CtrlNorm(kappa(i),e_4n(i),dR_4(i));
    e_4t(i)=alpha_4(i)-alpha_3(i)-pi/3;
    if e_4t(i)>pi
        e_4t(i)=e_4t(i)-2*pi;
    elseif e_4t(i)<-pi
        e_4t(i)=e_4t(i)+2*pi;
    else
        e_4t(i)=e_4t(i);
    end
    [omega_4(i),u_4t(i)]=CtrlTan(kappa(i),e_4t(i),r_4(i));
    u_4x(i)=u_4n(i)*cos(alpha_4(i))-u_4t(i)*sin(alpha_4(i));
    u_4y(i)=u_4n(i)*sin(alpha_4(i))+u_4t(i)*cos(alpha_4(i));
%----------  UAV 5  ----------%
    r_5(i)=((x_5(i)-x0(i))^2+(y_5(i)-y0(i))^2)^0.5;
    R_5(i)=(a(i)*b(i))/((b(i)^2*cos(alpha_5(i))^2+a(i)^2*sin(alpha_5(i))^2)^0.5);
    dR_5(i)=(-a(i)*b(i)*(a(i)^2-b(i)^2)*sin(alpha_5(i))*cos(alpha_5(i))*omega_d)/(((b(i)^2*cos(alpha_5(i))^2+a(i)^2*sin(alpha_5(i))^2)^0.5)^3);
    e_5n(i)=r_5(i)-R_5(i);
    [u_5n(i)]=CtrlNorm(kappa(i),e_5n(i),dR_5(i));
    e_5t(i)=alpha_5(i)-alpha_4(i)-pi/3;
    if e_5t(i)>pi
        e_5t(i)=e_5t(i)-2*pi;
    elseif e_5t(i)<-pi
        e_5t(i)=e_5t(i)+2*pi;
    else
        e_5t(i)=e_5t(i);
    end
    [omega_5(i),u_5t(i)]=CtrlTan(kappa(i),e_5t(i),r_5(i));
    u_5x(i)=u_5n(i)*cos(alpha_5(i))-u_5t(i)*sin(alpha_5(i));
    u_5y(i)=u_5n(i)*sin(alpha_5(i))+u_5t(i)*cos(alpha_5(i));
%----------  UAV 6  ----------%
    r_6(i)=((x_6(i)-x0(i))^2+(y_6(i)-y0(i))^2)^0.5;
    R_6(i)=(a(i)*b(i))/((b(i)^2*cos(alpha_6(i))^2+a(i)^2*sin(alpha_6(i))^2)^0.5);
    dR_6(i)=(-a(i)*b(i)*(a(i)^2-b(i)^2)*sin(alpha_6(i))*cos(alpha_6(i))*omega_d)/(((b(i)^2*cos(alpha_6(i))^2+a(i)^2*sin(alpha_6(i))^2)^0.5)^3);
    e_6n(i)=r_6(i)-R_6(i);
    [u_6n(i)]=CtrlNorm(kappa(i),e_6n(i),dR_6(i));
    e_6t(i)=alpha_6(i)-alpha_5(i)-pi/3;
    if e_6t(i)>pi
        e_6t(i)=e_6t(i)-2*pi;
    elseif e_6t(i)<-pi
        e_6t(i)=e_6t(i)+2*pi;
    else
        e_6t(i)=e_6t(i);
    end
    [omega_6(i),u_6t(i)]=CtrlTan(kappa(i),e_6t(i),r_6(i));
    u_6x(i)=u_6n(i)*cos(alpha_6(i))-u_6t(i)*sin(alpha_6(i));
    u_6y(i)=u_6n(i)*sin(alpha_6(i))+u_6t(i)*cos(alpha_6(i));
    
    i=i+1;
%     [v_1(i),psi_1(i),x_1(i),y_1(i)] = UAV_model(dt,u_1x(i-1),u_1y(i-1),omega_1(i-1),psi_1(i-1),x_1(i-1),y_1(i-1));
%     [v_2(i),psi_2(i),x_2(i),y_2(i)] = UAV_model(dt,u_2x(i-1),u_2y(i-1),omega_2(i-1),psi_2(i-1),x_2(i-1),y_2(i-1));
    x_1(i)=x_1(i-1)+u_1x(i-1)*dt;y_1(i)=y_1(i-1)+u_1y(i-1)*dt;
    x_2(i)=x_2(i-1)+u_2x(i-1)*dt;y_2(i)=y_2(i-1)+u_2y(i-1)*dt;
    x_3(i)=x_3(i-1)+u_3x(i-1)*dt;y_3(i)=y_3(i-1)+u_3y(i-1)*dt;
    x_4(i)=x_4(i-1)+u_4x(i-1)*dt;y_4(i)=y_4(i-1)+u_4y(i-1)*dt;
    x_5(i)=x_5(i-1)+u_5x(i-1)*dt;y_5(i)=y_5(i-1)+u_5y(i-1)*dt;
    x_6(i)=x_6(i-1)+u_6x(i-1)*dt;y_6(i)=y_6(i-1)+u_6y(i-1)*dt;
    psi_1(i)=psi_1(i-1)+omega_1(i-1)*dt;
    if psi_1(i)>pi
        psi_1(i)=psi_1(i)-2*pi;
    elseif psi_1(i)<-pi
        psi_1(i)=psi_1(i)+2*pi;
    else
        psi_1(i)=psi_1(i);
    end
    psi_2(i)=psi_2(i-1)+omega_2(i-1)*dt;
    if psi_2(i)>pi
        psi_2(i)=psi_2(i)-2*pi;
    elseif psi_2(i)<-pi
        psi_2(i)=psi_2(i)+2*pi;
    else
        psi_2(i)=psi_2(i);
    end
    psi_3(i)=psi_3(i-1)+omega_3(i-1)*dt;
    if psi_3(i)>pi
        psi_3(i)=psi_3(i)-2*pi;
    elseif psi_3(i)<-pi
        psi_3(i)=psi_3(i)+2*pi;
    else
        psi_3(i)=psi_3(i);
    end
    psi_4(i)=psi_4(i-1)+omega_4(i-1)*dt;
    if psi_4(i)>pi
        psi_4(i)=psi_4(i)-2*pi;
    elseif psi_4(i)<-pi
        psi_4(i)=psi_4(i)+2*pi;
    else
        psi_4(i)=psi_4(i);
    end
    psi_5(i)=psi_5(i-1)+omega_5(i-1)*dt;
    if psi_5(i)>pi
        psi_5(i)=psi_5(i)-2*pi;
    elseif psi_5(i)<-pi
        psi_5(i)=psi_5(i)+2*pi;
    else
        psi_5(i)=psi_5(i);
    end
    psi_6(i)=psi_6(i-1)+omega_6(i-1)*dt;
    if psi_6(i)>pi
        psi_6(i)=psi_6(i)-2*pi;
    elseif psi_6(i)<-pi
        psi_6(i)=psi_6(i)+2*pi;
    else
        psi_6(i)=psi_6(i);
    end
    xt_1(i)=xt_1(i-1)+vt_1(i-1)*cos(psit_1(i-1))*dt;yt_1(i)=yt_1(i-1)+vt_1(i-1)*sin(psit_1(i-1))*dt;
    xt_2(i)=xt_2(i-1)+vt_2(i-1)*cos(psit_2(i-1))*dt;yt_2(i)=yt_2(i-1)+vt_2(i-1)*sin(psit_2(i-1))*dt;
    xt_3(i)=xt_3(i-1)+vt_3(i-1)*cos(psit_3(i-1))*dt;yt_3(i)=yt_3(i-1)+vt_3(i-1)*sin(psit_3(i-1))*dt;
    
    t(i)=t(i-1)+dt;
    
end