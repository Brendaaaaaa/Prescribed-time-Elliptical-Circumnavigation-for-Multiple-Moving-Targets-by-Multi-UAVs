function [a,b,x0,y0,theta]=EllipticSelect(xt,yt)
global mo s
% mo: the number of the polygon convex
% xt,yt: position matrix, e.g., x=[x_1,...,x_mo]
% s: safe distance

if mo==2
    c=0.5*((xt(1)-xt(2))^2+(yt(1)-yt(2))^2)^0.5; 
    b=(s^2+c*s+(s^4+2*c^2*s^2+2*c*s^3)^0.5)^0.5;
    a=(b^2+c^2)^0.5;
    x0=0.5*(xt(1)+xt(2));
    y0=0.5*(yt(1)+yt(2));
    theta=atan2(yt(1)-yt(2),xt(1)-xt(2));
end

if mo==3
    x0=1/3*(xt(1)+xt(2)+xt(3));
    y0=1/3*(yt(1)+yt(2)+yt(3));
    for j=1:1:mo
        th(j)=atan2(yt(j)-y0,xt(j)-x0);
        xt(j)=xt(j)+s*cos(th(j));
        yt(j)=yt(j)+s*sin(th(j));
    end
    l1=((xt(1)-xt(2))^2+(yt(1)-yt(2))^2)^0.5;
    l2=((xt(2)-xt(3))^2+(yt(2)-yt(3))^2)^0.5;
    l3=((xt(3)-xt(1))^2+(yt(3)-yt(1))^2)^0.5;
    a=1/3*(l1^2+l2^2+l3^2+2*(l1^4+l2^4+l3^4-l1^2*l2^2-l2^2*l3^2-l3^2*l1^2)^0.5)^0.5;
    b=1/3*(l1^2+l2^2+l3^2-2*(l1^4+l2^4+l3^4-l1^2*l2^2-l2^2*l3^2-l3^2*l1^2)^0.5)^0.5;
    at=[atan2(yt(1)-yt(2),xt(1)-xt(2)),atan2(yt(2)-yt(3),xt(2)-xt(3)),atan2(yt(3)-yt(1),xt(3)-xt(1))];
    theta=min(at);
end

% if mo==4
%     j=1;theta(j)=-pi;dtheta=0.01*pi;
%     
%     while theta(j)<pi
%         AA=[]
%         j=j+1;
%         theta(j)=theta(j-1)+dtheta;
%         
%     end
%     
% end