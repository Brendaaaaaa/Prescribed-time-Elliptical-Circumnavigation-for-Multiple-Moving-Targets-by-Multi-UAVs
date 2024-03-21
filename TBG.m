function var = TBG(t,t_p)
global ep1;
a=10;
b=-24;
c=15;
if t<=t_p
    xi=a/(t_p^6)*t^6+b/(t_p^5)*t^5+c/(t_p^4)*t^4;
    dxi=6*a/(t_p^6)*t^5+5*b/(t_p^5)*t^4+4*c/(t_p^4)*t^3;
else
    xi=1;
    dxi=0;
end
var=dxi/(2*(1-xi+ep1));