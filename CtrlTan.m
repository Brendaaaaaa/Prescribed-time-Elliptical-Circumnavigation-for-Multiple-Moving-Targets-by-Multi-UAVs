function [omega,u_it]=CtrlTan(kappa,e_it,r)
global beta_it gamma_it omega_d
omega=-kappa*e_it-beta_it*abs(e_it)^(0.2)*sign(e_it)-gamma_it*e_it^3+omega_d;
u_it=omega*r;