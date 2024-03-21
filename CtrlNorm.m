function [u_in]=CtrlNorm(kappa,e_in,dR)
global beta_in gamma_in
u_in=-kappa*e_in-beta_in*abs(e_in)^(0.2)*sign(e_in)-gamma_in*e_in^3+dR;