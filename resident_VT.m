function dres_vt = resident_VT(res, b, v, mu, beta_ff, beta_fm, ...
    beta_mf, beta_mm, gamma_f, gamma_m, alpha_f, alpha_m)
% A script to compute the resident system of ODEs at a given state; order
% of the variables in res is S_f, S_m, I_f, I_m.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define total population size
N = sum(res);
% setup empty vector to store values of ODEs
dres_vt = zeros(1, 4);
% evaluate ODEs
dres_vt(1) = b*(res(1) + (1 - v)*res(3))*(res(2) + res(4))/N + gamma_f*res(3) ...
    - res(1)*beta_ff*res(3) - res(1)*beta_fm*res(4) - mu*res(1);
dres_vt(2) = b*(res(1) + (1 - v)*res(3))*(res(2) + res(4))/N + gamma_m*res(4) ...
    - res(2)*beta_mf*res(3) - res(2)*beta_mm*res(4) - mu*res(2);
dres_vt(3) = b*v*res(3)*(res(2) + res(4))/N + res(1)*beta_ff*res(3) + res(1)*beta_fm*res(4) ...
    - (mu + alpha_f + gamma_f)*res(3);
dres_vt(4) = b*v*res(3)*(res(2) + res(4))/N + res(2)*beta_mf*res(3) + res(2)*beta_mm*res(4) ...
    - (mu + alpha_m + gamma_m)*res(4);