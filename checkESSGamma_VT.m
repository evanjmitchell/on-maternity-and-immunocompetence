function check_gamma = checkESSGamma_VT(css, res, b_max, c_f, c_m, v, ...
    beta_ff, beta_fm, beta_mf, beta_mm, mu)
% A script to numerically evaluate the ESS condition on gamma_f and gamma_m, and return 
% either a 0 if it fails to satisfy the condition or a 1 otherwise. css is 
% a four-dimensional vector containing the candidate CSS values of 
% gamma_f, gamma_m, alpha_f, and alpha_m (in that order); res is a 
% four-dimensional vector giving the equilibrium values of S_f, S_m, I_f, 
% and I_m (in that order).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract values from css
gamma_f = css(1);
gamma_m = css(2);
alpha_f = css(3);
alpha_m = css(4);
% numerically approximate the Hessian matrix using a second-order
% centred difference formula
h = 0.01;
f1 = fitnessGamma_VT(res, b_max, c_f, c_m, v, beta_ff, beta_fm, beta_mf, beta_mm, mu, alpha_f, alpha_m, gamma_f, gamma_m, gamma_f, gamma_m);
f2 = fitnessGamma_VT(res, b_max, c_f, c_m, v, beta_ff, beta_fm, beta_mf, beta_mm, mu, alpha_f, alpha_m, gamma_f, gamma_m, gamma_f + h, gamma_m);
f3 = fitnessGamma_VT(res, b_max, c_f, c_m, v, beta_ff, beta_fm, beta_mf, beta_mm, mu, alpha_f, alpha_m, gamma_f, gamma_m, gamma_f - h, gamma_m);
f4 = fitnessGamma_VT(res, b_max, c_f, c_m, v, beta_ff, beta_fm, beta_mf, beta_mm, mu, alpha_f, alpha_m, gamma_f, gamma_m, gamma_f, gamma_m + h);
f5 = fitnessGamma_VT(res, b_max, c_f, c_m, v, beta_ff, beta_fm, beta_mf, beta_mm, mu, alpha_f, alpha_m, gamma_f, gamma_m, gamma_f, gamma_m - h);
f6 = fitnessGamma_VT(res, b_max, c_f, c_m, v, beta_ff, beta_fm, beta_mf, beta_mm, mu, alpha_f, alpha_m, gamma_f, gamma_m, gamma_f + h, gamma_m + h);
f7 = fitnessGamma_VT(res, b_max, c_f, c_m, v, beta_ff, beta_fm, beta_mf, beta_mm, mu, alpha_f, alpha_m, gamma_f, gamma_m, gamma_f + h, gamma_m - h);
f8 = fitnessGamma_VT(res, b_max, c_f, c_m, v, beta_ff, beta_fm, beta_mf, beta_mm, mu, alpha_f, alpha_m, gamma_f, gamma_m, gamma_f - h, gamma_m + h);
f9 = fitnessGamma_VT(res, b_max, c_f, c_m, v, beta_ff, beta_fm, beta_mf, beta_mm, mu, alpha_f, alpha_m, gamma_f, gamma_m, gamma_f - h, gamma_m - h);
H11 = (f2 - 2*f1 + f3)/h^2;
H12 = (f6 - f7 - f8 + f9)/(4*h^2);
H21 = H12;
H22 = (f4 - 2*f1 + f5)/h^2;
% evaluate the trace of the Hessian
Htr = H11 + H22;
% evaluate the determinant of the Hessian
Hdt = H11*H22 - H12*H21;
% check that the trace is negative and the determinant is positive, and
% return either a 0 (fail) or 1 (pass), as well as the values of the trace
% and determinant
if ((Htr < 0) && (Hdt > 0))
    check_gamma = [1, Htr, Hdt];
else
    check_gamma = [0, Htr, Hdt];
end