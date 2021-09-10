function css = findCSS_VT(init_vals, b_max, c_f, c_m, v, mu, beta_max, ...
    d, stepEE, tolEE, maxitEE, stepCSS, tolCSS)
% A script to numerically approximate the CSS values of gamma_f, gamma_m, 
% alpha_f, and alpha_m. init_vals is a four-element vector giving the starting 
% point for the mutant host trajectory (the order is gamma_f, gamma_m). The 
% function returns a 14-element vector with elements (in order): 
% CSS values of gamma_f, gamma_m, alpha_f, alpha_m, 0/1 (fail/pass)
% for ESS check on gamma, values of the trace and determinant used
% to check the ESS condition for gamma, 0/1 (fail/pass) for ESS check on
% alpha, values of the trace and determinant used to check the ESS 
% condition for alpha, equilibrium values S_f, S_m, I_f, I_m of resident
% system at CSS.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract values of gamma_f, gamma_m, alpha_f, and alpha_m
gamma_f = init_vals(1);
gamma_m = init_vals(2);
alpha_f = init_vals(3);
alpha_m = init_vals(4);
if gamma_f < 0
    gamma_f = 1e-5;
end
if gamma_m < 0
    gamma_m = 1e-5;
end
if alpha_f < 0
    alpha_f = 1e-5;
end
if alpha_m < 0
    alpha_m = 1e-5;
end
% compute the value of b
b = b_max - (c_f*gamma_f + c_m*gamma_m);
% end the process if the birth rate becomes too small
if b < 2*mu
    css = repmat(2000, 1, 14);
    return;
end
% compute the values of beta_ij
beta_ff = beta_max*alpha_f/(alpha_f + d);
beta_fm = beta_max*alpha_m/(alpha_m + d);
beta_mf = beta_max*alpha_f/(alpha_f + d);
beta_mm = beta_max*alpha_m/(alpha_m + d);
% approximate the equilibrium values of the resident system
res = findEE_VT(b, v, mu, beta_ff, beta_fm, beta_mf, beta_mm, ...
    gamma_f, gamma_m, alpha_f, alpha_m, stepEE, tolEE, maxitEE);
% if EE couldn't be found, return corresponding initial condition for
% closer inspection
if isnan(res)
    css = [init_vals(1), init_vals(2), init_vals(3), init_vals(4), ...
        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN];
    return;
end
% if system broke down when finding EE, end the process
if sum(res) == 4000
    css = repmat(1000, 1, 14);
    return;
end
% evaluate selection gradients
sgrad_gamma_f = sgradGammaF_VT(res, b_max, c_f, c_m, v, mu, beta_ff, ...
    beta_fm, beta_mf, beta_mm, gamma_f, gamma_m, alpha_f, alpha_m);
sgrad_gamma_m = sgradGammaM_VT(res, b_max, c_f, c_m, v, mu, beta_ff, ...
    beta_fm, beta_mf, beta_mm, gamma_f, gamma_m, alpha_f, alpha_m);
sgrad_alpha_f = sgradAlphaF_VT(res, b, v, beta_max, p_f, p_m, d, mu, ...
    gamma_f, gamma_m, alpha_f, alpha_m);
sgrad_alpha_m = sgradAlphaM_VT(res, b, v, beta_max, p_f, p_m, d, mu, ...
    gamma_f, gamma_m, alpha_f, alpha_m);
% set flag to indicate whether or not the selection gradients are
% within the tolerance of zero
if max([abs(sgrad_gamma_f), abs(sgrad_gamma_m), ...
        abs(sgrad_alpha_f), abs(sgrad_alpha_m)]) < tolCSS
    flag = 0;
else
    flag = 1;
end
% iterate the selection gradients until they are within the tolerance
% of zero
while flag == 1
    % update parameters
    gamma_f = gamma_f + stepCSS*sgrad_gamma_f;
    gamma_m = gamma_m + stepCSS*sgrad_gamma_m;
    alpha_f = alpha_f + stepCSS*sgrad_alpha_f;
    alpha_m = alpha_m + stepCSS*sgrad_alpha_m;
    % check to make sure parameters are still positive
    if gamma_f < 0
        gamma_f = 1e-5;
    end
    if gamma_m < 0
        gamma_m = 1e-5;
    end
    if alpha_f < 0
        alpha_f = 1e-5;
    end
    if alpha_m < 0
        alpha_m = 1e-5;
    end
    % update value of b
    b = b_max - (c_f*gamma_f + c_m*gamma_m);
    % end the process if the birth rate becomes too small
    if b < 2*mu
        css = repmat(2000, 1, 14);
        return;
    end
    % update values of beta_ij
    beta_ff = beta_max*alpha_f/(alpha_f + d);
    beta_fm = beta_max*alpha_m/(alpha_m + d);
    beta_mf = beta_max*alpha_f/(alpha_f + d);
    beta_mm = beta_max*alpha_m/(alpha_m + d);
    % approximate the equilibrium values of the resident system
    res = findEE_VT(b, v, mu, beta_ff, beta_fm, beta_mf, beta_mm, ...
        gamma_f, gamma_m, alpha_f, alpha_m, stepEE, tolEE, maxitEE);
    % if EE couldn't be found, return corresponding initial condition for
    % closer inspection
    if isnan(res)
        css = [init_vals(1), init_vals(2), init_vals(3), init_vals(4), ...
            NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN];
        return;
    end
    % if system broke down when finding EE, end the process
    if sum(res) == 4000
        css = repmat(1000, 1, 14);
        return;
    end
    % update selection gradient values
    sgrad_gamma_f = sgradGammaF_VT(res, b_max, c_f, c_m, v, mu, beta_ff, ...
       beta_fm, beta_mf, beta_mm, gamma_f, gamma_m, alpha_f, alpha_m);
    sgrad_gamma_m = sgradGammaM_VT(res, b_max, c_f, c_m, v, mu, beta_ff, ...
       beta_fm, beta_mf, beta_mm, gamma_f, gamma_m, alpha_f, alpha_m);
   sgrad_alpha_f = sgradAlphaF_VT(res, b, v, beta_max, d, mu, ...
       gamma_f, gamma_m, alpha_f, alpha_m);
   sgrad_alpha_m = sgradAlphaM_VT(res, b, v, beta_max, d, mu, ...
       gamma_f, gamma_m, alpha_f, alpha_m);
    disp('.');
    % update flag if necessary
    if max([abs(sgrad_gamma_f), abs(sgrad_gamma_m), ...
            abs(sgrad_alpha_f), abs(sgrad_alpha_m)]) < tolCSS
        flag = 0;
    elseif (gamma_f <= 1e-5 && gamma_m <= 1e-5) || (gamma_f <= 1e-5 ...
            && abs(sgrad_gamma_m) < tolCSS) || (gamma_m <= 1e-5 ...
            && abs(sgrad_gamma_f) < tolCSS)
        flag = 0;
    end
end
% store equilibrium values
result = [gamma_f, gamma_m, alpha_f, alpha_m];
% check ESS conditions to confirm that this is the CSS and return result
ess_check_gamma = checkESSGamma_VT(result, res, b_max, c_f, c_m, v, ...
    beta_ff, beta_fm, beta_mf, beta_mm, mu);
ess_check_alpha = checkESSAlpha_VT(result, res, b, v, beta_max, d, mu);
res = res.';
css = [result, ess_check_gamma, ess_check_alpha, res];