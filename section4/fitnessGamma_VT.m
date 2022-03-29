function fitnessGamma = fitnessGamma_VT(res, b_max, c_f, c_m, v, ...
    beta_ff, beta_fm, beta_mf, beta_mm, mu, alpha_f, alpha_m, ...
    gamma_f, gamma_m, gamma_f_mut, gamma_m_mut)
% A function to compute the fitness function for the host recovery
% based on output from Maple. Input res is a four-element vector containing
% the equilibrium resident population values S_f, S_m, I_f, I_m.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evaluate the fitness function
fitnessGamma = (-(res(3) * (((-beta_mf * res(3) - beta_mm * res(4) - mu - alpha_m - gamma_m_mut) * res(1) + res(3) * (v * alpha_m - beta_mf * res(3) - beta_mm * res(4) - mu - alpha_m - gamma_m_mut)) * (-gamma_f * c_f - c_m * gamma_m_mut + b_max) * alpha_f - mu * (-gamma_f * c_f - c_m * gamma_m_mut + b_max) * (beta_mf * res(3) + beta_mm * res(4) + mu + alpha_m + gamma_m_mut) * res(1) - ((res(2) + res(4)) * (-c_f * gamma_f_mut - gamma_m * c_m + b_max) * alpha_m + (res(3) * (-gamma_f * c_f - c_m * gamma_m_mut + b_max) + (res(2) + res(4)) * (-c_f * gamma_f_mut - gamma_m * c_m + b_max)) * mu) * res(3) * beta_mf - ((res(2) + res(4)) * (-c_f * gamma_f_mut - gamma_m * c_m + b_max) * alpha_m + (res(3) * (-gamma_f * c_f - c_m * gamma_m_mut + b_max) + (res(2) + res(4)) * (-c_f * gamma_f_mut - gamma_m * c_m + b_max)) * mu) * res(4) * beta_mm + (v * alpha_m * (-gamma_f * c_f - c_m * gamma_m_mut + b_max) * res(3) - (res(3) * (-gamma_f * c_f - c_m * gamma_m_mut + b_max) + (res(2) + res(4)) * (-c_f * gamma_f_mut - gamma_m * c_m + b_max)) * (mu + alpha_m + gamma_m_mut)) * mu) * beta_ff) + ((-((-beta_mf * res(3) - beta_mm * res(4) - mu - alpha_m - gamma_m_mut) * res(1) + res(3) * (v * alpha_m - beta_mf * res(3) - beta_mm * res(4) - mu - alpha_m - gamma_m_mut)) * (-gamma_f * c_f - c_m * gamma_m_mut + b_max) * res(4) * beta_fm + mu * (-gamma_f * c_f - c_m * gamma_m_mut + b_max) * (beta_mf * res(3) + beta_mm * res(4) + mu + alpha_m + gamma_m_mut) * res(1) + ((res(2) + res(4)) * (-c_f * gamma_f_mut - gamma_m * c_m + b_max) * alpha_m + (res(3) * (-gamma_f * c_f - c_m * gamma_m_mut + b_max) + (res(2) + res(4)) * (-c_f * gamma_f_mut - gamma_m * c_m + b_max)) * mu) * res(3) * beta_mf + ((res(2) + res(4)) * (-c_f * gamma_f_mut - gamma_m * c_m + b_max) * alpha_m + (res(3) * (-gamma_f * c_f - c_m * gamma_m_mut + b_max) + (res(2) + res(4)) * (-c_f * gamma_f_mut - gamma_m * c_m + b_max)) * mu) * res(4) * beta_mm - (v * alpha_m * (-gamma_f * c_f - c_m * gamma_m_mut + b_max) * res(3) - (res(3) * (-gamma_f * c_f - c_m * gamma_m_mut + b_max) + (res(2) + res(4)) * (-c_f * gamma_f_mut - gamma_m * c_m + b_max)) * (mu + alpha_m + gamma_m_mut)) * mu) * alpha_f) - ((-mu * (-gamma_f * c_f - c_m * gamma_m_mut + b_max) * (beta_mf * res(3) + beta_mm * res(4) + mu + alpha_m + gamma_m_mut) * res(1) - ((res(2) + res(4)) * (-c_f * gamma_f_mut - gamma_m * c_m + b_max) * alpha_m + (res(3) * (-gamma_f * c_f - c_m * gamma_m_mut + b_max) + (res(2) + res(4)) * (-c_f * gamma_f_mut - gamma_m * c_m + b_max)) * mu) * res(3) * beta_mf - ((res(2) + res(4)) * (-c_f * gamma_f_mut - gamma_m * c_m + b_max) * alpha_m + (res(3) * (-gamma_f * c_f - c_m * gamma_m_mut + b_max) + (res(2) + res(4)) * (-c_f * gamma_f_mut - gamma_m * c_m + b_max)) * mu) * res(4) * beta_mm + (v * alpha_m * (-gamma_f * c_f - c_m * gamma_m_mut + b_max) * res(3) - (res(3) * (-gamma_f * c_f - c_m * gamma_m_mut + b_max) + (res(2) + res(4)) * (-c_f * gamma_f_mut - gamma_m * c_m + b_max)) * (mu + alpha_m + gamma_m_mut)) * mu) * res(4) * beta_fm) + (mu * (-gamma_f * c_f - c_m * gamma_m_mut + b_max) * (mu + gamma_f_mut) * (beta_mf * res(3) + beta_mm * res(4) + mu + alpha_m + gamma_m_mut) * res(1)) + (res(3) * (mu * (-c_f * gamma_f_mut - gamma_m * c_m + b_max) * (mu + alpha_m) * (res(2) + res(4)) * v + ((res(2) + res(4)) * (-c_f * gamma_f_mut - gamma_m * c_m + b_max) * alpha_m + (res(3) * (-gamma_f * c_f - c_m * gamma_m_mut + b_max) + (res(2) + res(4)) * (-c_f * gamma_f_mut - gamma_m * c_m + b_max)) * mu) * (mu + gamma_f_mut)) * beta_mf) + ((mu * (-c_f * gamma_f_mut - gamma_m * c_m + b_max) * (mu + alpha_m) * (res(2) + res(4)) * v + ((res(2) + res(4)) * (-c_f * gamma_f_mut - gamma_m * c_m + b_max) * alpha_m + (res(3) * (-gamma_f * c_f - c_m * gamma_m_mut + b_max) + (res(2) + res(4)) * (-c_f * gamma_f_mut - gamma_m * c_m + b_max)) * mu) * (mu + gamma_f_mut)) * res(4) * beta_mm) + (mu * ((-res(3) * (-gamma_f * c_f - c_m * gamma_m_mut + b_max) * gamma_f_mut + (-res(3) * (-gamma_f * c_f - c_m * gamma_m_mut + b_max) + (res(2) + res(4)) * (-c_f * gamma_f_mut - gamma_m * c_m + b_max)) * mu) * alpha_m + mu * (-c_f * gamma_f_mut - gamma_m * c_m + b_max) * (mu + gamma_m_mut) * (res(2) + res(4))) * v) + ((res(3) * (-gamma_f * c_f - c_m * gamma_m_mut + b_max) + (res(2) + res(4)) * (-c_f * gamma_f_mut - gamma_m * c_m + b_max)) * (mu + gamma_f_mut) * mu * alpha_m) + ((res(3) * (-gamma_f * c_f - c_m * gamma_m_mut + b_max) + (res(2) + res(4)) * (-c_f * gamma_f_mut - gamma_m * c_m + b_max)) * mu * (mu + gamma_m_mut) * gamma_f_mut) + (mu ^ 2 * (-gamma_f * c_f - c_m * gamma_m_mut + b_max) * (mu + gamma_m_mut) * res(3)) + (mu ^ 2 * (-c_f * gamma_f_mut - gamma_m * c_m + b_max) * (res(2) + res(4)) * gamma_m_mut) + (mu ^ 3 * (-c_f * gamma_f_mut - gamma_m * c_m + b_max) * res(2)) + (mu ^ 3 * (-c_f * gamma_f_mut - gamma_m * c_m + b_max) * res(4)) + sqrt(((mu ^ 2 + (beta_mf * res(3) + beta_mm * res(4) + alpha_m + gamma_m_mut) * mu + alpha_m * (beta_mf * res(3) + beta_mm * res(4))) ^ 2) * (((-1 + v) ^ 2 * mu ^ 2) - (2 * (beta_ff * res(3) + beta_fm * res(4) + alpha_f + gamma_f_mut) * (-1 + v) * mu) + (beta_fm ^ 2 * res(4) ^ 2) - 0.4e1 * beta_fm * (-(beta_ff * res(3)) / 0.2e1 + (v - 0.1e1 / 0.2e1) * alpha_f - gamma_f_mut / 0.2e1) * res(4) + (beta_ff ^ 2 * res(3) ^ 2) - 0.4e1 * ((v - 0.1e1 / 0.2e1) * alpha_f - gamma_f_mut / 0.2e1) * beta_ff * res(3) + ((alpha_f + gamma_f_mut) ^ 2)) * ((res(2) + res(4)) ^ 2) * ((-c_f * gamma_f_mut - gamma_m * c_m + b_max) ^ 2) - 0.2e1 * (mu ^ 2 + (beta_mf * res(3) + beta_mm * res(4) + alpha_m + gamma_m_mut) * mu + alpha_m * (beta_mf * res(3) + beta_mm * res(4))) * (((res(3) + res(1)) * (-1 + v) * mu ^ 2) + (((res(3) + res(1)) * ((-1 + v) * beta_mm - beta_fm) * res(4) + ((-1 + v) * beta_mf - beta_ff) * res(3) ^ 2 + ((v ^ 2 - 1) * alpha_m + (2 * v - 1) * alpha_f + ((-1 + v) * beta_mf - beta_ff) * res(1) - gamma_f_mut + (-1 + v) * gamma_m_mut) * res(3) + ((-1 + v) * alpha_m - alpha_f - gamma_f_mut + (-1 + v) * gamma_m_mut) * res(1)) * mu) - (beta_fm * beta_mm * (res(3) + res(1)) * res(4) ^ 2) + (((-beta_ff * beta_mm - beta_fm * beta_mf) * res(3) ^ 2) + ((beta_fm * (-1 + v) * alpha_m + beta_mm * (2 * v - 1) * alpha_f + (-beta_ff * beta_mm - beta_fm * beta_mf) * res(1) - beta_fm * gamma_m_mut - beta_mm * gamma_f_mut) * res(3)) + 0.2e1 * (beta_fm * (v - 0.1e1 / 0.2e1) * alpha_m - (alpha_f * beta_mm) / 0.2e1 - (beta_fm * gamma_m_mut) / 0.2e1 - (beta_mm * gamma_f_mut) / 0.2e1) * res(1)) * res(4) - (beta_ff * beta_mf * res(3) ^ 3) + ((beta_ff * (-1 + v) * alpha_m + beta_mf * (2 * v - 1) * alpha_f - beta_ff * beta_mf * res(1) - beta_ff * gamma_m_mut - beta_mf * gamma_f_mut) * res(3) ^ 2) + (((alpha_f * (-1 + v) + res(1) * (2 * v - 1) * beta_ff - gamma_f_mut * (v + 1)) * alpha_m + (-res(1) * beta_mf + gamma_m_mut * (2 * v - 1)) * alpha_f + (-beta_ff * gamma_m_mut - beta_mf * gamma_f_mut) * res(1) - gamma_f_mut * gamma_m_mut) * res(3)) - (res(1) * (alpha_m + gamma_m_mut) * (alpha_f + gamma_f_mut))) * (mu ^ 2 + (beta_ff * res(3) + beta_fm * res(4) + alpha_f + gamma_f_mut) * mu + alpha_f * (beta_ff * res(3) + beta_fm * res(4))) * (res(2) + res(4)) * (-gamma_f * c_f - c_m * gamma_m_mut + b_max) * (-c_f * gamma_f_mut - gamma_m * c_m + b_max) + (((-res(3) - res(1)) * mu - beta_mm * (res(3) + res(1)) * res(4) - res(3) ^ 2 * beta_mf + ((-1 + v) * alpha_m - res(1) * beta_mf - gamma_m_mut) * res(3) - res(1) * (alpha_m + gamma_m_mut)) ^ 2 * (mu ^ 2 + (beta_ff * res(3) + beta_fm * res(4) + alpha_f + gamma_f_mut) * mu + alpha_f * (beta_ff * res(3) + beta_fm * res(4))) ^ 2 * (-gamma_f * c_f - c_m * gamma_m_mut + b_max) ^ 2))) / (beta_mf * (mu + alpha_m) * res(3) + beta_mm * (mu + alpha_m) * res(4) + mu * (mu + alpha_m + gamma_m_mut)) / (res(1) + res(3) + res(2) + res(4)) / (res(3) * (mu + alpha_f) * beta_ff + (beta_fm * res(4) + mu) * alpha_f + mu * (beta_fm * res(4) + mu + gamma_f_mut)) / 0.4e1;