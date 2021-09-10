function endemic = findEE_VT(b, v, mu, beta_ff, beta_fm, beta_mf, beta_mm, ...
    gamma_f, gamma_m, alpha_f, alpha_m, stepEE, tol, maxit)
% A script to approximate the endemic equilibrium of the resident system
% for a given set of parameter values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a set of random initial conditions
init_vals = 10*rand(50, 4);
% setup an empty array to store the equilibrium values for all initial
% conditions
endemic = zeros(4, size(init_vals, 1));
% setup a loop counter
itn = 0;
% run through each set of initial conditions
for i = 1:size(init_vals, 1)
    % check the loop counter
    if itn >= maxit
        endemic = NaN;
        disp('Maximum iterations exceeded');
        return
    end
    % run Euler's method to iterate the resident system
    x0 = init_vals(i, :);
    % if we get overflow in x0, return indicator that system breaks down
    if isnan(sum(x0))
        endemic = repmat(1000, 1, 4);
        return;
    end
    dxdt = resident_VT(x0, b, v, mu, beta_ff, beta_fm, beta_mf, ...
        beta_mm, gamma_f, gamma_m, alpha_f, alpha_m);
    % set a flag to indicate the resident system is away from equilibrium
    if max(abs(dxdt)) < tol
        flag = 0;
    else
        flag = 1;
    end
    % continue to iterate as long as the system is away from equilibrium
    while flag == 1
        x0 = x0 + stepEE*dxdt;
        % if we get overflow in x0, return indicator that system breaks down
        if isnan(sum(x0))
            endemic = repmat(1000, 1, 4);
            return;
        end
        % check that everything is still positive
        for k = 1:4
            if x0(k) < 0
                x0(k) = stepEE^2;
            end
        end
        dxdt = resident_VT(x0, b, v, mu, beta_ff, beta_fm, beta_mf, ...
            beta_mm, gamma_f, gamma_m, alpha_f, alpha_m);   
        % check if the derivatives are sufficiently close to zero
        if max(abs(dxdt)) < stepEE^2
            flag = 0;
        end
    end
    % store result for the equilibrium values
    endemic(:, i) = x0;
    disp('>');
    itn = itn + 1;
end
% check the distance between each pair of columns and track the maximum
% distance found
maxDist = 0;
for i = 1:(size(endemic, 2) - 1)
    for j = (i + 1):size(endemic, 2)
        % check distance and update maxDist if necessary
        colDist = norm(endemic(:, i) - endemic(:, j));
        if colDist > maxDist
            maxDist = colDist;
        end
    end
end
if maxDist < tol
    % return average across all columns
    endemic = mean(endemic, 2);
else
    % return message to indicate possibility of multiple stable equilibria
    endemic = NaN;
    disp('Possible instance of multiple endemic equilibria');
end