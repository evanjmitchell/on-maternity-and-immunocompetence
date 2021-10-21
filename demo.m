tic;
% A script to demo the use of the Matlab code files and generate a small example data set .
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% choose a set of cost differences
cdiff_vals = linspace(-0.03, 0.03, 5);
% choose a set of vertical transmission values
v_vals = linspace(0, 0.9, 4);
% setup empty array to store data
points = zeros(11, length(cdiff_vals)*length(v_vals));
points(1, :) = repelem(cdiff_vals, length(v_vals));
points(2, :) = repmat(v_vals, 1, length(cdiff_vals));
cdiff_vals = points(1, :);
v_vals = points(2, :);
% for each cost difference/vertical transmission pair, compute the coCSS
% and record an indicator: 1 for high female virulence/high ...
% female recovery, 2 for high female virulence/low female recovery, ...
% 3 for low female virulence/low female recovery, 4 for low female ...
% virulence/high female recovery, 5 for cases where the model breaks ...
% down, 6 for cases where the birth rate falls below 2*mu and the ...
% population crashes, or 7 for cases where there are potentially ...
% multiple stable endemic equilibria
tmp = zeros(8, size(points, 2));
parfor i = 1:size(points, 2)
    cdiff = cdiff_vals(i);
    v = v_vals(i);
    css = findCSS_VT([5, 5.5, 3, 3.2], 10, 0.15 + cdiff, 0.15 - cdiff, ...
        v, 4, 10, 1, 0.001, 0.001, 10000, 10, 1e-5);
    if sum (css)/14 == 1000
        points(3, i) = 5;
    elseif sum (css)/14 == 2000
        points(3, i) = 6;
    elseif isnan(sum(css))
        points(3, i) = 7;
    else
        gamma_diff = css(1) - css(2);
        alpha_diff = css(3) - css(4);
        if alpha_diff > 0 && gamma_diff > 0
            points(3, i) = 1;
        elseif alpha_diff > 0 && gamma_diff < 0
            points(3, i) = 2;
        elseif alpha_diff < 0 && gamma_diff < 0
            points(3, i) = 3;
        else
            points(3, i) = 4;
        end
    end
    eqmVals = [css(1:4), css(11:14)];
    tmp(:, i) = eqmVals.';
end
points(4:11, :) = tmp;
% save data as xlsx file
writematrix(points, 'data_demo.xlsx');
toc;