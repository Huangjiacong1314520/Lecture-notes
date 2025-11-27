%% Example usage and test function
% function test_skew_mu()
    % Test case: simple 4x4 system
    % 2 fixed dimensions, 2 varying dimensions
clc
clear
close
    M = [0.5+0.1i, 0.2, 0.1, 0.05;
         0.2, 0.6-0.1i, 0.15, 0.1;
         0.1, 0.15, 0.4+0.2i, 0.2;
         0.05, 0.1, 0.2, 0.5-0.1i];
    
    % Fixed: 1 real scalar (size 1), 1 complex scalar (size 1)
    % Varying: 1 real scalar (size 1), 1 complex scalar (size 1)
    k_f = [1, 1, 0];
    k_v = [1, 1, 0];
    
    [mu_s_lower, v_total] = skew_mu_lower_bound(M, k_f, k_v);
    
    fprintf('Skew mu lower bound: %.6f\n', mu_s_lower);
    fprintf('Total scaling factor: %.6f\n', v_total);
% end

