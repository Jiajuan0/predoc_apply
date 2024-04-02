function [grid, P, stationary_dist] = tauchen(rho, sigma_eps, n, m)
    % Inputs:
    % rho - AR(1) process persistence parameter
    % sigma_eps - standard deviation of the error term
    % n - number of grid points
    % m - width of the state space, in terms of standard deviations from the mean
    
    % Step 1: Construct state space
    Z_min = -(m * sigma_eps) / sqrt(1 - rho^2); % Minimum state
    Z_max = (m * sigma_eps) / sqrt(1 - rho^2);  % Maximum state
    step = (Z_max - Z_min) / (n - 1);           % Distance between grid points
    Z = linspace(Z_min, Z_max, n);              % State space grid
    
    % Step 2: Construct transition matrix
    P = zeros(n,n);
    for j = 1:n
        for k = 1:n
            if k == 1
                P(j,k) = normcdf((Z(k) - rho*Z(j) + step/2) / sigma_eps);
            elseif k == n
                P(j,k) = 1 - normcdf((Z(k) - rho*Z(j) - step/2) / sigma_eps);
            else
                P(j,k) = normcdf((Z(k) - rho*Z(j) + step/2) / sigma_eps) - ...
                         normcdf((Z(k) - rho*Z(j) - step/2) / sigma_eps);
            end
        end
    end
    
    % Step 3: Find stationary distribution
    [V,D] = eig(P');
    [~,ind] = max(diag(D));
    stationary_dist = V(:,ind)' / sum(V(:,ind));
    
    % Return output
    grid = Z;
end
