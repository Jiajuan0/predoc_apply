% Define parameters
rho = 0.7; % AR(1) persistence parameter
sigma_eps = 0.3; % Standard deviation of AR(1) shocks
n = 50 ; % Number of grid points
m = 1; % Number of standard deviations around the mean
alpham=0.5;
alphak=0.2;
alphal=0.3;
tol=0.01;
maxits=1000;
dif=tol+1000;
its=0;
beta=0.9;
delta=0.9;

% Call the tauchen function
[grid, P, stationary_dist] = tauchen(rho, sigma_eps, n, m);

% Parameters
tol = 1e-4; % Convergence tolerance
max_iter = 1000; % Maximum number of iterations
% State space for X




grid_min_x = 0; % Minimum state for X
grid_max_x = 0.3; % Maximum state for X
n_grid_x = 50; % Number of grid points for X
grid_x = linspace(grid_min_x, grid_max_x, n_grid_x)';


% State space for Z
grid_z=exp(grid);
n_grid_z = 50; % Number of states for Z


profit_function = @(x, z, x_next) z*x.^alphak...
    *(((alpham.^(1-alphal)*alphal.^alphal*z*x.^alphak).^(1/(1-alpham-alphal)))).^alpham ...
    *(((alphal.^(1-alpham)*alpham.^alpham*z*x.^alphak).^(1/(1-alpham-alphal)))).^alphal ...
    -(alpham.^(1-alphal)*alphal.^alphal*z*x.^alphak).^(1/(1-alpham-alphal))...
    -(alphal.^(1-alpham)*alpham.^alpham*z*x.^alphak).^(1/(1-alpham-alphal))...
    -(x_next-delta*x).^2;


% Value function initialization
V = zeros(n_grid_x, n_grid_z);
V_new = zeros(n_grid_x, n_grid_z);
% Value function iteration
for iter = 1:max_iter
for ix = 1:n_grid_x
for iz = 1:n_grid_z
V_temp = zeros(n_grid_x, 1);
for iz_next = 1:n_grid_z
    for ix_next=1:n_grid_x
% Compute the value for each action and each possible next state of Z
temp = profit_function(grid_x(ix), grid_z(iz), grid_x(ix_next)) + ...
beta * V(:, ix_next);
    end
V_temp = V_temp + P(iz, iz_next) * temp;
end
V_new(ix, iz) = max(V_temp);
end
end
% Check for convergence
if max(max(abs(V_new - V))) < tol
fprintf('Converged after %d iterations\n', iter);
break;
end
% Update value function
V = V_new;

end

% Optionally, you can plot the value function for a specific state of Z
z_plot = 1; % Choose which Z state to plot
plot(grid_x, V(:, z_plot));
title(['Value Function for Productivity']);
xlabel('Capital');
ylabel('Value');
hold on;
plot(grid_x, V(:, 5));
hold on;
plot(grid_x, V(:, 10));
hold on;
plot(grid_x, V(:, 15));
hold on;
plot(grid_x, V(:, 20));
hold on;
plot(grid_x, V(:, 25));
hold on;
plot(grid_x, V(:, 30));
hold on;
plot(grid_x, V(:, 35));
hold on;
plot(grid_x, V(:, 40));
hold on;
plot(grid_x, V(:, 45));
hold on;
plot(grid_x, V(:, 50));