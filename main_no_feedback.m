%% Esophagus Neuromechanical Model (No Distension Feedback)
% -------------------------------------------------------------------------
% Author: Guy Elishsa
%
% Description:
%   This script simulates the neuromechanics of the esophagus using a
%   coupled Wilson-Cowan oscillator framework and a 1-D mechanical model
%   of a fluid-filled, closed tube. It integrates coupled PDEs describing
%   excitatory (E) and inhibitory (I) neural populations, muscle
%   contraction (θ), and tube area and velocity.
%
%   In this simplified version, the sensory input from distension-sensitive
%   mechanoreceptors is replaced by constant terms:
%       S_E = C.P (excitatory input)
%       S_I = C.Q (inhibitory input)
%
% Requirements:
%   - MATLAB
%   - simulation_constants.m (defines model parameters)
%
% Outputs:
%   - Simulation data saved in Results1D_W_C/
%   - Figures showing spatiotemporal evolution of A, E, I, and θ
%
% References:
%   Gjorgjieva, J., et al. (2013). [Add full paper citation if applicable]
% -------------------------------------------------------------------------

clearvars; close all; clc
tic

%% Simulation Setup
casenum = 1;
C = simulation_constants();         % Load parameters
nxs = 70;                           % Spatial discretization points
nts = 1000;                         % Temporal discretization points
xspan = linspace(0, C.L, nxs);
tspan = linspace(0, C.duration, nts);

%% Initial Conditions
% u = [A, velocity, E, I, θ]
u0 = zeros(1, 5*nxs);
u0(1:nxs) = C.S_IC;                 % Initial area
u0(4*nxs+1:5*nxs) = 1;              % Initial θ
u0(2*nxs+1) = C.E_init;             % Initial E at oral end

%% Solver Options
reltol = 1e-4; 
abstol = 1e-4;
options = odeset('RelTol', reltol, 'AbsTol', abstol);

%% Integrate PDE System
[t2, u] = ode45(@(t,u) pde_1(t,u), tspan, u0, options);
timeToSolve = toc;

%% Unpack Variables
A   = u(:, 1:nxs);
v   = u(:, nxs+1:2*nxs);
E   = u(:, 2*nxs+1:3*nxs);
I   = u(:, 3*nxs+1:4*nxs);
tht = u(:, 4*nxs+1:5*nxs);

%% Pressure Calculation
P = zeros(nts, nxs);
for i = 1:nts
    P(i,:) = A(i,:) ./ tht(i,:) - 1;
end

%% Transpose for Visualization
A = A'; v = v'; E = E'; I = I'; tht = tht';

%% Visualization
figure('Color','w');
cmap = flipud(jet);

subplot(2,2,1);
imagesc(tspan, xspan, A); colormap(cmap); colorbar;
xlabel('Time (\tau)','Interpreter','latex'); ylabel('Length (\chi)','Interpreter','latex');
title('Cross-sectional Area (A)');

subplot(2,2,2);
imagesc(tspan, xspan, E); colormap(cmap); colorbar;
xlabel('Time (\tau)','Interpreter','latex'); ylabel('Length (\chi)','Interpreter','latex');
title('Excitatory Activity (E)');

subplot(2,2,3);
imagesc(tspan, xspan, I); colormap(cmap); colorbar;
xlabel('Time (\tau)','Interpreter','latex'); ylabel('Length (\chi)','Interpreter','latex');
title('Inhibitory Activity (I)');

subplot(2,2,4);
imagesc(tspan, xspan, tht); colormap(cmap); colorbar;
xlabel('Time (\tau)','Interpreter','latex'); ylabel('Length (\chi)','Interpreter','latex');
title('Muscle Activation (\theta)');

%% Save Results
dirname = fullfile('Results1D_W_C', sprintf('case_%02d', casenum));
if ~exist(dirname, 'dir'); mkdir(dirname); end
save(fullfile(dirname, 'results.mat'), 'A', 'v', 'E', 'I', 'tht', 'P', 'tspan', 'xspan');

%% ------------------------------------------------------------------------
% PDE System Definition
% ------------------------------------------------------------------------
function ut = pde_1(t, u)
    C = simulation_constants();
    n = length(u) / 5;
    dx = C.L / n;
    ut = zeros(size(u));

    %% (1) Area dynamics
    for i = 1:n
        j = i + n;
        if i == 1
            ut(i) = (-1/(2*dx)) * (u(i) * u(j) + u(i+1) * u(j+1));
        elseif i == n
            ut(i) = (1/(2*dx)) * (u(i) * u(j) + u(i-1) * u(j-1));
        else
            ut(i) = (1/(2*dx)) * (u(i-1) * u(j-1) - u(i+1) * u(j+1));
        end
    end

    %% (2) Velocity dynamics
    for i = n+1:2*n
        j = i - n;
        l = i + 3*n;
        if i == n+1
            ut(i) = -C.BETA * u(i)/u(j) + C.PSI/(2*dx) * (u(j)/u(l) - u(j+1)/u(l+1)) ...
                   - 1/(8*dx) * (u(i+1)^2 + 2*u(i+1)*u(i) + u(i)^2);
        elseif i == 2*n
            ut(i) = -C.BETA * u(i)/u(j) + C.PSI/(2*dx) * (u(j-1)/u(l-1) - u(j)/u(l)) ...
                   + 1/(8*dx) * (u(i-1)^2 + 2*u(i-1)*u(i) + u(i)^2);
        else
            ut(i) = -C.BETA * u(i)/u(j) + C.PSI/(2*dx) * (u(j-1)/u(l-1) - u(j+1)/u(l+1)) ...
                   + 1/(8*dx) * (u(i-1)^2 - u(i+1)^2 + 2*u(i-1)*u(i) - 2*u(i+1)*u(i));
        end
    end

    %% (3) Excitatory/Inhibitory and Muscle Dynamics
    for i = 2*n+1:3*n
        j = i + n;       % inhibitory neuron index
        l = j + n;       % theta index

        % Constant sensory inputs
        S_E = C.P;
        S_I = C.Q;

        % Neural coupling terms
        if i == 2*n+1
            x_e = C.a*u(i) - C.e*u(j) + S_E;
        else
            x_e = C.b_p*u(i-1) + C.a*u(i) - C.d_p*u(j-1) - C.e*u(j) + S_E;
        end
        x_i = C.c*u(i) - C.f*u(j) + S_I;

        % E and I dynamics
        ut(i) = (-u(i) + (C.k_e - u(i)) * S(C.a_e, C.tht_e, x_e)) / C.tau_e;
        ut(j) = (-u(j) + (C.k_i - u(j)) * S(C.a_i, C.tht_i, x_i)) / C.tau_i;

        % Muscle (θ) dynamics
        x_tht = u(i) - C.E_h;
        ut(l) = (1 - u(l) - sigma(x_tht, 0, 0.5*(1 - C.tht_0), C.g_tht, 0)) / C.tau_tht;
    end
end

%% ------------------------------------------------------------------------
% Supporting Functions
% ------------------------------------------------------------------------
function sig = S(a, theta, x)
    sig = 1 / (1 + exp(-a * (x - theta))) - 1 / (1 + exp(a * theta));
end

function sig = sigma(xk, xi, A, g, xs)
    sig = A * (1 + tanh(g * (xk - xi + xs)));
end
