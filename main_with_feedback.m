%% Esophagus Neuromechanical Model
% -------------------------------------------------------------------------
% Author: Guy Elishsa
% Description:
%   This script simulates the neuromechanics of the esophagus using a
%   coupled Wilson-Cowan oscillator framework and a 1-D mechanical model
%   of a fluid filled tube closed on both ends. 
%   It numerically integrates a set of coupled PDEs representing excitatory (E) and inhibitory (I)
%   neural populations, muscle contraction pattern (θ), and tube area and fluid velocity.
%
% Requirements:
%   - MATLAB 
%   - Function: simulation_constants.m (defines model constants)
%
% Outputs:
%   - Simulation data saved in Results1D_W_C/
%   - Figures showing spatiotemporal evolution of variables
%
% -------------------------------------------------------------------------
clearvars; close all; clc
timerightnow = clock;
tic

%% Simulation Case
casenum = 1;

%% Parameters
C = simulation_constants();        % Loading simulation parameters
nxs = 70;                % Number of spatial discretization points
nts = 1000;              % Number of temporal discretization points
xspan = linspace(0, C.L, nxs);
tspan = linspace(0, C.duration, nts);

%% Initial Conditions
% u = [u1 (area), u2 (velocity), E, I, θ]
u0 = zeros(1,5*nxs);
u0(1:nxs) = C.S_IC;           % Initial area (cross-sectional)
u0(4*nxs+1:5*nxs) = 1;        % Initial θ (normalized to 1)
u0(2*nxs+1) = C.E_init;       % Initial excitation at to oral end (set to zero)

%% ODE Solver Setup
reltol=1.0e-04; 
abstol=1.0e-04;
options=odeset('RelTol',reltol,'AbsTol',abstol);

%% Integrate Solver
[t2,u]=ode45(@(t,u) pde_1(t,u),tspan,u0,options);
timeToSolve = toc

%% Unpack Results
u1  = u(:, 1:nxs);                 % Area
u2  = u(:, nxs+1:2*nxs);           % Velocity
E   = u(:, 2*nxs+1:3*nxs);         % Excitatory neurons
I   = u(:, 3*nxs+1:4*nxs);         % Inhibitory neurons
tht = u(:, 4*nxs+1:5*nxs);         % Muscle contraction pattern

%% Calculate pressure
for i = 1:1:nts
    for j = 1:1:nxs
        p_tn(j)   = u1(i,j)/tht(i,j) - 1;
    end
    P(i,:) = p_tn';
end

%% Transpose for plotting
u1 = u1'; u2 = u2'; E = E'; I = I'; tht = tht';

%% Visualization
mainFigHandle = figure(6);
set(gcf,'color','w');
c = flipud(jet);

subplot(2,2,[1]);
colormap(c); colorbar
imagesc(tspan,xspan,u1)
ylabel('Length \chi','interpreter','latex')
xlabel('Time \tau','interpreter','latex')
title('CSA')

%
subplot(2,2,[2]);
colormap(c); colorbar
imagesc(tspan,xspan,E)
ylabel('Length \chi','interpreter','latex')
xlabel('Time \tau','interpreter','latex')
title('E')
%
subplot(2,2,[3]);
colormap(c); colorbar
imagesc(tspan,xspan,I)
ylabel('Length \chi','interpreter','latex')
xlabel('Time \tau','interpreter','latex')
title('I')
%
subplot(2,2,[4]);
colormap(c); colorbar
imagesc(tspan,xspan,tht)
ylabel('Length \chi','interpreter','latex')
xlabel('Time \tau','interpreter','latex')
title('\theta')

%% Save Results
dirname = fullfile('Results1D_W_C', sprintf('case_%02d', casenum));
if ~exist(dirname, 'dir')
    mkdir(dirname);
end
save(fullfile(dirname, 'results.mat'));


%% ------------------------------------------------------------------------
% PDE System Definition
% ------------------------------------------------------------------------
function ut=pde_1(t,u)
    C = simulation_constants();
    n = length(u) / 5;
    dx = C.L / n;

    %% Area (continuity-like equation)
    for i=1:n
        j = i+n;
        if(i==1)     ut(i)=(-1/(2*dx)) * (u(i) * u(j) + u(i+1) * (u(j+1)));
        elseif(i==n) ut(i)=(1/(2*dx)) * (u(i) * u(j) + u(i-1) * (u(j-1)));
        else         ut(i)=(1/(2*dx)) * (u(i-1) * u(j-1) - u(i+1) * (u(j+1)));
        end
    end

    %% velocity (momentum-like equation)
    for i=1+n:2*n
        j = i-n;
        l = i+3*n;
        if(i==n+1)     ut(i)=-C.BETA * u(i)/u(j) + C.PSI/(2*dx) * (u(j)/u(l) - u(j+1)/u(l+1))...
                - 1/(8*dx) * (u(i+1)*u(i+1) + 2 * u(i+1)*u(i) + u(i)*u(i));
        elseif(i==2 * n) ut(i)=-C.BETA * u(i)/u(j) + C.PSI/(2*dx) * (u(j-1)/u(l-1) - u(j)/u(l))...
                + 1/(8*dx) * (u(i-1)*u(i-1) + 2 * u(i-1)*u(i) + u(i)*u(i));
        else ut(i) = -C.BETA * u(i)/u(j) + C.PSI/(2*dx) * (u(j-1)/u(l-1) - u(j+1)/u(l+1))...
                + 1/(8*dx) * (u(i-1)*u(i-1) - u(i+1)*u(i+1) + 2 * u(i-1)*u(i) - 2*u(i+1)*u(i));
        end
    end

%% E and I, Mechanical coupling
    for i = 2*n+1:3*n
        j = i + n;
        l = j + n;
        x = dx/2 + (i-1-2*n)*dx;

        % Compute strain integrals for E and I
        StrainI = 0; StrainE = 0;
        for k = 1:i
            StrainI = StrainI + max(0,(u(k) / u(4*n + k)-1-C.alpha_h))*dx;
        end
        
        for k = i:n
            x_k = dx/2 + (k-1)*dx;
            StrainE = StrainE + max(0,(u(k) / u(4*n + k)-1-C.alpha_h))*sigma(x_k,x,0.5,-C.g,-C.xs)*dx;
        end
        % Input from mechanoreceptors to excitatory and inhibitory neuronal populations
        S_E = C.P * sigma_s(StrainE,0,1,C.g,0);
        S_I = C.Q .* sigma_s(StrainI,0,1,C.g,0);

        % E and I dynamics
        if i == 2*n+1
            x_e = C.a * u(i) - C.e * u(j) + S_E;
        else
            x_e = C.b_p * u(i-1) + C.a * u(i) - C.d_p * u(j-1) - C.e * u(j) + ...
                S_E;
        end
        x_i = C.c * u(i) - C.f * u(j) + S_I;

        ut(i) = (-u(i) + (C.k_e - u(i)) * S(C.a_e, C.tht_e, x_e))/C.tau_e;
        ut(j) = (-u(j) + (C.k_i - u(j)) * S(C.a_i, C.tht_i, x_i))/C.tau_i;

        % Muscle contraction pattern θ dynamics
        x_tht = u(i)-C.E_h;
        ut(l) = (1 - u(l)-sigma(x_tht,0,0.5*(1-C.tht_0),C.g_tht,0))/C.tau_tht;
    end

    ut=ut';
end

%% ------------------------------------------------------------------------
% Supporting Functions
% ------------------------------------------------------------------------
function sig = S(a,theta,x)
    sig = 1/(1+exp(-a * (x-theta))) - 1/(1+exp(a * (theta)));
end

%% sigmoid functions
function sig = sigma(xk,xi,A,g,xs)
    sig = A * (1 + tanh(g * (xk-xi+xs)));
end

function sig = sigma_s(xk,xi,A,g,xs)
    sig =  tanh(g * (xk-xi+xs));
end
