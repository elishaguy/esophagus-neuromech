function C = simulation_constants()

% -------------------------------------------------------------------------
% Defines all model parameters for the coupled Wilson–Cowan neuromechanical
% oscillator representing esophageal peristalsis.
%
% The constants define:
%   - Neural network coupling parameters
%   - Mechanoreceptor input strengths
%   - Sigmoid activation parameters
%   - Time constants for dynamic equations
%   - Mechanical and geometric properties
%   - Simulation grid and time duration
%
% Output:
%   C : Structure containing all model parameters
%
% Usage:
%   C = simulation_constants();
% -------------------------------------------------------------------------

%% --------------------------- Neural Parameters ---------------------------
% Coupling coefficients in the Wilson–Cowan equations
% These define how excitatory (E) and inhibitory (I) populations interact.
C.a   = 16;   % Self-excitation of excitatory neurons
C.b_p = 20;   % Excitatory input from neighboring neurons (from proximal)
C.c   = 12;   % Excitatory influence on inhibitory population
C.d_p = 40;   % Inhibitory input from neighboring neurons (from proximal)
C.e   = 15;   % Inhibitory influence on excitatory population
C.f   = 3;    % Self-inhibition of inhibitory neurons

%% ------------------- Strength of Mechanoreceptor Input -------------------
% Strength of sensory input from mechanical stretch receptors
C.P = 1.6;    % Input weight to excitatory population (w_E)
C.Q = 1.35;   % Input weight to inhibitory population (w_I)

%% -------------------------- Equation Constants ---------------------------
% Baseline scaling and time constants for dynamic equations
C.k_e     = 1;     % Scaling factor for excitatory dynamics
C.k_i     = 1;     % Scaling factor for inhibitory dynamics
C.tau_e   = 1;     % Excitatory neuron time constant
C.tau_i   = 4;     % Inhibitory neuron time constant
C.tau_tht = 0.2;   % Muscle deformation (θ) time constant

% Thresholds and steady-state parameters
C.E_h    = 0.3;    % Excitatory activity threshold for θ dynamics
C.tht_0  = 0.05;   % Resting θ (baseline deformation)
C.alpha_h  = 0.3;   % Threshold offset for distension activation

%% ---------------------------- Sigmoid Parameters -------------------------
% Parameters controlling the sigmoid function which characterizes the switching threshold
C.a_e      = 1.3;   % Gamma_E - slope of sigmoid for excitatory equation
C.a_i      = 2.0;   % Gamma_I - slope of sigmoid for inhibitory equation
C.tht_e  = 4.0;    % Threshold in sigmoid for excitatory neurons (location of sigmiod maximum slope)
C.tht_i  = 3.7;    % Threshold in sigmoid for inhibitory neurons (location of sigmiod maximum slope)
C.g        = 1000;  % Gain of spatial weight function
C.xs       = 0.1;   % Horizontal shift of spatial weight function
C.g_tht    = 5;     % Gain for θ sigmoid

%% -------------------------- Initial Conditions ---------------------------
C.E_init = 0.0;     % Initial excitatory activation
C.I_init = 0.0;     % Initial inhibitory activation

%% ------------------------- Simulation Parameters -------------------------
C.L        = 1.0;    % Total length of spatial domain (normalized)
C.duration = 200;    % Total simulation time
C.nxs      = 70;     % Number of spatial grid points
C.nts      = 1000;   % Number of temporal grid points

%% -------------------------- Mechanical Parameters ------------------------
C.S_IC = 2;          % Initial cross-sectional area (baseline)
C.BETA = 100;        % Non-dim fluid resistance parameter
C.PSI  = 1000;        % Non-dim stiffness parameter


end