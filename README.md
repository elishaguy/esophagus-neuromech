# Esophagus Neuromechanical Model (1D Wilson–Cowan)

**Author:** Guy Elishsa  
**Short description:** MATLAB implementation of a 1-D neuromechanical model of peristalsis using coupled Wilson–Cowan oscillators and a fluid-filled flexible tube model. Includes versions with and without sensory feedback.

---

## Contents
- `main_with_feedback.m` — Main script including distension → mechanoreceptor feedback (if present).
- `main_no_feedback.m` — Main script (simplified version where constant input value is applied `S_E = C.P`, `S_I = C.Q`).
- `simulation_constants.m` — All model parameters (constants).
- `Results1D_W_C/` — (gitignored) folder where simulation outputs are saved. 
- `README.md` — This file.

---
## Output format
Results includes a madlab data file with the solution of E, I, muscle activation function, tube area, and fluid velocity.
Simulation output is saved in the following path: Results1D_W_C/case_##/results.mat

## Requirements
- MATLAB (R2018b or later recommended)  
- No additional toolboxes required for basic execution (but check `simulation_constants.m` and scripts for dependencies).

---
## How the model is organized

- u vector layout: [A, velocity, E, I, theta] stacked in space.
- Neural coupling uses Wilson–Cowan equations; mechanical part discretizes area/velocity.

Two modes:
- Constant sensory input (S_E = C.P, S_I = C.Q) — simplified.
- Distension-dependent sensory input — computes h(x) = max(A/θ - α̂, 0) and integrates with kernels to form S_E, S_I.

Parameters:
See simulation_constants.m for parameter definitions and short comments on biophysical meaning (e.g., C.P, C.Q, C.alpha_h, ...).

---

## See additional details and cite the following
(1) Elisha, G., Halder, S., Liu, X., Carlson, D.A., Kahrilas, P.J., Pandolfino, J.E. and Patankar, N.A., 2024. Neurological disorders leading to mechanical dysfunction of the esophagus: an emergent behavior of a neuromechanical dynamical system. arXiv preprint arXiv:2402.18103.

(2) Elisha, G., Gast, R., Halder, S., Solla, S.A., Kahrilas, P.J., Pandolfino, J.E. and Patankar, N.A., 2025. Direct and Retrograde Wave Propagation in Unidirectionally Coupled Wilson-Cowan Oscillators. Physical review letters, 134(5), p.058401.

### For input on mechanical parameters and simulations
(1) Elisha, G., Acharya, S., Halder, S., Carlson, D.A., Kou, W., Kahrilas, P.J., Pandolfino, J.E. and Patankar, N.A., 2023. Peristaltic regimes in esophageal transport. Biomechanics and modeling in mechanobiology, 22(1), pp.23-41.

(2) Acharya, S., Kou, W., Halder, S., Carlson, D.A., Kahrilas, P.J., Pandolfino, J.E. and Patankar, N.A., 2021. Pumping patterns and work done during peristalsis in finite-length elastic tubes. Journal of Biomechanical Engineering, 143(7), p.071001.