# Li-Na Selectivity Model Using an LFP Electrode

This repository contains MATLAB codes used to investigate Li/Na selectivity in LFP-based intercalation electrodes. The codes are organized into four modules that progressively introduce thermodynamic, kinetic, solid-state diffusion, and electrode-scale transport effects.

The **Frumkin fitting** module uses a single-salt Frumkin model to fit galvanostatic charge-discharge data and determine the thermodynamic parameters for Li or Na intercalation.

The **Thermodynamic analysis** module applies a competitive Frumkin model to mixed Li/Na electrolytes. It calculates competitive Li and Na uptake and quantifies the intrinsic thermodynamic contribution to selectivity under equilibrium conditions.

The **Flow-through simulation** implements a zero-dimensional flow-through model in which the electrolyte composition within the electrode macropores is assumed to be spatially uniform. The model couples competitive Frumkin thermodynamics, Butler-Volmer interfacial kinetics, and radial solid-state diffusion within spherical particles. This model is used to evaluate interfacial kinetic and solid-state diffusion effects without explicitly resolving macropore concentration gradients.

The **Flow-by simulation** implements a one-dimensional model along the electrode thickness direction. It explicitly resolves ion transport and electrolyte-potential distributions within the electrode macropores while coupling competitive Frumkin thermodynamics, Butler-Volmer interfacial kinetics, electronic conduction, and radial solid-state diffusion. This model is used to evaluate the combined thermodynamic, kinetic, solid-state diffusion, and electrode-scale transport origins of Li/Na selectivity.

## Repository structure

```text
.
├── Frumkin fitting/
│   └── Fit_frumkin.m
├── Thermodynamic analysis/
│   └── Thermo_analysis.m
├── Flow-through simulation/
│   ├── 0D_main.m
│   ├── OD_electrochemical_equations_reduced.m
│   └── solve_theta_diffusion_all.m
└── Flow-by simulation/
    ├── main.m
    ├── solve_step_linear_explicit.m
    ├── solve_local_simple.m
    ├── compute_flux.m
    └── solve_theta_diffusion_all.m
```

> **Important:** The source files in the provided archive are stored with the `.txt` extension. Rename all MATLAB source files from `.txt` to `.m` before running the codes.

> MATLAB function names must begin with a letter, and the filename of a primary function must match the function name. The file currently named `0D_electrochemical_equations.txt` contains the function `OD_electrochemical_equations_reduced`, where `OD` begins with the letter O rather than the numeral zero. The file should therefore be renamed as `OD_electrochemical_equations_reduced.m`.

## Requirements

- MATLAB
- Optimization Toolbox, used by `fsolve`
- Curve Fitting Toolbox, used by `fit` in the Frumkin fitting module
- Symbolic Math Toolbox, used by `syms` and `vpasolve` during initialization of the flow-through model

Some scripts export simulation results to Excel-compatible `.xlsx` files.

## Quick start

1. Clone or download this repository.
2. Rename all MATLAB source files from `.txt` to `.m`.
3. Open MATLAB.
4. Set the MATLAB current folder to the module that you want to run.
5. Modify the parameter definitions near the beginning of the corresponding main script, if needed.
6. Run the main script.

Each module should remain in its own folder so that MATLAB can locate the associated helper functions.

---

## 1. Frumkin fitting

### Purpose

`Fit_frumkin.m` fits single-salt electrochemical potential-composition data using a modified Frumkin expression. This module is applied separately to Li-containing and Na-containing single-salt electrolytes to determine the thermodynamic parameters for each ion.

The fitted expression is:

```text
E = E_ref - V_T * ln[(k*theta)/(1-k*theta)] - g_i*k*theta
```

where:

- `E_ref` is the reference intercalation potential.
- `V_T = RT/F` is the thermal voltage.
- `g_i` is the self-interaction parameter for the intercalating ion.
- `k` is the site-occupation coefficient.
- `theta` is the normalized degree of intercalation.

The fitted single-ion thermodynamic parameters are subsequently used in the competitive Li/Na models.

### Required input

The script currently expects an input file named:

```text
Capacity_LFP_NaCl.txt
```

The data file should be placed in the same folder as `Fit_frumkin.m`. Alternatively, the variable `fname` near the beginning of the script can be changed to specify another filename.

The input table must contain:

- a potential column, such as `E/V`, `Ewe`, `E`, `Voltage`, or `Ewe(V)`; and
- either an occupancy column, such as `theta`, or a capacity column, such as `capacity`, `Capacity`, `Q`, or `Capacity_mAh_g`.

When capacity data are supplied, the script normalizes the capacity by `C_max_manual` to calculate `theta`. The default value of `C_max_manual` is 135 mAh/g.

### Main user-defined parameters

- `fname`: input-data filename
- `C_max_manual`: capacity used to normalize the experimental data
- `theta_min`, `theta_max`: lower and upper limits of the occupancy range included in the fitting
- `n`: number of retained data points
- `k_fixed`: fixed site-occupation coefficient

### Output

The script plots the experimental and fitted potential-composition curves and writes the fitting results to:

```text
Frumkin_Langmuir_Results.xlsx
```

The workbook contains:

- processed fitting data;
- fitted thermodynamic parameters;
- goodness-of-fit statistics; and
- a description of the fitted model.

---

## 2. Thermodynamic analysis

### Purpose

`Thermo_analysis.m` applies a competitive Frumkin thermodynamic model to mixed Li/Na electrolytes.

Unlike the single-salt fitting module, this model accounts for simultaneous Li and Na occupation of the intercalation sites. Competitive uptake is calculated by combining:

- Li and Na mass balances;
- the prescribed final electrode occupation;
- equality of the Li and Na equilibrium potentials; and
- self- and cross-interaction terms in the competitive Frumkin model.

The model is used to quantify the intrinsic equilibrium selectivity arising from the thermodynamic preference of the electrode for Li over Na.

### Parameter sweep

The script performs a parameter sweep over:

- the initial Na/Li concentration ratio, defined by `ratio_list`; and
- the Li-Na cross-interaction parameter, defined by `g_cross_list`.

For each condition, the nonlinear equilibrium equations are solved using `fsolve`.

### Main user-defined parameters

- `E_Li_ref`, `E_Na_ref`: reference intercalation potentials for Li and Na
- `g_Li`, `g_Na`: Li-Li and Na-Na self-interaction parameters
- `g_cross_list`: Li-Na cross-interaction parameters included in the parameter sweep
- `k`: site-occupation coefficient for Na
- `ratio_list`: initial Na/Li concentration ratios
- `target_theta`: prescribed final total electrode occupation
- `C_sp_Na_feed`: initial Na concentration in the feed solution
- `V_elec`, `V_sp_feed`: electrode and solution volumes used in the mass balances

### Output

The script plots Li/Na selectivity as a function of the initial Na/Li concentration ratio and writes the calculated results to:

```text
selectivity_g_cross_ratio_analysis.xlsx
```

The workbook contains:

- calculated Li/Na selectivity;
- final Li occupation;
- final Na occupation; and
- results for each combination of Na/Li ratio and cross-interaction parameter.

---

## 3. Flow-through simulation

### Purpose

The flow-through module implements a zero-dimensional electrochemical model for a flow-through electrode.

The electrolyte composition within the electrode macropores is assumed to be spatially uniform. Therefore, concentration and potential gradients are not resolved along the electrode thickness direction. This approximation represents rapid convective replenishment in the flow-through configuration and removes explicit electrode-scale macropore transport limitations from the model.

The model couples:

- competitive Frumkin thermodynamics for Li/Na co-intercalation;
- Butler-Volmer interfacial kinetics for Li and Na;
- galvanostatic current partitioning between Li and Na reactions;
- spatially uniform electrolyte mass balances; and
- radial solid-state diffusion inside spherical active-material particles.

The flow-through model is used to isolate the effects of interfacial reaction kinetics and solid-state diffusion on Li/Na selectivity while assuming negligible concentration gradients within the electrode macropores.

### Files

- `0D_main.m`: defines the operating conditions and model parameters, initializes the system, advances the solution in time, stores the calculated variables, and plots the results.
- `OD_electrochemical_equations_reduced.m`: defines the nonlinear algebraic electrochemical equations solved at each time step, including current partitioning, Butler-Volmer kinetics, equilibrium potentials, and overpotentials.
- `solve_theta_diffusion_all.m`: solves radial Li and Na diffusion within spherical particles and updates the particle concentration profiles.

### Main user-defined parameters

- `NLR`: initial Na/Li concentration ratio
- `C_Na_0`: initial Na concentration
- `I`: applied current density; negative values represent adsorption or intercalation under the sign convention used in the code
- `dt`, `t_max`: time-step size and total simulation time
- `C_max`: maximum intercalation-site concentration
- `E_Li_ref`, `E_Na_ref`: reference intercalation potentials
- `g_Li`, `g_Na`, `g_cross`, `k`: competitive Frumkin thermodynamic parameters
- `K_c_Li`, `K_a_Li`, `K_c_Na`, `K_a_Na`: cathodic and anodic kinetic prefactors for Li and Na
- `r_p`: active-material particle radius
- `D_in_Li`, `D_in_Na`: Li and Na solid-state diffusion coefficients
- `V_sp`, `A_sp`, `L`, `P_IHC`: solution volume and electrode geometry parameters

### Termination criterion

The simulation terminates when the weighted electrode occupation approaches the prescribed capacity limit:

```text
C_s,Li + k*C_s,Na > 0.95*C_max
```

### Output

The script generates figures showing quantities including:

- spatially uniform electrolyte Li and Na concentrations;
- particle-surface Li and Na concentrations;
- Li and Na equilibrium potentials;
- Li and Na overpotentials;
- Li/Na selectivity;
- Li and Na reaction rates; and
- radial Li and Na occupation profiles within the particles.

The current flow-through script displays the calculated results but does not automatically export all simulation variables. Permanent outputs can be added using MATLAB functions such as `save`, `writetable`, or `exportgraphics`.

---

## 4. Flow-by simulation

### Purpose

The flow-by module implements a one-dimensional, spatially resolved model along the electrode thickness direction.

In contrast to the zero-dimensional flow-through model, the flow-by model explicitly resolves electrolyte transport within the electrode macropores. It therefore captures the development of concentration and potential gradients from the electrode-spacer interface into the electrode interior.

The model couples:

- one-dimensional Nernst-Planck transport of Li, Na, and Cl within the electrode macropores;
- electrolyte-phase potential distribution;
- electronic conduction and solid-phase potential distribution;
- local competitive Frumkin thermodynamics;
- local Butler-Volmer interfacial kinetics;
- local current partitioning between Li and Na reactions;
- exchange with a finite, well-mixed spacer reservoir; and
- radial solid-state diffusion within spherical active-material particles at each axial location.

This model is used to evaluate the combined influence of thermodynamics, interfacial kinetics, solid-state diffusion, and electrode-scale macropore transport on Li/Na selectivity.

### Boundary conditions

The electrode is resolved from the interior current-collector side to the electrode-spacer interface.

The model applies:

- a zero ionic-flux condition at the closed electrode boundary; and
- ion exchange between the electrode macropores and the well-mixed spacer reservoir at the electrode-spacer boundary.

### Files

- `main.m`: defines model parameters, initializes the concentration and occupation fields, performs time integration, stores simulation results, generates figures, and exports results to Excel.
- `solve_step_linear_explicit.m`: advances the coupled macropore transport and reaction model by one time step.
- `solve_local_simple.m`: solves the local electrochemical problem and reconstructs the solid-phase and electrolyte-phase potential and current distributions.
- `compute_flux.m`: calculates Li, Na, and Cl fluxes within the electrode macropores.
- `solve_theta_diffusion_all.m`: updates radial Li and Na occupation profiles at every axial position in the electrode.

### Main user-defined parameters

- `NLR`, `C_Na_0`: initial electrolyte composition
- `I`: applied current density; negative values represent adsorption or intercalation under the sign convention used in the code
- `dt`, `t_max`, `n_save`: time integration and output intervals
- `N_x`, `L`: number of axial grid points and electrode thickness
- `P_mA`, `P_IHC`: macropore and active-material volume fractions
- `D_Li`, `D_Na`, `D_Cl`: aqueous diffusion coefficients of Li, Na, and Cl
- `sigma_s`: effective solid-phase electronic conductivity
- `C_max`: maximum intercalation-site concentration
- `E_Li_ref`, `E_Na_ref`: reference intercalation potentials
- `g_Li`, `g_Na`, `g_cross`, `k`: competitive Frumkin thermodynamic parameters
- `K_c_Li`, `K_a_Li`, `K_c_Na`, `K_a_Na`: Li and Na kinetic prefactors
- `r_p`, `D_in_Li`, `D_in_Na`, `n_sub`: particle radius, solid-state diffusion coefficients, and particle-diffusion substep parameters
- `A_sp`, `V_sp`, `L_sp`: spacer or external reservoir geometry

### Output

The script writes an Excel workbook whose filename is determined by the simulated Na/Li ratio and current density:

```text
simulation_NLR<value>_I<value>.xlsx
```

The workbook contains quantities including:

- spacer-reservoir Li and Na concentrations;
- axial macropore concentration profiles;
- particle-surface Li and Na occupations;
- particle-averaged Li and Na occupations;
- local Li and Na current densities;
- local Li and Na reaction rates;
- overpotential profiles;
- solid-phase potential profiles;
- electrolyte-phase potential profiles;
- Li and Na equilibrium-potential profiles;
- global electrode-average quantities;
- solution-depletion-based selectivity;
- electrode-uptake-based selectivity; and
- final radial occupation profiles at selected axial positions.

---

## Selectivity definitions

For selectivity calculated from electrolyte depletion, the implemented expression is:

```text
S_Li/Na = (1 - C_Li/C_Li,0) / (1 - C_Na/C_Na,0)
```

This definition compares the fractional removal of Li with the fractional removal of Na.

The flow-by model also calculates an electrode-based selectivity using changes in the average Li and Na occupations of the active material.

At very early simulation times, Na uptake may be close to zero. The denominator of the selectivity expression can therefore become very small, producing numerically large or noisy selectivity values. Selectivity should consequently be interpreted together with the absolute Li and Na uptake and the simulation time.

## Units and sign convention

Unless otherwise specified, the codes use SI units:

- concentration: mol/m^3
- length: m
- time: s
- diffusivity: m^2/s
- current density: A/m^2
- potential: V

The applied current variable `I` follows the sign convention used in the source codes. Negative current corresponds to ion adsorption or intercalation.

## Reproducing parameter sweeps

The current flow-through and flow-by main scripts simulate one operating condition at a time.

To reproduce a parameter sweep over quantities such as:

- initial Na/Li concentration ratio;
- applied current density;
- electrode thickness;
- macropore volume fraction;
- particle radius; or
- solid-state diffusion coefficient,

place the relevant main script inside an outer loop, update the desired parameter before each simulation, and save each result using a unique filename.

Because the one-dimensional flow-by model may require a large number of time steps, a new parameter set should first be tested using:

- a shorter `t_max`;
- a larger output interval `n_save`; or
- a coarser spatial grid.

The numerical solution should be checked for convergence with respect to:

- time-step size `dt`;
- axial grid size `N_x`;
- radial particle grid size `N_r`; and
- particle-diffusion substep number `n_sub`.

## Notes and model limitations

- The parameter values supplied in the scripts correspond to example or manuscript-specific simulation conditions and should be reviewed before applying the model to another electrode material or operating condition.
- The Frumkin fitting module uses single-salt experimental data and does not account for direct Li-Na competition.
- The thermodynamic analysis module assumes equilibrium and does not include Butler-Volmer kinetics or transport limitations.
- The flow-through model assumes that the electrolyte composition is spatially uniform within the electrode macropores. Electrode-scale macropore concentration and potential gradients are therefore not explicitly represented.
- The flow-by model uses a one-dimensional description along the electrode thickness direction.
- Active-material particles are represented as spheres with radial solid-state diffusion.
- Activity-coefficient corrections, parasitic reactions, mechanical deformation, and detailed moving phase-boundary dynamics are not explicitly included.
- Numerical safeguards, including concentration floors, capacity limits, and limits on selected potential updates, are included to improve numerical robustness.
- The experimental data file required by the Frumkin fitting module is not included in the current code archive and must be supplied separately.

## Citation

When using these codes, please cite the associated publication:

```text
[Authors]. "Thermodynamic, kinetic, and transport origins of selectivity in lithium electrosorption with intercalation electrodes." [Journal, year, DOI].
```

Replace this placeholder with the final bibliographic information once the associated publication is available.

## License

No license file is currently included in the repository. An appropriate open-source license, such as the MIT License or BSD 3-Clause License, should be added before permitting redistribution or reuse beyond the terms specified by the authors.

## Contact

For questions regarding the model or source codes, please open an issue in this repository or contact the corresponding author of the associated publication.
