# Li-Na Selectivity Model Using an LFP Electrode

This repository contains MATLAB codes used to investigate Li/Na selectivity in LFP-based intercalation electrodes. The codes are organized into four modules that progressively introduce thermodynamic, kinetic, solid-state diffusion, and electrode-scale transport effects.

The **Frumkin fitting** module uses a single-salt Frumkin model to fit galvanostatic charge-discharge data and determine the thermodynamic parameters for Li or Na intercalation.

The **Thermodynamic analysis** module applies a competitive Frumkin model to mixed Li/Na electrolytes. It calculates competitive Li and Na uptake and quantifies the intrinsic thermodynamic contribution to selectivity under equilibrium conditions.

The **Flow-through simulation** implements a zero-dimensional flow-through model in which the electrolyte composition within the electrode macropores is assumed to be spatially uniform. The model couples competitive Frumkin thermodynamics, Butler-Volmer interfacial kinetics, and radial solid-state diffusion within spherical particles. This model is used to evaluate interfacial kinetic and solid-state diffusion effects without explicitly resolving macropore concentration gradients.

The **Flow-by simulation** implements a one-dimensional model along the electrode thickness direction. It explicitly resolves ion transport within the electrode macropores while coupling competitive Frumkin thermodynamics, Butler-Volmer interfacial kinetics, and radial solid-state diffusion. This model is used to evaluate the combined thermodynamic, kinetic, solid-state diffusion, and electrode-scale transport origins of Li/Na selectivity.

## Repository structure

```text
.
├── Frumkin fitting/
│   └── Fit_frumkin.m
├── Thermodynamic analysis/
│   └── Thermo_analysis.m
├── Flow-through simulation/
│   ├── zeroD_main.m
│   ├── zeroD_electrochemical_equations.m
│   └── solve_theta_diffusion_all.m
└── Flow-by simulation/
    ├── main.m
    ├── solve_step_linear_explicit.m
    ├── solve_local_simple.m
    ├── compute_flux.m
    └── solve_theta_diffusion_all.m
```

## Requirements

- MATLAB
- Curve Fitting Toolbox, used by `fit` in the Frumkin fitting module 
- Optimization Toolbox, used by `lsqnonlin` in the Thermodynamic analysis module and used by `fsolve` in the flow-through simulation
- Symbolic Math Toolbox, used by `syms` and `vpasolve` during initialization of the flow-through model

The Flow-by simulation uses standard MATLAB functions and does not require an additional toolbox beyond MATLAB.

Some scripts export simulation results to Excel-compatible `.xlsx` files.

## Quick start

1. Clone or download this repository.
2. Open MATLAB.
3. Set the MATLAB current folder to the module that you want to run.
4. Modify the parameter definitions near the beginning of the corresponding main script, if needed.
5. Run the main script.

Each module should remain in its own folder so that MATLAB can locate the associated helper functions.

---

## 1. Frumkin fitting

`Fit_frumkin.m` fits single-salt potential-composition data using a modified Frumkin thermodynamic model to determine `E_ref` and `g_i` for Li or Na intercalation.

The fitted expression is:

    E = E_ref
        - V_T ln[(k_capacity theta)/(1 - k_capacity theta)]
        + V_T ln(c_ion/c_0)
        - g_i (k_capacity theta)

Use:

- `k_capacity = 1.00` for Li fitting
- `k_capacity = 1.63` for Na fitting

The input file must contain a potential column and either a `theta` column or a capacity column. When capacity data are provided, the script calculates:

    theta = capacity/C_max

The input filename is specified by `fname`. The script plots the fitted curve and exports the results to:

    Frumkin_Fitting_Results.xlsx

---

## 2. Thermodynamic analysis

`Thermo_analysis.m` calculates equilibrium Li/Na selectivity in mixed-salt electrolytes using a competitive Frumkin thermodynamic model.

For each initial Na/Li concentration ratio and cross-interaction parameter, the script solves the Li and Na mass balances together with equality of their equilibrium potentials at the prescribed final site occupancy:

    theta_Li + k_capacity theta_Na = target_theta

where:

    theta_Li = C_s_Li/C_max
    theta_Na = C_s_Na/C_max

The parameter sweep is defined by `NLR_list` and `g_cross_list`. The equilibrium problem is solved using `lsqnonlin`.

The script plots Li/Na selectivity as a function of the initial Na/Li concentration ratio and exports the results to:

    selectivity_g_cross_NLR_analysis.xlsx

---

## 3. Flow-through simulation

The flow-through module implements a zero-dimensional electrochemical model for a flow-through electrode. The electrolyte composition in the electrode macropores is assumed to be spatially uniform, so electrode-scale concentration and potential gradients are not explicitly resolved.

The model couples:

- competitive Frumkin thermodynamics;
- Butler-Volmer interfacial kinetics;
- Li/Na current partitioning; and
- radial solid-state diffusion in spherical particles.

Main files:

- `0D_main.m`: defines operating conditions, runs the simulation, and plots the results
- `OD_electrochemical_equations.m`: solves the electrochemical equations at each time step
- `solve_theta_diffusion_all.m`: solves radial Li and Na diffusion in the particles

The simulation stops when the weighted electrode occupation reaches the prescribed capacity limit:

    C_s,Li + k C_s,Na >= 0.95 C_max

The module reports electrolyte concentrations, particle occupations, equilibrium potentials, overpotentials, reaction rates, and Li/Na selectivity.

---

## 4. Flow-by simulation

The flow-by module implements a one-dimensional electrochemical model along the electrode thickness direction. Unlike the flow-through model, it explicitly resolves ion transport and potential distributions within the electrode macropores.

The model couples:

- one-dimensional Nernst-Planck transport of Li, Na, and Cl;
- electrolyte and solid-phase potential distributions;
- competitive Frumkin thermodynamics;
- Butler-Volmer interfacial kinetics; and
- radial solid-state diffusion at each axial location.

Main files:

- `main.m`: defines parameters, runs the simulation, plots results, and exports data
- `solve_step_linear_explicit.m`: advances the coupled transport-reaction model
- `solve_local_simple.m`: solves the local electrochemical problem
- `compute_flux.m`: calculates macropore ion fluxes
- `solve_theta_diffusion_all.m`: solves radial particle diffusion

The script exports the simulation results to:

    simulation_NLR<value>_I<value>.xlsx

The output includes electrolyte concentration profiles, particle occupations, reaction rates, potential distributions, and Li/Na selectivity.

---

## Selectivity definitions

The electrolyte-depletion-based selectivity is calculated as:

    S_Li/Na = (1 - C_Li/C_Li,0) / (1 - C_Na/C_Na,0)

The flow-by model also calculates an electrode-uptake-based selectivity from the changes in average Li and Na occupation.

Very early-time selectivity values may be unstable when Na uptake is close to zero.

## Units and sign convention

Unless otherwise specified, SI units are used:

- concentration: mol/m^3
- length: m
- time: s
- diffusivity: m^2/s
- current density: A/m^2
- potential: V

Negative applied current corresponds to ion adsorption or intercalation.

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
