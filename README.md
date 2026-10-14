
# Li-Na Selectivity Model Using an LFP Electrode

This repository contains MATLAB codes used to investigate Li/Na selectivity in LFP-based intercalation electrodes. The codes are organized into four modules that progressively introduce thermodynamic, kinetic, solid-state diffusion, and electrode-scale transport effects.

The **Frumkin fitting** module uses a single-salt Frumkin model to fit galvanostatic charge-discharge data and determine the thermodynamic parameters for Li or Na intercalation.

The **Thermodynamic analysis** module applies a competitive Frumkin model to mixed Li/Na electrolytes. It calculates competitive Li and Na uptake and quantifies the intrinsic thermodynamic contribution to selectivity under equilibrium conditions.

The **Flow-through simulation** implements a zero-dimensional flow-through model in which the electrolyte composition within the electrode macropores is assumed to be spatially uniform. The model couples competitive Frumkin thermodynamics, Butler-Volmer interfacial kinetics, and radial solid-state diffusion within spherical particles. This model is used to evaluate interfacial kinetic and solid-state diffusion effects without explicitly resolving macropore concentration gradients. Li/Na selectivity is calculated from the net amounts intercalated into the solid phase.

The **Flow-by simulation** implements a one-dimensional model along the electrode thickness direction. It explicitly resolves ion transport within the electrode macropores while coupling competitive Frumkin thermodynamics, Butler-Volmer interfacial kinetics, and radial solid-state diffusion. This model is used to evaluate the combined thermodynamic, kinetic, solid-state diffusion, and electrode-scale transport origins of Li/Na selectivity. Selectivity is evaluated from the electrode-wide solid-phase Li and Na inventories.

## Repository structure

```text
.
└── Li-Na-selectivity-model-using-LFP-electrode/
    ├── Frumkin fitting/
    │   ├── Fit_frumkin.m
    │   └── Capacity_LFP_NaCl.txt
    ├── Thermodynamic analysis/
    │   └── Thermo_analysis.m
    ├── Flow-through simulation/
    │   ├── zeroD_main.m
    │   ├── zeroD_electrochemical_equations.m
    │   └── solve_theta_diffusion_all.m
    ├── Flow-by simulation/
    │   ├── main.m
    │   ├── solve_step_linear_explicit.m
    │   ├── solve_local_simple.m
    │   ├── compute_flux.m
    │   └── solve_theta_diffusion_all.m
    └── README.md
```

## Requirements

- MATLAB
- Curve Fitting Toolbox, used by `fit` in the Frumkin fitting module
- Optimization Toolbox, used by `lsqnonlin` in the Thermodynamic analysis module and by `fsolve` in the flow-through simulation
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
        - V_T ln[(k_capacity * theta)/(1 - k_capacity * theta)]
        + V_T ln(c_ion/c_0)
        - g_i (k_capacity * theta)

Use:

- `k_capacity = 1.00` for Li fitting
- `k_capacity = 1.63` for Na fitting

The input file must contain a potential column and either a `theta` column or a capacity column. When capacity data are provided, the script calculates:

    theta = capacity/C_max

The script currently expects an input file named:

    Capacity_LFP_NaCl.txt

The script plots the fitted curve and exports the results to:

    Frumkin_Fitting_Results.xlsx

---

## 2. Thermodynamic analysis

`Thermo_analysis.m` calculates equilibrium Li/Na selectivity in mixed-salt electrolytes using a competitive Frumkin thermodynamic model.

For each initial Na/Li concentration ratio and cross-interaction parameter, the script solves the Li and Na mass balances together with equality of their equilibrium potentials at the prescribed final site occupancy:

    theta_Li + k_capacity * theta_Na = target_theta

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

- `zeroD_main.m`: defines operating conditions, runs the simulation, plots the results, and exports the calculated data
- `zeroD_electrochemical_equations.m`: solves the electrochemical equations at each time step
- `solve_theta_diffusion_all.m`: solves radial Li and Na diffusion in the particles using a conservative spherical finite-volume discretization

Li/Na selectivity is calculated from the changes in particle-volume-averaged Li and Na occupancies relative to their initial values, rather than from the decrease in bulk electrolyte concentrations.

The simulation stops when the weighted electrode occupation reaches the prescribed capacity limit:

    C_s,Li + k C_s,Na >= 0.95 C_max

The script generates figures and exports the results to an Excel file named:

    Results_Li<value>_Na<value>_<timestamp>.xlsx

The workbook contains time-series results, concentration histories, Li and Na particle profiles, and model parameters.

---

## 4. Flow-by simulation

The flow-by module implements a one-dimensional electrochemical model along the electrode thickness direction. Unlike the flow-through model, it explicitly resolves ion transport and potential distributions within the electrode macropores.

The model couples:

- one-dimensional Nernst-Planck transport of Li, Na, and Cl;
- electrolyte- and solid-phase potential distributions;
- competitive Frumkin thermodynamics;
- Butler-Volmer interfacial kinetics and Li/Na current partitioning; and
- radial solid-state diffusion at each axial location.

The electrochemical solver determines the local reaction currents and potential distributions under a prescribed total applied current. Local reaction availability is constrained by Li depletion and weighted site occupancy. These constraints are incorporated into the electrochemical solution rather than applied as corrections after the current distribution has been calculated.

The model uses explicit time integration for macropore concentrations, with adaptive time-step subdivision when a trial update produces negative Li or Na concentrations. Nernst-Planck fluxes are evaluated at control-volume faces, with the electrolyte-potential gradient reconstructed from the local face current to satisfy ionic current closure. Radial solid-state diffusion uses a conservative spherical finite-volume discretization.

Local electroneutrality is imposed for Li+/Na+/Cl- in the electrode macropores and the well-mixed spacer. The spacer Cl concentration is set by electroneutrality; an external charge-compensation process is assumed but not explicitly resolved. Consequently, the spacer Cl inventory is not obtained solely from the electrode-spacer Cl flux.

Main files:

- `main.m`: defines parameters, runs the simulation, and exports data
- `solve_step_linear_explicit.m`: advances the coupled transport-reaction model
- `solve_local_simple.m`: solves the coupled electrochemical potential and current distributions under the applied-current constraint
- `compute_flux.m`: calculates Nernst-Planck ion fluxes and enforces ionic current closure at computational faces
- `solve_theta_diffusion_all.m`: solves radial particle diffusion using a conservative spherical finite-volume discretization

The simulation terminates when at least 99% of the axial locations reach the prescribed surface-occupancy threshold:

    theta_Li,surface + k * theta_Na,surface >= 0.95

The simulation is performed under continuous constant-current operation.

The script exports the simulation results to:

    simulation_NLR<value>_I<value>.xlsx

The output includes electrolyte concentration profiles, particle occupations, local reaction currents, reaction rates, electrolyte- and solid-phase potential distributions, and Li/Na selectivity calculated from net solid-phase intercalation. Particle occupancies are volume-averaged radially and then averaged across the electrode thickness to determine the total solid-phase inventories.

---

## Selectivity definition

Li/Na selectivity is defined using the net amounts of Li and Na intercalated into the electrode active material:

    S_Li/Na(t) = [Delta_n_Li,solid(t) / Delta_n_Na,solid(t)] * NLR
    NLR = C_Na,0 / C_Li,0

Here, `Delta_n_i,solid(t)` is the change in the total solid-phase inventory of ion `i` relative to its initial value. Both Li and Na inventories are calculated using the same active-material volume and maximum site concentration.

For spherical particles, the volume-averaged occupancy is:

    theta_i,avg(t) = 3 * integral_0^1 [theta_i(r,t) * r^2] dr

The radial diffusion solver uses finite-volume shell weights to evaluate this average numerically. In the flow-through model, the selectivity is therefore:

    S_Li/Na(t) = [(theta_Li,avg(t) - theta_Li,0) /
                  (theta_Na,avg(t) - theta_Na,0)] * NLR

In the flow-by model, the particle-volume-averaged occupancies are additionally averaged across the electrode thickness to determine the net Li and Na inventories of the entire electrode. The Na capacity-weighting factor `k` enters the thermodynamic and site-occupancy relations; it is **not** an extra multiplier in the ratio of intercalated Li and Na moles.

This definition represents **net solid-phase intercalation selectivity**. It excludes Li and Na remaining in the electrode macropore electrolyte, which can contribute to apparent selectivity measured without complete rinsing. It is not necessarily equal to the equilibrium thermodynamic selectivity under transport-limited operation.

Very early-time selectivity values may be undefined or numerically unstable when net Na intercalation is close to zero.

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
- The thermodynamic analysis module assumes equilibrium and does not include Butler-Volmer kinetics or transport limitations.
- The flow-through model assumes that the electrolyte composition is spatially uniform within the electrode macropores. Electrode-scale macropore concentration and potential gradients are therefore not explicitly represented.
- The flow-by model uses a one-dimensional description along the electrode thickness direction.
- In the flow-by model, reaction availability is constrained by prescribed Li-depletion and site-occupancy thresholds. These thresholds are numerical modeling assumptions rather than independently derived kinetic laws.
- Active-material particles are represented as spheres with radial solid-state diffusion, discretized using a conservative finite-volume method.
- Activity-coefficient corrections, parasitic reactions, mechanical deformation, and detailed moving phase-boundary dynamics are not explicitly included.
- The flow-by spacer model imposes electroneutrality through an assumed external charge-compensation process; it does not explicitly model the counter-electrode or establish an independent Cl mass balance for the spacer.
- The flow-by current-closure calculation does not by itself verify consistency of all electrolyte potential and concentration fields; coupled numerical validation remains necessary.
- The flow-through code retains post-solve reaction-rate limiting and a fallback after unsuccessful nonlinear solves. These safeguards can affect mass-balance consistency and should be checked when interpreting quantitative results.
- Numerical safeguards include concentration floors, occupancy limits, bounded potential updates, and adaptive time-step subdivision following negative-concentration errors. Numerical convergence and Li/Na inventory balances should be verified for the selected operating conditions.
- The experimental data file required by the Frumkin fitting module is not included in the current code archive and must be supplied separately.

## Citation

When using these codes, please cite the associated publication:

```text
[Authors]. "Thermodynamic, kinetic, and transport origins of selectivity in lithium electrosorption with intercalation electrodes." [Journal, year, DOI].
```

Replace this placeholder with the final bibliographic information once the associated publication is available.

## License

No license file is currently included in this repository. Therefore, the code may be viewed, but reuse, modification, or redistribution is not permitted unless explicit permission is granted by the authors.

## Contact

For questions regarding the model or source codes, please open an issue in this repository or contact the corresponding author of the associated publication.
