# Li–Na Selectivity Model Using an LFP Electrode

This repository contains MATLAB codes used to analyze competitive Li/Na intercalation in lithium iron phosphate (LFP)/iron phosphate electrodes. The codes combine a competitive Frumkin thermodynamic model, Butler–Volmer interfacial kinetics, ion transport, and solid-state diffusion to examine the thermodynamic, kinetic, and transport origins of Li/Na selectivity.

The repository includes four components:

1. **Frumkin fitting** – fits single-salt electrochemical data to obtain thermodynamic parameters.
2. **Thermodynamic analysis** – predicts equilibrium Li/Na selectivity as a function of the Na/Li ratio and cross-interaction parameter.
3. **Flow-through simulation** – simulates competitive Li/Na intercalation under a spatially uniform, well-mixed flow-through approximation.
4. **Flow-by simulation** – resolves concentration, potential, reaction-rate, and solid-occupancy profiles through the electrode thickness.

## Repository structure

```text
.
├── Frumkin fitting/
│   └── Fit_frumkin.m
├── Thermodynamic analysis/
│   └── Thermo_analysis.m
├── Flow-through simulation/
│   ├── 0D_main.m
│   ├── 0D_electrochemical_equations.m
│   └── solve_theta_diffusion_all.m
└── Flow-by simulation/
    ├── main.m
    ├── solve_step_linear_explicit.m
    ├── solve_local_simple.m
    ├── compute_flux.m
    └── solve_theta_diffusion_all.m
```

> **Important:** The files in the uploaded archive are stored with the `.txt` extension. Rename all MATLAB source files from `.txt` to `.m` before running them. Function filenames must match the corresponding function names. In particular, rename `0D_electrochemical_equations.txt` to `0D_electrochemical_equations_reduced.m`, because the function defined inside the file is named `OD_electrochemical_equations_reduced` (letter **O**, not zero). Alternatively, rename both the function and its call in `0D_main.m` consistently.

## Requirements

- MATLAB
- Optimization Toolbox, used by `fsolve`
- Curve Fitting Toolbox, used by `fit` in the Frumkin-fitting code
- Symbolic Math Toolbox, used by `syms` and `vpasolve` in the flow-through initialization

The codes also write results to Microsoft Excel-compatible `.xlsx` files.

## Quick start

1. Clone or download this repository.
2. Rename the source files from `.txt` to `.m`.
3. Open MATLAB and set the current folder to the subfolder corresponding to the desired analysis.
4. Modify the parameter block at the beginning of the main script, if needed.
5. Run the main script for that subfolder.

Each simulation folder should remain self-contained so that MATLAB can locate its helper functions.

## 1. Frumkin fitting

### Purpose

`Fit_frumkin.m` fits electrochemical potential–composition data using a modified Frumkin-type expression:

\[
E = E_{\mathrm{ref}}
- V_T \ln\left(\frac{k\theta}{1-k\theta}\right)
- g_i k\theta,
\]

where \(E_{\mathrm{ref}}\) is the reference potential, \(g_i\) is the interaction parameter, \(k\) is a fixed site-competition coefficient, and \(\theta\) is the degree of intercalation.

### Required input

The script currently expects an input file named:

```text
Capacity_LFP_NaCl.txt
```

Place this file in the same folder as `Fit_frumkin.m`, or change the variable `fname` near the beginning of the script.

The input table must contain:

- a potential column, such as `E/V`, `Ewe`, `E`, `Voltage`, or `Ewe(V)`; and
- either an occupancy column such as `theta`, or a capacity column such as `capacity`, `Capacity`, `Q`, or `Capacity_mAh_g`.

If capacity is supplied, it is normalized by `C_max_manual` to obtain \(\theta\). The default value is 135 mAh g\(^{-1}\).

### Main user-defined parameters

- `fname` – input-data filename
- `C_max_manual` – capacity used to normalize the data
- `theta_min`, `theta_max` – fitting range
- `n` – number of retained data points
- `k_fixed` – fixed site-competition coefficient

### Output

The script displays the fitted curve and writes:

```text
Frumkin_Langmuir_Results.xlsx
```

The workbook contains the fitted data, fitted parameters, goodness-of-fit statistics, and model description.

## 2. Thermodynamic analysis

### Purpose

`Thermo_analysis.m` calculates equilibrium Li/Na uptake and selectivity using mass balances, a prescribed final electrode occupation, and equality of the Li and Na equilibrium potentials.

The script performs a parameter sweep over:

- the feed Na/Li concentration ratio (`ratio_list`); and
- the Li–Na cross-interaction parameter (`g_cross_list`).

For each condition, the nonlinear equilibrium equations are solved using `fsolve`.

### Main user-defined parameters

- `E_Li_ref`, `E_Na_ref` – reference potentials
- `g_Li`, `g_Na` – self-interaction parameters
- `g_cross_list` – cross-interaction parameters
- `k` – site-competition coefficient
- `ratio_list` – feed Na/Li ratios
- `target_theta` – prescribed final total occupation
- `C_sp_Na_feed` – initial Na concentration
- `V_elec`, `V_sp_feed` – electrode and solution volumes

### Output

The script plots selectivity against the feed concentration ratio and writes:

```text
selectivity_g_cross_ratio_analysis.xlsx
```

The workbook contains calculated selectivity and final Li and Na occupancies for each parameter combination.

## 3. Flow-through simulation

### Purpose

The flow-through model uses a spatially uniform, well-mixed approximation for the electrolyte within the cell. It couples:

- competitive Frumkin thermodynamics;
- Butler–Volmer kinetics for Li and Na;
- galvanostatic current partitioning;
- bulk ion mass balances; and
- radial solid-state diffusion inside spherical particles.

### Files

- `0D_main.m` – defines parameters, initializes the system, advances the solution in time, and plots the results
- `OD_electrochemical_equations_reduced.m` – defines the 15-equation nonlinear electrochemical system solved at each time step
- `solve_theta_diffusion_all.m` – updates radial Li and Na occupancy profiles in the particles

### Main user-defined parameters

- `NLR` – initial Na/Li concentration ratio
- `C_Na_0` – initial Na concentration
- `I` – applied current density; negative values represent adsorption/intercalation in the present sign convention
- `dt`, `t_max` – time-step size and total simulation time
- `C_max` – maximum intercalation-site concentration
- `E_Li_ref`, `E_Na_ref`, `g_Li`, `g_Na`, `g_cross`, `k` – thermodynamic parameters
- `K_c_Li`, `K_a_Li`, `K_c_Na`, `K_a_Na` – kinetic prefactors
- `r_p`, `D_in_Li`, `D_in_Na` – particle radius and solid-state diffusivities
- `V_sp`, `A_sp`, `L`, `P_IHC` – cell and electrode geometry

### Termination criterion

The simulation stops when

\[
C_{s,\mathrm{Li}} + kC_{s,\mathrm{Na}} > 0.95C_{\max}.
\]

### Output

The script generates figures showing:

- bulk Li and Na concentrations;
- particle-surface concentrations;
- equilibrium potentials;
- overpotentials;
- Li/Na selectivity;
- reaction rates; and
- radial Li and Na occupancy profiles.

The current version displays results but does not export them automatically. Add a `save`, `writetable`, or `exportgraphics` command if permanent output files are needed.

## 4. Flow-by simulation

### Purpose

The flow-by model resolves variations through the electrode thickness. It couples:

- Nernst–Planck transport of Li\(^+\), Na\(^+\), and Cl\(^-\) in electrode macropores;
- ionic and electronic potential distributions;
- local competitive Frumkin thermodynamics;
- local Butler–Volmer kinetics;
- a finite external solution/spacer reservoir; and
- radial solid-state diffusion within particles.

The left electrode boundary uses a zero-flux condition, while the right boundary exchanges ions with the well-mixed spacer reservoir.

### Files

- `main.m` – parameter definition, initialization, time integration, data storage, and Excel export
- `solve_step_linear_explicit.m` – advances the coupled transport and reaction model by one time step
- `solve_local_simple.m` – solves the local galvanostatic electrochemical problem and reconstructs phase-potential and current distributions
- `compute_flux.m` – calculates Li, Na, and Cl macropore fluxes
- `solve_theta_diffusion_all.m` – updates radial Li and Na occupancy profiles at each axial location

### Main user-defined parameters

- `NLR`, `C_Na_0` – feed composition
- `I` – applied current density; negative values represent adsorption/intercalation
- `dt`, `t_max`, `n_save` – integration and output intervals
- `N_x`, `L` – axial grid and electrode thickness
- `P_mA`, `P_IHC` – macropore and active-material volume fractions
- `D_Li`, `D_Na`, `D_Cl` – aqueous diffusion coefficients
- `sigma_s` – solid-phase conductivity
- `C_max`, `E_Li_ref`, `E_Na_ref`, `g_Li`, `g_Na`, `g_cross`, `k` – thermodynamic parameters
- `K_c_Li`, `K_a_Li`, `K_c_Na`, `K_a_Na` – kinetic prefactors
- `r_p`, `D_in_Li`, `D_in_Na`, `n_sub` – particle-scale transport parameters
- `A_sp`, `V_sp`, `L_sp` – spacer/reservoir geometry

### Output

The script writes an Excel workbook named according to the simulated conditions:

```text
simulation_NLR<value>_I<value>.xlsx
```

The workbook contains:

- spacer concentrations;
- axial macropore concentration profiles;
- particle-surface and particle-averaged occupancies;
- local Li and Na currents;
- local reaction rates;
- overpotential and phase-potential profiles;
- equilibrium-potential profiles;
- global summary quantities;
- chamber-based and electrode-based selectivity; and
- final radial occupancy profiles at the left, middle, and right sides of the electrode.

## Selectivity definitions

The codes report concentration-normalized Li/Na selectivity. For solution depletion, the implemented form is

\[
S_{\mathrm{Li/Na}}
=
\frac{1-C_{\mathrm{Li}}/C_{\mathrm{Li},0}}
{1-C_{\mathrm{Na}}/C_{\mathrm{Na},0}}.
\]

The flow-by code also reports an electrode-based quantity calculated from changes in the mean Li and Na occupancies.

At very early times, Na removal or uptake can be close to zero, making the calculated selectivity numerically large or noisy. Interpret such values together with the absolute Li and Na uptake.

## Units and sign convention

Unless otherwise noted, the codes use SI units:

- concentration: mol m\(^{-3}\)
- length: m
- time: s
- diffusivity: m\(^2\) s\(^{-1}\)
- current density: A m\(^{-2}\)
- potential: V

The applied current `I` follows the convention used in the source code: negative current corresponds to ion adsorption/intercalation.

## Reproducing a parameter sweep

The current scripts run one flow-through or flow-by condition at a time. To reproduce a sweep over Na/Li ratio, current density, electrode thickness, or particle radius, place the relevant main script inside an outer loop, update the desired parameter before each run, and save each output using a unique filename.

Because the flow-by model may require many time steps, first test a new parameter set using a shorter `t_max`, a larger output interval `n_save`, or a coarser spatial grid. Numerical results should be checked for convergence with respect to `dt`, `N_x`, `N_r`, and `n_sub`.

## Notes and known limitations

- The parameter values in the scripts are example or manuscript-specific values and should be reviewed before applying the model to another electrode or operating condition.
- The flow-through model assumes a spatially uniform electrolyte concentration within the cell.
- The flow-by model uses a one-dimensional electrode description and spherical-particle diffusion.
- Activity coefficients, side reactions, mechanical effects, and detailed phase-boundary dynamics are not explicitly represented.
- Some numerical safeguards, such as concentration floors, potential-step clipping, and capacity limits, are included to improve robustness.
- The Frumkin fitting input data are not included in the current archive and must be supplied separately.

## Citation

When using this code, please cite the associated publication:

```text
[Authors]. “Thermodynamic, kinetic, and transport origins of selectivity in
lithium electrosorption with intercalation electrodes.” [Journal, year, DOI].
```

Replace the placeholder above with the final bibliographic information once available.

## License

No license file is currently included. Add an open-source license, such as MIT or BSD-3-Clause, before redistributing or reusing the code beyond the terms permitted by the authors.

## Contact

For questions about the model or code, please open an issue in this repository or contact the corresponding author of the associated publication.
