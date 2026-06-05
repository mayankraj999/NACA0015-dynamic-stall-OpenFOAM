# NACA 0015 Dynamic Stall — OpenFOAM Comparative Study

Comparative study of three mesh motion methods for simulating dynamic stall
on a NACA 0015 oscillating airfoil using OpenFOAM.

## Methods Compared
- **interpolatingSolidBody** — whole-mesh rotation
- **displacementSBRStress** — mesh deformation with SBR stress solver
- **solidBody + AMI** — sliding mesh with Arbitrary Mesh Interface

## Flow Regimes
| Regime | Mach | Solver | ω (rad/s) | k* |
|--------|------|--------|-----------|-----|
| Incompressible | 0.05 | pimpleFoam | 6.28 | 0.18 |
| Compressible (k* matched) | 0.5 | rhoPimpleFoam | 62.8 | 0.18 |
| Compressible (ω same) | 0.5 | rhoPimpleFoam | 6.28 | 0.018 |

## Results

### Incompressible — 3 Method Comparison
![Cl vs alpha](plots/incompressible/01_Cl_vs_alpha.png)
![Cd vs alpha](plots/incompressible/02_Cd_vs_alpha.png)
![Cm vs alpha](plots/incompressible/03_Cm_vs_alpha.png)
![Time history](plots/incompressible/04_Cl_Cd_time_history.png)
![Drag polar](plots/incompressible/05_drag_polar.png)
![Single cycle](plots/incompressible/06_single_cycle_comparison.png)

### Compressible k* Matched — 3 Method Comparison
![Cl vs alpha](plots/compressible_k_matched/01_Cl_vs_alpha.png)
![Cd vs alpha](plots/compressible_k_matched/02_Cd_vs_alpha.png)
![Cm vs alpha](plots/compressible_k_matched/03_Cm_vs_alpha.png)
![Time history](plots/compressible_k_matched/04_Cl_Cd_time_history.png)
![Drag polar](plots/compressible_k_matched/05_drag_polar.png)
![Single cycle](plots/compressible_k_matched/06_single_cycle_comparison.png)

### Compressible ω Same — 3 Method Comparison
![Cl vs alpha](plots/compressible_omega_same/01_Cl_vs_alpha.png)
![Cd vs alpha](plots/compressible_omega_same/02_Cd_vs_alpha.png)
![Cm vs alpha](plots/compressible_omega_same/03_Cm_vs_alpha.png)
![Time history](plots/compressible_omega_same/04_Cl_Cd_time_history.png)
![Drag polar](plots/compressible_omega_same/05_drag_polar.png)
![Single cycle](plots/compressible_omega_same/06_single_cycle_comparison.png)

### Compressible vs Incompressible — k* Matched
![Per method](plots/comp_vs_incomp_k_matched/09_comp_vs_incomp_per_method.png)
![Overlay](plots/comp_vs_incomp_k_matched/10_comp_vs_incomp_overlay.png)

### Compressible vs Incompressible — ω Same
![Per method](plots/comp_vs_incomp_omega_same/09_comp_vs_incomp_per_method.png)
![Overlay](plots/comp_vs_incomp_omega_same/10_comp_vs_incomp_overlay.png)

### Pressure Field at α = 15°

**Incompressible**
| interpolatingSolidBody | displacementSBRStress | solidBody + AMI |
|---|---|---|
| ![](images/p/incompressible/incomp_interpolating_pressure_alpha15.png) | ![](images/p/incompressible/incomp_displacementSBRStress_pressure_alpha15.png) | ![](images/p/incompressible/incomp_solidBody_pressure_alpha15.png) |

**Compressible — k* Matched**
| interpolatingSolidBody | displacementSBRStress | solidBody + AMI |
|---|---|---|
| ![](images/p/compressible/comp_k_matched/comp_interpolating_k_matched_pressure_alpha15.png) | ![](images/p/compressible/comp_k_matched/comp_displacementSBRStress_k_matched_pressure_alpha15.png) | ![](images/p/compressible/comp_k_matched/comp_solidBody_k_matched_pressure_alpha15.png) |

**Compressible — ω Same**
| interpolatingSolidBody | displacementSBRStress | solidBody + AMI |
|---|---|---|
| ![](images/p/compressible/comp_omega_same/comp_interpolating_omega_same_pressure_alpha15.png) | ![](images/p/compressible/comp_omega_same/comp_displacementSBRStress_omega_same_pressure_alpha15.png) | ![](images/p/compressible/comp_omega_same/comp_solidBody_omega_same_pressure_alpha15.png) |

### Velocity Field at α = 15°

**Incompressible**
| interpolatingSolidBody | displacementSBRStress | solidBody + AMI |
|---|---|---|
| ![](images/U/incompressible/incompressible_interpolating_velocity_alpha15.png) | ![](images/U/incompressible/incompressible_displacementSBRStress_velocity_alpha15.png) | ![](images/U/incompressible/incompressible_solidBody_velocity_alpha15.png) |

**Compressible — k* Matched**
| interpolatingSolidBody | displacementSBRStress | solidBody + AMI |
|---|---|---|
| ![](images/U/compressible/comp_k_matched/compressible_interpolating_k_matched_velocity_alpha15.png) | ![](images/U/compressible/comp_k_matched/compressible_displacementSBRStress_k_matched_velocity_alpha15.png) | ![](images/U/compressible/comp_k_matched/compressible_solidBody_k_matched_velocity_alpha15.png) |

**Compressible — ω Same**
| interpolatingSolidBody | displacementSBRStress | solidBody + AMI |
|---|---|---|
| ![](images/U/compressible/comp_omega_same/compressible_interpolating_omega_same_velocity_alpha15.png) | ![](images/U/compressible/comp_omega_same/compressible_displacementSBRStress_omega_same_velocity_alpha15.png) | ![](images/U/compressible/comp_omega_same/compressible_solidBody_omega_same_velocity_alpha15.png) |

### Mesh
| interpolatingSolidBody | displacementSBRStress | solidBody + AMI |
|---|---|---|
| ![](images/meshes/interpolatingSolidBody_alpha15.png) | ![](images/meshes/displacementSBRStress_alpha15.png) | ![](images/meshes/solidBody_AMI_alpha15.png) |

## Summary Statistics

| Regime | Cl_max | Cd_mean | Cm_max | Loop Area |
|--------|--------|---------|--------|-----------|
| Incompressible | 1.79 – 1.83 | 0.24 – 0.25 | 0.42 – 0.46 | 13 – 22 |
| Compressible (k* matched) | 1.49 – 1.61 | 0.24 – 0.25 | 0.33 – 0.40 | 27 – 30 |
| Compressible (ω same) | 0.41 – 0.42 | 0.12 | 0.04 – 0.06 | 6.1 – 6.3 |

Detailed statistics for each case are in the [stats/](stats/) folder.

## Key Findings
- Compressibility reduces peak lift by **10–18%** at matched reduced frequency (k* = 0.18)
- Pitching moment drops by **10–24%** under compressible conditions
- Mesh motion method sensitivity **increases with Mach number** — incompressible Cl_max spread is 2%, compressible is 8%
- At same ω but different k* (omega_same case), the compressible flow is **quasi-steady** with no dynamic stall (k* = 0.018)
- All three methods show excellent symmetry (Cl_mean ≈ 0) confirming mesh quality

## Setup
- **Airfoil**: NACA 0015, chord = 1 m
- **Turbulence model**: k-ω SST
- **Mesh**: 6-block O-grid, ~35 cells wall-normal, grading ratio 50
- **Oscillation**: α = 0° ± 15°
- **OpenFOAM versions**: v2412 (ESI) and v6 (Foundation)

## Repository Structure
├── cases/
│   ├── incompressible/
│   │   ├── interpolatingSolidBody/
│   │   ├── displacementSBRStress/
│   │   └── solidBodyAMI/
│   └── compressible/
│       ├── k_matched/
│       │   ├── interpolatingSolidBody/
│       │   ├── displacementSBRStress/
│       │   └── solidBodyAMI/
│       └── omega_same/
│           ├── interpolatingSolidBody/
│           ├── displacementSBRStress/
│           └── solidBodyAMI/
├── plots/
│   ├── incompressible/
│   ├── compressible_k_matched/
│   ├── compressible_omega_same/
│   ├── comp_vs_incomp_k_matched/
│   └── comp_vs_incomp_omega_same/
├── stats/
├── images/
│   ├── p/
│   ├── U/
│   └── meshes/
├── scripts/
│   └── postprocess_unified.py
└── README.md

## Post-Processing Usage
```bash
# Incompressible 3-method comparison
python3 scripts/postprocess_unified.py ./case1 ./case2 ./case3

# Compressible k_matched (ω=62.8)
python3 scripts/postprocess_unified.py --omega 62.8 --discard 2 ./case1 ./case2 ./case3

# Compressible omega_same (ω=6.28)
python3 scripts/postprocess_unified.py ./case1 ./case2 ./case3

# Compressible vs Incompressible comparison
python3 scripts/postprocess_unified.py --compare \
    --incomp ./inc1 ./inc2 ./inc3 \
    --comp ./comp1 ./comp2 ./comp3 \
    --omega-inc 6.28 --omega-comp 62.8
```

## Author
Mayank Raj — Independent CFD study
READMEEOF
