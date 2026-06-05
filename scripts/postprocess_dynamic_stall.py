#!/usr/bin/env python3
"""
NACA 0015 Dynamic Stall — Post-Processing Suite (Unified)
==========================================================
Generates all key plots for the comparative study of mesh motion methods.
Supports incompressible, compressible, and cross-comparison modes.

Usage:
    # Incompressible (ω=6.28, default)
    python3 postprocess_dynamic_stall.py ./interp ./disp ./ami

    # Compressible k_matched (ω=62.8)
    python3 postprocess_dynamic_stall.py --omega 62.8 ./interp_km ./disp_km ./ami_km

    # Compressible omega_same (ω=6.28, same as default)
    python3 postprocess_dynamic_stall.py ./interp_os ./disp_os ./ami_os

    # Compressible vs Incompressible (per-method comparison)
    python3 postprocess_dynamic_stall.py --compare \\
        --incomp ./interp_inc ./disp_inc ./ami_inc \\
        --comp ./interp_comp ./disp_comp ./ami_comp \\
        --omega-inc 6.28 --omega-comp 62.8

    # Grid convergence
    python3 postprocess_dynamic_stall.py --grid ./coarse ./medium ./fine

    # Time scheme comparison
    python3 postprocess_dynamic_stall.py --timescheme ./euler ./backward

Each case directory must contain:
    postProcessing/forceCoeffs/0/coefficient.dat  (or forceCoeffs.dat)
"""

import sys
import os
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

# ══════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════

AMPLITUDE = 15.0
ALPHA_MEAN = 0.0
AREF_CORRECTION = 1.0 / 0.2  # = 5.0

DEFAULT_LABELS = [
    "interpolatingSolidBody",
    "displacementSBRStress",
    "solidBody + AMI"
]

COLORS = ['#2563eb', '#dc2626', '#16a34a']
LINEWIDTHS = [2.0, 1.8, 1.8]
LINESTYLES = ['-', '--', '-.']

OUTPUT_DIR = "./plots"

# ══════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ══════════════════════════════════════════════════════════════════

def compute_alpha(t, omega):
    return ALPHA_MEAN + AMPLITUDE * np.sin(omega * t)


def get_period(omega):
    return 2.0 * np.pi / omega


def read_force_coeffs(case_dir):
    """
    Read OpenFOAM forceCoeffs output.
    Auto-detects OF6 (6-col), OF7+ (7-col), and v2412 ESI (13-col) formats.
    """
    candidates = [
        os.path.join(case_dir, "postProcessing", "forceCoeffs", "0", "coefficient.dat"),
        os.path.join(case_dir, "postProcessing", "forceCoeffs", "0", "forceCoeffs.dat"),
    ]

    pp_dir = os.path.join(case_dir, "postProcessing", "forceCoeffs")
    if os.path.isdir(pp_dir):
        for d in sorted(os.listdir(pp_dir)):
            for fname in ["coefficient.dat", "forceCoeffs.dat"]:
                candidates.append(os.path.join(pp_dir, d, fname))

    filepath = None
    for c in candidates:
        if os.path.isfile(c):
            filepath = c
            break

    if filepath is None:
        raise FileNotFoundError(
            f"No forceCoeffs data found in {case_dir}.\n"
            f"Searched: {candidates[:2]}\n"
            f"Make sure you ran with the forceCoeffs function object."
        )

    print(f"  Reading: {filepath}")

    data_lines = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#') and not line.startswith('Time'):
                data_lines.append(line)

    if len(data_lines) == 0:
        raise ValueError(f"No data rows found in {filepath}")

    rows = []
    for line in data_lines:
        parts = line.split()
        parts = [p.strip('()') for p in parts]
        try:
            rows.append([float(x) for x in parts])
        except ValueError:
            continue

    data = np.array(rows)
    ncols = data.shape[1]
    result = {'t': data[:, 0]}

    if ncols >= 13:
        # v2412 ESI: Time Cd Cd(f) Cd(r) Cl Cl(f) Cl(r) CmPitch CmRoll CmYaw Cs Cs(f) Cs(r)
        result['Cd']      = data[:, 1] * AREF_CORRECTION
        result['Cl']      = data[:, 4] * AREF_CORRECTION
        result['CmPitch'] = data[:, 7] * AREF_CORRECTION
        result['CmRoll']  = data[:, 8] * AREF_CORRECTION
        result['CmYaw']   = data[:, 9] * AREF_CORRECTION
        result['Cs']      = data[:, 10] * AREF_CORRECTION
    elif ncols == 6:
        # OF6 foundation: Time Cm Cd Cl Cl(f) Cl(r)
        result['CmPitch'] = data[:, 1] * AREF_CORRECTION
        result['Cd']      = data[:, 2] * AREF_CORRECTION
        result['Cl']      = data[:, 3] * AREF_CORRECTION
        result['Cs']      = np.zeros_like(data[:, 0])
        result['CmRoll']  = np.zeros_like(data[:, 0])
        result['CmYaw']   = np.zeros_like(data[:, 0])
    elif ncols >= 7:
        # OF7+ foundation: Time Cd Cs Cl CmRoll CmPitch CmYaw
        result['Cd']      = data[:, 1] * AREF_CORRECTION
        result['Cs']      = data[:, 2] * AREF_CORRECTION
        result['Cl']      = data[:, 3] * AREF_CORRECTION
        result['CmRoll']  = data[:, 4] * AREF_CORRECTION
        result['CmPitch'] = data[:, 5] * AREF_CORRECTION
        result['CmYaw']   = data[:, 6] * AREF_CORRECTION
    elif ncols >= 4:
        result['Cd']      = data[:, 1] * AREF_CORRECTION
        result['Cl']      = data[:, 2] * AREF_CORRECTION
        result['CmPitch'] = data[:, 3] * AREF_CORRECTION
        result['Cs']      = np.zeros_like(data[:, 0])
        result['CmRoll']  = np.zeros_like(data[:, 0])
        result['CmYaw']   = np.zeros_like(data[:, 0])
    else:
        raise ValueError(f"Unexpected number of columns ({ncols}) in {filepath}")

    print(f"    Detected {ncols} columns")
    return result


def trim_transient(data, omega, discard_cycles):
    t_start = discard_cycles * get_period(omega)
    mask = data['t'] >= t_start
    return {k: v[mask] for k, v in data.items()}


def setup_plot_style():
    plt.rcParams.update({
        'figure.figsize': (10, 7),
        'figure.dpi': 150,
        'font.size': 13,
        'font.family': 'serif',
        'axes.linewidth': 1.2,
        'axes.grid': True,
        'grid.alpha': 0.3,
        'grid.linewidth': 0.8,
        'lines.linewidth': 2.0,
        'legend.fontsize': 11,
        'legend.framealpha': 0.9,
        'legend.edgecolor': '0.8',
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'xtick.major.size': 5,
        'ytick.major.size': 5,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.15,
    })


# ══════════════════════════════════════════════════════════════════
# STANDARD PLOT FUNCTIONS
# ══════════════════════════════════════════════════════════════════

def plot_cl_vs_alpha(all_data, labels, omega, output_dir):
    fig, ax = plt.subplots()
    for i, (data, label) in enumerate(zip(all_data, labels)):
        alpha = compute_alpha(data['t'], omega)
        ax.plot(alpha, data['Cl'],
                color=COLORS[i % len(COLORS)],
                linewidth=LINEWIDTHS[i % len(LINEWIDTHS)],
                linestyle=LINESTYLES[i % len(LINESTYLES)],
                label=label, alpha=0.85)
    ax.set_xlabel(r'Angle of Attack $\alpha$ [°]')
    ax.set_ylabel(r'Lift Coefficient $C_L$')
    ax.set_title(f'Dynamic Stall — Lift Coefficient Hysteresis Loop\n'
                 f'NACA 0015, $\\alpha = {ALPHA_MEAN}° \\pm {AMPLITUDE}°$, '
                 f'$\\omega = {omega}$ rad/s')
    ax.legend(loc='best')
    ax.axhline(y=0, color='k', linewidth=0.5, alpha=0.3)
    ax.axvline(x=0, color='k', linewidth=0.5, alpha=0.3)
    filepath = os.path.join(output_dir, "01_Cl_vs_alpha.png")
    fig.savefig(filepath)
    plt.close(fig)
    print(f"  Saved: {filepath}")


def plot_cd_vs_alpha(all_data, labels, omega, output_dir):
    fig, ax = plt.subplots()
    for i, (data, label) in enumerate(zip(all_data, labels)):
        alpha = compute_alpha(data['t'], omega)
        ax.plot(alpha, data['Cd'],
                color=COLORS[i % len(COLORS)],
                linewidth=LINEWIDTHS[i % len(LINEWIDTHS)],
                linestyle=LINESTYLES[i % len(LINESTYLES)],
                label=label, alpha=0.85)
    ax.set_xlabel(r'Angle of Attack $\alpha$ [°]')
    ax.set_ylabel(r'Drag Coefficient $C_D$')
    ax.set_title(f'Dynamic Stall — Drag Coefficient Hysteresis Loop\n'
                 f'NACA 0015, $\\alpha = {ALPHA_MEAN}° \\pm {AMPLITUDE}°$, '
                 f'$\\omega = {omega}$ rad/s')
    ax.legend(loc='best')
    filepath = os.path.join(output_dir, "02_Cd_vs_alpha.png")
    fig.savefig(filepath)
    plt.close(fig)
    print(f"  Saved: {filepath}")


def plot_cm_vs_alpha(all_data, labels, omega, output_dir):
    fig, ax = plt.subplots()
    for i, (data, label) in enumerate(zip(all_data, labels)):
        alpha = compute_alpha(data['t'], omega)
        ax.plot(alpha, data['CmPitch'],
                color=COLORS[i % len(COLORS)],
                linewidth=LINEWIDTHS[i % len(LINEWIDTHS)],
                linestyle=LINESTYLES[i % len(LINESTYLES)],
                label=label, alpha=0.85)
    ax.set_xlabel(r'Angle of Attack $\alpha$ [°]')
    ax.set_ylabel(r'Pitching Moment Coefficient $C_M$')
    ax.set_title(f'Dynamic Stall — Moment Coefficient Hysteresis Loop\n'
                 f'NACA 0015, $\\alpha = {ALPHA_MEAN}° \\pm {AMPLITUDE}°$, '
                 f'$\\omega = {omega}$ rad/s')
    ax.legend(loc='best')
    ax.axhline(y=0, color='k', linewidth=0.5, alpha=0.3)
    filepath = os.path.join(output_dir, "03_Cm_vs_alpha.png")
    fig.savefig(filepath)
    plt.close(fig)
    print(f"  Saved: {filepath}")


def plot_cl_cd_vs_time(all_data, labels, omega, discard_cycles, output_dir):
    period = get_period(omega)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    for i, (data, label) in enumerate(zip(all_data, labels)):
        ax1.plot(data['t'], data['Cl'],
                 color=COLORS[i % len(COLORS)], linewidth=1.5,
                 linestyle=LINESTYLES[i % len(LINESTYLES)],
                 label=label, alpha=0.8)
        ax2.plot(data['t'], data['Cd'],
                 color=COLORS[i % len(COLORS)], linewidth=1.5,
                 linestyle=LINESTYLES[i % len(LINESTYLES)],
                 label=label, alpha=0.8)
    t_max = max(d['t'][-1] for d in all_data)
    for cycle in range(int(t_max / period) + 2):
        t_cycle = cycle * period
        ax1.axvline(x=t_cycle, color='gray', linewidth=0.5, alpha=0.3)
        ax2.axvline(x=t_cycle, color='gray', linewidth=0.5, alpha=0.3)
    t_discard = discard_cycles * period
    ax1.axvspan(0, t_discard, alpha=0.08, color='red', label='Transient (discarded)')
    ax2.axvspan(0, t_discard, alpha=0.08, color='red')
    ax1.set_ylabel(r'$C_L$')
    ax1.set_title('Time History of Aerodynamic Coefficients\nNACA 0015 Oscillating Airfoil')
    ax1.legend(loc='upper right', fontsize=9)
    ax2.set_xlabel('Time [s]')
    ax2.set_ylabel(r'$C_D$')
    ax2.legend(loc='upper right', fontsize=9)
    plt.tight_layout()
    filepath = os.path.join(output_dir, "04_Cl_Cd_time_history.png")
    fig.savefig(filepath)
    plt.close(fig)
    print(f"  Saved: {filepath}")


def plot_cl_cd_polar(all_data, labels, output_dir):
    fig, ax = plt.subplots()
    for i, (data, label) in enumerate(zip(all_data, labels)):
        ax.plot(data['Cd'], data['Cl'],
                color=COLORS[i % len(COLORS)],
                linewidth=LINEWIDTHS[i % len(LINEWIDTHS)],
                linestyle=LINESTYLES[i % len(LINESTYLES)],
                label=label, alpha=0.85)
    ax.set_xlabel(r'Drag Coefficient $C_D$')
    ax.set_ylabel(r'Lift Coefficient $C_L$')
    ax.set_title('Drag Polar — Dynamic Stall\nNACA 0015')
    ax.legend(loc='best')
    filepath = os.path.join(output_dir, "05_drag_polar.png")
    fig.savefig(filepath)
    plt.close(fig)
    print(f"  Saved: {filepath}")


def plot_single_cycle_comparison(all_data, labels, omega, output_dir):
    period = get_period(omega)
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    coeffs = ['Cl', 'Cd', 'CmPitch']
    ylabels = [r'$C_L$', r'$C_D$', r'$C_M$']
    titles = ['Lift', 'Drag', 'Pitching Moment']
    for ax, coeff, ylabel, title in zip(axes, coeffs, ylabels, titles):
        for i, (data, label) in enumerate(zip(all_data, labels)):
            t_end = data['t'][-1]
            mask = data['t'] >= (t_end - period)
            alpha_cycle = compute_alpha(data['t'][mask], omega)
            ax.plot(alpha_cycle, data[coeff][mask],
                    color=COLORS[i % len(COLORS)],
                    linewidth=LINEWIDTHS[i % len(LINEWIDTHS)],
                    linestyle=LINESTYLES[i % len(LINESTYLES)],
                    label=label, alpha=0.85)
        ax.set_xlabel(r'$\alpha$ [°]')
        ax.set_ylabel(ylabel)
        ax.set_title(f'{title} Coefficient')
        ax.legend(fontsize=9)
        ax.axhline(y=0, color='k', linewidth=0.5, alpha=0.3)
        ax.grid(True, alpha=0.3)
    fig.suptitle(f'Last Complete Cycle — Three Mesh Motion Methods\n'
                 f'NACA 0015, $\\alpha = {ALPHA_MEAN}° \\pm {AMPLITUDE}°$, '
                 f'$\\omega = {omega}$ rad/s',
                 fontsize=14, y=1.02)
    plt.tight_layout()
    filepath = os.path.join(output_dir, "06_single_cycle_comparison.png")
    fig.savefig(filepath)
    plt.close(fig)
    print(f"  Saved: {filepath}")


def compute_and_print_stats(all_data, labels, omega, discard_cycles, output_dir):
    period = get_period(omega)
    filepath = os.path.join(output_dir, "stats_summary.txt")
    lines = []
    lines.append("=" * 72)
    lines.append("NACA 0015 Dynamic Stall — Summary Statistics")
    lines.append(f"Oscillation: α = {ALPHA_MEAN}° ± {AMPLITUDE}°, ω = {omega} rad/s")
    lines.append(f"Period = {period:.4f} s, Transient discarded = {discard_cycles} cycle(s)")
    lines.append("=" * 72)
    lines.append("")
    for i, (data, label) in enumerate(zip(all_data, labels)):
        alpha = compute_alpha(data['t'], omega)
        lines.append(f"{'─' * 72}")
        lines.append(f"  Method: {label}")
        lines.append(f"{'─' * 72}")
        lines.append(f"  Time range:     {data['t'][0]:.4f} – {data['t'][-1]:.4f} s")
        lines.append(f"  Data points:    {len(data['t'])}")
        lines.append(f"  Cycles covered: {(data['t'][-1] - data['t'][0]) / period:.1f}")
        lines.append(f"")
        lines.append(f"  Cl  max = {data['Cl'].max():+.6f}   at α = {alpha[np.argmax(data['Cl'])]:.2f}°")
        lines.append(f"  Cl  min = {data['Cl'].min():+.6f}   at α = {alpha[np.argmin(data['Cl'])]:.2f}°")
        lines.append(f"  Cl  mean = {data['Cl'].mean():+.6f}")
        lines.append(f"")
        lines.append(f"  Cd  max = {data['Cd'].max():+.6f}   at α = {alpha[np.argmax(data['Cd'])]:.2f}°")
        lines.append(f"  Cd  min = {data['Cd'].min():+.6f}")
        lines.append(f"  Cd  mean = {data['Cd'].mean():+.6f}")
        lines.append(f"")
        lines.append(f"  Cm  max = {data['CmPitch'].max():+.6f}")
        lines.append(f"  Cm  min = {data['CmPitch'].min():+.6f}")
        lines.append(f"  Cm  mean = {data['CmPitch'].mean():+.6f}")
        t_end = data['t'][-1]
        mask = data['t'] >= (t_end - period)
        alpha_cycle = compute_alpha(data['t'][mask], omega)
        cl_cycle = data['Cl'][mask]
        area = 0.5 * np.abs(np.sum(
            alpha_cycle[:-1] * cl_cycle[1:] - alpha_cycle[1:] * cl_cycle[:-1]
        ))
        lines.append(f"")
        lines.append(f"  Cl-α loop area (last cycle) = {area:.4f} [°]")
        lines.append(f"")
    text = "\n".join(lines)
    with open(filepath, 'w') as f:
        f.write(text)
    print(f"\n{text}")
    print(f"\n  Saved: {filepath}")


# ══════════════════════════════════════════════════════════════════
# COMPRESSIBLE vs INCOMPRESSIBLE COMPARISON
# ══════════════════════════════════════════════════════════════════

def plot_comp_vs_incomp(inc_data, comp_data, method_labels, omega_inc, omega_comp, output_dir):
    """
    Per-method comparison: incompressible (solid) vs compressible (dashed).
    Generates 3×3 panel: rows = Cl, Cd, Cm; columns = each method.
    Uses last complete cycle from each.
    """
    n_methods = min(len(inc_data), len(comp_data))
    period_inc = get_period(omega_inc)
    period_comp = get_period(omega_comp)

    coeffs = ['Cl', 'Cd', 'CmPitch']
    ylabels = [r'$C_L$', r'$C_D$', r'$C_M$']
    row_titles = ['Lift', 'Drag', 'Pitching Moment']

    fig, axes = plt.subplots(3, n_methods, figsize=(6 * n_methods, 15), squeeze=False)

    for col in range(n_methods):
        # Incompressible last cycle
        t_end_i = inc_data[col]['t'][-1]
        mask_i = inc_data[col]['t'] >= (t_end_i - period_inc)
        alpha_i = compute_alpha(inc_data[col]['t'][mask_i], omega_inc)

        # Compressible last cycle
        t_end_c = comp_data[col]['t'][-1]
        mask_c = comp_data[col]['t'] >= (t_end_c - period_comp)
        alpha_c = compute_alpha(comp_data[col]['t'][mask_c], omega_comp)

        for row, (coeff, ylabel, rtitle) in enumerate(zip(coeffs, ylabels, row_titles)):
            ax = axes[row, col]

            ax.plot(alpha_i, inc_data[col][coeff][mask_i],
                    color='#2563eb', linewidth=2.0, linestyle='-',
                    label=f'Incomp (ω={omega_inc})', alpha=0.85)
            ax.plot(alpha_c, comp_data[col][coeff][mask_c],
                    color='#dc2626', linewidth=1.8, linestyle='--',
                    label=f'Comp (ω={omega_comp})', alpha=0.85)

            ax.set_xlabel(r'$\alpha$ [°]')
            ax.set_ylabel(ylabel)
            ax.axhline(y=0, color='k', linewidth=0.5, alpha=0.3)
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize=8)

            if row == 0:
                ax.set_title(method_labels[col], fontsize=13, fontweight='bold')

    fig.suptitle('Compressible vs Incompressible — Last Cycle Comparison\n'
                 f'NACA 0015, $\\alpha = {ALPHA_MEAN}° \\pm {AMPLITUDE}°$',
                 fontsize=15, y=1.01)
    plt.tight_layout()
    filepath = os.path.join(output_dir, "09_comp_vs_incomp_per_method.png")
    fig.savefig(filepath)
    plt.close(fig)
    print(f"  Saved: {filepath}")

    # Also make overlay plots: all 6 curves on one axis per coefficient
    fig2, axes2 = plt.subplots(1, 3, figsize=(20, 7))
    method_colors = ['#2563eb', '#dc2626', '#16a34a']

    for ax, coeff, ylabel, rtitle in zip(axes2, coeffs, ylabels, row_titles):
        for col in range(n_methods):
            t_end_i = inc_data[col]['t'][-1]
            mask_i = inc_data[col]['t'] >= (t_end_i - period_inc)
            alpha_i = compute_alpha(inc_data[col]['t'][mask_i], omega_inc)

            t_end_c = comp_data[col]['t'][-1]
            mask_c = comp_data[col]['t'] >= (t_end_c - period_comp)
            alpha_c = compute_alpha(comp_data[col]['t'][mask_c], omega_comp)

            ax.plot(alpha_i, inc_data[col][coeff][mask_i],
                    color=method_colors[col], linewidth=2.0, linestyle='-',
                    label=f'{method_labels[col]} (Incomp)', alpha=0.85)
            ax.plot(alpha_c, comp_data[col][coeff][mask_c],
                    color=method_colors[col], linewidth=1.8, linestyle='--',
                    label=f'{method_labels[col]} (Comp)', alpha=0.7)

        ax.set_xlabel(r'$\alpha$ [°]')
        ax.set_ylabel(ylabel)
        ax.set_title(f'{rtitle} Coefficient')
        ax.axhline(y=0, color='k', linewidth=0.5, alpha=0.3)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, ncol=2)

    fig2.suptitle('All Methods — Compressible (dashed) vs Incompressible (solid)\n'
                  f'NACA 0015, last cycle',
                  fontsize=15, y=1.01)
    plt.tight_layout()
    filepath2 = os.path.join(output_dir, "10_comp_vs_incomp_overlay.png")
    fig2.savefig(filepath2)
    plt.close(fig2)
    print(f"  Saved: {filepath2}")

    # Stats comparison table
    filepath3 = os.path.join(output_dir, "comp_vs_incomp_stats.txt")
    lines = []
    lines.append("=" * 80)
    lines.append("Compressible vs Incompressible — Summary (Last Cycle)")
    lines.append("=" * 80)
    lines.append(f"{'Method':<28} {'Type':<10} {'Cl_max':>10} {'Cl_min':>10} "
                 f"{'Cd_max':>10} {'Cd_mean':>10} {'Cm_max':>10} {'Cm_min':>10}")
    lines.append("─" * 80)

    for col in range(n_methods):
        for tag, data_list, omega_val, period_val in [
            ("Incomp", inc_data, omega_inc, period_inc),
            ("Comp", comp_data, omega_comp, period_comp),
        ]:
            d = data_list[col]
            t_end = d['t'][-1]
            mask = d['t'] >= (t_end - period_val)
            cl_c = d['Cl'][mask]
            cd_c = d['Cd'][mask]
            cm_c = d['CmPitch'][mask]
            lines.append(
                f"{method_labels[col]:<28} {tag:<10} "
                f"{cl_c.max():>+10.5f} {cl_c.min():>+10.5f} "
                f"{cd_c.max():>+10.5f} {cd_c.mean():>+10.5f} "
                f"{cm_c.max():>+10.5f} {cm_c.min():>+10.5f}"
            )
        lines.append("")

    text = "\n".join(lines)
    with open(filepath3, 'w') as f:
        f.write(text)
    print(f"\n{text}")
    print(f"\n  Saved: {filepath3}")


# ══════════════════════════════════════════════════════════════════
# GRID CONVERGENCE
# ══════════════════════════════════════════════════════════════════

def plot_grid_convergence(case_dirs, grid_labels, omega, discard_cycles, output_dir):
    period = get_period(omega)
    fig, ax = plt.subplots()
    peak_cls = []
    grid_colors = ['#f59e0b', '#2563eb', '#16a34a']
    for i, (case_dir, label) in enumerate(zip(case_dirs, grid_labels)):
        data = read_force_coeffs(case_dir)
        data = trim_transient(data, omega, discard_cycles)
        t_end = data['t'][-1]
        mask = data['t'] >= (t_end - period)
        alpha = compute_alpha(data['t'][mask], omega)
        ax.plot(alpha, data['Cl'][mask],
                color=grid_colors[i], linewidth=2.0 - 0.3*i,
                label=label, alpha=0.85)
        peak_cls.append(data['Cl'][mask].max())
    ax.set_xlabel(r'$\alpha$ [°]')
    ax.set_ylabel(r'$C_L$')
    ax.set_title('Grid Convergence Study — Cl Hysteresis\nNACA 0015')
    ax.legend()
    filepath = os.path.join(output_dir, "07_grid_convergence.png")
    fig.savefig(filepath)
    plt.close(fig)
    print(f"  Saved: {filepath}")
    if len(peak_cls) == 3:
        f1, f2, f3 = peak_cls[2], peak_cls[1], peak_cls[0]
        r = 2.0
        p_order = np.log(abs((f3 - f2) / (f2 - f1 + 1e-15))) / np.log(r)
        Fs = 1.25
        GCI_fine = Fs * abs((f1 - f2) / (f1 + 1e-15)) / (r**p_order - 1) * 100
        print(f"\n  Grid Convergence Index (GCI):")
        print(f"    Peak Cl: coarse={f3:.6f}, medium={f2:.6f}, fine={f1:.6f}")
        print(f"    Observed order of convergence p = {p_order:.2f}")
        print(f"    GCI (fine grid) = {GCI_fine:.2f}%")


# ══════════════════════════════════════════════════════════════════
# TIME SCHEME COMPARISON
# ══════════════════════════════════════════════════════════════════

def plot_time_scheme_comparison(euler_dir, backward_dir, omega, discard_cycles, output_dir):
    period = get_period(omega)
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for case_dir, label, color, ls in [
        (euler_dir, 'Euler (1st order)', '#dc2626', '--'),
        (backward_dir, 'backward (2nd order)', '#2563eb', '-'),
    ]:
        data = read_force_coeffs(case_dir)
        data = trim_transient(data, omega, discard_cycles)
        t_end = data['t'][-1]
        mask = data['t'] >= (t_end - period)
        alpha = compute_alpha(data['t'][mask], omega)
        axes[0].plot(alpha, data['Cl'][mask], color=color, linestyle=ls, linewidth=2, label=label)
        axes[1].plot(alpha, data['CmPitch'][mask], color=color, linestyle=ls, linewidth=2, label=label)
    axes[0].set_xlabel(r'$\alpha$ [°]'); axes[0].set_ylabel(r'$C_L$')
    axes[0].set_title('Lift Coefficient'); axes[0].legend()
    axes[1].set_xlabel(r'$\alpha$ [°]'); axes[1].set_ylabel(r'$C_M$')
    axes[1].set_title('Pitching Moment'); axes[1].legend()
    fig.suptitle('Effect of Time Discretization on Dynamic Stall\nNACA 0015', fontsize=14)
    plt.tight_layout()
    filepath = os.path.join(output_dir, "08_time_scheme_comparison.png")
    fig.savefig(filepath)
    plt.close(fig)
    print(f"  Saved: {filepath}")


# ══════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description='NACA 0015 Dynamic Stall Post-Processing Suite',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Incompressible (default ω=6.28)
  %(prog)s ./interp ./disp ./ami

  # Compressible k_matched (ω=62.8)
  %(prog)s --omega 62.8 --discard 2 ./interp_km ./disp_km ./ami_km

  # Compressible vs Incompressible
  %(prog)s --compare \\
      --incomp ./interp_inc ./disp_inc ./ami_inc \\
      --comp ./interp_comp ./disp_comp ./ami_comp \\
      --omega-inc 6.28 --omega-comp 62.8

  # Grid convergence
  %(prog)s --grid ./coarse ./medium ./fine

  # Time scheme comparison
  %(prog)s --timescheme ./euler ./backward
        """
    )

    parser.add_argument('cases', nargs='*', help='Case directories')
    parser.add_argument('--omega', type=float, default=6.28,
                        help='Oscillation frequency in rad/s (default: 6.28)')
    parser.add_argument('--discard', type=int, default=1,
                        help='Number of transient cycles to discard (default: 1)')
    parser.add_argument('--labels', type=str, default='',
                        help='Comma-separated labels for cases')
    parser.add_argument('--output', type=str, default='./plots',
                        help='Output directory (default: ./plots)')

    # Modes
    parser.add_argument('--grid', action='store_true', help='Grid convergence mode')
    parser.add_argument('--timescheme', action='store_true', help='Time scheme comparison mode')
    parser.add_argument('--compare', action='store_true',
                        help='Compressible vs Incompressible comparison mode')

    # Compare mode arguments
    parser.add_argument('--incomp', nargs='*', default=[],
                        help='Incompressible case directories (for --compare)')
    parser.add_argument('--comp', nargs='*', default=[],
                        help='Compressible case directories (for --compare)')
    parser.add_argument('--omega-inc', type=float, default=6.28,
                        help='ω for incompressible cases (default: 6.28)')
    parser.add_argument('--omega-comp', type=float, default=62.8,
                        help='ω for compressible cases (default: 62.8)')
    parser.add_argument('--discard-inc', type=int, default=1,
                        help='Transient cycles to discard for incompressible (default: 1)')
    parser.add_argument('--discard-comp', type=int, default=2,
                        help='Transient cycles to discard for compressible (default: 2)')
    parser.add_argument('--method-labels', type=str, default='',
                        help='Comma-separated method labels for compare mode')

    args = parser.parse_args()

    setup_plot_style()
    os.makedirs(args.output, exist_ok=True)

    # ── Compare mode ─────────────────────────────────────────────
    if args.compare:
        if not args.incomp or not args.comp:
            print("Error: --compare requires --incomp and --comp with case directories")
            print("Example:")
            print("  python3 postprocess_dynamic_stall.py --compare \\")
            print("      --incomp ./interp_inc ./disp_inc ./ami_inc \\")
            print("      --comp ./interp_comp ./disp_comp ./ami_comp")
            sys.exit(1)

        method_labels = (args.method_labels.split(',') if args.method_labels
                         else DEFAULT_LABELS[:min(len(args.incomp), len(args.comp))])

        print(f"\n{'='*60}")
        print("Compressible vs Incompressible Comparison")
        print(f"{'='*60}")
        print(f"Incompressible (ω={args.omega_inc}): {args.incomp}")
        print(f"Compressible   (ω={args.omega_comp}): {args.comp}")
        print()

        print("Reading incompressible cases...")
        inc_data = []
        for cd in args.incomp:
            d = read_force_coeffs(cd)
            inc_data.append(trim_transient(d, args.omega_inc, args.discard_inc))

        print("Reading compressible cases...")
        comp_data = []
        for cd in args.comp:
            d = read_force_coeffs(cd)
            comp_data.append(trim_transient(d, args.omega_comp, args.discard_comp))

        print(f"\nGenerating comparison plots → {args.output}/")
        plot_comp_vs_incomp(inc_data, comp_data, method_labels,
                            args.omega_inc, args.omega_comp, args.output)

        print(f"\n{'='*60}")
        print(f"Done! Outputs in {args.output}/")
        print(f"{'='*60}")
        return

    # ── Grid convergence mode ────────────────────────────────────
    if args.grid:
        if len(args.cases) < 2:
            print("Need at least 2 grid levels")
            sys.exit(1)
        grid_labels = [f"Grid {i+1}" for i in range(len(args.cases))]
        plot_grid_convergence(args.cases, grid_labels, args.omega, args.discard, args.output)
        return

    # ── Time scheme comparison ───────────────────────────────────
    if args.timescheme:
        if len(args.cases) < 2:
            print("Need: --timescheme <euler_case> <backward_case>")
            sys.exit(1)
        plot_time_scheme_comparison(args.cases[0], args.cases[1],
                                    args.omega, args.discard, args.output)
        return

    # ── Standard mode ────────────────────────────────────────────
    if not args.cases:
        parser.print_help()
        sys.exit(1)

    labels = args.labels.split(',') if args.labels else DEFAULT_LABELS[:len(args.cases)]

    print(f"\n{'='*60}")
    print("NACA 0015 Dynamic Stall Post-Processing")
    print(f"{'='*60}")
    print(f"ω = {args.omega} rad/s, Period = {get_period(args.omega):.4f} s")
    print(f"Discard {args.discard} cycle(s)")
    print(f"Cases: {len(args.cases)}")
    for cd, lb in zip(args.cases, labels):
        print(f"  [{lb}] → {cd}")
    print()

    print("Reading force coefficients...")
    all_data_raw = []
    for cd in args.cases:
        all_data_raw.append(read_force_coeffs(cd))

    print(f"\nDiscarding first {args.discard} cycle(s) as transient...")
    all_data = [trim_transient(d, args.omega, args.discard) for d in all_data_raw]

    print(f"\nGenerating plots → {args.output}/")
    print("─" * 40)

    plot_cl_vs_alpha(all_data, labels, args.omega, args.output)
    plot_cd_vs_alpha(all_data, labels, args.omega, args.output)
    plot_cm_vs_alpha(all_data, labels, args.omega, args.output)
    plot_cl_cd_vs_time(all_data_raw, labels, args.omega, args.discard, args.output)
    plot_cl_cd_polar(all_data, labels, args.output)
    plot_single_cycle_comparison(all_data, labels, args.omega, args.output)

    print(f"\n{'─' * 40}")
    print("Computing statistics...")
    compute_and_print_stats(all_data, labels, args.omega, args.discard, args.output)

    print(f"\n{'='*60}")
    print(f"Done! All outputs in {args.output}/")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
