import os
import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import streamlit as st

# ON/OFF: write PNG, SVG, and CSV copies of the current figure.
ENABLE_SAVE_PLOT_FILES = True
ENABLE_SAVE_MWD_FILES = True


def _overlay_series(ax, compare_runs, x_key, y_key):
    if not compare_runs:
        return
    for run in compare_runs:
        ax.plot(run[x_key], run[y_key], linestyle="--", label=run.get("label", "saved"))


def plot_results(
    t_values,
    R_pol,
    ef,
    R_data,
    thiele_data,
    polymer_mass,
    selected_Ca,
    radial_positions,
    AC_Conc,
    Cas,
    axis_locked=False,
    save_plots=False,
    save_csv=False,
    plots_dir="plots",
    csv_dir="plot_data",
    compare_runs=None,
    selected_Ca2=None,
    Cas2=None,
):
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))

    axs[0, 0].plot(t_values, R_pol, label="current")
    _overlay_series(axs[0, 0], compare_runs, "t", "R_pol")
    axs[0, 0].set_xlabel('Time (s)')
    axs[0, 0].set_ylabel('Polymerization Rate (gr pol/gr cat.hr)')
    axs[0, 0].set_title('Polymerization Rate')

    axs[0, 1].plot(t_values, ef, label="current")
    _overlay_series(axs[0, 1], compare_runs, "t", "ef")
    axs[0, 1].set_xlabel('Time (s)')
    axs[0, 1].set_ylabel('Efficiency')
    axs[0, 1].set_title('Efficiency')

    axs[0, 2].plot(t_values, R_data, label="current")
    _overlay_series(axs[0, 2], compare_runs, "t", "R_data")
    axs[0, 2].set_xlabel('Time (s)')
    axs[0, 2].set_ylabel('Particle Radius (m)')
    axs[0, 2].set_title('Particle Radius')

    axs[1, 0].plot(t_values, thiele_data, label="current")
    _overlay_series(axs[1, 0], compare_runs, "t", "thiele_data")
    axs[1, 0].set_xlabel('Time (s)')
    axs[1, 0].set_ylabel('Thiele Modulus')
    axs[1, 0].set_title('Thiele Modulus')

    axs[1, 1].plot(t_values, polymer_mass, label="current")
    _overlay_series(axs[1, 1], compare_runs, "t", "polymer_mass")
    axs[1, 1].set_xlabel('Time (s)')
    axs[1, 1].set_ylabel('Cumulative polymer mass (grams)')
    axs[1, 1].set_title('Polymerization yield')

    axs[1, 2].plot(radial_positions, selected_Ca, label="Monomer 1")
    if selected_Ca2 is not None:
        axs[1, 2].plot(radial_positions, selected_Ca2, label="Monomer 2")
    if compare_runs:
        for run in compare_runs:
            axs[1, 2].plot(run["radial_positions"], run["selected_Ca"], linestyle="--", label=run.get("label", "saved"))
    axs[1, 2].set_xlabel('Radial Position (m)')
    axs[1, 2].set_ylabel('Monomer Concentration (mol/m³)')
    axs[1, 2].set_title('Monomer Concentration')
    cas_hi = Cas if np.isfinite(Cas) and Cas > 0 else 1.0
    if Cas2 is not None:
        cas_hi = max(cas_hi, Cas2)
    axs[1, 2].set_ylim([0, cas_hi])
    if compare_runs or selected_Ca2 is not None:
        axs[0, 0].legend(fontsize=8)
        axs[1, 2].legend(fontsize=8)

    if axis_locked and 'y_limits' in st.session_state:
        y_limits = st.session_state.y_limits
        axs[0, 0].set_ylim(y_limits[0])
        axs[0, 1].set_ylim(y_limits[1])
        axs[0, 2].set_ylim(y_limits[2])
        axs[1, 0].set_ylim(y_limits[3])
        axs[1, 1].set_ylim(y_limits[4])

    plt.tight_layout()
    st.pyplot(fig)

    png_path = os.path.join(plots_dir, "simulation_results.png")
    svg_path = os.path.join(plots_dir, "simulation_results.svg")
    timeseries_csv_path = os.path.join(csv_dir, "timeseries.csv")
    concentration_csv_path = os.path.join(csv_dir, "concentration_profile.csv")

    if save_plots or ENABLE_SAVE_PLOT_FILES:
        os.makedirs(plots_dir, exist_ok=True)
        fig.savefig(png_path, format="png", dpi=200, bbox_inches="tight")
        fig.savefig(svg_path, format="svg", bbox_inches="tight")

    if save_csv or ENABLE_SAVE_PLOT_FILES:
        os.makedirs(csv_dir, exist_ok=True)
        pd.DataFrame(
            {
                "Time (s)": t_values,
                "Polymerization Rate": R_pol,
                "Efficiency": ef,
                "Particle Radius (m)": R_data,
                "Thiele Modulus": thiele_data,
                "Cumulative polymer mass (grams)": polymer_mass,
                "Active Sites Conc. (mol/m3)": AC_Conc,
            }
        ).to_csv(timeseries_csv_path, index=False)
        conc_out = {
            "Radial Position (m)": radial_positions,
            "Monomer Concentration (mol/m3)": selected_Ca,
        }
        if selected_Ca2 is not None:
            conc_out["Monomer 2 Concentration (mol/m3)"] = selected_Ca2
        pd.DataFrame(conc_out).to_csv(concentration_csv_path, index=False)

    return fig, png_path, svg_path


def plot_mwd(n, M, w_frac, Mn, Mw, PDI):
    fig, axs = plt.subplots(1, 2, figsize=(12, 4))
    axs[0].plot(n, w_frac)
    axs[0].set_xlabel('Chain length n')
    axs[0].set_ylabel('Weight fraction')
    axs[0].set_title('MWD (chain length)')
    axs[1].plot(np.log10(np.maximum(M, 1e-30)), w_frac)
    axs[1].set_xlabel('log10 M (g/mol)')
    axs[1].set_ylabel('Weight fraction')
    axs[1].set_title('MWD (log M)')
    fig.suptitle(f'Mn = {Mn:.4e}   Mw = {Mw:.4e}   PDI = {PDI:.3f}')
    plt.tight_layout()
    if ENABLE_SAVE_MWD_FILES:
        os.makedirs('plot_data', exist_ok=True)
        fig.savefig('plot_data/mwd.png', dpi=150)
        fig.savefig('plot_data/mwd.svg')
        with open('plot_data/mwd.csv', 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['n', 'M_g_per_mol', 'weight_fraction', 'Mn', 'Mw', 'PDI'])
            for i in range(len(n)):
                writer.writerow([n[i], M[i], w_frac[i], Mn, Mw, PDI])
    st.pyplot(fig)
