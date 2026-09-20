import os
import csv
import matplotlib.pyplot as plt
import streamlit as st

# ON/OFF: write PNG, SVG, and CSV copies of the current figure.
ENABLE_SAVE_PLOT_FILES = True
PLOT_DATA_DIR = "plot_data"


def _overlay_series(ax, compare_runs, x_key, y_key):
    if not compare_runs:
        return
    for run in compare_runs:
        ax.plot(run[x_key], run[y_key], linestyle="--", label=run.get("label", "saved"))


def _write_plot_files(fig, t_values, R_pol, ef, R_data, thiele_data, polymer_mass, selected_Ca, radial_positions, compare_runs):
    os.makedirs(PLOT_DATA_DIR, exist_ok=True)
    fig.savefig(os.path.join(PLOT_DATA_DIR, "current_plots.png"), dpi=150)
    fig.savefig(os.path.join(PLOT_DATA_DIR, "current_plots.svg"))
    csv_path = os.path.join(PLOT_DATA_DIR, "current_plots.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["time_s", "polymerization_rate", "efficiency", "particle_radius_m", "thiele_modulus", "polymer_mass_g"])
        for i in range(len(t_values)):
            writer.writerow([t_values[i], R_pol[i], ef[i], R_data[i], thiele_data[i], polymer_mass[i]])
    conc_path = os.path.join(PLOT_DATA_DIR, "current_concentration.csv")
    with open(conc_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["radial_position", "monomer_concentration"])
        for x, y in zip(radial_positions, selected_Ca):
            writer.writerow([x, y])
    if compare_runs:
        cmp_path = os.path.join(PLOT_DATA_DIR, "compare_runs.csv")
        with open(cmp_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["label", "time_s", "polymerization_rate", "efficiency", "particle_radius_m", "thiele_modulus", "polymer_mass_g"])
            for run in compare_runs:
                n = len(run["t"])
                for i in range(n):
                    writer.writerow([
                        run.get("label", "saved"),
                        run["t"][i],
                        run["R_pol"][i],
                        run["ef"][i],
                        run["R_data"][i],
                        run["thiele_data"][i],
                        run["polymer_mass"][i],
                    ])


def plot_results(t_values, R_pol, ef, R_data, thiele_data, polymer_mass, selected_Ca, radial_positions, AC_Conc, Cas, axis_locked=False, compare_runs=None):
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))

    # Plot Polymerization Rate
    axs[0, 0].plot(t_values, R_pol, label="current")
    _overlay_series(axs[0, 0], compare_runs, "t", "R_pol")
    axs[0, 0].set_xlabel('Time (s)')
    axs[0, 0].set_ylabel('Polymerization Rate (gr pol/gr cat.hr)')
    axs[0, 0].set_title('Polymerization Rate')

    # Plot Efficiency
    axs[0, 1].plot(t_values, ef, label="current")
    _overlay_series(axs[0, 1], compare_runs, "t", "ef")
    axs[0, 1].set_xlabel('Time (s)')
    axs[0, 1].set_ylabel('Efficiency')
    axs[0, 1].set_title('Efficiency')

    # Plot Particle Radius
    axs[0, 2].plot(t_values, R_data, label="current")
    _overlay_series(axs[0, 2], compare_runs, "t", "R_data")
    axs[0, 2].set_xlabel('Time (s)')
    axs[0, 2].set_ylabel('Particle Radius (m)')
    axs[0, 2].set_title('Particle Radius')

    # Plot Thiele Modulus
    axs[1, 0].plot(t_values, thiele_data, label="current")
    _overlay_series(axs[1, 0], compare_runs, "t", "thiele_data")
    axs[1, 0].set_xlabel('Time (s)')
    axs[1, 0].set_ylabel('Thiele Modulus')
    axs[1, 0].set_title('Thiele Modulus')

    # Plot Polymer Mass
    axs[1, 1].plot(t_values, polymer_mass, label="current")
    _overlay_series(axs[1, 1], compare_runs, "t", "polymer_mass")
    axs[1, 1].set_xlabel('Time (s)')
    axs[1, 1].set_ylabel('Cumulative polymer mass (grams)')
    axs[1, 1].set_title('Polymerization yield')

    # Plot Monomer Concentration
    axs[1, 2].plot(radial_positions, selected_Ca, label="current")
    if compare_runs:
        for run in compare_runs:
            axs[1, 2].plot(run["radial_positions"], run["selected_Ca"], linestyle="--", label=run.get("label", "saved"))
    axs[1, 2].set_xlabel('Radial Position (m)')
    axs[1, 2].set_ylabel('Monomer Concentration (mol/m³)')
    axs[1, 2].set_title('Monomer Concentration')
    axs[1, 2].set_ylim([0, Cas])

    if compare_runs:
        axs[0, 0].legend(fontsize=8)

    # Set axis limits if locked
    if axis_locked and 'y_limits' in st.session_state:
        y_limits = st.session_state.y_limits
        
        axs[0, 0].set_ylim(y_limits[0])
        axs[0, 1].set_ylim(y_limits[1])
        axs[0, 2].set_ylim(y_limits[2])
        axs[1, 0].set_ylim(y_limits[3])
        axs[1, 1].set_ylim(y_limits[4])

    plt.tight_layout()
    if ENABLE_SAVE_PLOT_FILES:
        _write_plot_files(fig, t_values, R_pol, ef, R_data, thiele_data, polymer_mass, selected_Ca, radial_positions, compare_runs)
    st.pyplot(fig)
