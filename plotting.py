import os
import pandas as pd
import matplotlib.pyplot as plt
import streamlit as st

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
):
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))

    # Plot Polymerization Rate
    axs[0, 0].plot(t_values, R_pol)
    axs[0, 0].set_xlabel('Time (s)')
    axs[0, 0].set_ylabel('Polymerization Rate (gr pol/gr cat.hr)')
    axs[0, 0].set_title('Polymerization Rate')

    # Plot Efficiency
    axs[0, 1].plot(t_values, ef)
    axs[0, 1].set_xlabel('Time (s)')
    axs[0, 1].set_ylabel('Efficiency')
    axs[0, 1].set_title('Efficiency')

    # Plot Particle Radius
    axs[0, 2].plot(t_values, R_data)
    axs[0, 2].set_xlabel('Time (s)')
    axs[0, 2].set_ylabel('Particle Radius (m)')
    axs[0, 2].set_title('Particle Radius')

    # Plot Thiele Modulus
    axs[1, 0].plot(t_values, thiele_data)
    axs[1, 0].set_xlabel('Time (s)')
    axs[1, 0].set_ylabel('Thiele Modulus')
    axs[1, 0].set_title('Thiele Modulus')

    # Plot Polymer Mass
    axs[1, 1].plot(t_values, polymer_mass)
    axs[1, 1].set_xlabel('Time (s)')
    axs[1, 1].set_ylabel('Cumulative polymer mass (grams)')
    axs[1, 1].set_title('Polymerization yield')

    # Plot Monomer Concentration
    axs[1, 2].plot(radial_positions, selected_Ca)
    axs[1, 2].set_xlabel('Radial Position (m)')
    axs[1, 2].set_ylabel('Monomer Concentration (mol/m³)')
    axs[1, 2].set_title('Monomer Concentration')
    axs[1, 2].set_ylim([0, Cas])

    # Set axis limits if locked
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

    if save_plots:
        os.makedirs(plots_dir, exist_ok=True)
        fig.savefig(png_path, format="png", dpi=200, bbox_inches="tight")
        fig.savefig(svg_path, format="svg", bbox_inches="tight")

    if save_csv:
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
        pd.DataFrame(
            {
                "Radial Position (m)": radial_positions,
                "Monomer Concentration (mol/m3)": selected_Ca,
            }
        ).to_csv(concentration_csv_path, index=False)

    return fig, png_path, svg_path
