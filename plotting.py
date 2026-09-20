import matplotlib.pyplot as plt
import streamlit as st
import numpy as np

def plot_results(t_values, R_pol, ef, R_data, thiele_data, polymer_mass, selected_Ca, radial_positions, AC_Conc, Cas, axis_locked=False):
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


# ON/OFF: write MWD PNG/SVG/CSV next to the on-screen figure.
ENABLE_SAVE_MWD_FILES = True


def plot_mwd(n, M, w_frac, Mn, Mw, PDI):
    import os
    import csv
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
