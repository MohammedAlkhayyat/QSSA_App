import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from calculations import SS_run, FS_run
from plotting import plot_results
import math

# Streamlit layout
st.set_page_config(layout="wide")
st.sidebar.title("Polymer Flow Model Simulation")
st.sidebar.markdown("""
Welcome to the Polymer Flow Model Simulation app! This tool allows you to explore and visualize the effects of various parameters on polymer particle growth.

Use the sliders in the sidebar to adjust parameters and see real-time updates of the simulation results.
""")
st.sidebar.markdown("**App Version:** 1.0.3")
st.sidebar.markdown("""
**Warning:** The numerical solver may take a significant amount of time to compute.
""")

# Sidebar header
st.sidebar.header("Simulation Settings")

# Define solver options and add a note for the disabled feature
options = ["QSSA", "Numerical (Under Development)"]
disabled_index = 1  # Index of the "under development" option

# Use st.radio to allow only one selection
solver_type = st.sidebar.radio(
    "Select Solver",
    options,
    index=0,  # Default selection to the first option
    help="Choose the solver to use for the simulation. Note that 'Numerical' is under development."
)

if solver_type == options[disabled_index]:
    # Handle the "under development" feature
    st.warning("The 'Numerical' solver is currently under development. Please select 'QSSA' for now.")
    st.stop()  # Stop further execution if 'Numerical' is selected

enable_thiele = st.sidebar.checkbox('Use presets', value=False)
show_tables = st.sidebar.checkbox('Show Tables', value=False)  # Checkbox to enable/disable tables

D_current = 0.5
d_current = 0.5
# Sidebar for sliders
with st.sidebar:
    if enable_thiele:
        st.sidebar.header("Presets")
        thiele_value = st.sidebar.selectbox("Select a preset", [f"ϕ = {0.50:.3f}",f"ϕ = {1.00:.3f}",f"ϕ = {3.00:.3f}", f"ϕ = {5.00:.3f}", f"ϕ = {10.00:.3f}"], help="Choose the Thiele modulus for the simulation.")
    if not enable_thiele:
        D = st.slider('Diffusivity (m²/s) × 10⁻¹⁰', min_value=D_current, max_value=50 * D_current, value=D_current, help="Diffusivity of the monomer in the polymer medium (multiplied by 10⁻¹⁰ for scale).")
    d = st.slider('Catalyst Diameter (m) × 10⁻⁵', min_value=d_current, max_value=30 * d_current, value=d_current, help="Diameter of the catalyst particle (multiplied by 10⁻⁵ for scale).")
    kd = st.slider('Deactivation Constant (s⁻¹) × 10⁻⁴', min_value=0.0, max_value=10.0, value=1.5, help="Rate constant for catalyst deactivation (multiplied by 10⁻⁴ for scale).")
    kp = st.slider('Propagation Constant (m³·mol/s)', min_value=1.0, max_value=10000.0, value=500.0, help="Rate constant for chain propagation.")
    Cas = st.slider('Monomer Concentration (mol/m³)', min_value=0.0, max_value=1000.0, value=100.0, help="Concentration of monomer at the surface of the particle.")
    eps = st.slider('Porosity', min_value=0.0, max_value=1.0, value=0.5, help="Porosity of the polymer.")
    t = st.slider('Simulation Time (s)', min_value=10, max_value=100000, value=15000, help="Total simulation time in seconds.")
    C1 = st.slider('Active Sites Concentration (mol/m³)', min_value=0.1, max_value=10.0, value=5.0, help="Concentration of active sites in the polymer.")
    time_idx = st.slider('Select Time Index', min_value=0, max_value=t, value=10, help="Index for selecting a specific time point for concentration profile.")

    if enable_thiele:
        if thiele_value == f"ϕ = {0.50:.3f}":
            D = 1e8 * (math.sqrt(kp * C1 * (1 - eps)) * d * 1e-6 / (2 * 3 * 0.5)) ** 2
        elif thiele_value == f"ϕ = {1.00:.3f}":
            D = 1e8 * (math.sqrt(kp * C1 * (1 - eps)) * d * 1e-6 / (2 * 3 * 1.0)) ** 2
        elif thiele_value == f"ϕ = {3.00:.3f}":
            D = 1e8 * (math.sqrt(kp * C1 * (1 - eps)) * d * 1e-6 / (2 * 3 * 3.0)) ** 2
        elif thiele_value == f"ϕ = {5.00:.3f}":
            D = 1e8 * (math.sqrt(kp * C1 * (1 - eps)) * d * 1e-6 / (2 * 3 * 5.0)) ** 2
        else:
            D = 1e8 * (math.sqrt(kp * C1 * (1 - eps)) * d * 1e-6 / (2 * 3 * 10.0)) ** 2

dencat = 900
denpol = 2300

# Lock Axes Checkbox
lock_axes = st.checkbox('Lock Y-Axis', value=False, help="Lock the y-axis limits for all plots to the current range to allow comparison across different simulations.")

# Compute data based on selected solver
if solver_type == "QSSA":
    t_values, R_pol, ef, R_data, thiele_data, polymer_mass, AC_Conc, Ca = SS_run(
        D * 1e-8, kp, Cas, d * 1e-6, t, C1, kd * 1e-4, dencat, denpol, eps)
else:  # FS_run
    R_data, R_pol, ef, AC_Conc, Ca, Avg_Conc_Profile, R_pol, polymer_mass, t_values = FS_run(
        D * 1e-8, kp, Cas, d * 1e-6, t, C1, kd * 1e-4, dencat, denpol, eps)

# Store axis limits in session state if checkbox is checked
if lock_axes:
    if 'y_limits' not in st.session_state:
        st.session_state.y_limits = [
            (min(R_pol), max(R_pol)),
            (min(ef), max(ef)),
            (min(R_data), max(R_data)),
            (min(thiele_data), max(thiele_data)),
            (min(polymer_mass), max(polymer_mass)),
            (min(AC_Conc), max(AC_Conc))
        ]
else:
    if 'y_limits' in st.session_state:
        del st.session_state.y_limits

# Extract the concentration data for the selected time
selected_time = t_values[time_idx]
selected_Ca = Ca[time_idx][1]
selected_R = R_data[time_idx]
radial_positions = np.linspace(1e-9 / selected_R, selected_R / selected_R, len(selected_Ca))  # Radial positions from rls to R

# Plot results
plot_results(
    t_values, R_pol, ef, R_data, thiele_data, polymer_mass, selected_Ca, radial_positions, AC_Conc, Cas, axis_locked=lock_axes
)

# Convert tuples to a DataFrame-friendly format
# Determine the maximum number of columns needed
max_concentration_length = max(len(item[1]) for item in Ca)

# Prepare data for DataFrame
conc_data = {
    'Time': [item[0] for item in Ca]
}
# Add each concentration value as a separate column
for i in range(max_concentration_length):
    conc_data[f'{i+1}'] = [item[1][i] if i < len(item[1]) else None for item in Ca]

# Create DataFrame
output_conc = pd.DataFrame(conc_data)

# Create DataFrames for the other outputs
output_data = {
    'Time (s)': t_values,
    'Polymerization Rate': R_pol,
    'Efficiency': ef,
    'Particle Radius (m)': R_data,
    'Thiele Modulus': thiele_data,
    'Cumulative polymer mass (grams)': polymer_mass,
    'Active Sites Conc. (mol/m³)': AC_Conc,
}

output_df = pd.DataFrame(output_data)

# Define a function to format numeric values
def format_value(x):
    try:
        return f'{float(x):.6e}'  # Format as scientific notation
    except ValueError:
        return x  # Return non-numeric values as is
if show_tables:
    # Create formatted DataFrame for display
    output_df_display = output_df.map(format_value)
    output_conc_display = output_conc.map(format_value)

    # Display the tables
    st.write("### Simulation Results")
    st.dataframe(output_df_display)
    st.write("### Concentration Profile")
    st.dataframe(output_conc_display)

    # Add download button for saving the table
    csv = output_df.to_csv(index=False)
    st.download_button(
        label="Download Simulation Results CSV",
        data=csv,
        file_name="simulation_results.csv",
        mime="text/csv"
    )

    # Add download button for saving the concentration data
    con = output_conc.to_csv(index=False)
    st.download_button(
        label="Download Conc. CSV",
        data=con,
        file_name="Conc_results.csv",
        mime="text/csv"
    )
