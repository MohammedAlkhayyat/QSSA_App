import streamlit as st
import numpy as np
import time
import matplotlib.pyplot as plt
from calculations import SS_run
from calculations import FS_run
from plotting import plot_results

# Streamlit layout
st.set_page_config(layout="wide")

st.title("Polymer Flow Model Simulation")
st.markdown("""
Welcome to the Polymer Flow Model Simulation app! This tool allows you to explore and visualize the effects of various parameters on polymer partilce growth.

Use the sliders in the sidebar to adjust parameters and see real-time updates of the simulation results.
""")

st.sidebar.markdown("""
**Warning:** The numerical solver may take a significant amount of time to compute
""")


st.sidebar.header("Simulation Settings")
solver_type = st.sidebar.selectbox("Select Solver", ["QSSA", "Numerical"])

# Define current values
D_current = 1
d_current = 3

# Sidebar for sliders
with st.sidebar:
    D = st.slider('Diffusivity (m²/s) × 10⁻¹⁰', min_value=D_current, max_value=10*D_current, value=D_current, help="Diffusivity of the monomer in the medium.")
    d = st.slider('Catalyst Diameter (m) × 10⁻⁵', min_value=d_current, max_value=2*d_current, value=d_current)
    kd = st.slider('Deactivation Constant (s⁻¹) × 10⁻⁴', min_value=0.0, max_value=5.0, value=1.5)
    kp = st.slider('Propagation Constant (m³·mol/s)', min_value=1.0, max_value=10000.0, value=500.0)
    Cas = st.slider('Monomer Concentration (mol/m³)', min_value=0.0, max_value=500.0, value=100.0)
    eps = st.slider('Porosity', min_value=0.0, max_value=1.0, value=0.1)
    t = st.slider('Simulation Time (s)', min_value=1000, max_value=100000, value=15000)
    C1 = st.slider('Active Sites Concentration (mol/m³)', min_value=0.1, max_value=10.0, value=5.0)

    time_idx = st.slider('Select Time Index', min_value=0, max_value=t, value=10)

dencat = 900
denpol = 2300

# Lock Axes Checkbox
lock_axes = st.checkbox('Lock Y-Axis', value=False)

# Compute data based on selected solver
if solver_type == "QSSA":
    t_values, R_pol, ef, R_data, thiele_data, polymer_mass, AC_Conc, Ca = SS_run(
        D * 1e-8, kp, Cas, d * 1e-6, t, C1, kd * 1e-4, dencat, denpol, eps)
else:  # FS_run
    t_values, R_pol, ef, R_data, thiele_data, polymer_mass, AC_Conc, Ca = FS_run(
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
