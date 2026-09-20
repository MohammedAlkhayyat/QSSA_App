import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from calculations import SS_run, FS_run, compute_mwd, ENABLE_MWD_CALC
from plotting import plot_results, plot_mwd
import math
import os

# Feature flags (ON/OFF)
ENABLE_PAGE_HEADER = True
ENABLE_ABOUT_SECTION = True
ENABLE_PARAM_SUMMARY = True
ENABLE_SAVE_PLOTS = True
ENABLE_SAVE_PLOT_CSV = True
ENABLE_PLOT_DOWNLOADS = True
ENABLE_SAFE_TIME_INDEX = True
ENABLE_BLANK_SCREEN_GUARD = True
ENABLE_NUMERICAL = True
ENABLE_COMPARE_PLOTS = True
ENABLE_MULTICOMPONENT = True
ENABLE_MWD = True

# Streamlit layout
st.set_page_config(
    page_title="QSSA Polymer Flow Model",
    page_icon="⚗️",
    layout="wide",
)
st.sidebar.title("Polymer Flow Model Simulation")
st.sidebar.markdown("""
Welcome to the Polymer Flow Model Simulation app! This tool allows you to explore and visualize the effects of various parameters on polymer particle growth.

Use the sliders in the sidebar to adjust parameters and see real-time updates of the simulation results.
""")
st.sidebar.markdown("**App Version:** 1.0.3")
st.sidebar.markdown("""
**Warning:** The numerical solver may take a significant amount of time to compute.
""")

if ENABLE_PAGE_HEADER:
    st.title("QSSA Polymer Flow Model")
    st.markdown(
        "Quasi-steady-state polymer particle growth: set the kinetic and transport inputs in the sidebar, then inspect rate, efficiency, radius, Thiele modulus, yield, and the radial monomer profile."
    )

if ENABLE_ABOUT_SECTION:
    with st.expander("About this page", expanded=False):
        st.markdown(
            """
This page runs a **QSSA** (quasi-steady-state approximation) polymer-flow model.

- **Polymerization rate** is the instantaneous yield of polymer per gram of catalyst per hour.
- **Efficiency** is the effectiveness factor for monomer use inside the growing particle.
- **Particle radius** grows as polymer accumulates around the catalyst.
- **Thiele modulus** compares reaction and diffusion; larger values mean stronger intraparticle gradients.
- **Presets** back-calculate diffusivity from a chosen Thiele modulus.

Python packages live in `requirements.txt`. Optional Linux packages for Streamlit Cloud / Codespaces live in `packages.txt`. Theme settings live in `.streamlit/config.toml`.
            """
        )

# Sidebar header
st.sidebar.header("Simulation Settings")

if ENABLE_NUMERICAL:
    options = ["QSSA", "Numerical"]
else:
    options = ["QSSA", "Numerical (Under Development)"]
disabled_index = 1

solver_type = st.sidebar.radio(
    "Select Solver",
    options,
    index=0,
    help="Choose the solver to use for the simulation."
)

if (not ENABLE_NUMERICAL) and solver_type == options[disabled_index]:
    st.warning("The 'Numerical' solver is currently under development. Please select 'QSSA' for now.")
    st.stop()

enable_thiele = st.sidebar.checkbox('Use Presets', value=False)
show_tables = st.sidebar.checkbox('Show Tables', value=False)

D_current = 0.5
d_current = 0.5
enable_second = False
D2 = kp2 = Cas2 = MW2 = None
ktr = 0.0
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
    if ENABLE_SAFE_TIME_INDEX:
        time_select_s = st.slider('Select Time (s)', min_value=0, max_value=int(t), value=min(10, int(t)), help="Simulation time at which the monomer concentration profile is shown.")
        time_idx = time_select_s
    else:
        time_idx = st.slider('Select Time Index', min_value=0, max_value=t, value=10, help="Index for selecting a specific time point for concentration profile.")
    if ENABLE_MWD:
        ktr = st.slider('Chain Transfer Constant (s⁻¹)', min_value=0.0, max_value=50.0, value=1.0, help="Chain transfer / chain-stop rate used for the Schulz-Flory MWD.")
    if ENABLE_MULTICOMPONENT:
        enable_second = st.checkbox('Multicomponent (second monomer)', value=False, help="Add a second monomer with its own diffusivity, kp, and surface concentration.")
        if enable_second:
            D2 = st.slider('Monomer 2 Diffusivity (m²/s) × 10⁻¹⁰', min_value=D_current, max_value=50 * D_current, value=D_current)
            kp2 = st.slider('Monomer 2 Propagation Constant (m³·mol/s)', min_value=1.0, max_value=10000.0, value=200.0)
            Cas2 = st.slider('Monomer 2 Concentration (mol/m³)', min_value=0.0, max_value=1000.0, value=50.0)
            MW2 = st.slider('Monomer 2 Molecular Weight (g/mol)', min_value=1.0, max_value=200.0, value=42.08)

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

lock_axes = st.checkbox('Lock Y-Axis', value=False, help="Lock the y-axis limits for all plots to the current range to allow comparison across different simulations.")


def run_solver():
    if solver_type == "QSSA":
        if ENABLE_MULTICOMPONENT and enable_second:
            return SS_run(
                D * 1e-8, kp, Cas, d * 1e-6, t, C1, kd * 1e-4, dencat, denpol, eps,
                D2=D2 * 1e-8, kp2=kp2, Cas2=Cas2, MW2=MW2)
        return SS_run(
            D * 1e-8, kp, Cas, d * 1e-6, t, C1, kd * 1e-4, dencat, denpol, eps)
    fs_out = FS_run(
        D * 1e-8, kp, Cas, d * 1e-6, t, C1, kd * 1e-4, dencat, denpol, eps)
    if len(fs_out) == 8:
        return fs_out
    R_data, R_pol, ef, AC_Conc, Ca, Avg_Conc_Profile, R_pol, polymer_mass, t_values = fs_out
    return t_values, R_pol, ef, R_data, None, polymer_mass, AC_Conc, Ca

if ENABLE_BLANK_SCREEN_GUARD:
    try:
        ss_out = run_solver()
    except Exception as exc:
        st.error("The solver could not finish with the current sliders. Try a different combination of Diffusivity, Catalyst Diameter, or Simulation Time.")
        st.exception(exc)
        st.stop()
else:
    ss_out = run_solver()

t_values, R_pol, ef, R_data, thiele_data, polymer_mass, AC_Conc, Ca = ss_out[:8]
Ca2_profiles = ss_out[8] if len(ss_out) > 8 else None
selected_Ca2 = None

if thiele_data is None:
    thiele_data = [0.0 for _ in t_values]


def _finite_limits(arr):
    a = np.asarray(arr, dtype=float).ravel()
    a = a[np.isfinite(a)]
    if a.size == 0:
        return (0.0, 1.0)
    lo = float(np.min(a))
    hi = float(np.max(a))
    if lo == hi:
        return (lo - 1.0, hi + 1.0)
    return (lo, hi)

if lock_axes:
    if 'y_limits' not in st.session_state:
        st.session_state.y_limits = [
            _finite_limits(R_pol),
            _finite_limits(ef),
            _finite_limits(R_data),
            _finite_limits(thiele_data),
            _finite_limits(polymer_mass),
            _finite_limits(AC_Conc)
        ]
else:
    if 'y_limits' in st.session_state:
        del st.session_state.y_limits

n_times = len(t_values)
if n_times == 0:
    st.error("The solver returned no time points.")
    st.stop()
if ENABLE_SAFE_TIME_INDEX:
    time_idx = int(np.argmin(np.abs(np.asarray(t_values) - float(time_select_s))))
time_idx = max(0, min(int(time_idx), n_times - 1))
selected_time = t_values[time_idx]
ca_entry = Ca[time_idx]
selected_Ca = ca_entry[1] if isinstance(ca_entry, tuple) else ca_entry
selected_R = R_data[time_idx]
if selected_R == 0 or not np.isfinite(selected_R):
    selected_R = 1e-9
radial_positions = np.linspace(1e-9 / selected_R, selected_R / selected_R, len(selected_Ca))
if Ca2_profiles is not None:
    selected_Ca2 = Ca2_profiles[time_idx][1]

if ENABLE_PARAM_SUMMARY:
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Time (s)", f"{selected_time:.4g}")
    m2.metric("Thiele modulus", f"{thiele_data[time_idx]:.4g}")
    m3.metric("Efficiency", f"{ef[time_idx]:.4g}")
    m4.metric("Particle radius (m)", f"{selected_R:.4g}")

compare_runs = None
if ENABLE_COMPARE_PLOTS:
    if "compare_runs" not in st.session_state:
        st.session_state.compare_runs = []
    st.markdown("### Compare plots")
    c1, c2 = st.columns([3, 1])
    run_label = c1.text_input("Run label", value=f"Run {len(st.session_state.compare_runs) + 1}")
    if c1.button("Save current run for comparison"):
        st.session_state.compare_runs.append({
            "label": run_label,
            "t": list(t_values),
            "R_pol": list(R_pol),
            "ef": list(ef),
            "R_data": list(R_data),
            "thiele_data": list(thiele_data),
            "polymer_mass": list(polymer_mass),
            "radial_positions": list(radial_positions),
            "selected_Ca": list(selected_Ca),
        })
    if c2.button("Clear saved runs"):
        st.session_state.compare_runs = []
    if st.session_state.compare_runs:
        st.caption("Saved runs: " + ", ".join(run["label"] for run in st.session_state.compare_runs))
        compare_runs = st.session_state.compare_runs

if ENABLE_BLANK_SCREEN_GUARD:
    try:
        fig, png_path, svg_path = plot_results(
            t_values, R_pol, ef, R_data, thiele_data, polymer_mass, selected_Ca, radial_positions, AC_Conc, Cas,
            axis_locked=lock_axes,
            save_plots=ENABLE_SAVE_PLOTS,
            save_csv=ENABLE_SAVE_PLOT_CSV,
            compare_runs=compare_runs,
            selected_Ca2=selected_Ca2,
            Cas2=Cas2,
        )
    except Exception as exc:
        st.error("Plotting failed for the current results. The simulation data is still available in the tables if enabled.")
        st.exception(exc)
        png_path = svg_path = ""
else:
    fig, png_path, svg_path = plot_results(
        t_values, R_pol, ef, R_data, thiele_data, polymer_mass, selected_Ca, radial_positions, AC_Conc, Cas,
        axis_locked=lock_axes,
        save_plots=ENABLE_SAVE_PLOTS,
        save_csv=ENABLE_SAVE_PLOT_CSV,
        compare_runs=compare_runs,
        selected_Ca2=selected_Ca2,
        Cas2=Cas2,
    )

if ENABLE_MWD and ENABLE_MWD_CALC:
    n_mwd, M_mwd, w_mwd, Mn, Mw, PDI = compute_mwd(Ca, polymer_mass, kp, kd * 1e-4, ktr, MW=28.05)
    st.write(f"**MWD:** Mn = {Mn:.4e} g/mol, Mw = {Mw:.4e} g/mol, PDI = {PDI:.3f}")
    plot_mwd(n_mwd, M_mwd, w_mwd, Mn, Mw, PDI)

if ENABLE_PLOT_DOWNLOADS and ENABLE_SAVE_PLOTS:
    if png_path and svg_path and os.path.exists(png_path) and os.path.exists(svg_path):
        d1, d2 = st.columns(2)
        with d1:
            with open(png_path, "rb") as png_file:
                st.download_button(
                    label="Download plots PNG",
                    data=png_file,
                    file_name="simulation_results.png",
                    mime="image/png",
                )
        with d2:
            with open(svg_path, "rb") as svg_file:
                st.download_button(
                    label="Download plots SVG",
                    data=svg_file,
                    file_name="simulation_results.svg",
                    mime="image/svg+xml",
                )

def _conc_row(item):
    if isinstance(item, tuple) and len(item) > 1:
        return item[0], list(item[1])
    return None, list(item)

conc_times = []
conc_rows = []
for item in Ca:
    t_i, vals = _conc_row(item)
    conc_times.append(t_i if t_i is not None else np.nan)
    conc_rows.append(vals)

max_concentration_length = max((len(row) for row in conc_rows), default=0)
conc_data = {'Time': conc_times}
for i in range(max_concentration_length):
    conc_data[f'{i+1}'] = [row[i] if i < len(row) else None for row in conc_rows]
output_conc = pd.DataFrame(conc_data)

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

def format_value(x):
    try:
        return f'{float(x):.6e}'
    except ValueError:
        return x
if show_tables:
    try:
        output_df_display = output_df.map(format_value)
        output_conc_display = output_conc.map(format_value)
    except Exception:
        output_df_display = output_df.applymap(format_value)
        output_conc_display = output_conc.applymap(format_value)

    st.write("### Simulation Results")
    st.dataframe(output_df_display)
    st.write("### Concentration Profile")
    st.dataframe(output_conc_display)

    csv = output_df.to_csv(index=False)
    st.download_button(
        label="Download Simulation Results CSV",
        data=csv,
        file_name="simulation_results.csv",
        mime="text/csv"
    )
    con = output_conc.to_csv(index=False)
    st.download_button(
        label="Download Conc. CSV",
        data=con,
        file_name="Conc_results.csv",
        mime="text/csv"
    )
