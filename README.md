# QSSA Polymer Flow Model Simulation

Streamlit app for exploring polymer particle growth with a quasi-steady-state approximation (QSSA) solver. Adjust diffusivity, catalyst size, kinetics, and porosity in the sidebar and watch rate, efficiency, radius, Thiele modulus, yield, and the radial monomer profile update.

**App version:** 1.0.3

## Run locally

```bash
pip install -r requirements.txt
streamlit run app.py
```

The app listens on port `8501` by default (`8080` in Docker).

## GitHub Codespaces

This repo includes a `.devcontainer` setup. After the container starts, Streamlit is launched with:

```bash
streamlit run app.py --server.enableCORS false --server.enableXsrfProtection false
```

## Docker

```bash
docker build -t qssa-app .
docker run -p 8080:8080 qssa-app
```

## Parameters

| Input | Meaning |
| --- | --- |
| Solver | QSSA (numerical solver is listed but not enabled yet) |
| Diffusivity | Monomer diffusivity in the polymer medium |
| Catalyst diameter | Initial catalyst particle size |
| Deactivation / propagation constants | Catalyst deactivation and chain growth rates |
| Monomer concentration | Surface monomer concentration |
| Porosity | Polymer porosity |
| Simulation time | Total time in seconds |
| Active sites concentration | Initial active-site loading |
| Presets | Optional Thiele modulus presets that set diffusivity |

## Project files

- `app.py` — Streamlit page and controls
- `calculations.py` — QSSA (`SS_run`) and numerical (`FS_run`) solvers
- `plotting.py` — figure layout
- `requirements.txt` — Python packages
- `packages.txt` — optional Linux packages for Streamlit Cloud / Codespaces
- `.streamlit/config.toml` — page theme
- `plot_data/` — CSV of plotted series (written when saving is on)
- `plots/` — PNG and SVG copies of the figure (written when saving is on)

## Requirements

See [`requirements.txt`](requirements.txt) for pinned versions of Streamlit, NumPy, pandas, Matplotlib, and Pillow.
