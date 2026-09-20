# QSSA Polymer Flow Model Simulation

Streamlit app for the **polymer flow model (PFM)** — also written **polymeric flow model** — of **single-particle growth** in **heterogeneous olefin polymerization**. This is a **polyolefin** problem: **polyethylene (PE)**, **polypropylene (PP)**, and **ethylene/1-olefin copolymers** made on a **supported solid catalyst** (Ziegler–Natta or metallocene). Polymer accumulates on the catalyst fragment, the particle grows, and **intraparticle mass transfer** can limit the rate. That morphology does not apply to typical homogeneous polymerizations of other resins.

The solver uses the **quasi-steady-state approximation (QSSA)** of the PFM so you can inspect **Thiele modulus**, effectiveness factor, polymerization rate, particle radius, yield, and the radial **monomer concentration** profile. A related single-particle description in the same literature is the **multigrain model (MGM)** / **multi-grain model**; this app implements the PFM/QSSA route, not a full MGM.

Adjust diffusivity, catalyst size, kinetics, and porosity in the sidebar and watch the plots update.

**App version:** 1.0.3

## Keywords

polymer flow model, polymeric flow model, PFM, QSSA, quasi-steady-state approximation, single-particle model, polyolefin, polyethylene, PE, polypropylene, PP, olefin polymerization, olefin homopolymerization, ethylene 1-olefin copolymerization, heterogeneous catalysis, supported catalyst, Ziegler–Natta, metallocene, catalyst particle growth, catalyst fragmentation, intraparticle diffusion, mass-transfer limitation, Thiele modulus, effectiveness factor, multigrain model, multi-grain model, MGM, random-pore polymeric flow model, slurry-phase polyethylene, gas-phase olefin polymerization, polymer reactor engineering

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

See [`requirements.txt`](requirements.txt) for pinned versions of Streamlit, NumPy, pandas, Matplotlib, and Pillow (chosen so the page installs on Python 3.9–3.12).
