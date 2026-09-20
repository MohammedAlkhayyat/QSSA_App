# Polymer Flow Model Simulation (QSSA)

Streamlit app for a quasi-steady-state approximation (QSSA) polymer particle growth model. Adjust diffusivity, catalyst size, kinetics, and other parameters, then plot polymerization rate, efficiency, particle radius, Thiele modulus, yield, and a radial monomer profile.

**App version:** 1.0.3

## Run locally

```bash
python -m pip install -r requirements.txt
streamlit run app.py
```

The app opens at [http://localhost:8501](http://localhost:8501). Use the **QSSA** solver. The **Numerical** option is still under development and will stop the run.

Default QSSA settings take a noticeable amount of time (tens of seconds) because the solver uses 15,000 time samples.

## Docker

```bash
docker build -t qssa-app .
docker run -p 8080:8080 qssa-app
```

Then open [http://localhost:8080](http://localhost:8080).

## Repository layout

| File | Role |
| --- | --- |
| `app.py` | Streamlit UI, sliders, tables, CSV download |
| `calculations.py` | QSSA (`SS_run`) and unfinished numerical (`FS_run`) solvers |
| `plotting.py` | Matplotlib figures shown in the app |
| `requirements.txt` | Python dependencies |
| `Dockerfile` | Container image (Python 3.9, Streamlit on port 8080) |

## Notes

- Turning on **Show Tables** also enables CSV downloads.
- **Use presets** sets diffusivity from a chosen Thiele modulus \(\phi\).
- No API keys or secrets are required.
