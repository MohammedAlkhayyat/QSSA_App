import numpy as np
import math

pi = math.pi

# ON/OFF: use overflow-safe QSSA concentration and efficiency formulas.
ENABLE_STABLE_QSSA = True


def stable_sinh_ratio(r, R, s):
    """Return sinh(r*s) / sinh(R*s) without overflowing for large Thiele moduli."""
    r = np.asarray(r, dtype=float)
    Rs = float(s) * float(R)
    rs = float(s) * r
    if abs(Rs) < 1e-12:
        return np.ones_like(rs, dtype=float) * (r / R if R != 0 else 1.0)
    if abs(Rs) > 40.0:
        return np.exp(rs - Rs) * (1.0 - np.exp(-2.0 * rs)) / (1.0 - np.exp(-2.0 * Rs))
    return np.sinh(rs) / np.sinh(Rs)


def qssa_ca(Cas, r, R, b_sqrt):
    r = max(float(r), 1e-30)
    R = max(float(R), 1e-30)
    if ENABLE_STABLE_QSSA:
        return Cas * (R / r) * float(stable_sinh_ratio(r, R, b_sqrt))
    return Cas * R / r * np.sinh(r * b_sqrt) / (np.sinh(R * b_sqrt))


def x_coth_x_minus_one(x):
    """Return x*coth(x) - 1 without overflow."""
    x = float(x)
    ax = abs(x)
    if ax < 1e-8:
        return x * x / 3.0
    if ax > 40.0:
        return ax - 1.0
    return x * np.cosh(x) / np.sinh(x) - 1.0


def qssa_R_over_sinh_times_integral(Cas, R, b_sqrt, b):
    """Cas * R / sinh(R s) * (R s cosh(R s) - sinh(R s)) / b  ==  Cas * R / b * (R s coth(R s) - 1)."""
    if abs(b) < 1e-30:
        return 0.0
    Rs = float(b_sqrt) * float(R)
    if ENABLE_STABLE_QSSA:
        return Cas * R / b * x_coth_x_minus_one(Rs)
    return Cas * R / (np.sinh(Rs)) * (Rs * np.cosh(Rs) - np.sinh(Rs)) / b

def A(r, j):
    if j == 0:
        return -1.0  # error
    
    A = (2.0 * r[j] - r[j + 1]) / (r[j] * (r[j - 1] - r[j]) * (r[j - 1] - r[j + 1]))
    return A

def B(r, j):
    if j == 0:
        return -1.0  # error
    
    B = (3.0 * r[j] - r[j + 1] - r[j - 1]) / (r[j] * (r[j] - r[j - 1]) * (r[j] - r[j + 1]))
    return B

def C(r, j):
    if j == 0:
        return -1.0  # error
    
    C = (2.0 * r[j] - r[j - 1]) / (r[j] * (r[j + 1] - r[j - 1]) * (r[j + 1] - r[j]))
    return C

def B0(r, r0):
    B0 = -6.0 / ((r - r0) ** 2)
    return B0

def rate(Conc, Y, kp, j):
    rate = kp * Y[j] * 0.5 * (Conc[j] + Conc[j + 1])
    return rate

def rate_t(Conc, Y, kp, j):
    rate = kp[j] * Y[j] * 0.5 * (Conc[j] + Conc[j + 1])
    return rate


def linspace(start, stop, num):
    return np.linspace(start, stop, num)

def calculate_C0(C1, C2, kd, t):
    return C1 * np.exp(-kd * t)

def beta(t, C0, kp, kd, Dae, alpha, eps):
    return (kp * C0 * np.exp(-kd * t)) / Dae * (1 - eps) / (alpha ** 3)

def calculate_Vcati(Vcat_init, vtot, r):
    return Vcat_init * (4.0 / 3.0 * pi * r ** 3) / vtot

def SS_run(D, kp, Cas, d, t, C1, kd, dencat, denpol, eps):
    Nt = 15000
    MW = 28.05
    denampol = 903.5 * (1 - eps)
    rls = 1e-9
    C2 = 0

    Cas = Cas
    iter = 50
    Dae = D
    rcat = d / 2
    Vcat_init = 4.0 / 3.0 * pi * rcat ** 3

    t_values = linspace(0, t, Nt)
    vtot = Vcat_init
    vpol = Vcat_init
    R_pol = []
    R = rcat
    ddtc = t / (Nt - 1.0)
    ddt = 0
    Ca = []
    R_data = []
    ef = []
    thiele_data = []
    polymer_mass = []
    AC_Conc = []
    acc_mass = []
    acmass = 0
    for t in t_values:
        delta_r = (R - rls) / (iter - 1.0)
        alpha = R / rcat
        C0 = calculate_C0(C1, C2, kd, t)
        b = beta(t, C0, kp, kd, Dae, alpha, eps)
        if (not np.isfinite(b)) or b < 0:
            b = 0.0
        b_sqrt = np.sqrt(b) if b > 0 else 0.0

        Ca_values = []
        for r in np.linspace(rls, R, iter):
            Ca_value = qssa_ca(Cas, r, R, b_sqrt)
            Ca_values.append(Ca_value)
        if len(Ca_values) != iter:
            Ca_values.append(Cas)

        if vpol == 0:
            vpol += 1e-9

        flux_term = qssa_R_over_sinh_times_integral(Cas, R, b_sqrt, b)
        pmass = Vcat_init / vpol * kp * C0 * np.exp(-kd * t) * 28 * ddt * 4 * pi * flux_term
        if not np.isfinite(pmass):
            pmass = 0.0
        acmass += pmass
        Rins_pol = pmass / Vcat_init / dencat / (ddt+0.00000001) * 3600.0 / 1000.0

        vpol_tot = pmass / denpol / 1000
        vpol += vpol_tot
        vtot += vpol_tot / (1 - eps)
        thiele = R / 3.0 * b_sqrt
        R = (3.0 / 4.0 / pi * vpol) ** (1.0 / 3.0)

        Ca.append((t, Ca_values))
        if abs(b) < 1e-30 or vpol == 0:
            eff = 0.0
        else:
            eff = 4 * pi * flux_term / vpol
        if not np.isfinite(eff):
            eff = 0.0

        ddt = ddtc
        R_data.append(R)
        R_pol.append(Rins_pol)
        ef.append(eff)
        AC_Conc.append(C0 * np.exp(-kd * t) / (alpha ** 3))
        thiele_data.append(thiele)
        polymer_mass.append(acmass)

    return t_values, R_pol, ef, R_data, thiele_data, polymer_mass, AC_Conc, Ca



def FS_run(D, kp, Surf_Conc, d, set_time, Y0 ,kd, dencat, denpol, eps):
    MW = 28.05
    denspol = 903.5  # g/L
    dencat = 2300.0  # g/L
    v0_cat = 4.0 / 3.0 * pi * (d / 2.0) ** 3
    Vpol = v0_cat
    pol_R = d / 2.0
    dt = 0.0
    total_mass = 0.0
    mass_i = 0.0
    vtot = v0_cat
    save = 0
    timer = 1
    ins_rate = 0
    eff = 0
    nd = 20
    # Initialize arrays
    r_values = np.linspace(0, d / 2.0, nd + 1)
    step = d / 2.0 / nd  # Calculate the step size
    r = [i * step for i in range(nd + 1)]  # Initial radius calculation
    Conc = np.zeros(nd + 1)
    rate_pol = np.zeros(nd)
    updatedConc = np.zeros(nd + 1)
    Y = np.full(nd, Y0)
    updatedY = np.zeros(nd)
    Shell_Conc = np.zeros(nd)
    inst_mass = np.zeros(nd)
    v = np.zeros(nd)
    v_cat = np.zeros(nd)
    v0 = np.zeros(nd)
    R_data = []
    R_pol = []
    AC_Conc = []
    ef = []
    save_time = []
    dt_frame = []
    mass = []
    Conc_Profile = []
    Avg_Conc_Profile = []
    rate_Profile = []
    AC_Profile = []
    mass_Profile = []

    # Initialize grid volumes and concentrations
    step = d / (2.0 * nd)
    for k in range(nd):
        v[k] = 4.0 / 3.0 * pi * (r[k + 1] ** 3 - r[k] ** 3)
        v0[k] = v[k]

    Conc[-1] = Surf_Conc

    R_data.append(pol_R)
    R_pol.append(ins_rate)
    ef.append(eff)
    save_time.append(0.0)
    Conc_Profile.append(Conc.copy())
    Avg_Conc_Profile.append(Shell_Conc.copy())
    AC_Profile.append(Y.copy())
    rate_Profile.append(rate_pol.copy())
    mass_Profile.append(inst_mass.copy())
    dt_frame.append(dt)

    time = 0.0
    while time < set_time:
        # Setting Boundary Conditions
        Conc[-1] = Surf_Conc
        r[0] = 0

        # Calculate B0 and rate for center monomer concentration
        B0_ = B0(r[1], r[0])
        rate_pol[0] = rate(Conc, Y, kp, 0)

        # Monomer concentration at the center
        Cent_Conc = Conc[0] + dt * (D * B0_ * (Conc[0] - Conc[1]) - rate_pol[0])
        if np.isnan(Cent_Conc):
            print("Parameter became NaN. Stopping the loop.")
            break

        # Monomer concentration at grid sites
        for k in range(1, nd):
            A_ = A(r, k)
            B_ = B(r, k)
            C_ = C(r, k)
            rate_pol[k] = rate(Conc, Y, kp, k)
            Conc_at_k = Conc[k] + dt * (2.0 * D * (A_ * Conc[k - 1] + B_ * Conc[k] + C_ * Conc[k + 1]) - rate_pol[k])
            updatedConc[k] = Conc_at_k

        # Saving results
        Conc[0] = Cent_Conc
        Conc[1:nd] = updatedConc[1:nd]
        Shell_Conc = 0.5 * (Conc[:-1] + Conc[1:])

        total_mass = 0
        for k in range(nd):
            mass_i = MW * v[k] * dt * Y[k] * kp * Shell_Conc[k]
            inst_mass[k] = mass_i
            total_mass += mass_i

        for k in range(nd):
            volume_k = v[k] * ((kp * Shell_Conc[k] * Y[k] * MW * dt / (1 - eps) / denspol / 1000) + 1)
            v[k] = volume_k

        for k in range(1, nd + 1):
            rad = (3.0 / 4.0 / pi * v[k - 1] + r[k - 1] ** 3) ** (1.0 / 3.0)
            r[k] = rad

        # Y0 calculation
        for k in range(nd):
            S = Y0 * np.exp(-kd * time)
            updatedY[k] = S * v0[k] / v[k]

        Y[:] = updatedY

        pol_R = r[nd - 1]
        ins_rate = total_mass / dencat / v0_cat / dt / 1000.0 * 3600.0
        vtot += total_mass / denspol / 1000

        pol_RR = (3.0 / 4.0 / pi * vtot) ** (1.0 / 3.0)
        Vpol = 4.0 / 3.0 * pi * pol_R ** 3
        AvgConc = Conc[0:nd]
        esum = np.sum(AvgConc * v)
        eff = esum / Vpol / Surf_Conc

        # Time step and time increment
        smallest_dt = float('inf')
        dt2 = r[1] ** 2 / (6.0 * D + 0.5 * r[1] ** 2 * kp * Y[0])

        for k in range(1, nd):
            dt1 = -1.0 / (2.0 * D * B(r, k) - 0.5 * kp * Y[k])
            dt = min(dt1, dt2)
            if dt < smallest_dt:
                dt2 = dt
                smallest_dt = dt

        time += smallest_dt
        dt = smallest_dt

        save += dt
        if save > timer:
            mass.append(total_mass)
            R_data.append(pol_RR)
            R_pol.append(ins_rate)
            ef.append(eff)
            save_time.append(time)
            Avg_Conc_Profile.append(Shell_Conc.copy())
            Conc_Profile.append(Conc.copy())
            rate_Profile.append(rate_pol.copy())
            AC_Profile.append(Y.copy())
            mass_Profile.append(total_mass.copy())
            dt_frame.append(dt)
            save = 0
            print(f"time : {time} dt: {smallest_dt}")

    return (R_data, R_pol, ef, AC_Conc, Conc_Profile, Avg_Conc_Profile, rate_Profile, mass_Profile, save_time)
