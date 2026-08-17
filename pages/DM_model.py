import streamlit as st
import numpy as np
from scipy.special import expi
from scipy.stats import norm

#st.set_page_config(page_title="DM Two-State Markov — Uncertainty Analysis", layout="centered")

st.title("DM Two-State Markov Process Uncertainty Analysis")
st.caption("Dye-Cycling Mode")
st.divider()

def expint(x):
    x = np.asarray(x, dtype=float)
    return -expi(-x)

def igamma(a, x):
    return expint(x)

def Porb1(tau, delta_t):
    return 1 + tau / delta_t * (np.exp(-delta_t / tau) - 1)

def Phit(tau, delta_t):
    return tau / delta_t * (1 - np.exp(-delta_t / tau))

def missN(tau, delta_t, cover=0.9):
    y = int(np.log(1 - cover) / np.log(1 + tau / delta_t * (np.exp(-delta_t / tau) - 1))) - 1
    return 0 if y < 0 else y

def bM(tau, delta_t):
    tau = np.asarray(tau, dtype=float)
    A = np.exp(-delta_t / tau)
    return delta_t * (A + 1) / (1 - A) - 2 * tau

def VarM(tau, delta_t):
    tau = np.asarray(tau, dtype=float)
    A = np.exp(-delta_t / tau)
    return delta_t**2 * A / (6 * tau / delta_t * (1 - A))

def bT(tau, T):
    tau = np.asarray(tau, dtype=float)
    T   = np.asarray(T,   dtype=float)
    tau, T = np.broadcast_arrays(tau, T)
    tau = tau.astype(float).copy()
    T   = T.astype(float).copy()
    with np.errstate(all='ignore'):
        y = (2 * T**2 / tau * np.exp(T / tau) *
             (expint(2 * T / tau) - expint(T / tau)) +
             (T - tau) * (2 - np.exp(-T / tau)) + tau)
        bad = ~np.isfinite(y)
        if np.any(bad):
            tb, Tb = tau[bad], T[bad]
            y[bad] = np.exp(-Tb / tb) * (tb / 2 + tb**2 / (2 * Tb)) - 4 * tb**2 / Tb + tb
        I1a = (2 * T / tau * np.exp(T / tau) *
               (igamma(0, T / tau) - igamma(0, 2 * T / tau)) +
               np.exp(-T / tau) - 1)
        bad_i = ~np.isfinite(I1a)
        if np.any(bad_i):
            I1a[bad_i] = 1.0
        y = y / I1a - tau

    return y

def VarT(tau, T):
    tau = np.asarray(tau, dtype=float)
    T   = np.asarray(T,   dtype=float)
    tau, T = np.broadcast_arrays(tau, T)
    tau = tau.astype(float).copy()
    T   = T.astype(float).copy()
    with np.errstate(all='ignore'):
        y = ((-2 * T**2 / tau * np.exp(T / tau) *
              (expint(2 * T / tau) - expint(T / tau)) * (T + 2 * tau) +
              (np.exp(-T / tau) - 1) * (T**2 + 2 * T * tau) -
              T**2 + tau**2) -
             (2 * T**2 / tau * np.exp(T / tau) *
              (expint(2 * T / tau) - expint(T / tau)) +
              (T - tau) * (2 - np.exp(-T / tau)))**2)
        bad = ~np.isfinite(y)
        if np.any(bad):
            y[bad] = tau[bad]**2
        I1a = (2 * T / tau * np.exp(T / tau) *
               (igamma(0, T / tau) - igamma(0, 2 * T / tau)) +
               np.exp(-T / tau) - 1)
        bad_i = ~np.isfinite(I1a)
        if np.any(bad_i):
            I1a[bad_i] = 1.0
        y = y / I1a
    return y

#count 
def computeNaNb(tau_a, tau_b, T_on, T_off):
    A = tau_a / (tau_a + tau_b)
    B = tau_b / (tau_a + tau_b)
    E = 0.0  # large-T limit
    P1 = A * (A + B * E);  P2 = A * (B - B * E)
    P3 = B * (A - A * E);  P4 = B * (B + A * E)
    TA = (2 * P1 + 1.5 * P2 + 1.5 * P3 + P4) * (tau_a + tau_b)
    TB = (2 * P4 + 1.5 * P2 + 1.5 * P3 + P1) * (tau_a + tau_b)
    Na = (TA + T_on) / (tau_a + tau_b) * np.exp(-TA / T_on)
    Nb = (TB + T_on) / (tau_a + tau_b) * np.exp(-TB / T_on)
    return Na, Nb

def Npts_from_countA(countA, tau_a, tau_b, T_on, T_off, f0):
    Na, _ = computeNaNb(tau_a, tau_b, T_on, T_off)
    M = countA / Na
    T = M * (T_on + T_off)
    return T * f0

# core computation

def computeDM(tau_a, tau_b, T_on, T_off, f0, N_pts):
    alpha   = 0.05
    z_alpha = -norm.ppf(alpha / 2)
    delta_t = 1.0 / f0
    T       = N_pts / f0
    cover   = 0.9

    M    = max(1.0, T / (T_on + T_off))
    Step = 500

    A = tau_a / (tau_a + tau_b)
    B = tau_b / (tau_a + tau_b)
    E = np.exp(-T * (tau_a + tau_b) / (tau_a * tau_b))
    P1 = A * (A + B * E);  P2 = A * (B - B * E)
    P3 = B * (A - A * E);  P4 = B * (B + A * E)
    TA = (2 * P1 + 1.5 * P2 + 1.5 * P3 + P4) * (tau_a + tau_b)
    TB = (2 * P4 + 1.5 * P2 + 1.5 * P3 + P1) * (tau_a + tau_b)

    xa_edges = (10.0**(0.01 * np.arange(0, Step + 1))) * 1e-2
    xb_edges = (10.0**(0.01 * np.arange(0, Step + 1))) * 1e-2
    xa = xa_edges[:Step];  dxa = np.diff(xa_edges)
    xb = xb_edges[:Step];  dxb = np.diff(xb_edges)

    Na = (TA + T_on) / (tau_a + tau_b) * np.exp(-TA / T_on)
    Nb = (TB + T_on) / (tau_a + tau_b) * np.exp(-TB / T_on)
    Ifactora = np.exp(-xa / T_on) / T_on * dxa
    Ifactorb = np.exp(-xb / T_on) / T_on * dxb

    def side(tau_x, tau_other, x, Ifactor, Nx):
        MissN = missN(tau_other, delta_t, cover)
        phit  = Phit(tau_other, delta_t)

        Bias = phit * (bM(tau_x, delta_t) + bT(tau_x, x))
        MSE  = phit * (VarM(tau_x, delta_t) + VarT(tau_x, x) +
                       (bM(tau_x, delta_t) + bT(tau_x, x))**2)

        if MissN > 0:
            k_vec = np.arange(1, MissN + 1).reshape(-1, 1)
            tau_k = (k_vec + 1) * tau_x
            TAU_mat, X_mat = np.meshgrid(tau_k.ravel(), x, indexing='ij')
            PORB_mat = np.tile(Porb1(tau_other, delta_t)**k_vec, (1, Step))
            BM_mat   = np.tile(bM(tau_k, delta_t), (1, Step))
            K_mat    = np.tile(k_vec, (1, Step))
            BT_mat   = bT(TAU_mat, X_mat)
            VT_mat   = VarT(TAU_mat, X_mat)
            bias_n = PORB_mat * phit * (BM_mat + BT_mat + K_mat * tau_x)
            MSE_n  = PORB_mat * phit * (VarM(TAU_mat, delta_t) + VT_mat +
                     (BM_mat + BT_mat)**2 + 2 * K_mat * tau_x * (BM_mat + BT_mat) +
                     (K_mat * tau_x)**2)
            Bias = Bias + np.sum(bias_n, axis=0)
            MSE  = MSE  + np.sum(MSE_n, axis=0)

        if MissN != 0:
            tau_last  = (MissN + 2) * tau_x
            PORB_last = Porb1(tau_other, delta_t)**(MissN + 1)
            BM_last   = bM(tau_last, delta_t)
            bt_last   = bT(tau_last, x)
            bias_n1 = PORB_last * (BM_last + bt_last + (MissN + 1) * tau_x)
            MSE_n1  = PORB_last * (VarM(tau_last, delta_t) + VarT(tau_last, x) +
                      (BM_last + bt_last)**2 +
                      2 * (MissN + 1) * tau_x * (BM_last + bt_last) +
                      ((MissN + 1) * tau_x)**2)
            Bias = Bias + bias_n1
            MSE  = MSE  + MSE_n1

        temp1 = (MSE - Bias**2) / Nx
        IVT = np.sum(temp1 * Ifactor)
        IBT = np.sum(Bias * Ifactor)
        return IVT, IBT

    IVTa, IBTa = side(tau_a, tau_b, xa, Ifactora, Na)
    IVTb, IBTb = side(tau_b, tau_a, xb, Ifactorb, Nb)

    Var_a = IVTa / M;  Bias_a = IBTa
    Var_b = IVTb / M;  Bias_b = IBTb
    Std_a = np.sqrt(max(0.0, Var_a))
    Std_b = np.sqrt(max(0.0, Var_b))

    Rerr_total   = z_alpha * Std_a / tau_a * 100 + z_alpha * Std_b / tau_b * 100
    Rebias_total = abs(Bias_a) / tau_a * 100 + abs(Bias_b) / tau_b * 100
    return Rebias_total, Rerr_total, Bias_a, Bias_b, Std_a, Std_b, M, Na, Nb

#  UI
col1, col2 = st.columns(2)
with col1:
    tau_a = st.number_input(r"$\tau_a$ (s)", min_value=1e-9, value=1.0, format="%.4f")
    T_on  = st.number_input(r"$T_{on}$ (s)", min_value=1e-9, value=4.0, format="%.4f")
    f0    = st.number_input(r"Sampling frequency $f_0$ (Hz)", min_value=1e-9, value=100.0, format="%.4f")
with col2:
    tau_b = st.number_input(r"$\tau_b$ (s)", min_value=1e-9, value=2.0, format="%.4f")
    T_off = st.number_input(r"$T_{off}$ (s)", min_value=1e-9, value=0.5, format="%.4f")

st.write("")
st.markdown("**Experiment scale — choose ONE input:**")
mode = st.radio(
    "input mode",
    ["Total sampling points (N)", "Dwell A average count"],
    horizontal=True,
    label_visibility="collapsed",
)


if mode == "Total sampling points (N)":
    N_pts_in = st.number_input("Total number of sampling points N pts", min_value=1, value=1000000, step=100000)
    N_pts = float(N_pts_in)
else:
    countA = st.number_input("Dwell A average count", min_value=1e-9, value=2180.0, format="%.4f")
    st.caption("Note: enter the Dwell A count. Dwell B is derived from the fixed Na:Nb ratio.")
    N_pts = None

st.write("")
if st.button("▶  Calculate", type="primary", use_container_width=True):
    try:
        if N_pts is None:
            N_pts_use = Npts_from_countA(countA, tau_a, tau_b, T_on, T_off, f0)
        else:
            N_pts_use = N_pts

        Rebias, Rerr, Bias_a, Bias_b, Std_a, Std_b, M, Na, Nb = \
            computeDM(tau_a, tau_b, T_on, T_off, f0, N_pts_use)

        st.divider()
        m1, m2 = st.columns(2)
        m1.metric("Total Relative Bias", f"{Rebias:.4f} %")
        m2.metric("Total Relative Variation", f"{Rerr:.4f} %")

        st.divider()
        c1, c2 = st.columns(2)
        with c1:
            st.markdown(r"$\tau_a$")
            st.write(f"Bias: `{Bias_a:.4e}`")
            st.write(f"Variation: `{Std_a:.4e}`")
        with c2:
            st.markdown(r"$\tau_b$")
            st.write(f"Bias: `{Bias_b:.4e}`")
            st.write(f"Variation: `{Std_b:.4e}`")


        st.divider()
        st.write(f"Dwell A average count: `{Na * M:.1f}`")
        st.write(f"Dwell B average count: `{Nb * M:.1f}`")

        if N_pts is None:
            st.caption(f"Derived total sampling points N = {N_pts_use:,.0f}")
        else:
            st.caption(f"Derived Dwell A count = {Na * M:.1f}   |   N pts = {N_pts_use:,.0f}")

    except Exception as e:
        st.error(f"Error: {e}")
