import streamlit as st
import numpy as np
from scipy.special import expi
from scipy.stats import norm

#st.set_page_config(page_title="CM Two-State Markov — Uncertainty Analysis", layout="centered")

st.title("CM Two-State Markov Process Uncertainty Analysis")
st.divider()

# special functions

def expint(x):
    x = np.asarray(x, dtype=float)
    return -expi(-x)

def igamma(a, x):
    return expint(x)

#  dependency functions (from bM.m, bT.m, VarM.m, VarT.m)

def bM(tau, delta_t):
    tau = np.asarray(tau, dtype=float)
    A = np.exp(-delta_t / tau)
    return delta_t * (A + 1) / (1 - A) - 2 * tau

def _bT_scalar(tau, T):
    y = (2 * T**2 / tau * np.exp(T / tau) *
         (expint(2 * T / tau) - expint(T / tau)) +
         (T - tau) * (2 - np.exp(-T / tau)) + tau)
    if not np.isfinite(y):
        y = np.exp(-T / tau) * (tau / 2 + tau**2 / (2 * T)) - 4 * tau**2 / T + tau
    I1a = (2 * T / tau * np.exp(T / tau) *
           (igamma(0, T / tau) - igamma(0, 2 * T / tau)) +
           np.exp(-T / tau) - 1)
    if not np.isfinite(I1a):
        I1a = 1.0
    return y / I1a - tau

def bT(tau, T):
    tau = np.atleast_1d(np.asarray(tau, dtype=float))
    T   = np.atleast_1d(np.asarray(T,   dtype=float))
    if tau.size == 1 and T.size > 1: tau = np.full(T.shape, tau[0])
    if T.size == 1 and tau.size > 1: T = np.full(tau.shape, T[0])
    return np.array([_bT_scalar(t, tt) for t, tt in zip(tau, T)])

def VarM(tau, delta_t):
    tau = np.asarray(tau, dtype=float)
    A = np.exp(-delta_t / tau)
    return delta_t**2 * A / (6 * tau / delta_t * (1 - A))

def _VarT_scalar(tau, T):
    I1a = (2 * T / tau * np.exp(T / tau) *
           (igamma(0, T / tau) - igamma(0, 2 * T / tau)) +
           np.exp(-T / tau) - 1)
    if not np.isfinite(I1a):
        I1a = 1.0
    y = ((-2 * T**3 / tau * np.exp(T / tau) *
          (expint(2 * T / tau) - expint(T / tau)) +
          (2 * tau**2 + T**2) * (np.exp(-T / tau) - 1) +
          T * (2 * tau - T)) -
         (2 * T**2 / tau * np.exp(T / tau) *
          (expint(2 * T / tau) - expint(T / tau)) +
          (T - tau) * (1 - np.exp(-T / tau)) + T)**2)
    y = y / I1a
    if not np.isfinite(y):
        y = tau**2
    return y

def VarT(tau, T):
    tau = np.atleast_1d(np.asarray(tau, dtype=float))
    T   = np.atleast_1d(np.asarray(T,   dtype=float))
    if tau.size == 1 and T.size > 1: tau = np.full(T.shape, tau[0])
    if T.size == 1 and tau.size > 1: T = np.full(tau.shape, T[0])
    return np.array([_VarT_scalar(t, tt) for t, tt in zip(tau, T)])

def Porb1(tau, delta_t):
    return 1 + tau / delta_t * (np.exp(-delta_t / tau) - 1)

def Phit(tau, delta_t):
    return tau / delta_t * (1 - np.exp(-delta_t / tau))

def missN(tau, delta_t, cover=0.9):
    val = np.log(1 - cover) / np.log(1 + tau / delta_t * (np.exp(-delta_t / tau) - 1))
    return max(0, int(val) - 1)


def computeDiscard(T, tau_a, tau_b):
    A = tau_a / (tau_a + tau_b)
    B = tau_b / (tau_a + tau_b)
    E = np.exp(-T * (tau_a + tau_b) / (tau_a * tau_b))
    Pa = A * (A + B * E);  Pb = A * (B - B * E)
    Pc = B * (A - A * E);  Pd = B * (B + A * E)
    return Pa * tau_a + Pd * tau_b + (Pb + Pc) * (tau_a + tau_b) / 2

def T_from_count(count_target, tau_a, tau_b):
    T = count_target * (tau_a + tau_b)
    for _ in range(200):
        discard_T = computeDiscard(T, tau_a, tau_b)
        T_new = count_target * (tau_a + tau_b) + discard_T
        if abs(T_new - T) < 1e-13:
            T = T_new
            break
        T = T_new
    return T

#  core computation

def computeMarkov(tau_a, tau_b, f0, T):
    alpha   = 0.05
    z_alpha = -norm.ppf(alpha / 2)
    delta_t = 1.0 / f0
    cover   = 0.9

    discard_T = computeDiscard(T, tau_a, tau_b)
    N = max(0.0, (T - discard_T) / (tau_a + tau_b))

    def one_side(tau_x, tau_other):
        MissN = missN(tau_other, delta_t, cover)
        bm   = float(bM(tau_x, delta_t))
        bt   = float(bT(tau_x, T)[0])
        phit = float(Phit(tau_other, delta_t))

        Bias = phit * (bm + bt)
        MSE  = phit * (float(VarM(tau_x, delta_t)) + float(VarT(tau_x, T)[0]) + (bm + bt)**2)

        if MissN > 0:
            k      = np.arange(1, MissN + 1)
            tau_k  = (k + 1) * tau_x
            porb_k = Porb1(tau_other, delta_t)**k
            bm_v   = bM(tau_k, delta_t)
            bt_v   = bT(tau_k, T)
            vm_v   = VarM(tau_k, delta_t)
            vt_v   = VarT(tau_k, T)
            Bias += np.sum(porb_k * phit * (bm_v + bt_v + k * tau_x))
            MSE  += np.sum(porb_k * phit * (vm_v + vt_v + (bm_v + bt_v)**2 +
                          2 * k * tau_x * (bm_v + bt_v) + (k * tau_x)**2))

        if MissN != 0:
            pl  = Porb1(tau_other, delta_t)**(MissN + 1)
            tl  = (MissN + 2) * tau_x
            bml = float(bM(tl, delta_t))
            btl = float(bT(tl, T)[0])
            Bias += pl * (bml + btl + (MissN + 1) * tau_x)
            MSE  += pl * (float(VarM(tl, delta_t)) + float(VarT(tl, T)[0]) +
                         (bml + btl)**2 + 2 * (MissN + 1) * tau_x * (bml + btl) +
                         ((MissN + 1) * tau_x)**2)

        Std = np.sqrt(max(0.0, (MSE - Bias**2) / N)) if N > 0 else 0.0
        return Bias, Std

    Bias_a, Std_a = one_side(tau_a, tau_b)
    Bias_b, Std_b = one_side(tau_b, tau_a)

    Rerr_total   = z_alpha * Std_a / tau_a * 100 + z_alpha * Std_b / tau_b * 100
    Rebias_total = abs(Bias_a) / tau_a * 100 + abs(Bias_b) / tau_b * 100
    return Rebias_total, Rerr_total, Bias_a, Bias_b, Std_a, Std_b, N

# UI

col1, col2 = st.columns(2)
with col1:
    tau_a = st.number_input(r"$\tau_a$ (s)", min_value=1e-9, value=10.0, format="%.4f")
    f0    = st.number_input(r"Sampling frequency f0 (Hz)", min_value=1e-9, value=100.0, format="%.4f")
with col2:
    tau_b = st.number_input(r"$\tau_b$ (s)", min_value=1e-9, value=10.0, format="%.4f")

st.write("")
st.markdown("**Experiment scale — choose ONE input:**")
mode = st.radio(
    "input mode",
    ["Total sampling points (N)", "Dwell average count"],
    horizontal=True,
    label_visibility="collapsed",
)

if mode == "Total sampling points (N)":
    N_pts = st.number_input("Total number of sampling points N pts", min_value=1, value=10000, step=1000)
    T_input = N_pts / f0
    count_input = None
else:
    count_input = st.number_input("Dwell average count (same for A & B)",
                                  min_value=1e-9, value=4.5, format="%.4f")
    st.caption("Note: one count value is used for both state A and state B.")
    T_input = None

st.write("")
if st.button("▶  Calculate", type="primary", use_container_width=True):
    try:
        if T_input is None:
            T = T_from_count(count_input, tau_a, tau_b)
        else:
            T = T_input

        Rebias, Rerr, Bias_a, Bias_b, Std_a, Std_b, N = computeMarkov(tau_a, tau_b, f0, T)

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
        if T_input is None:
            st.caption(f"Derived total sampling points N = {T * f0:,.0f}   (dwell count = {N:.4f})")
        else:
            st.caption(f"Dwell A/B average count = {N:.4f}")

    except Exception as e:
        st.error(f"Error: {e}")
