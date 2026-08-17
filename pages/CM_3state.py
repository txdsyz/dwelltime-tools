#%%
import streamlit as st
import numpy as np
from scipy.stats import norm

#st.set_page_config(page_title="CM Three-State Markov — Simulation Analysis", layout="centered")

st.title("CM Three-State Markov Process Uncertainty")
st.caption("simulation-based")
st.divider()


def hmmgenerate(L, TRANS, EMIS, rng):
    L = int(L)
    trc = np.cumsum(TRANS, axis=1)
    emc = np.cumsum(EMIS, axis=1)
    states = np.zeros(L, dtype=np.int64)
    seq = np.zeros(L, dtype=np.int64)
    r1 = rng.random(L)
    r2 = rng.random(L)
    cur = 0                                   # the original MATLAB starts in state 1
    for i in range(L):
        st = int(np.searchsorted(trc[cur], r1[i], side='right'))
        if st >= trc.shape[1]:
            st = trc.shape[1] - 1
        em = int(np.searchsorted(emc[st], r2[i], side='right'))
        if em >= emc.shape[1]:
            em = emc.shape[1] - 1
        states[i] = st + 1
        seq[i] = em + 1
        cur = st
    return seq, states


def _forward_backward(seq0, tr, e):
    numStates = tr.shape[0]
    L = len(seq0)
    fs = np.zeros((numStates, L + 1))
    fs[0, 0] = 1.0                             # assume we start in state 1
    s = np.ones(L + 1)
    for t in range(1, L + 1):
        fs[:, t] = e[:, seq0[t - 1]] * (fs[:, t - 1] @ tr)
        s[t] = fs[:, t].sum()
        if s[t] > 0:
            fs[:, t] /= s[t]
    bs = np.ones((numStates, L + 1))
    for t in range(L - 1, -1, -1):
        if s[t + 1] > 0:
            bs[:, t] = (tr @ (bs[:, t + 1] * e[:, seq0[t]])) / s[t + 1]
    return fs, bs, s


def hmmtrain(seq, guessTR, guessE, maxiter=500, tol=1e-6):
    seq0 = np.asarray(seq, dtype=np.int64) - 1
    tr = np.array(guessTR, dtype=float)
    e = np.array(guessE, dtype=float)
    numStates, numEmissions = e.shape
    L = len(seq0)
    loglik = 1.0
    for _ in range(int(maxiter)):
        oldLL, oldTR, oldE = loglik, tr.copy(), e.copy()
        fs, bs, s = _forward_backward(seq0, tr, e)
        loglik = np.sum(np.log(s[1:][s[1:] > 0]))

        # expected transition counts
        pseudoTR = np.zeros_like(tr)
        for t in range(L - 1):
            pseudoTR += (fs[:, t + 1][:, None] * tr
                         * (e[:, seq0[t + 1]] * bs[:, t + 2])[None, :]) / s[t + 2]
        # expected emission counts
        gamma = fs[:, 1:] * bs[:, 1:]
        pseudoE = np.zeros_like(e)
        for k in range(numEmissions):
            pos = (seq0 == k)
            if pos.any():
                pseudoE[:, k] = gamma[:, pos].sum(axis=1)

        rs = pseudoTR.sum(axis=1); rs[rs == 0] = 1.0
        tr = pseudoTR / rs[:, None]
        rs = pseudoE.sum(axis=1); rs[rs == 0] = 1.0
        e = pseudoE / rs[:, None]

        if (abs(loglik - oldLL) / (1 + abs(oldLL)) < tol
                and np.abs(tr - oldTR).max() / numStates < tol
                and np.abs(e - oldE).max() / numEmissions < tol):
            break
    return tr, e


def hmmviterbi(seq, tr, e):
    seq0 = np.asarray(seq, dtype=np.int64) - 1
    numStates = tr.shape[0]
    L = len(seq0)
    if L == 0:
        return np.zeros(0, dtype=np.int64)
    with np.errstate(divide='ignore'):
        logTR = np.log(tr)
        logE = np.log(e)
    pTR = np.zeros((numStates, L), dtype=np.int64)
    v = np.full(numStates, -np.inf)
    v[0] = 0.0                                 # start in state 1
    for t in range(L):
        cand = v[:, None] + logTR              # (from, to)
        best = np.argmax(cand, axis=0)
        pTR[:, t] = best
        v = logE[:, seq0[t]] + cand[best, np.arange(numStates)]
    path = np.zeros(L, dtype=np.int64)
    path[L - 1] = int(np.argmax(v))
    for t in range(L - 2, -1, -1):
        path[t] = pTR[path[t + 1], t + 1]
    return path + 1



# 3-State Simulation Core
def calcStats(Tau_data, tau_true, z_alpha):
    if Tau_data.size == 0:
        return np.nan, np.nan, np.nan, np.nan
    tau_mean = np.mean(Tau_data)
    std_val = np.std(Tau_data, ddof=1) if Tau_data.size > 1 else 0.0
    bias = tau_mean - tau_true
    Restd = z_alpha * std_val / tau_true * 100
    Rebi = bias / tau_true * 100
    return bias, std_val, Rebi, Restd


def ifempty(val, default):
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return default
    return val


def run3StateSim(tau_a0, tau_b0, tau_c0, T, f0, M, rng):
    from scipy.linalg import expm

    alpha = 0.05
    z_alpha = -norm.ppf(alpha / 2)
    dt = 1e-1
    N = int(np.fix(T / dt))
    sigma0 = 0
    dtt = 1 / f0
    Ra = int(np.fix(dtt / dt))
    NN = int(np.fix(T / dtt))
    indx = np.arange(NN) * Ra                 

    sta_num = 10
    sta_level = (np.arange(1, sta_num + 1)) * ((2 - 0) / sta_num)
    pro = 0.5      # symmetric exit split

    a0 = 1 / tau_a0; a1 = a0 * pro
    b0 = 1 / tau_b0; b1 = b0 * pro
    c0 = 1 / tau_c0; c1 = c0 * pro

    Q = np.array([[-a0, a1, a0 - a1],
                  [b1, -b0, b0 - b1],
                  [c1, c0 - c1, -c0]])
    CQ = Q * dt
    TRANS = expm(CQ)
    EMIS = np.eye(3)

    A_GUESS = np.full((3, 3), 1 / 3)
    B_GUESS = np.array([[1, 0, 0], [1, 0, 0], [1, 0, 0],
                        [0, 1, 0], [0, 1, 0], [0, 1, 0], [0, 1, 0],
                        [0, 0, 1], [0, 0, 1], [0, 0, 1]], dtype=float).T

    Na_sum = 0
    Nb_sum = 0
    Nc_sum = 0

    Tau_a = np.zeros(M); Tau_b = np.zeros(M); Tau_c = np.zeros(M)

    for k in range(M):
        _, states = hmmgenerate(N, TRANS, EMIS, rng)
        # M is the Monte Carlo count

        states = np.concatenate([states, [states[M - 1]]])
        states = states[1:]
        RawI = states + rng.normal(0, sigma0, N) - 1 if sigma0 > 0 else states - 1

        Sa_RawI = RawI[indx]
        seq = np.zeros(NN, dtype=np.int64)
        for i in range(NN):
            diffs = np.abs(Sa_RawI[i] - sta_level)
            seq[i] = int(np.argmin(diffs)) + 1

        A_EST, B_EST = hmmtrain(seq, A_GUESS, B_GUESS)
        STATES = hmmviterbi(seq, A_EST, B_EST) - 1

        # Measure Events
        swpt = [0]
        stmark = []
        for i in range(1, NN):
            if STATES[i] != STATES[i - 1]:
                swpt.append(i)
                stmark.append(STATES[i])
        swpt = np.array(swpt); stmark = np.array(stmark, dtype=np.int64)

        Events = np.diff(swpt) * dtt
        No = Events.size

        if No <= 1:
            Tau_a[k] = -1; Tau_b[k] = -1; Tau_c[k] = -1
        elif No == 2:
            Tau_a[k] = -1; Tau_b[k] = -1; Tau_c[k] = -1
            if stmark[0] == 0: Tau_a[k] = Events[1]
            if stmark[0] == 1: Tau_b[k] = Events[1]
            if stmark[0] == 2: Tau_c[k] = Events[1]
        else:
            Events = Events[1:]
            stmark = stmark[:-1]

            Event_a = Events[stmark == 0]
            Event_b = Events[stmark == 1]
            Event_c = Events[stmark == 2]

            Na_sum += Event_a.size
            Nb_sum += Event_b.size
            Nc_sum += Event_c.size

            Tau_a[k] = ifempty(np.mean(Event_a) if Event_a.size else np.nan, -1)
            Tau_b[k] = ifempty(np.mean(Event_b) if Event_b.size else np.nan, -1)
            Tau_c[k] = ifempty(np.mean(Event_c) if Event_c.size else np.nan, -1)

    Na = Na_sum / M
    Nb = Nb_sum / M
    Nc = Nc_sum / M

    Tau_a = Tau_a[Tau_a > 0]
    Tau_b = Tau_b[Tau_b > 0]
    Tau_c = Tau_c[Tau_c > 0]

    biasa, stda, Rebia, Restda = calcStats(Tau_a, tau_a0, z_alpha)
    biasb, stdb, Rebib, Restdb = calcStats(Tau_b, tau_b0, z_alpha)
    biasc, stdc, Rebic, Restdc = calcStats(Tau_c, tau_c0, z_alpha)

    Restd_total = np.nansum([Restda, Restdb, Restdc])
    Rebi_total = np.nansum([abs(Rebia), abs(Rebib), abs(Rebic)])

    return (Rebi_total, Restd_total, biasa, biasb, biasc,
            stda, stdb, stdc, Na, Nb, Nc)



# UI

col1, col2 = st.columns(2)
with col1:
    tau_a0 = st.number_input(r"$\tau_a$ (s)", min_value=1e-9, value=10.0, format="%.4f")
    tau_c0 = st.number_input(r"$\tau_c$ (s)", min_value=1e-9, value=2.0, format="%.4f")
    M = st.number_input("Monte Carlo Iterations (M)", min_value=1, value=50, step=10)
with col2:
    tau_b0 = st.number_input(r"$\tau_b$ (s)", min_value=1e-9, value=5.0, format="%.4f")
    f0 = st.number_input(r"Sampling frequency $f_0$ (Hz)", min_value=1e-9, value=1.0, format="%.4f")

st.write("")
st.markdown("**Experiment scale — choose ONE input:**")
mode = st.radio(
    "input mode",
    ["Time T (s)", "Dwell count"],
    horizontal=True,
    label_visibility="collapsed",
)

tau_sum = tau_a0 + tau_b0 + tau_c0

if mode == "Time T (s)":
    T_in = st.number_input("Time T (s)", min_value=1e-9, value=500.0, format="%.4f")
    count_in = None
else:
    count_in = st.number_input("Expected dwell average count", min_value=1e-9,
                               value=29.4, format="%.4f")
    st.caption("Expected (true) number of complete dwells per state. Because p = 0.5 "
               "the three states share the same event rate, so one value applies to "
               "A, B and C. The counts reported below are the events actually "
               "recovered by Viterbi, which are fewer.")
    T_in = None

st.write("")
if st.button("▶  Calculate", type="primary", use_container_width=True):
    try:
        T = T_in if T_in is not None else count_in * tau_sum

        with st.spinner(f"Running ..."):
            rng = np.random.default_rng()
            (Rebias, Rerr, BiasA, BiasB, BiasC,
             StdA, StdB, StdC, Na, Nb, Nc) = run3StateSim(
                tau_a0, tau_b0, tau_c0, T, f0, int(M), rng)

        st.divider()
        m1, m2 = st.columns(2)
        m1.metric("Total Relative Bias", f"{Rebias:.4f} %")
        m2.metric("Total Relative Variation", f"{Rerr:.4f} %")

        st.divider()
        c1, c2, c3 = st.columns(3)
        for col, lab, bias, sd in ((c1, r"$\tau_a$", BiasA, StdA),
                                   (c2, r"$\tau_b$", BiasB, StdB),
                                   (c3, r"$\tau_c$", BiasC, StdC)):
            with col:
                st.markdown(lab)
                st.write(f"Bias: `{bias:.6e}`")
                st.write(f"Variation: `{sd:.6e}`")

        st.divider()
        st.write(f"Dwell A average count : `{Na:.1f}`")
        st.write(f"Dwell B average count : `{Nb:.1f}`")
        st.write(f"Dwell C average count : `{Nc:.1f}`")

    except Exception as ex:
        st.error(f"Error: {ex}")
