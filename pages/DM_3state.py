import streamlit as st
import numpy as np
from scipy.stats import norm
from scipy.linalg import expm

#st.set_page_config(page_title="DM Three-State Markov — Simulation Analysis", layout="centered")

st.title("DM Three-State Markov Process Uncertainty Analysis")
st.divider()


def hmmgenerate(L, TRANS, EMIS, rng):
    L = int(L)
    trc = np.cumsum(TRANS, axis=1)
    emc = np.cumsum(EMIS, axis=1)
    nS, nE = trc.shape[1], emc.shape[1]
    states = np.zeros(L, dtype=np.int64)
    seq = np.zeros(L, dtype=np.int64)
    r1 = rng.random(L)
    r2 = rng.random(L)
    cur = 0                                   # MATLAB starts in state 1
    for i in range(L):
        st = int(np.searchsorted(trc[cur], r1[i], side='right'))
        if st >= nS:
            st = nS - 1
        em = int(np.searchsorted(emc[st], r2[i], side='right'))
        if em >= nE:
            em = nE - 1
        states[i] = st + 1
        seq[i] = em + 1
        cur = st
    return seq, states


def _forward_backward(seq0, tr, e):
    numStates = tr.shape[0]
    L = len(seq0)
    fs = np.zeros((numStates, L + 1))
    fs[0, 0] = 1.0                             # assume start in state 1
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
        pos_s = s[1:][s[1:] > 0]
        loglik = float(np.sum(np.log(pos_s))) if pos_s.size else 0.0

        # expected transition counts
        pseudoTR = np.zeros_like(tr)
        for t in range(L - 1):
            if s[t + 2] > 0:
                pseudoTR += (fs[:, t + 1][:, None] * tr
                             * (e[:, seq0[t + 1]] * bs[:, t + 2])[None, :]) / s[t + 2]
        # expected emission counts
        gamma = fs[:, 1:] * bs[:, 1:]
        pseudoE = np.zeros_like(e)
        for k in range(numEmissions):
            posk = (seq0 == k)
            if posk.any():
                pseudoE[:, k] = gamma[:, posk].sum(axis=1)

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
        cand = v[:, None] + logTR            
        best = np.argmax(cand, axis=0)
        pTR[:, t] = best
        v = logE[:, seq0[t]] + cand[best, np.arange(numStates)]
    path = np.zeros(L, dtype=np.int64)
    path[L - 1] = int(np.argmax(v))
    for t in range(L - 2, -1, -1):
        path[t] = pTR[path[t + 1], t + 1]
    return path + 1


# 3-State DM Simulation Core

def calcStats(Tau_data, tau_meas, z_alpha):
    if Tau_data.size == 0:
        return np.nan, np.nan, np.nan, np.nan
    tau_sim_mean = np.mean(Tau_data)
    std_val = np.std(Tau_data, ddof=1) if Tau_data.size > 1 else 0.0
    bias = tau_sim_mean - tau_meas
    Restd = z_alpha * std_val / tau_meas * 100
    Rebi = bias / tau_meas * 100
    return bias, std_val, Rebi, Restd


def eval3StateDMError(tau_a, tau_b, tau_c, T_on, T_off, T, f0, M, rng, progress=None):
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
    sta_level = (np.arange(1, sta_num + 1)) * ((3 - 0) / sta_num)
    pro = 0.5      # symmetric exit split

    TB_win = (1 / T_on) / (1 / T_on + 1 / T_off)
    TA_win = (1 / T_off) / (1 / T_on + 1 / T_off)
    TE_win = np.exp(-(1 / T_on + 1 / T_off) * dt)
    T_TRANS = np.array([[TB_win + TA_win * TE_win, TA_win - TA_win * TE_win],
                        [TB_win - TB_win * TE_win, TA_win + TB_win * TE_win]])
    T_EMIS = np.array([[1.0, 0.0], [0.0, 1.0]])

    a0 = 1 / tau_a; a1 = a0 * pro
    b0 = 1 / tau_b; b1 = b0 * pro
    c0 = 1 / tau_c; c1 = c0 * pro
    Q = np.array([[-a0, a1, a0 - a1],
                  [b1, -b0, b0 - b1],
                  [c1, c0 - c1, -c0]])
    CQ = Q * dt
    TRANS = expm(CQ)
    EMIS = np.eye(3)

    Tau_a_sim = np.zeros(M)
    Tau_b_sim = np.zeros(M)
    Tau_c_sim = np.zeros(M)

    A_GUESS = np.full((4, 4), 1 / 4)
    B_GUESS = np.array([[1, 0, 0, 0], [1, 0, 0, 0],
                        [0, 1, 0, 0], [0, 1, 0, 0], [0, 1, 0, 0],
                        [0, 0, 1, 0], [0, 0, 1, 0], [0, 0, 1, 0],
                        [0, 0, 0, 1], [0, 0, 0, 1]], dtype=float).T

    Na_sum = 0
    Nb_sum = 0
    Nc_sum = 0

    for k in range(M):
        if progress is not None:
            progress(k, M)

        _, states = hmmgenerate(N, TRANS, EMIS, rng)
        states = np.concatenate([states, [states[-1]]])[1:]

        _, Tstates = hmmgenerate(N, T_TRANS, T_EMIS, rng)
        Tstates = np.concatenate([Tstates, [Tstates[-1]]])[1:]
        WindowI = Tstates - 1                  # 0 = off, 1 = on

        RawI = states * WindowI                # 0 = off, 1 = a, 2 = b, 3 = c
        if sigma0 > 0:
            RawI = RawI + rng.normal(0, sigma0, N)

        Sa_RawI = RawI[indx]
        seq = np.argmin(np.abs(Sa_RawI[:, None] - sta_level[None, :]), axis=1) + 1

        try:
            A_EST, B_EST = hmmtrain(seq, A_GUESS, B_GUESS)
            STATES = hmmviterbi(seq, A_EST, B_EST) - 1     # 0=dark, 1=a, 2=b, 3=c
        except Exception:
            Tau_a_sim[k] = np.nan; Tau_b_sim[k] = np.nan; Tau_c_sim[k] = np.nan
            continue

        swpt = [0]
        Inista = [int(STATES[0])]
        for i in range(1, NN):
            if STATES[i] != STATES[i - 1]:
                swpt.append(i)
                Inista.append(int(STATES[i]))
        swpt = np.array(swpt)

        if swpt.size < 3:
            Tau_a_sim[k] = np.nan; Tau_b_sim[k] = np.nan; Tau_c_sim[k] = np.nan
            continue

        Events = np.diff(swpt) * dtt           # length = numel(swpt) - 1
        Inista = np.array(Inista[:-1])

        if Inista[0] == 0:
            blank = [1, 2]
        elif Inista.size >= 2 and Inista[1] == 0:
            blank = [1, 2, 3]
        else:
            blank = [1]
        for i in range(3, swpt.size + 1):    
            if i - 1 < Inista.size and Inista[i - 1] == 0:
                blank += [i - 1, i, i + 1]
        blank = np.unique([b for b in blank if b <= Events.size])
        keep = np.setdiff1d(np.arange(1, Events.size + 1), blank) - 1
        Events = Events[keep]
        Inista = Inista[keep]

        if Events.size == 0:
            Tau_a_sim[k] = np.nan; Tau_b_sim[k] = np.nan; Tau_c_sim[k] = np.nan
            continue

        Event_a = Events[Inista == 1]
        Event_b = Events[Inista == 2]
        Event_c = Events[Inista == 3]

        Na_sum += Event_a.size
        Nb_sum += Event_b.size
        Nc_sum += Event_c.size

        Tau_a_sim[k] = np.mean(Event_a) if Event_a.size else np.nan
        Tau_b_sim[k] = np.mean(Event_b) if Event_b.size else np.nan
        Tau_c_sim[k] = np.mean(Event_c) if Event_c.size else np.nan

    Na = Na_sum / M
    Nb = Nb_sum / M
    Nc = Nc_sum / M

    Tau_a_sim = Tau_a_sim[~(np.isnan(Tau_a_sim) | (Tau_a_sim <= 0))]
    Tau_b_sim = Tau_b_sim[~(np.isnan(Tau_b_sim) | (Tau_b_sim <= 0))]
    Tau_c_sim = Tau_c_sim[~(np.isnan(Tau_c_sim) | (Tau_c_sim <= 0))]

    biasa, stda, Rebia, Restda = calcStats(Tau_a_sim, tau_a, z_alpha)
    biasb, stdb, Rebib, Restdb = calcStats(Tau_b_sim, tau_b, z_alpha)
    biasc, stdc, Rebic, Restdc = calcStats(Tau_c_sim, tau_c, z_alpha)

    Restd_total = np.nansum([Restda, Restdb, Restdc])
    Rebi_total = np.nansum([abs(Rebia), abs(Rebib), abs(Rebic)])

    return (Rebi_total, Restd_total, biasa, biasb, biasc,
            stda, stdb, stdc, Na, Nb, Nc)


def formatOutput(val, unit):
    if not np.isfinite(val):
        return 'N/A (Insufficient events)'
    if unit == '%':
        return f'{val:.4f} %'
    return f'{val:.6e} {unit}'



# UI

col1, col2 = st.columns(2)
with col1:
    tau_m_a = st.number_input(r"Measured $\tau_a$ (s)", min_value=1e-9, value=0.1, format="%.4f")
    tau_m_c = st.number_input(r"Measured $\tau_c$ (s)", min_value=1e-9, value=1.0, format="%.4f")
    T_off = st.number_input(r"Time $T_{off}$ (s)", min_value=1e-9, value=4.0, format="%.4f")
    M = st.number_input("Monte Carlo Iterations (M)", min_value=1, value=5, step=5)
with col2:
    tau_m_b = st.number_input(r"Measured $\tau_b$ (s)", min_value=1e-9, value=0.1, format="%.4f")
    T_on = st.number_input(r"Time $T_{on}$ (s)", min_value=1e-9, value=42.0, format="%.4f")
    f0 = st.number_input(r"Sampling freq $f_0$ (Hz)", min_value=1e-9, value=1.0, format="%.4f")

st.write("")
st.markdown("**Experiment scale — choose ONE input:**")
mode = st.radio(
    "input mode",
    ["Time T (s)", "Dwell count"],
    horizontal=True,
    label_visibility="collapsed",
)

tau_sum = tau_m_a + tau_m_b + tau_m_c
duty = T_on / (T_on + T_off)

if mode == "Time T (s)":
    T_in = st.number_input("Time T (s)", min_value=1e-9, value=1100.0, format="%.4f")
    count_in = None
else:
    count_in = st.number_input("Dwell average count", min_value=1e-9,
                               value=1100.0 * duty / (tau_sum), format="%.4f")
    st.caption("Expected number of complete dwells per state (same for A, B and C).")
    T_in = None


st.write("")
if st.button("▶  Calculate", type="primary", use_container_width=True):
    try:
        T = T_in if T_in is not None else count_in * tau_sum / duty
        N_pts = T * f0

        bar = st.progress(0.0, text="Running...")

        def _prog(k, tot):
            bar.progress(k / tot, text=f"run {k + 1} / {tot}")

        rng = np.random.default_rng()
        (Rebias, Rerr, BiasA, BiasB, BiasC,
         StdA, StdB, StdC, Na, Nb, Nc) = eval3StateDMError(
            tau_m_a, tau_m_b, tau_m_c, T_on, T_off, T, f0, int(M), rng, _prog)
        bar.empty()

        st.divider()
        m1, m2 = st.columns(2)
        m1.metric("Total Expected Relative Bias", f"{Rebias:.4f} %")
        m2.metric("Total Expected Variation", formatOutput(Rerr, '%'))

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
