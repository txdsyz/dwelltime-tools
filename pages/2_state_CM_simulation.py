import streamlit as st
import numpy as np

#st.set_page_config(page_title="CM Two-State Markov — Simulation", layout="centered")

st.title("CM Two-State Markov Process Simulation")
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


def _round_half_away(x):
    return int(np.floor(x + 0.5)) if x >= 0 else -int(np.floor(-x + 0.5))


# 2-State Simulation Core (ported from CM_simulation_GUI.m / runCMSim)
def runCMSim(ta0, tb0, N, M, L, dt, Ra, dtt, z_alpha, sigma0,
             sta_level, A_GUESS, B_GUESS, EMIS, minValid, rng):

    alph = 1 / ta0; beta = 1 / tb0
    Bc = beta / (alph + beta); Ac = alph / (alph + beta); Ec = np.exp(-(alph + beta) * dt)
    TRANS = np.array([[Bc + Ac * Ec, Ac - Ac * Ec],
                       [Bc - Bc * Ec, Ac + Bc * Ec]])

    Tau_a = np.full(M, np.nan)
    Tau_b = np.full(M, np.nan)
    cntA = np.full(M, np.nan)

    for k in range(M):
        # pool all dwells from the N traces of this run into one set, then take one mean
        allEa = []
        allEb = []
        for _tr in range(N):
            Lon = rng.exponential(L)
            N_tr = _round_half_away(Lon / dt)
            if N_tr < 20 * Ra:
                continue

            _, states = hmmgenerate(N_tr, TRANS, EMIS, rng)
            states = np.concatenate([states[1:], [states[-1]]])
            RawI = states.astype(float) - 1.0
            if sigma0 > 0:
                RawI = RawI + rng.normal(0, sigma0, N_tr)

            idx = np.arange(0, N_tr, Ra)
            Sa = RawI[idx]
            NN = Sa.size
            seq = np.zeros(NN, dtype=np.int64)
            for i in range(NN):
                seq[i] = int(np.argmin(np.abs(Sa[i] - sta_level))) + 1

            A_EST, B_EST = hmmtrain(seq, A_GUESS, B_GUESS)
            STATES = hmmviterbi(seq, A_EST, B_EST) - 1

            swpt = [0]
            for i in range(1, NN):
                if STATES[i] != STATES[i - 1]:
                    swpt.append(i)
            if len(swpt) < 3:
                continue

            Events = np.diff(swpt) * dtt
            Events = Events[1:]

            if STATES[swpt[1]] == 0:
                Ea = Events[0::2]; Eb = Events[1::2]
            else:
                Eb = Events[0::2]; Ea = Events[1::2]

            allEa.extend(Ea.tolist())
            allEb.extend(Eb.tolist())

        # one simulation instance k: mean is taken once over all traces' pooled events
        if allEa:
            Tau_a[k] = np.mean(allEa)
        if allEb:
            Tau_b[k] = np.mean(allEb)
        cntA[k] = len(allEa)

    nva = int(np.sum(~np.isnan(Tau_a)))
    nEvt = np.nanmean(cntA) if np.any(~np.isnan(cntA)) else np.nan

    ta = tb = sa = sb = ba = bb = Rsa = Rsb = Rba = Rbb = np.nan
    if nva >= minValid:
        ta = np.nanmean(Tau_a); tb = np.nanmean(Tau_b)
        sa = np.nanstd(Tau_a);  sb = np.nanstd(Tau_b)
        ba = ta - ta0; bb = tb - tb0
        Rsa = z_alpha * sa / ta0 * 100; Rsb = z_alpha * sb / tb0 * 100
        Rba = ba / ta0 * 100;           Rbb = bb / tb0 * 100

    return ta, tb, sa, sb, ba, bb, Rsa, Rsb, Rba, Rbb, nva, nEvt


# UI

col1, col2 = st.columns(2)
with col1:
    tau_a = st.number_input(r"$\tau_a$ (s)", min_value=1e-9, value=10.0, format="%.4f")
    L = st.number_input("Mean trace length L (s)", min_value=1e-9, value=100.0, format="%.4f")
    N_traces = st.number_input("Number of traces N", min_value=1, value=5, step=1)
with col2:
    tau_b = st.number_input(r"$\tau_b$ (s)", min_value=1e-9, value=10.0, format="%.4f")
    f0 = st.number_input(r"Sampling frequency $f_0$ (Hz)", min_value=1e-9, value=1.0, format="%.4f")
    M = st.number_input("Monte Carlo runs M", min_value=1, value=10, step=1)

st.caption("Traces are terminated at random, so individual lengths are drawn from an "
           "exponential distribution with mean L. Number of traces N is per simulation "
           "run; Monte Carlo runs M is the number of independent runs bias/variation are "
           "computed across.")

st.write("")
if st.button("▶  Run simulation", type="primary", use_container_width=True):
    try:
        dt = 1e-3
        dtt = 1 / f0
        Ra = _round_half_away(dtt / dt)
        if Ra < 1:
            raise ValueError("Sampling interval 1/f0 is shorter than the simulation step "
                              "dt = 1e-3 s. Reduce f0 below 1000 Hz.")

        alpha = 0.05
        from scipy.stats import norm
        z_alpha = -norm.ppf(alpha / 2)
        sigma0 = 0
        minValid = 5
        sta_num = 10
        sta_level = (np.arange(1, sta_num + 1)) * (1 / sta_num)
        A_GUESS = np.full((2, 2), 1 / 2)
        B_GUESS = np.array([[1, 0], [1, 0], [1, 0], [1, 0], [1, 0],
                            [0, 1], [0, 1], [0, 1], [0, 1], [0, 1]], dtype=float).T
        EMIS = np.eye(2)

        with st.spinner(f"Simulating {int(M)} Monte Carlo runs of {int(N_traces)} traces each..."):
            rng = np.random.default_rng()
            (ta, tb, sa, sb, ba, bb, Rsa, Rsb, Rba, Rbb, nva, nEvt) = runCMSim(
                tau_a, tau_b, int(N_traces), int(M), L, dt, Ra, dtt, z_alpha, sigma0,
                sta_level, A_GUESS, B_GUESS, EMIS, minValid, rng)

        st.divider()
        if np.isnan(ta):
            st.error(f"Not enough valid runs: {nva} of {int(M)} (minimum {minValid}). "
                      "The traces are too short, or the dwells are unresolvable at this f0. "
                      "Try a longer L, a higher f0, or larger N.")
        else:
            m1, m2 = st.columns(2)
            m1.metric("Total Relative Bias", f"{abs(Rba) + abs(Rbb):.4f} %")
            m2.metric("Total Relative Variation", f"{Rsa + Rsb:.4f} %")

            st.divider()
            c1, c2 = st.columns(2)
            with c1:
                st.markdown(r"$\tau_a$")
                st.write(f"true : `{tau_a:.4g}` s   measured : `{ta:.4g}` s")
                st.write(f"Bias: `{ba:.4e}`")
                st.write(f"Variation: `{sa:.4e}`")
            with c2:
                st.markdown(r"$\tau_b$")
                st.write(f"true : `{tau_b:.4g}` s   measured : `{tb:.4g}` s")
                st.write(f"Bias: `{bb:.4e}`")
                st.write(f"Variation: `{sb:.4e}`")

            st.divider()
            st.caption(f"Valid runs: {nva} / {int(M)}   |   mean pooled dwell-A events per run: {nEvt:.1f}")

    except Exception as ex:
        st.error(f"Error: {ex}")
