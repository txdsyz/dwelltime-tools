import streamlit as st

st.set_page_config(page_title="Dwell-time analysis tools", layout="centered")

st.title("Accuracy and precision limits of dwell-time analysis")
st.caption("Two-state and three-state Markov kinetics, Continuous Mode and DyeCycling Mode")
st.divider()

st.markdown("""
Pick a tool from the sidebar.

**Analytical error models** — closed-form bias and confidence interval, results are
instant.

- *CM model* — a single trace of fixed length
- *DM model* — alternating bright and dark windows, bias averaged over the
  exponential window-length distribution

**Monte Carlo simulations** — trajectories are generated and re-analysed, so these
take from seconds to minutes depending on the settings.

- *CM / DM simulation* — two-state
- *CM / DM 3-state* — three-state, symmetric branching (p = 0.5)
""")

st.info(
    "Simulation runtime scales with trace length x number of traces x Monte Carlo "
    "runs. The defaults here are reduced for interactive use; the production "
    "settings used in the paper are best run locally.",
    icon=None,
)