import streamlit as st

st.set_page_config(page_title="Dwell-time analysis tools", layout="centered")

st.title("Accuracy and precision limits of dwell-time analysis")
st.caption("Two-state and three-state Markov kinetics, "
           "Continuous Mode and DyeCycling Mode")
st.divider()

st.markdown("""
Pick a tool from the sidebar.

**Analytical error models** — closed-form bias and confidence interval,
results are instant.

- *CM model* — a single trace of fixed length
- *DM model* — alternating bright and dark windows

**Monte Carlo simulations** — trajectories are generated and re-analysed,
so these take from seconds to minutes.

- *CM / DM simulation* — two-state
- *CM / DM 3-state* — three-state, symmetric branching (p = 0.5)
""")
