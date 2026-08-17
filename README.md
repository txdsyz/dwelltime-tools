# Dwell-time analysis tools

Interactive tools for the accuracy and precision limits of dwell-time analysis of
two-state and three-state Markov processes, in Continuous Mode (CM) and
DyeCycling Mode (DM).

## Tools

**Analytical error models** — closed-form bias and 95 % confidence interval.
Results are instant.

| Page | Model |
|---|---|
| CM model | single trace of fixed length |
| DM model | alternating bright/dark windows, bias averaged over the exponential window-length distribution |

**Monte Carlo simulations** — trajectories are generated, sampled, decoded with
`hmmtrain`/`hmmviterbi`-equivalent routines, and cut into dwells.

| Page | Model |
|---|---|
| CM simulation | two-state, traces terminated by photobleaching |
| DM simulation | two-state, bright/dark windows |
| CM 3-state | three-state, symmetric branching (p = 0.5) |
| DM 3-state | three-state with dye cycling |

## Run locally

```
pip install -r requirements.txt
streamlit run streamlit_app.py
```

Requires Python 3.9 or newer.

## Notes

- Simulation runtime scales as (trace length / dt) x (number of traces) x
  (Monte Carlo runs). The defaults in the app are reduced so that a run finishes
  interactively; the production settings used in the paper are considerably
  larger and are best run locally or on a cluster.
- The three-state models assume symmetric branching: every state splits its exit
  rate equally between the two other states. This satisfies detailed balance and
  gives all three states the same event rate.
- The simulation pages report the number of dwells actually recovered by the
  decoder, which is smaller than the expected number of true dwells. The
  difference is the missed-event effect these tools are meant to quantify.

## License

MIT
