# Laminar flame speed and NOx formation in NH₃–H₂–air blends

Computational Methods in Combustion (MKWS) project — Cantera study of the
trade-off between laminar burning velocity and NOx emission for ammonia
flames enriched with hydrogen.

## Motivation

Ammonia (NH₃) is a carbon-free fuel and a practical hydrogen carrier, but it
burns slowly (S_L ≈ 6–7 cm/s at stoichiometric vs. ~37 cm/s for methane) and
tends to produce high NOx because of its fuel-bound nitrogen. Blending in
hydrogen accelerates the flame dramatically, but can raise NOx. This project
maps that trade-off across three control variables and locates favourable
operating points.

## Control variables

| Variable | Symbol | Range |
|----------|--------|-------|
| H₂ fraction in the NH₃/H₂ fuel | `X_H2` | 0.0 – 0.6 (mole fraction) |
| Equivalence ratio | `phi` | 0.7 – 1.3 |
| Preheat (unburned) temperature | `T_in` | 300, 450, 600 K |

All cases at p = 1 atm, air = O₂:1, N₂:3.76 (molar), GRI-Mech 3.0 mechanism
(`gri30.yaml`, includes the full N-chemistry: NO, NO₂, N₂O).

## Method

* 1-D freely-propagating premixed flame (`cantera.FreeFlame`).
* Mixture-averaged transport.
* **Continuation**: along each φ-sweep the previous converged flame seeds the
  next case, which speeds up the sweep and stabilises it near the lean/rich
  limits.
* NOx is read as (NO + NO₂) mole fraction in the burned gas, reported in ppm.

## Files

| File | Purpose |
|------|---------|
| `flame_nh3_h2.py` | Core model: `solve_flame()` and the full `run_sweep()`. |
| `run_chunk.py` | Runs one `T_in` level, appending each row to a CSV (resumable). |
| `plots.py` | Builds all six figures from `results/results.csv`. |
| `results/results.csv` | Tidy table of every converged case. |
| `figures/*.png` | Output figures (see below). |

## How to run

```bash
pip install cantera numpy pandas matplotlib

# Option A — single run of the whole grid
python3 flame_nh3_h2.py

# Option B — one preheat level at a time (resumable; good for slow machines)
python3 run_chunk.py 300
python3 run_chunk.py 450
python3 run_chunk.py 600
# then merge the chunk_*.csv files into results/results.csv

# Figures
python3 plots.py
```

## Figures

* `fig1_SL_contour.png` — S_L over the (φ, X_H₂) plane at T_in = 300 K.
* `fig2_NOx_contour.png` — NOx over the same plane.
* `fig3_SL_vs_XH2.png` — flame-speed gain from H₂ enrichment, per φ.
* `fig4_pareto.png` — NOx vs S_L trade-off, coloured by φ, with the Pareto front.
* `fig5_preheat.png` — effect of preheat temperature on S_L and NOx.
* `fig6_flame_profile.png` — resolved temperature and species structure of one flame.

## Headline findings (T_in = 300 K)

* Pure NH₃ peaks at S_L ≈ 6.5 cm/s near φ = 1.1; adding 60 % H₂ raises the peak
  to ≈ 52 cm/s — an ~8× increase in burning velocity.
* NOx peaks slightly lean (φ ≈ 0.8–0.9) and falls sharply on the rich side as
  O-atom availability drops.
* The Pareto front shows that mildly rich operation (φ ≈ 1.2–1.3) gives the best
  compromise: high flame speed at substantially reduced NOx.

## Notes / possible extensions

* Swap `gri30.yaml` for a dedicated NH₃/H₂ mechanism (e.g. Okafor or
  Stagni) for higher-fidelity NOx and compare — a natural model-sensitivity study.
* Add a multicomponent-transport final run for the best operating points.
* Extend φ to 1.4–1.5 and add intermediate X_H₂ steps for smoother contours.
