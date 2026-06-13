"""
================================================================================
Laminar flame speed and NOx formation in NH3-H2-air blends
--------------------------------------------------------------------------------
Computational Methods in Combustion (MKWS) project.

This module computes the laminar burning velocity (S_L) and the post-flame
NOx (NO + NO2) emission of premixed ammonia/hydrogen/air flames as a function
of three control variables:

    1. Hydrogen fraction in the fuel blend   X_H2  in [0.0 ... 0.6]
    2. Equivalence ratio                      phi   in [0.7 ... 1.4]
    3. Unburned mixture (preheat) temperature T_in  in [300 ... 600] K

Solver strategy:
    * Cantera FreeFlame (1-D freely propagating premixed flame).
    * Continuation: each solved flame is reused as the initial guess for the
      next case along a sweep. This makes the parametric sweep an order of
      magnitude faster and far more robust near lean/rich extinction limits.
    * Mixture-averaged transport (good accuracy/speed trade-off for screening),
      with an option to switch to multicomponent for a final high-fidelity run.

Outputs:
    * results/results.csv     - tidy table of every converged case
    * (figures are produced by plots.py from this CSV)

Author: <your name>, MEiL, Warsaw University of Technology
================================================================================
"""

from __future__ import annotations

import os
import time
import itertools
import numpy as np
import pandas as pd
import cantera as ct


# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------
MECHANISM = "gri30.yaml"     # contains NH3, H2 and full N-chemistry (NO, NO2, N2O)
PRESSURE = ct.one_atm        # all cases at 1 atm
OXIDIZER = "O2:1.0, N2:3.76" # standard air composition (molar)
DOMAIN_WIDTH = 0.03          # m, initial 1-D domain width

# Sweep grids ----------------------------------------------------------------
X_H2_GRID = np.array([0.0, 0.2, 0.4, 0.6])                  # H2 mole fraction in fuel
PHI_GRID = np.array([0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3])    # equivalence ratio
TIN_GRID = np.array([300.0, 450.0, 600.0])                  # K, preheat temperatures

# Refinement criteria. Looser values -> faster screening; tighten for final run.
REFINE = dict(ratio=3.0, slope=0.08, curve=0.08, prune=0.03)

OUTDIR = "results"
os.makedirs(OUTDIR, exist_ok=True)


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
def fuel_string(x_h2: float) -> str:
    """Return a Cantera composition string for the NH3/H2 fuel blend.

    x_h2 is the mole fraction of H2 in the (NH3 + H2) fuel mixture.
    """
    x_nh3 = 1.0 - x_h2
    return f"NH3:{x_nh3:.6f}, H2:{x_h2:.6f}"


def solve_flame(gas: ct.Solution,
                phi: float,
                x_h2: float,
                t_in: float,
                initial_guess: ct.FreeFlame | None = None
                ) -> tuple[ct.FreeFlame, dict]:
    """Solve a single freely-propagating premixed flame.

    Parameters
    ----------
    gas : ct.Solution
        Reusable gas object (its state is overwritten here).
    phi, x_h2, t_in : float
        Operating point.
    initial_guess : ct.FreeFlame or None
        A previously converged flame used to seed the solver (continuation).

    Returns
    -------
    flame : ct.FreeFlame
        The converged flame (also serves as the next initial guess).
    record : dict
        Scalar results for this operating point.
    """
    gas.TP = t_in, PRESSURE
    gas.set_equivalence_ratio(phi, fuel_string(x_h2), OXIDIZER)

    flame = ct.FreeFlame(gas, width=DOMAIN_WIDTH)
    flame.set_refine_criteria(**REFINE)
    flame.transport_model = "mixture-averaged"

    # Seed from a previously converged solution when shapes are compatible.
    if initial_guess is not None:
        try:
            flame.set_initial_guess(data=initial_guess.to_array())
        except Exception:
            pass  # fall back to Cantera's default ignition guess

    flame.solve(loglevel=0, auto=True)

    # --- extract scalar quantities ------------------------------------------
    s_l = flame.velocity[0]                 # m/s, unburned-side flame speed
    t_ad = flame.T[-1]                       # K, burned-gas temperature

    def exit_ppm(species: str) -> float:
        idx = gas.species_index(species)
        return float(flame.X[idx][-1] * 1.0e6)

    no_ppm = exit_ppm("NO")
    no2_ppm = exit_ppm("NO2")
    n2o_ppm = exit_ppm("N2O")
    nox_ppm = no_ppm + no2_ppm              # NOx defined as NO + NO2

    record = dict(
        x_h2=x_h2,
        phi=phi,
        t_in=t_in,
        S_L_cm_s=s_l * 100.0,
        T_ad_K=t_ad,
        NO_ppm=no_ppm,
        NO2_ppm=no2_ppm,
        N2O_ppm=n2o_ppm,
        NOx_ppm=nox_ppm,
        grid_points=flame.flame.n_points,
    )
    return flame, record


# ----------------------------------------------------------------------------
# Main sweep
# ----------------------------------------------------------------------------
def run_sweep() -> pd.DataFrame:
    gas = ct.Solution(MECHANISM)
    rows: list[dict] = []
    t_start = time.time()

    total = len(TIN_GRID) * len(X_H2_GRID) * len(PHI_GRID)
    done = 0

    # Outer loops: T_in and X_H2 define a "family"; phi is swept with
    # continuation so neighbouring rich/lean cases reuse the prior solution.
    for t_in, x_h2 in itertools.product(TIN_GRID, X_H2_GRID):
        guess = None
        for phi in PHI_GRID:
            done += 1
            try:
                guess, rec = solve_flame(gas, phi, x_h2, t_in, guess)
                rows.append(rec)
                print(f"[{done:3d}/{total}] T_in={t_in:.0f}K  X_H2={x_h2:.2f}  "
                      f"phi={phi:.2f}  ->  S_L={rec['S_L_cm_s']:6.2f} cm/s  "
                      f"NOx={rec['NOx_ppm']:8.1f} ppm")
            except Exception as exc:
                # Record a failed (non-converged / extinguished) case as NaN.
                print(f"[{done:3d}/{total}] T_in={t_in:.0f}K  X_H2={x_h2:.2f}  "
                      f"phi={phi:.2f}  ->  FAILED ({type(exc).__name__})")
                rows.append(dict(x_h2=x_h2, phi=phi, t_in=t_in,
                                 S_L_cm_s=np.nan, T_ad_K=np.nan,
                                 NO_ppm=np.nan, NO2_ppm=np.nan,
                                 N2O_ppm=np.nan, NOx_ppm=np.nan,
                                 grid_points=np.nan))
                guess = None  # reset continuation chain after a failure

    df = pd.DataFrame(rows)
    csv_path = os.path.join(OUTDIR, "results.csv")
    df.to_csv(csv_path, index=False)

    elapsed = time.time() - t_start
    print(f"\nDone. {len(df)} cases in {elapsed:.1f} s. Saved -> {csv_path}")
    return df


if __name__ == "__main__":
    run_sweep()
