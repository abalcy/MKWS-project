"""Run the sweep for a single T_in value and save to a per-temperature CSV.

Usage: python3 run_chunk.py <T_in>
"""
import sys, os, time, itertools
import numpy as np, pandas as pd, cantera as ct
import flame_nh3_h2 as m

t_in = float(sys.argv[1])
gas = ct.Solution(m.MECHANISM)
out = f"results/chunk_{int(t_in)}.csv"
cols = ["x_h2", "phi", "t_in", "S_L_cm_s", "T_ad_K",
        "NO_ppm", "NO2_ppm", "N2O_ppm", "NOx_ppm", "grid_points"]

# Resume support: skip (x_h2, phi) pairs already present in the CSV.
done = set()
if os.path.exists(out):
    prev = pd.read_csv(out)
    done = {(round(r.x_h2, 3), round(r.phi, 3)) for r in prev.itertuples()}
else:
    pd.DataFrame(columns=cols).to_csv(out, index=False)

t0 = time.time()
for x_h2 in m.X_H2_GRID:
    guess = None
    for phi in m.PHI_GRID:
        if (round(float(x_h2), 3), round(float(phi), 3)) in done:
            continue
        try:
            guess, rec = m.solve_flame(gas, phi, x_h2, t_in, guess)
            print(f"T={t_in:.0f} X_H2={x_h2:.2f} phi={phi:.2f} "
                  f"S_L={rec['S_L_cm_s']:6.2f} NOx={rec['NOx_ppm']:8.1f}", flush=True)
        except Exception as e:
            print(f"T={t_in:.0f} X_H2={x_h2:.2f} phi={phi:.2f} FAILED {type(e).__name__}", flush=True)
            rec = dict(x_h2=x_h2, phi=phi, t_in=t_in, S_L_cm_s=np.nan,
                       T_ad_K=np.nan, NO_ppm=np.nan, NO2_ppm=np.nan,
                       N2O_ppm=np.nan, NOx_ppm=np.nan, grid_points=np.nan)
            guess = None
        # append this single row immediately
        pd.DataFrame([rec])[cols].to_csv(out, mode="a", header=False, index=False)
print(f"\nDone T_in={t_in:.0f} in {time.time()-t0:.0f}s -> {out}", flush=True)
