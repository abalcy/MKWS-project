"""
================================================================================
Post-processing and visualisation for the NH3-H2-air flame study.
--------------------------------------------------------------------------------
Reads results/results.csv (produced by flame_nh3_h2.py) and generates:

    fig1_SL_contour.png      S_L(phi, X_H2) contour map at a chosen T_in
    fig2_NOx_contour.png     NOx(phi, X_H2) contour map at a chosen T_in
    fig3_SL_vs_XH2.png       S_L vs X_H2 for several phi (effect of H2 enrichment)
    fig4_pareto.png          NOx vs S_L trade-off (Pareto view) coloured by phi
    fig5_preheat.png         Effect of preheat temperature on S_L and NOx
    fig6_flame_profile.png   T and species profiles for one representative flame

All figures use a clean, report-ready style.
================================================================================
"""

from __future__ import annotations

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import cantera as ct

import flame_nh3_h2 as model  # reuse mechanism / solver settings


OUTDIR = "results"
FIGDIR = "figures"
os.makedirs(FIGDIR, exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 130,
    "font.size": 11,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "axes.axisbelow": True,
})


def load() -> pd.DataFrame:
    df = pd.read_csv(os.path.join(OUTDIR, "results.csv"))
    return df.dropna(subset=["S_L_cm_s", "NOx_ppm"])


# ----------------------------------------------------------------------------
# Figure 1 & 2 : contour maps on the (phi, X_H2) plane at fixed T_in
# ----------------------------------------------------------------------------
def contour(df: pd.DataFrame, t_in: float, z_col: str, label: str,
            fname: str, cmap: str) -> None:
    sub = df[df["t_in"] == t_in]
    if sub.empty:
        print(f"  (skip {fname}: no data at T_in={t_in})")
        return

    x = sub["phi"].to_numpy()
    y = sub["x_h2"].to_numpy()
    z = sub[z_col].to_numpy()

    tri = mtri.Triangulation(x, y)

    fig, ax = plt.subplots(figsize=(7, 5))
    cf = ax.tricontourf(tri, z, levels=14, cmap=cmap)
    cs = ax.tricontour(tri, z, levels=8, colors="k", linewidths=0.4, alpha=0.5)
    ax.clabel(cs, inline=True, fontsize=7, fmt="%.0f")
    ax.scatter(x, y, s=8, c="k", alpha=0.4)

    cbar = fig.colorbar(cf, ax=ax)
    cbar.set_label(label)
    ax.set_xlabel(r"Equivalence ratio $\varphi$ (-)")
    ax.set_ylabel(r"H$_2$ fraction in fuel $X_{H_2}$ (-)")
    ax.set_title(f"{label}   |   $T_{{in}}$ = {t_in:.0f} K, p = 1 atm")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, fname))
    plt.close(fig)
    print(f"  wrote {fname}")


# ----------------------------------------------------------------------------
# Figure 3 : S_L vs X_H2 for several phi (the headline "enrichment" effect)
# ----------------------------------------------------------------------------
def sl_vs_xh2(df: pd.DataFrame, t_in: float) -> None:
    sub = df[df["t_in"] == t_in]
    fig, ax = plt.subplots(figsize=(7, 5))
    for phi in sorted(sub["phi"].unique()):
        s = sub[sub["phi"] == phi].sort_values("x_h2")
        ax.plot(s["x_h2"], s["S_L_cm_s"], "o-", label=fr"$\varphi$={phi:.1f}")
    ax.set_xlabel(r"H$_2$ fraction in fuel $X_{H_2}$ (-)")
    ax.set_ylabel(r"Laminar flame speed $S_L$ (cm/s)")
    ax.set_title(fr"Effect of H$_2$ enrichment   |   $T_{{in}}$ = {t_in:.0f} K")
    ax.legend(ncol=2, fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig3_SL_vs_XH2.png"))
    plt.close(fig)
    print("  wrote fig3_SL_vs_XH2.png")


# ----------------------------------------------------------------------------
# Figure 4 : the trade-off / Pareto view  NOx vs S_L
# ----------------------------------------------------------------------------
def pareto(df: pd.DataFrame, t_in: float) -> None:
    sub = df[df["t_in"] == t_in]
    fig, ax = plt.subplots(figsize=(7, 5))
    sc = ax.scatter(sub["S_L_cm_s"], sub["NOx_ppm"], c=sub["phi"],
                    cmap="viridis", s=45, edgecolor="k", linewidth=0.3)
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label(r"Equivalence ratio $\varphi$ (-)")

    # Pareto front: minimise NOx, maximise S_L.
    pts = sub[["S_L_cm_s", "NOx_ppm"]].to_numpy()
    order = np.argsort(-pts[:, 0])          # fastest first
    front_idx, best_nox = [], np.inf
    for i in order:
        if pts[i, 1] < best_nox:
            best_nox = pts[i, 1]
            front_idx.append(i)
    front = pts[front_idx]
    front = front[np.argsort(front[:, 0])]
    ax.plot(front[:, 0], front[:, 1], "r--", lw=1.5,
            label="Pareto front (max $S_L$, min NOx)")

    ax.set_xlabel(r"Laminar flame speed $S_L$ (cm/s)")
    ax.set_ylabel("Exhaust NOx (ppm)")
    ax.set_title(fr"Flame-speed / NOx trade-off   |   $T_{{in}}$ = {t_in:.0f} K")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig4_pareto.png"))
    plt.close(fig)
    print("  wrote fig4_pareto.png")


# ----------------------------------------------------------------------------
# Figure 5 : effect of preheat temperature (at stoichiometric, fixed X_H2)
# ----------------------------------------------------------------------------
def preheat(df: pd.DataFrame, x_h2: float = 0.4, phi: float = 1.0) -> None:
    sub = df[(np.isclose(df["x_h2"], x_h2)) & (np.isclose(df["phi"], phi))]
    sub = sub.sort_values("t_in")
    if sub.empty:
        print("  (skip fig5: no matching data)")
        return
    fig, ax1 = plt.subplots(figsize=(7, 5))
    ax2 = ax1.twinx()
    l1, = ax1.plot(sub["t_in"], sub["S_L_cm_s"], "o-", color="tab:blue",
                   label=r"$S_L$")
    l2, = ax2.plot(sub["t_in"], sub["NOx_ppm"], "s--", color="tab:red",
                   label="NOx")
    ax1.set_xlabel(r"Preheat temperature $T_{in}$ (K)")
    ax1.set_ylabel(r"$S_L$ (cm/s)", color="tab:blue")
    ax2.set_ylabel("NOx (ppm)", color="tab:red")
    ax1.tick_params(axis="y", labelcolor="tab:blue")
    ax2.tick_params(axis="y", labelcolor="tab:red")
    ax1.set_title(fr"Preheat effect   |   $X_{{H_2}}$={x_h2:.1f}, $\varphi$={phi:.1f}")
    ax1.legend(handles=[l1, l2], fontsize=9, loc="upper left")
    ax2.grid(False)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig5_preheat.png"))
    plt.close(fig)
    print("  wrote fig5_preheat.png")


# ----------------------------------------------------------------------------
# Figure 6 : a single resolved flame structure (T and key species profiles)
# ----------------------------------------------------------------------------
def flame_profile(x_h2: float = 0.4, phi: float = 1.0,
                  t_in: float = 300.0) -> None:
    gas = ct.Solution(model.MECHANISM)
    flame, _ = model.solve_flame(gas, phi, x_h2, t_in)
    z_mm = flame.grid * 1000.0

    fig, ax1 = plt.subplots(figsize=(7, 5))
    ax2 = ax1.twinx()
    ax1.plot(z_mm, flame.T, color="k", lw=2, label="Temperature")
    for sp, col in [("NH3", "tab:green"), ("H2", "tab:blue"),
                    ("NO", "tab:red"), ("OH", "tab:orange")]:
        idx = gas.species_index(sp)
        ax2.plot(z_mm, flame.X[idx], "--", color=col, label=sp)

    ax1.set_xlabel("Distance through flame (mm)")
    ax1.set_ylabel("Temperature (K)")
    ax2.set_ylabel("Mole fraction (-)")
    ax1.set_title(fr"Flame structure   |   $X_{{H_2}}$={x_h2:.1f}, "
                  fr"$\varphi$={phi:.1f}, $T_{{in}}$={t_in:.0f} K")
    lines = ax1.get_lines() + ax2.get_lines()
    ax1.legend(lines, [l.get_label() for l in lines], fontsize=9, loc="center right")
    ax2.grid(False)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig6_flame_profile.png"))
    plt.close(fig)
    print("  wrote fig6_flame_profile.png")


# ----------------------------------------------------------------------------
def main() -> None:
    df = load()
    t_ref = 300.0
    print("Generating figures:")
    contour(df, t_ref, "S_L_cm_s", r"$S_L$ (cm/s)",
            "fig1_SL_contour.png", "plasma")
    contour(df, t_ref, "NOx_ppm", "NOx (ppm)",
            "fig2_NOx_contour.png", "inferno")
    sl_vs_xh2(df, t_ref)
    pareto(df, t_ref)
    preheat(df)
    flame_profile()
    print("All figures saved to ./figures/")


if __name__ == "__main__":
    main()
