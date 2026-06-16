"""
run_demo.py — Quantumograph v14 desktop demo
© 2025–2026 Sergej Materov <sergejmaterov2@gmail.com>
=============================================
Full pipeline: Z⁴ torus → spectral dimension → plaquette Hamiltonian → visualisation.
Supports gauge groups U(1), SU(2), SU(3).

Usage
-----
    python run_demo.py                         # U(1), L=6, sparse spectrum
    python run_demo.py --gauge SU2 --L 4       # SU(2), L=4
    python run_demo.py --gauge SU3 --L 4       # SU(3), L=4
    python run_demo.py --gauge all --L 4       # compare all three groups
    python run_demo.py --L 8 --mode exact      # larger torus, exact spectrum
    python run_demo.py --mode walk --n_walks 100000  # Monte Carlo walks
    python run_demo.py --fss                   # finite-size scaling L=3..7
    python run_demo.py --animate               # animated t-slice scan
    python run_demo.py --gauge SU2 --wilson    # Wilson loop vs R,T

Author  : © 2025–2026 Sergej Materov <sergejmaterov2@gmail.com>
License : CC BY-NC 4.0
"""

import argparse
import time
import numpy as np

from quantumograph import (Z4Torus, SpectralAnalysis, PlaquetteAnalysis,
                            plot_full_analysis, plot_fss,
                            fss_extrapolate, check_4d_criterion,
                            wilson_loop_su2, polyakov_loop_su3)


# ──────────────────────────────────────────────────────────────────────────────
# Single-gauge run
# ──────────────────────────────────────────────────────────────────────────────

def run_single(L: int, gauge: str, mode: str, k: int,
               n_walks: int, t_slice: int, animate: bool,
               show_wilson: bool) -> None:

    print(f"\n{'='*62}")
    print(f"  Quantumograph v14 — Z⁴ Torus Demo")
    print(f"  L={L}  N={L**4}  gauge={gauge}  mode={mode}")
    print(f"{'='*62}\n")

    # 1. Build torus
    t0    = time.time()
    torus = Z4Torus(L=L, gauge_group=gauge, seed=42)
    print(f"[1/3] Torus built in {time.time()-t0:.2f}s  "
          f"— {torus.N} nodes, {torus.graph.number_of_edges()} edges")

    # 2. Spectral dimension (always on the graph, independent of gauge group)
    t0 = time.time()
    sa = SpectralAnalysis(torus, mode=mode, k=k, n_walks=n_walks).run()
    print(f"[2/3] Spectral analysis done in {time.time()-t0:.2f}s")
    print(sa.summary())

    # 3. Plaquette / gauge analysis
    t0 = time.time()
    pa = PlaquetteAnalysis(torus, J=1.0).run()
    print(f"[3/3] Plaquette analysis done in {time.time()-t0:.2f}s")
    print(pa.summary())

    # 4. Wilson loop table for SU(2)
    if show_wilson and gauge == 'SU2':
        _print_wilson_table(torus)

    # 5. Polyakov loop vs temperature for SU(3)
    if gauge == 'SU3':
        _print_polyakov(torus)

    # 6. Visualise
    import matplotlib.pyplot as plt
    fig = plot_full_analysis(torus, sa, pa, t_slice=t_slice, show=False)

    if animate:
        from quantumograph import animate_torus
        ani = animate_torus(torus, pa, field='energy', show=False)

    plt.tight_layout()
    plt.show()


# ──────────────────────────────────────────────────────────────────────────────
# Compare all three gauge groups side by side
# ──────────────────────────────────────────────────────────────────────────────

def run_all_gauges(L: int, mode: str, k: int) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    print(f"\n{'='*62}")
    print(f"  Quantumograph v14 — Gauge Group Comparison")
    print(f"  L={L}  N={L**4}  mode={mode}")
    print(f"{'='*62}\n")

    gauges  = ['U1', 'SU2', 'SU3']
    results = {}

    for gauge in gauges:
        print(f"\n── {gauge} ──")
        t0    = time.time()
        torus = Z4Torus(L=L, gauge_group=gauge, seed=42)
        sa    = SpectralAnalysis(torus, mode=mode, k=k).run()
        pa    = PlaquetteAnalysis(torus, J=1.0).run()
        print(f"   Built+analysed in {time.time()-t0:.2f}s")
        print(f"   dₛ(∞) = {sa.ds_inf:.4f}  {'✓' if sa.passes_criterion else '✗'}")
        print(f"   H/site = {pa.H_total/torus.N:.6f}")
        results[gauge] = (torus, sa, pa)

    # ── Comparison figure ──────────────────────────────────────────────────
    fig = plt.figure(figsize=(15, 10), facecolor='#080818')
    gs  = GridSpec(2, 3, figure=fig, hspace=0.38, wspace=0.32)

    colours = {'U1': '#00d4ff', 'SU2': '#ff6600', 'SU3': '#00ff88'}

    # Row 0: dₛ(σ) curves
    ax_ds = fig.add_subplot(gs[0, :])
    ax_ds.set_facecolor('#0a0a1a')
    ax_ds.tick_params(colors='white')
    ax_ds.spines[:].set_color('#222244')

    for gauge, (torus, sa, pa) in results.items():
        ax_ds.semilogx(sa.sigma_mid, sa.ds_arr,
                       color=colours[gauge], lw=1.6,
                       label=fr'{gauge}  $d_s(\infty)={sa.ds_inf:.3f}$')

    ax_ds.axhline(4.0, color='white', lw=0.8, ls='--', alpha=0.5,
                  label=r'$d_s=4$ target')
    ax_ds.axhspan(3.9, 4.1, alpha=0.08, color='white')
    ax_ds.set_xlabel(r'$\sigma$', color='white')
    ax_ds.set_ylabel(r'$d_s(\sigma)$', color='white')
    ax_ds.set_title('Spectral dimension — gauge group comparison',
                    color='white', fontsize=11)
    ax_ds.legend(facecolor='#111133', edgecolor='#333355', labelcolor='white')
    ax_ds.set_ylim(bottom=0)

    # Row 1: energy density histograms per gauge
    for col, (gauge, (torus, sa, pa)) in enumerate(results.items()):
        ax = fig.add_subplot(gs[1, col])
        ax.set_facecolor('#0a0a1a')
        ax.tick_params(colors='white')
        ax.spines[:].set_color('#222244')

        eps = pa.energy_density
        ax.hist(eps, bins=40, color=colours[gauge], alpha=0.75, edgecolor='none')
        ax.axvline(eps.mean(), color='white', lw=1, ls='--',
                   label=f'mean={eps.mean():.3f}')
        ax.set_title(fr'{gauge}  energy density $\varepsilon(x)$',
                     color='white', fontsize=9)
        ax.set_xlabel(r'$\varepsilon$', color='white')
        ax.set_ylabel('count', color='white')
        ax.legend(facecolor='#111133', labelcolor='white', fontsize=8)

    fig.suptitle(
        f'Quantumograph v14 — Z⁴ Torus  L={L}  N={L**4}  '
        f'Gauge group comparison',
        color='white', fontsize=12)

    plt.tight_layout()
    plt.show()


# ──────────────────────────────────────────────────────────────────────────────
# SU(2) Wilson loop table
# ──────────────────────────────────────────────────────────────────────────────

def _print_wilson_table(torus) -> None:
    print(f"\n  SU(2) Wilson loops <Re Tr W(R,T)>/2")
    print(f"  {'R\\T':>4}", end='')
    T_vals = [1, 2, 3, 4]
    R_vals = [1, 2, 3, 4]
    for T in T_vals:
        print(f"  T={T:2d}   ", end='')
    print()
    for R in R_vals:
        print(f"  R={R:2d}  ", end='')
        for T in T_vals:
            W = wilson_loop_su2(torus, mu=0, nu=3, R=R, T=T)
            print(f"  {W:+.4f}", end='')
        print()
    print()

    # Check area vs perimeter law
    # W(R,T) ~ exp(-σ·R·T)  area law (confinement)
    # W(R,T) ~ exp(-μ·(2R+2T))  perimeter law (deconfinement)
    w11 = wilson_loop_su2(torus, mu=0, nu=3, R=1, T=1)
    w22 = wilson_loop_su2(torus, mu=0, nu=3, R=2, T=2)
    if w11 > 1e-6 and w22 > 1e-6:
        ratio = np.log(abs(w22)) / np.log(abs(w11))
        print(f"  ln W(2,2)/ln W(1,1) = {ratio:.3f}")
        print(f"  Area law predicts 4.0, perimeter law predicts 2.0")


# ──────────────────────────────────────────────────────────────────────────────
# SU(3) Polyakov loop
# ──────────────────────────────────────────────────────────────────────────────

def _print_polyakov(torus) -> None:
    P = polyakov_loop_su3(torus, mu=3)
    print(f"\n  SU(3) Polyakov loop (temporal direction μ=3):")
    print(f"  Re<P> = {P.real:.6f}")
    if abs(P.real) < 0.1:
        print(f"  → Confined phase (Re<P> ≈ 0)")
    else:
        print(f"  → Deconfined phase (Re<P> ≠ 0)")


# ──────────────────────────────────────────────────────────────────────────────
# FSS
# ──────────────────────────────────────────────────────────────────────────────

def run_fss(L_values: list, gauge: str, mode: str, k: int) -> None:
    print(f"\n{'='*62}")
    print(f"  Quantumograph v14 — Finite-Size Scaling  gauge={gauge}")
    print(f"  L values: {L_values}")
    print(f"{'='*62}\n")

    ds_results = []
    for L in L_values:
        print(f"\n── L={L}  N={L**4} ──")
        torus = Z4Torus(L=L, gauge_group=gauge, seed=42)
        sa    = SpectralAnalysis(torus, mode=mode,
                                  k=min(k, torus.N - 2)).run()
        print(f"   dₛ = {sa.ds_inf:.4f}  {'✓' if sa.passes_criterion else '✗'}")
        ds_results.append(sa.ds_inf)

    L_arr  = np.array(L_values, dtype=float)
    ds_arr = np.array(ds_results)
    ds_inf, coeffs = fss_extrapolate(L_arr, ds_arr, order=2)

    print(f"\n{'─'*40}")
    print(f"FSS extrapolation  dₛ(∞) = {ds_inf:.4f}")
    print(f"4D criterion: {'✓ PASSES' if check_4d_criterion(ds_inf) else '✗ FAILS'}")

    import matplotlib.pyplot as plt
    plot_fss(L_arr, ds_arr, ds_inf, coeffs, show=True)


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Quantumograph v14 — Z⁴ torus analysis (U1/SU2/SU3)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_demo.py                        # U(1), L=6, sparse
  python run_demo.py --gauge SU2 --L 4     # SU(2), L=4
  python run_demo.py --gauge SU3 --L 4     # SU(3), L=4
  python run_demo.py --gauge all --L 4     # compare all three
  python run_demo.py --gauge SU2 --wilson  # Wilson loop table
  python run_demo.py --fss --gauge U1      # FSS extrapolation
  python run_demo.py --animate             # animated t-slice
        """)

    parser.add_argument('--L',       type=int,   default=6,
                        help='Torus side L (N=L^4, default 6)')
    parser.add_argument('--gauge',   type=str,   default='U1',
                        choices=['U1', 'SU2', 'SU3', 'all'],
                        help='Gauge group (default U1)')
    parser.add_argument('--mode',    type=str,   default='exact',
                        choices=['exact', 'sparse', 'walk'],
                        help='Spectrum mode (default exact)')
    parser.add_argument('--k',       type=int,   default=256,
                        help='Eigenvalues to compute in sparse mode')
    parser.add_argument('--n_walks', type=int,   default=50_000,
                        help='MC walks (walk mode)')
    parser.add_argument('--t_slice', type=int,   default=0,
                        help='Time slice index for visualisation')
    parser.add_argument('--animate', action='store_true',
                        help='Show animated t-slice scan')
    parser.add_argument('--wilson',  action='store_true',
                        help='Print Wilson loop table (SU2 only)')
    parser.add_argument('--fss',     action='store_true',
                        help='Run finite-size scaling')
    parser.add_argument('--fss_L',   type=int,   nargs='+',
                        default=[3, 4, 5, 6, 7],
                        help='L values for FSS (default: 3 4 5 6 7)')
    args = parser.parse_args()

    if args.gauge == 'all':
        run_all_gauges(args.L, args.mode, args.k)
    elif args.fss:
        run_fss(args.fss_L, args.gauge, args.mode, args.k)
    else:
        run_single(args.L, args.gauge, args.mode, args.k,
                   args.n_walks, args.t_slice, args.animate, args.wilson)


if __name__ == '__main__':
    main()
