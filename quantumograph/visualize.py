"""
quantumograph.visualize
=======================
Physically correct visualisation of the Z⁴ torus and its observables.

All layout is derived from the actual graph/spectral data — no spring_layout,
no imported continuum metrics. What you see is what the theory predicts.

Panels
------
1. Spectral dimension curve  dₛ(σ) with 4D criterion line
2. Plaquette energy density  ε(x) projected onto 2D time-slice
3. Topological charge density q(x) on the same slice
4. Animated time-slice of the torus graph coloured by ε(x)

Author  : © 2025–2026 Sergej Materov <sergejmaterov2@gmail.com>
License : CC BY-NC 4.0
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.widgets as widgets
from matplotlib.animation import FuncAnimation
from matplotlib.collections import LineCollection
from matplotlib.gridspec import GridSpec
from typing import Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from .torus import Z4Torus
    from .spectral import SpectralAnalysis
    from .plaquette import PlaquetteAnalysis


# ──────────────────────────────────────────────────────────────────────────────
# Colour helpers
# ──────────────────────────────────────────────────────────────────────────────

def _normalise(arr: np.ndarray) -> np.ndarray:
    lo, hi = arr.min(), arr.max()
    if hi == lo:
        return np.zeros_like(arr)
    return (arr - lo) / (hi - lo)


# ──────────────────────────────────────────────────────────────────────────────
# 1. Spectral dimension plot
# ──────────────────────────────────────────────────────────────────────────────

def plot_spectral_dimension(sa: 'SpectralAnalysis',
                            ax: Optional[plt.Axes] = None,
                            show: bool = True) -> plt.Figure:
    """
    Plot dₛ(σ) curve with the 4D criterion band |dₛ - 4| < 0.1.

    Parameters
    ----------
    sa   : SpectralAnalysis (must have been .run())
    ax   : existing Axes or None (creates new figure)
    show : call plt.show()
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(9, 5), facecolor='#0a0a1a')
    else:
        fig = ax.figure

    ax.set_facecolor('#0a0a1a')

    # dₛ(σ) curve
    ax.semilogx(sa.sigma_mid, sa.ds_arr,
                color='#00d4ff', lw=1.8, label=r'$d_s(\sigma)$')

    # 4D criterion band
    ax.axhline(4.0, color='#ffaa00', lw=1.2, ls='--',
               label=r'$d_s = 4$ (target)')
    ax.axhspan(3.9, 4.1, alpha=0.12, color='#ffaa00',
               label=r'$|d_s - 4| < 0.1$ criterion')

    # dₛ(∞) estimate
    if sa.ds_inf is not None:
        color_inf = '#00ff88' if sa.passes_criterion else '#ff4444'
        ax.axhline(sa.ds_inf, color=color_inf, lw=1.0, ls=':',
                   label=fr'$d_s(\infty) = {sa.ds_inf:.3f}$')

    ax.set_xlabel(r'Diffusion time $\sigma$', color='white')
    ax.set_ylabel(r'Spectral dimension $d_s(\sigma)$', color='white')
    ax.set_title(
        fr'Quantumograph v14 — Spectral Dimension  '
        fr'$Z_{{{sa.torus.L}}}^4$ torus  '
        fr'($N = {sa.torus.N}$)',
        color='white', fontsize=11)
    ax.tick_params(colors='white')
    ax.spines[:].set_color('#333355')
    leg = ax.legend(facecolor='#111133', edgecolor='#333355', labelcolor='white')
    ax.set_ylim(bottom=0)

    status = '✓ PASSES' if sa.passes_criterion else '✗ FAILS'
    color  = '#00ff88'  if sa.passes_criterion else '#ff4444'
    ax.text(0.98, 0.05, f'4D criterion: {status}',
            transform=ax.transAxes, ha='right', va='bottom',
            color=color, fontsize=10, fontweight='bold')

    if show:
        plt.tight_layout()
        plt.show()
    return fig


# ──────────────────────────────────────────────────────────────────────────────
# 2. Torus slice projection
# ──────────────────────────────────────────────────────────────────────────────

def _project_slice(torus: 'Z4Torus',
                   t_slice: int,
                   scalar_field: Optional[np.ndarray] = None
                   ):
    """
    Project nodes at fixed t-coordinate onto (x₀, x₁) plane.

    Returns
    -------
    node_ids : (M,) array — node IDs in this slice
    xs, ys   : (M,) arrays — 2D coordinates (x₀, x₁)
    values   : (M,) normalised scalar field values (or None)
    edges    : list of ((x0,y0),(x1,y1)) for edges within slice
    """
    L      = torus.L
    coords = torus.coords           # (N, 4)
    G      = torus.graph

    # Nodes in this time slice (coordinate 3 = t)
    mask     = coords[:, 3] == (t_slice % L)
    node_ids = np.where(mask)[0]
    xs       = coords[node_ids, 0].astype(float)
    ys       = coords[node_ids, 1].astype(float)

    values = None
    if scalar_field is not None:
        vals   = scalar_field[node_ids]
        values = _normalise(vals)

    node_set = set(node_ids.tolist())
    edges = []
    for u, v in G.edges():
        if u in node_set and v in node_set:
            cu = torus.coord(u); cv = torus.coord(v)
            edges.append(((cu[0], cu[1]), (cv[0], cv[1])))

    return node_ids, xs, ys, values, edges


# ──────────────────────────────────────────────────────────────────────────────
# 3. Static four-panel figure
# ──────────────────────────────────────────────────────────────────────────────

def plot_full_analysis(torus: 'Z4Torus',
                       sa:    'SpectralAnalysis',
                       pa:    'PlaquetteAnalysis',
                       t_slice: int = 0,
                       show: bool = True) -> plt.Figure:
    """
    Four-panel figure — works for U(1), SU(2), SU(3):

      [0,0] dₛ(σ) curve
      [0,1] Plaquette energy density ε(x) on t-slice
      [1,0] Group-specific observable:
              U(1)  → topological charge q(x)
              SU(2) → Wilson loop W(R,T) heatmap (mean over nodes)
              SU(3) → plaquette Re Tr[P]/3 distribution
      [1,1] Torus graph coloured by ε(x)
    """
    gauge = torus.gauge_group
    fig   = plt.figure(figsize=(14, 10), facecolor='#080818')
    gs    = GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.3)

    ax_ds  = fig.add_subplot(gs[0, 0])
    ax_eps = fig.add_subplot(gs[0, 1])
    ax_obs = fig.add_subplot(gs[1, 0])
    ax_g   = fig.add_subplot(gs[1, 1])

    for ax in (ax_ds, ax_eps, ax_obs, ax_g):
        ax.set_facecolor('#0a0a1a')
        ax.tick_params(colors='white')
        ax.spines[:].set_color('#222244')

    # ── Panel 0: dₛ(σ) ────────────────────────────────────────────────────
    plot_spectral_dimension(sa, ax=ax_ds, show=False)

    # ── Panel 1: energy density on t-slice ────────────────────────────────
    node_ids, xs, ys, _, _ = _project_slice(torus, t_slice, pa.energy_density)
    if len(xs) > 0:
        jitter = np.random.default_rng(0).uniform(-0.1, 0.1, (len(xs), 2))
        sc = ax_eps.scatter(xs + jitter[:,0], ys + jitter[:,1],
                            c=pa.energy_density[node_ids],
                            cmap='inferno', s=60, alpha=0.85, linewidths=0)
        cb = plt.colorbar(sc, ax=ax_eps, pad=0.02)
        cb.ax.yaxis.set_tick_params(color='white')
        cb.set_label(r'$\varepsilon$', color='white')
        plt.setp(cb.ax.yaxis.get_ticklabels(), color='white')
    ax_eps.set_title(fr'Energy density $\varepsilon(x)$, $t={t_slice}$  [{gauge}]',
                     color='white', fontsize=9)
    ax_eps.set_xlabel(r'$x_0$', color='white')
    ax_eps.set_ylabel(r'$x_1$', color='white')

    # ── Panel 2: group-specific observable ────────────────────────────────
    if gauge == 'U1' and pa.q_density is not None:
        # Topological charge density
        node_ids2, xs2, ys2, _, _ = _project_slice(torus, t_slice, pa.q_density)
        if len(xs2) > 0:
            jitter2 = np.random.default_rng(1).uniform(-0.1, 0.1, (len(xs2), 2))
            sc2 = ax_obs.scatter(xs2 + jitter2[:,0], ys2 + jitter2[:,1],
                                 c=pa.q_density[node_ids2],
                                 cmap='coolwarm', s=60, alpha=0.85, linewidths=0)
            cb2 = plt.colorbar(sc2, ax=ax_obs, pad=0.02)
            cb2.ax.yaxis.set_tick_params(color='white')
            cb2.set_label(r'$q(x)$', color='white')
            plt.setp(cb2.ax.yaxis.get_ticklabels(), color='white')
        ax_obs.set_title(fr'Topological charge $q(x)$, $t={t_slice}$',
                         color='white', fontsize=9)
        ax_obs.set_xlabel(r'$x_0$', color='white')
        ax_obs.set_ylabel(r'$x_1$', color='white')

    elif gauge == 'SU2' and pa.wilson_loops is not None:
        # Wilson loop W(R,T) table as colour grid
        wl   = pa.wilson_loops
        R_max = max(r for r, t in wl)
        T_max = max(t for r, t in wl)
        grid  = np.full((R_max, T_max), np.nan)
        for (R, T), val in wl.items():
            grid[R-1, T-1] = val
        im = ax_obs.imshow(grid, cmap='RdYlGn', vmin=-1, vmax=1,
                           origin='lower', aspect='auto')
        cb3 = plt.colorbar(im, ax=ax_obs, pad=0.02)
        cb3.set_label(r'$\langle \mathrm{Re\,Tr}\,W(R,T)\rangle/2$', color='white')
        cb3.ax.yaxis.set_tick_params(color='white')
        plt.setp(cb3.ax.yaxis.get_ticklabels(), color='white')
        ax_obs.set_xticks(range(T_max)); ax_obs.set_xticklabels([f'T={i+1}' for i in range(T_max)])
        ax_obs.set_yticks(range(R_max)); ax_obs.set_yticklabels([f'R={i+1}' for i in range(R_max)])
        ax_obs.set_title(r'Wilson loops $W(R,T)$ [SU(2)]', color='white', fontsize=9)

        # Annotate cells
        for (R, T), val in wl.items():
            ax_obs.text(T-1, R-1, f'{val:.3f}', ha='center', va='center',
                        color='white', fontsize=8, fontweight='bold')

    elif gauge == 'SU3' and pa.plaq_values is not None:
        # Plaquette value distribution (all 6 orientations)
        plaq_flat = pa.plaq_values.flatten()
        ax_obs.hist(plaq_flat, bins=60, color='#00ff88', alpha=0.75, edgecolor='none')
        ax_obs.axvline(plaq_flat.mean(), color='white', lw=1.2, ls='--',
                       label=fr'mean={plaq_flat.mean():.4f}')
        if pa.polyakov is not None:
            ax_obs.set_title(
                fr'SU(3) plaquette $\mathrm{{Re\,Tr}}[P]/3$   '
                fr'Polyakov: {pa.polyakov.real:.4f}',
                color='white', fontsize=9)
        else:
            ax_obs.set_title(r'SU(3) plaquette $\mathrm{Re\,Tr}[P]/3$ distribution',
                             color='white', fontsize=9)
        ax_obs.set_xlabel(r'$\mathrm{Re\,Tr}[P_{x,\mu\nu}]/3$', color='white')
        ax_obs.set_ylabel('count', color='white')
        ax_obs.legend(facecolor='#111133', labelcolor='white', fontsize=8)
    else:
        ax_obs.text(0.5, 0.5, 'No group-specific\nobservable available',
                    color='white', ha='center', va='center',
                    transform=ax_obs.transAxes, fontsize=10)

    # ── Panel 3: graph coloured by energy density ──────────────────────────
    # ── Panel 3: torus graph coloured by ε(x) ─────────────────────────────
    node_ids_g, xs_g, ys_g, _, edges_g = _project_slice(
        torus, t_slice, pa.energy_density)

    if len(edges_g) > 0:
        segs = [[(x0, y0), (x1, y1)] for (x0,y0),(x1,y1) in edges_g]
        lc   = LineCollection(segs, colors='#1a3a5a', linewidths=0.6,
                              alpha=0.5, zorder=1)
        ax_g.add_collection(lc)

    if len(xs_g) > 0:
        sc_g = ax_g.scatter(xs_g, ys_g, c=pa.energy_density[node_ids_g],
                            cmap='plasma', s=80, alpha=0.9,
                            linewidths=0, zorder=2)
        cb_g = plt.colorbar(sc_g, ax=ax_g, pad=0.02)
        cb_g.set_label(r'$\varepsilon(x)$', color='white')
        plt.setp(cb_g.ax.yaxis.get_ticklabels(), color='white')
        cb_g.ax.yaxis.set_tick_params(color='white')

    ax_g.set_xlim(-0.5, torus.L - 0.5)
    ax_g.set_ylim(-0.5, torus.L - 0.5)
    ax_g.set_title(fr'Torus graph [{gauge}] t={t_slice}, $\varepsilon(x)$',
                   color='white', fontsize=9)
    ax_g.set_xlabel(r'$x_0$', color='white')
    ax_g.set_ylabel(r'$x_1$', color='white')

    # Global title
    status = '✓' if sa.passes_criterion else '✗'
    fig.suptitle(
        f'Quantumograph v14 — Z⁴ Torus Analysis\n'
        f'L={torus.L}  N={torus.N}  gauge={torus.gauge_group}  '
        f'dₛ(∞)={sa.ds_inf:.3f} {status}',
        color='white', fontsize=12, y=1.01)

    if show:
        plt.tight_layout()
        plt.show()
    return fig


# ──────────────────────────────────────────────────────────────────────────────
# 4. Animated time-slice scan
# ──────────────────────────────────────────────────────────────────────────────

def animate_torus(torus: 'Z4Torus',
                  pa:    'PlaquetteAnalysis',
                  field: str = 'energy',
                  interval: int = 150,
                  show: bool = True) -> FuncAnimation:
    """
    Animate scan through t-slices of the torus, coloured by a scalar field.

    Parameters
    ----------
    torus    : Z4Torus
    pa       : PlaquetteAnalysis (run())
    field    : 'energy' | 'charge' — which scalar field to colour by
    interval : ms between frames
    show     : call plt.show()
    """
    scalar = pa.energy_density if field == 'energy' else pa.q_density
    cmap   = 'plasma'          if field == 'energy' else 'coolwarm'
    label  = r'$\varepsilon(x)$' if field == 'energy' else r'$q(x)$'

    vmin, vmax = scalar.min(), scalar.max()

    fig, ax = plt.subplots(figsize=(7, 7), facecolor='#080818')
    ax.set_facecolor('#0a0a1a')
    ax.set_xlim(-0.5, torus.L - 0.5)
    ax.set_ylim(-0.5, torus.L - 0.5)
    ax.tick_params(colors='white'); ax.spines[:].set_color('#222244')
    ax.set_xlabel(r'$x_0$', color='white'); ax.set_ylabel(r'$x_1$', color='white')

    sc   = ax.scatter([], [], c=[], cmap=cmap, vmin=vmin, vmax=vmax,
                      s=100, alpha=0.9, linewidths=0)
    cb   = plt.colorbar(sc, ax=ax, pad=0.02)
    cb.set_label(label, color='white')
    plt.setp(cb.ax.yaxis.get_ticklabels(), color='white')
    cb.ax.yaxis.set_tick_params(color='white')

    lc_container = [None]
    title_obj    = ax.set_title('', color='white', fontsize=11)

    # Precompute all slices / Предвычисляем все срезы
    slices = []
    for t in range(torus.L):
        node_ids, xs, ys, _, edges = _project_slice(torus, t, scalar)
        slices.append((node_ids, xs, ys, edges))

    def update(frame):
        t = frame % torus.L
        node_ids, xs, ys, edges = slices[t]

        sc.set_offsets(np.column_stack([xs, ys]) if len(xs) > 0
                       else np.empty((0, 2)))
        sc.set_array(scalar[node_ids] if len(node_ids) > 0
                     else np.array([]))

        if lc_container[0] is not None:
            lc_container[0].remove()
        if edges:
            segs = [[(x0, y0), (x1, y1)] for (x0,y0),(x1,y1) in edges]
            lc   = LineCollection(segs, colors='#1a3a5a',
                                  linewidths=0.7, alpha=0.5, zorder=1)
            ax.add_collection(lc)
            lc_container[0] = lc
        else:
            lc_container[0] = None

        title_obj.set_text(
            f'Quantumograph v14 — Z⁴ Torus  '
            f'(L={torus.L}, t-slice = {t}/{torus.L-1})\n'
            f'Nodes in slice: {len(node_ids)}   '
            f'Edges: {len(edges)}')
        return sc,

    ani = FuncAnimation(fig, update, frames=torus.L * 3,
                        interval=interval, blit=False, repeat=True)
    if show:
        plt.tight_layout()
        plt.show()
    return ani


# ──────────────────────────────────────────────────────────────────────────────
# 5. FSS extrapolation plot
# ──────────────────────────────────────────────────────────────────────────────

def plot_fss(L_arr: np.ndarray,
             ds_arr: np.ndarray,
             ds_inf: float,
             coeffs: np.ndarray,
             show: bool = True) -> plt.Figure:
    """
    Plot finite-size scaling fit dₛ(L) with extrapolation to L→∞.
    """
    fig, ax = plt.subplots(figsize=(7, 5), facecolor='#0a0a1a')
    ax.set_facecolor('#0a0a1a'); ax.tick_params(colors='white')
    ax.spines[:].set_color('#333355')

    inv_L = 1.0 / np.asarray(L_arr, dtype=float)
    L_fit = np.linspace(0, inv_L.max() * 1.05, 200)
    A     = np.vander(L_fit, len(coeffs), increasing=True)
    ds_fit = A @ coeffs

    ax.plot(inv_L, ds_arr, 'o', color='#00d4ff', ms=8,
            label=r'$d_s(L)$ measured')
    ax.plot(L_fit, ds_fit, '-', color='#ffaa00', lw=1.5,
            label=fr'FSS fit, $d_s(\infty)={ds_inf:.3f}$')
    ax.axhline(4.0, ls='--', color='white', lw=0.8, alpha=0.4,
               label=r'$d_s = 4$')
    ax.axvline(0,   ls=':',  color='#888888', lw=0.8)
    ax.scatter([0], [ds_inf], color='#00ff88', s=100, zorder=5,
               label=fr'$d_s(\infty) = {ds_inf:.3f}$')

    ax.set_xlabel(r'$1/L$', color='white')
    ax.set_ylabel(r'$d_s$', color='white')
    ax.set_title('Finite-Size Scaling Extrapolation', color='white', fontsize=11)
    leg = ax.legend(facecolor='#111133', edgecolor='#333355', labelcolor='white')

    if show:
        plt.tight_layout(); plt.show()
    return fig
