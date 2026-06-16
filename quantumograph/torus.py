"""
quantumograph.torus
===================
Z⁴ hypercubic torus — the fundamental spacetime substrate of Quantumograph v14.

The torus T = Z_L⁴ is the discrete analogue of the 4-torus T⁴ with side L.
Nodes are labeled by 4-tuples (n₀, n₁, n₂, n₃) ∈ {0,…,L-1}⁴.
Edges connect nearest neighbours with periodic identification n ~ n + L.

Mathematical reference:
  Quantumograph v14 monograph, §2: "Discrete Spacetime on the Z⁴ Torus"

Author  : © 2025–2026 Sergej Materov <sergejmaterov2@gmail.com>
License : CC BY-NC 4.0
"""

# Z⁴ hypercubic torus — фундаментальный субстрат пространства-времени Quantumograph v14.
# Тор T = Z_L⁴ — дискретный аналог 4-тора T⁴ со стороной L.

from __future__ import annotations

import numpy as np
import networkx as nx
from itertools import product
from typing import Tuple, Dict, Optional

Coord4 = Tuple[int, int, int, int]


class Z4Torus:
    """
    Hypercubic torus Z_L⁴ with periodic boundary conditions.

    Parameters
    ----------
    L : int
        Side length. Total nodes = L⁴.
        Desktop recommended: L ≤ 8  (~4 096 nodes)
        HPC recommended    : L ≥ 16 (~65 536 nodes)
    gauge_group : str
        'U1' | 'SU2' | 'SU3' — determines link variable dimension.
        Currently used for phase assignment on edges.
    seed : int
        RNG seed for reproducible gauge configurations.
    """

    def __init__(self, L: int = 6, gauge_group: str = 'U1', seed: int = 42):
        if L < 2:
            raise ValueError("L must be ≥ 2")
        self.L           = L
        self.gauge_group = gauge_group
        self.seed        = seed
        self.N           = L ** 4          # total number of nodes / кол-во узлов

        self._rng        = np.random.default_rng(seed)
        self._graph: Optional[nx.Graph] = None
        self._coord_to_id: Dict[Coord4, int] = {}
        self._id_to_coord: Dict[int, Coord4] = {}
        self._link_phases: Optional[np.ndarray] = None  # (E, d_gauge)

        self._build()

    # ──────────────────────────────────────────────────────────────────────
    # Internal construction / Внутренняя сборка
    # ──────────────────────────────────────────────────────────────────────

    def _build(self) -> None:
        """Build the torus graph with nearest-neighbour edges."""
        L = self.L
        G = nx.Graph()

        # Assign integer IDs to 4-tuples / Присваиваем целочисленные ID четвёркам
        for idx, coord in enumerate(product(range(L), repeat=4)):
            self._coord_to_id[coord] = idx
            self._id_to_coord[idx]   = coord
            G.add_node(idx, coord=coord)

        # Nearest-neighbour edges with PBC / Рёбра ближайших соседей с ПГУ
        # We add only the +μ direction edge per node (forward step).
        # This correctly handles PBC: node at coord L-1 wraps to 0.
        # Each undirected edge is added exactly once (forward direction only).
        seen  = set()
        edges = []
        for coord, idx in self._coord_to_id.items():
            for mu in range(4):                        # 4 directions / 4 направления
                nb      = list(coord)
                nb[mu]  = (nb[mu] + 1) % L            # periodic / периодично
                nb_idx  = self._coord_to_id[tuple(nb)]
                key     = (min(idx, nb_idx), max(idx, nb_idx), mu)
                if key not in seen:
                    seen.add(key)
                    edges.append((idx, nb_idx, mu))    # (u, v, direction)

        # Store direction as edge attribute / Направление ребра как атрибут
        for u, v, mu in edges:
            G.add_edge(u, v, mu=mu, weight=1.0)

        self._graph = G
        self._assign_link_phases()

    def _assign_link_phases(self) -> None:
        """
        Assign U(1)/SU(2)/SU(3) link variables θ_{x,μ} to each directed edge.

        U(1)  : θ ∈ [0, 2π)                           — 1 real phase per link
        SU(2) : (α, β, γ) ∈ [0,π) × [0,2π) × [0,4π) — 3 Euler angles
        SU(3) : 8 real parameters (Gell-Mann basis)

        These are stored as raw parameters; actual group elements are
        computed on demand by plaquette.py.
        """
        E = self._graph.number_of_edges()
        if self.gauge_group == 'U1':
            d = 1
            phases = self._rng.uniform(0, 2 * np.pi, (E, d))
        elif self.gauge_group == 'SU2':
            d = 3
            phases = self._rng.uniform(0, np.pi, (E, d))
            phases[:, 1] *= 2   # β ∈ [0, 2π)
            phases[:, 2] *= 4   # γ ∈ [0, 4π)
        elif self.gauge_group == 'SU3':
            d = 8
            phases = self._rng.uniform(0, 2 * np.pi, (E, d))
        else:
            raise ValueError(f"Unknown gauge group: {self.gauge_group}")

        self._link_phases = phases

        # Store per-edge / Сохраняем на каждом ребре
        for eidx, (u, v) in enumerate(self._graph.edges()):
            self._graph[u][v]['phase'] = phases[eidx]

    # ──────────────────────────────────────────────────────────────────────
    # Public API
    # ──────────────────────────────────────────────────────────────────────

    @property
    def graph(self) -> nx.Graph:
        """NetworkX graph of the torus / Граф тора."""
        return self._graph

    @property
    def coords(self) -> np.ndarray:
        """
        Node coordinates as (N, 4) array.
        coords[i] = (n₀, n₁, n₂, n₃) for node i.
        """
        arr = np.zeros((self.N, 4), dtype=np.int32)
        for idx, coord in self._id_to_coord.items():
            arr[idx] = coord
        return arr

    def coord(self, node_id: int) -> Coord4:
        """4-tuple coordinate of a node / 4-кортеж координат узла."""
        return self._id_to_coord[node_id]

    def node_id(self, coord: Coord4) -> int:
        """Node ID from 4-tuple coordinate / ID узла по координатам."""
        coord = tuple(c % self.L for c in coord)
        return self._coord_to_id[coord]

    def link_phase(self, u: int, v: int) -> np.ndarray:
        """
        Gauge link variable θ_{u→v}.
        Returns array of shape (d_gauge,).
        """
        return self._graph[u][v]['phase']

    def laplacian_matrix(self, normalized: bool = False) -> np.ndarray:
        """
        Graph Laplacian L = D - A.

        On Z_L⁴: each node has degree 2·4 = 8 (two neighbours per direction).
        Returns dense (N, N) array — use sparse version for large L.

        Parameters
        ----------
        normalized : bool
            If True, returns L_sym = D^{-1/2} L D^{-1/2}.
        """
        if normalized:
            return nx.normalized_laplacian_matrix(self._graph).toarray()
        return nx.laplacian_matrix(self._graph).toarray().astype(float)

    def laplacian_sparse(self):
        """
        Sparse Laplacian (scipy.sparse.csr_matrix).
        Recommended for L ≥ 6 / Рекомендуется для L ≥ 6.
        """
        from scipy.sparse import csr_matrix
        return csr_matrix(nx.laplacian_matrix(self._graph))

    def __repr__(self) -> str:
        return (f"Z4Torus(L={self.L}, N={self.N}, "
                f"edges={self._graph.number_of_edges()}, "
                f"gauge={self.gauge_group})")
