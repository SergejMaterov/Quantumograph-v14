"""
tests/test_core.py — smoke tests for Quantumograph v14
Run with: pytest tests/ -v
© 2025–2026 Sergej Materov <sergejmaterov2@gmail.com>
"""

import numpy as np
import pytest

from quantumograph.torus    import Z4Torus
from quantumograph.spectral import (heat_kernel_exact, spectral_dimension,
                                     sigma_midpoints, check_4d_criterion,
                                     fss_extrapolate, compute_spectrum_sparse)
from quantumograph.plaquette import (compute_plaquettes_u1, hamiltonian_u1,
                                      plaquette_energy_density,
                                      total_topological_charge, MU_NU_PAIRS)


# ── Torus ─────────────────────────────────────────────────────────────────────

class TestZ4Torus:
    def test_node_count(self):
        for L in [3, 4, 5]:
            t = Z4Torus(L=L)
            assert t.N == L**4

    def test_edge_count(self):
        """Each node has degree 8 (2 neighbours × 4 directions)."""
        t = Z4Torus(L=4)
        assert t.graph.number_of_edges() == t.N * 4   # 8 edges/2 (undirected)

    def test_periodic_boundary(self):
        """Node (0,0,0,0) and (L-1,0,0,0) are neighbours in direction 0."""
        L = 4
        t = Z4Torus(L=L)
        id0  = t.node_id((0,   0, 0, 0))
        id_n = t.node_id((L-1, 0, 0, 0))
        assert t.graph.has_edge(id0, id_n)

    def test_coords_roundtrip(self):
        t = Z4Torus(L=5)
        for node in range(0, t.N, 100):
            coord = t.coord(node)
            assert t.node_id(coord) == node

    def test_link_phases_u1(self):
        t = Z4Torus(L=4, gauge_group='U1')
        for u, v, data in t.graph.edges(data=True):
            phase = data['phase']
            assert 0.0 <= float(phase[0]) < 2*np.pi + 1e-9

    def test_laplacian_positive_semidefinite(self):
        t    = Z4Torus(L=3)
        L_m  = t.laplacian_matrix()
        eigv = np.linalg.eigvalsh(L_m)
        assert eigv[0] >= -1e-10               # smallest ~ 0

    def test_laplacian_zero_eigenvalue(self):
        """Connected graph → exactly one zero eigenvalue."""
        t    = Z4Torus(L=3)
        L_m  = t.laplacian_matrix()
        eigv = np.linalg.eigvalsh(L_m)
        n_zero = np.sum(np.abs(eigv) < 1e-8)
        assert n_zero == 1


# ── Spectral dimension ────────────────────────────────────────────────────────

class TestSpectral:
    @pytest.fixture(scope='class')
    def torus_l4(self):
        return Z4Torus(L=4)

    def test_heat_kernel_at_zero(self, torus_l4):
        """K(σ→0) → 1."""
        t   = torus_l4
        eig = np.linalg.eigvalsh(t.laplacian_matrix())
        eig = np.sort(np.maximum(eig, 0))
        sig = np.array([1e-6])
        K   = heat_kernel_exact(eig, sig)
        assert abs(K[0] - 1.0) < 0.01

    def test_heat_kernel_decreasing(self, torus_l4):
        """K(σ) must be monotonically decreasing."""
        t   = torus_l4
        eig = np.linalg.eigvalsh(t.laplacian_matrix())
        eig = np.sort(np.maximum(eig, 0))
        sig = np.geomspace(0.01, 20, 50)
        K   = heat_kernel_exact(eig, sig)
        assert np.all(np.diff(K) <= 1e-10)

    def test_spectral_dimension_positive(self, torus_l4):
        t   = torus_l4
        eig = np.linalg.eigvalsh(t.laplacian_matrix())
        eig = np.sort(np.maximum(eig, 0))
        sig = np.geomspace(0.01, 10, 100)
        K   = heat_kernel_exact(eig, sig)
        ds  = spectral_dimension(sig, K)
        # dₛ should be positive in the bulk
        assert np.all(ds[5:-5] >= 0)

    def test_fss_extrapolate(self):
        """Synthetic dₛ(L) = 4 + 1/L → dₛ(∞) = 4."""
        L_arr  = np.array([3, 4, 5, 6, 7], dtype=float)
        ds_arr = 4.0 + 1.0 / L_arr
        ds_inf, _ = fss_extrapolate(L_arr, ds_arr, order=1)
        assert abs(ds_inf - 4.0) < 0.05

    def test_check_4d_criterion(self):
        assert check_4d_criterion(4.05)  is True
        assert check_4d_criterion(3.85)  is False
        assert check_4d_criterion(4.099) is True
        assert check_4d_criterion(4.101) is False

    def test_sparse_spectrum_shape(self, torus_l4):
        k    = 64
        eigs = compute_spectrum_sparse(torus_l4, k=k)
        assert eigs.shape == (k,)
        assert np.all(eigs >= -1e-9)


# ── Plaquette ─────────────────────────────────────────────────────────────────

class TestPlaquette:
    @pytest.fixture(scope='class')
    def torus_u1(self):
        return Z4Torus(L=4, gauge_group='U1', seed=0)

    def test_plaquette_shape(self, torus_u1):
        plaq = compute_plaquettes_u1(torus_u1)
        assert plaq.shape == (torus_u1.N, 6)

    def test_plaquette_range(self, torus_u1):
        """W_{x,μν} is a sum of 4 phases ∈ [0,2π), so |W| ≤ 8π."""
        plaq = compute_plaquettes_u1(torus_u1)
        assert np.all(np.abs(plaq) <= 8 * np.pi + 1e-9)

    def test_hamiltonian_finite(self, torus_u1):
        H = hamiltonian_u1(torus_u1)
        assert np.isfinite(H)

    def test_energy_density_shape(self, torus_u1):
        eps = plaquette_energy_density(torus_u1)
        assert eps.shape == (torus_u1.N,)
        assert np.all(np.isfinite(eps))

    def test_hamiltonian_equals_sum_density(self, torus_u1):
        H   = hamiltonian_u1(torus_u1)
        eps = plaquette_energy_density(torus_u1)
        assert abs(H - eps.sum()) < 1e-6

    def test_ordered_config_min_energy(self):
        """Ordered config (all θ=0) should give minimum energy H = -J·6·N."""
        import networkx as nx
        t = Z4Torus(L=3, gauge_group='U1', seed=0)
        # Set all phases to 0
        for u, v in t.graph.edges():
            t.graph[u][v]['phase'] = np.array([0.0])
        H     = hamiltonian_u1(t, J=1.0)
        H_min = -1.0 * 6 * t.N
        assert abs(H - H_min) < 1e-6

    def test_mu_nu_pairs_count(self):
        """C(4,2) = 6 pairs."""
        assert len(MU_NU_PAIRS) == 6

    def test_topological_charge_finite(self, torus_u1):
        Q = total_topological_charge(torus_u1)
        assert np.isfinite(Q)


# ── SU(2) plaquette ───────────────────────────────────────────────────────────

class TestSU2:
    @pytest.fixture(scope='class')
    def torus_su2(self):
        return Z4Torus(L=4, gauge_group='SU2', seed=0)

    def test_su2_matrix_unitary(self):
        """SU(2) matrices must be unitary: U†U = I."""
        from quantumograph import su2_matrix
        rng = np.random.default_rng(0)
        for _ in range(20):
            params = rng.uniform(0, np.pi, 3)
            U = su2_matrix(params)
            assert np.allclose(U.conj().T @ U, np.eye(2), atol=1e-12)

    def test_su2_matrix_det1(self):
        """det(U) = 1 for SU(2)."""
        from quantumograph import su2_matrix
        rng = np.random.default_rng(1)
        for _ in range(20):
            params = rng.uniform(0, np.pi, 3)
            U = su2_matrix(params)
            assert abs(np.linalg.det(U) - 1.0) < 1e-12

    def test_su2_plaquette_shape(self, torus_su2):
        from quantumograph import compute_plaquettes_su2
        plaq = compute_plaquettes_su2(torus_su2)
        assert plaq.shape == (torus_su2.N, 6)

    def test_su2_plaquette_range(self, torus_su2):
        """Re Tr[P]/2 ∈ [-1, 1] for SU(2) fundamental."""
        from quantumograph import compute_plaquettes_su2
        plaq = compute_plaquettes_su2(torus_su2)
        assert np.all(plaq >= -1.0 - 1e-9)
        assert np.all(plaq <=  1.0 + 1e-9)

    def test_su2_hamiltonian_finite(self, torus_su2):
        from quantumograph import hamiltonian_su2
        H = hamiltonian_su2(torus_su2)
        assert np.isfinite(H)

    def test_su2_ordered_config(self):
        """Ordered config (U=I on all links) → H = -J·6·N."""
        t = Z4Torus(L=3, gauge_group='SU2', seed=0)
        # Set all link phases to zero → U = I
        for u, v in t.graph.edges():
            t.graph[u][v]['phase'] = np.zeros(3)
        from quantumograph import hamiltonian_su2
        H     = hamiltonian_su2(t, J=1.0)
        H_min = -1.0 * 6 * t.N           # -J * 6 plaquettes/site * N * Re Tr[I]/2 = -J*6*N
        assert abs(H - H_min) < 1e-6

    def test_wilson_loop_finite(self, torus_su2):
        from quantumograph import wilson_loop_su2
        W = wilson_loop_su2(torus_su2, mu=0, nu=1, R=1, T=1)
        assert np.isfinite(W)
        assert -1.0 - 1e-9 <= W <= 1.0 + 1e-9

    def test_su2_analysis_runs(self, torus_su2):
        from quantumograph import PlaquetteAnalysis
        pa = PlaquetteAnalysis(torus_su2).run()
        assert pa.H_total is not None
        assert pa.wilson_loops is not None
        assert (1,1) in pa.wilson_loops


# ── SU(3) plaquette ───────────────────────────────────────────────────────────

class TestSU3:
    @pytest.fixture(scope='class')
    def torus_su3(self):
        return Z4Torus(L=4, gauge_group='SU3', seed=0)

    def test_su3_matrix_unitary(self):
        """SU(3) matrices must be unitary."""
        from quantumograph import su3_matrix
        rng = np.random.default_rng(0)
        for _ in range(10):
            params = rng.uniform(0, np.pi, 8)
            U = su3_matrix(params)
            assert np.allclose(U.conj().T @ U, np.eye(3), atol=1e-10)

    def test_su3_matrix_det1(self):
        """det(U) = 1 for SU(3)."""
        from quantumograph import su3_matrix
        rng = np.random.default_rng(2)
        for _ in range(10):
            params = rng.uniform(0, np.pi, 8)
            U = su3_matrix(params)
            assert abs(np.linalg.det(U) - 1.0) < 1e-10

    def test_su3_plaquette_shape(self, torus_su3):
        from quantumograph import compute_plaquettes_su3
        plaq = compute_plaquettes_su3(torus_su3)
        assert plaq.shape == (torus_su3.N, 6)

    def test_su3_plaquette_range(self, torus_su3):
        """Re Tr[P]/3 ∈ [-1, 1] for SU(3) fundamental."""
        from quantumograph import compute_plaquettes_su3
        plaq = compute_plaquettes_su3(torus_su3)
        assert np.all(plaq >= -1.0 - 1e-9)
        assert np.all(plaq <=  1.0 + 1e-9)

    def test_su3_ordered_config(self):
        """Ordered config (all params=0 → U=I) → H = -J·6·N."""
        t = Z4Torus(L=3, gauge_group='SU3', seed=0)
        for u, v in t.graph.edges():
            t.graph[u][v]['phase'] = np.zeros(8)
        from quantumograph import hamiltonian_su3
        H     = hamiltonian_su3(t, J=1.0)
        H_min = -1.0 * 6 * t.N
        assert abs(H - H_min) < 1e-6

    def test_polyakov_loop_finite(self, torus_su3):
        from quantumograph import polyakov_loop_su3
        P = polyakov_loop_su3(torus_su3, mu=3)
        assert np.isfinite(P.real)
        assert abs(P.real) <= 1.0 + 1e-9

    def test_su3_analysis_runs(self, torus_su3):
        from quantumograph import PlaquetteAnalysis
        pa = PlaquetteAnalysis(torus_su3).run()
        assert pa.H_total is not None
        assert pa.polyakov is not None

    def test_gell_mann_traceless(self):
        """All Gell-Mann matrices must be traceless."""
        from quantumograph.plaquette import _gell_mann_matrices
        lam = _gell_mann_matrices()
        for k in range(8):
            assert abs(np.trace(lam[k])) < 1e-12

    def test_gell_mann_hermitian(self):
        """All Gell-Mann matrices must be Hermitian."""
        from quantumograph.plaquette import _gell_mann_matrices
        lam = _gell_mann_matrices()
        for k in range(8):
            assert np.allclose(lam[k], lam[k].conj().T, atol=1e-12)
