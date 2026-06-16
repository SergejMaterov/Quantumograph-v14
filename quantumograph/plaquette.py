"""
quantumograph.plaquette
=======================
Plaquette Hamiltonian on the Z⁴ torus for U(1), SU(2), SU(3) gauge groups.

Theory (Quantumograph v14 monograph, §5):

  A plaquette □_{x,μν} is the minimal closed loop on the torus:
      x → x+μ̂ → x+μ̂+ν̂ → x+ν̂ → x

  Wilson action (all groups):
      S_W = -J · Re Tr[ U_{x,μ} U_{x+μ̂,ν} U†_{x+ν̂+μ̂,μ} U†_{x+ν̂,ν} ] / d_r

  where d_r is the dimension of the representation:
      U(1)  : d_r = 1,  U = exp(iθ)                    scalar phase
      SU(2) : d_r = 2,  U = exp(i α·σ/2)               2×2 unitary
      SU(3) : d_r = 3,  U = exp(i α·λ/2)               3×3 unitary

  U(1) reduces to: W_{x,μν} = θ_{x,μ}+θ_{x+μ̂,ν}-θ_{x+μ̂+ν̂,μ}-θ_{x+ν̂,ν}
                   H = -J Σ cos(W_{x,μν})

  SU(2): H = -J Σ Re Tr[P_{x,μν}] / 2
  SU(3): H = -J Σ Re Tr[P_{x,μν}] / 3

  Z⁴ has C(4,2) = 6 plaquette orientations per node.

Author  : © 2025–2026 Sergej Materov <sergejmaterov2@gmail.com>
License : CC BY-NC 4.0
"""

from __future__ import annotations

import numpy as np
from itertools import combinations, permutations
from typing import List, Tuple, Optional

# 6 independent μ<ν pairs on Z⁴ / 6 независимых пар μ<ν на Z⁴
MU_NU_PAIRS: List[Tuple[int, int]] = list(combinations(range(4), 2))


# ══════════════════════════════════════════════════════════════════════════════
# GROUP ELEMENT CONSTRUCTORS
# Конструкторы групповых элементов
# ══════════════════════════════════════════════════════════════════════════════

def u1_matrix(theta: float) -> complex:
    """
    U(1) link variable: U = exp(iθ)
    θ ∈ [0, 2π)
    """
    return complex(np.cos(theta), np.sin(theta))


def su2_matrix(params: np.ndarray) -> np.ndarray:
    """
    SU(2) link variable via exponential map:
      U = exp(i · αₖ · σₖ / 2)  where σₖ are Pauli matrices.

    Parameters stored in torus as Euler angles (α, β, γ):
      U = exp(i γ σ₃/2) · exp(i β σ₂/2) · exp(i α σ₃/2)

    Returns (2,2) complex array.

    SU(2) через отображение экспоненты и углы Эйлера.
    """
    alpha, beta, gamma = params[0], params[1], params[2]

    # Euler angle decomposition: ZYZ convention
    # Разложение Эйлера: соглашение ZYZ
    ca, sa = np.cos(alpha/2), np.sin(alpha/2)
    cb, sb = np.cos(beta/2),  np.sin(beta/2)
    cg, sg = np.cos(gamma/2), np.sin(gamma/2)

    # U = Rz(γ) · Ry(β) · Rz(α)
    U = np.array([
        [ cg*cb*ca - sg*sa*1j*0 + (cg*cb*sa*1j + sg*ca*1j)*0,
          -cg*sb*np.exp(1j*alpha/2) ],   # placeholder row
        [  sg*sb*np.exp(-1j*alpha/2),
           sg*cb*ca + cg*sa*1j*0    ]
    ], dtype=complex)

    # Direct construction is cleaner:
    # Прямое построение через параметры
    U = np.array([
        [ (cg*cb + 1j*0)*ca - 1j*(cg*cb)*sa,  0 ],
        [ 0, 0 ]
    ], dtype=complex)

    # Most reliable: explicit ZYZ Euler product
    Rz_a = np.array([[np.exp(-1j*alpha/2), 0],
                     [0, np.exp(1j*alpha/2)]], dtype=complex)
    Ry_b = np.array([[cb,  -sb],
                     [sb,   cb]], dtype=complex)
    Rz_g = np.array([[np.exp(-1j*gamma/2), 0],
                     [0, np.exp(1j*gamma/2)]], dtype=complex)

    return Rz_g @ Ry_b @ Rz_a


def su2_matrix_vec(params_arr: np.ndarray) -> np.ndarray:
    """
    Vectorised SU(2): params_arr (E, 3) → matrices (E, 2, 2).
    Векторизованный вариант для массива параметров.
    """
    alpha = params_arr[:, 0]
    beta  = params_arr[:, 1]
    gamma = params_arr[:, 2]

    ca2, sa2 = np.cos(alpha/2), np.sin(alpha/2)
    cb2, sb2 = np.cos(beta/2),  np.sin(beta/2)
    cg2, sg2 = np.cos(gamma/2), np.sin(gamma/2)

    # Rz(γ) @ Ry(β) @ Rz(α) vectorised
    # U[0,0] = exp(-iγ/2)(cos(β/2)exp(-iα/2))
    # U[0,1] = exp(-iγ/2)(-sin(β/2)exp(iα/2))
    # U[1,0] = exp(iγ/2)(sin(β/2)exp(-iα/2))
    # U[1,1] = exp(iγ/2)(cos(β/2)exp(iα/2))
    E = len(alpha)
    U = np.zeros((E, 2, 2), dtype=complex)

    exp_mg2 = np.exp(-1j * gamma / 2)
    exp_pg2 = np.exp( 1j * gamma / 2)
    exp_ma2 = np.exp(-1j * alpha / 2)
    exp_pa2 = np.exp( 1j * alpha / 2)

    U[:, 0, 0] = exp_mg2 *  cb2 * exp_ma2
    U[:, 0, 1] = exp_mg2 * -sb2 * exp_pa2
    U[:, 1, 0] = exp_pg2 *  sb2 * exp_ma2
    U[:, 1, 1] = exp_pg2 *  cb2 * exp_pa2

    return U


def su3_matrix(params: np.ndarray) -> np.ndarray:
    """
    SU(3) link variable via exponential map:
      U = exp(i · αₖ · λₖ / 2)

    where λₖ (k=1..8) are the Gell-Mann matrices.
    params: (8,) array of real coefficients.

    Returns (3,3) complex array.

    SU(3) через матрицы Гелл-Манна.
    """
    # Gell-Mann matrices / Матрицы Гелл-Манна
    lam = _gell_mann_matrices()                # (8, 3, 3) complex
    # Generator: H = i/2 · Σₖ αₖ λₖ
    H = 0.5j * np.einsum('k,kij->ij', params, lam)  # (3,3)
    # Matrix exponential via diagonalisation
    # Матричная экспонента через диагонализацию
    return _matrix_exp_3x3(H)


def su3_matrix_vec(params_arr: np.ndarray) -> np.ndarray:
    """
    Vectorised SU(3): params_arr (E, 8) → matrices (E, 3, 3).
    Uses scipy.linalg.expm per element (no closed-form for SU(3)).
    Векторизованный вариант для массива параметров.
    """
    from scipy.linalg import expm
    lam = _gell_mann_matrices()                # (8, 3, 3)
    E   = len(params_arr)
    U   = np.zeros((E, 3, 3), dtype=complex)
    for e in range(E):
        H      = 0.5j * np.einsum('k,kij->ij', params_arr[e], lam)
        U[e]   = expm(H)
    return U


def _gell_mann_matrices() -> np.ndarray:
    """
    The 8 Gell-Mann matrices λ₁ … λ₈ (traceless Hermitian 3×3).
    8 матриц Гелл-Манна.
    """
    lam = np.zeros((8, 3, 3), dtype=complex)

    lam[0] = [[0,1,0],[1,0,0],[0,0,0]]          # λ₁
    lam[1] = [[0,-1j,0],[1j,0,0],[0,0,0]]        # λ₂
    lam[2] = [[1,0,0],[0,-1,0],[0,0,0]]          # λ₃
    lam[3] = [[0,0,1],[0,0,0],[1,0,0]]           # λ₄
    lam[4] = [[0,0,-1j],[0,0,0],[1j,0,0]]        # λ₅
    lam[5] = [[0,0,0],[0,0,1],[0,1,0]]           # λ₆
    lam[6] = [[0,0,0],[0,0,-1j],[0,1j,0]]        # λ₇
    lam[7] = np.array([[1,0,0],[0,1,0],[0,0,-2]],
                       dtype=complex) / np.sqrt(3)  # λ₈
    return lam


def _matrix_exp_3x3(H: np.ndarray) -> np.ndarray:
    """Matrix exponential for general 3x3 complex matrix via scipy."""
    from scipy.linalg import expm
    return expm(H)
def _build_link_table_u1(torus) -> np.ndarray:
    """
    phase_mat[node, mu] = θ_{node,μ} ∈ [0,2π)
    (N, 4) float array.
    """
    L         = torus.L
    N         = torus.N
    phase_mat = np.zeros((N, 4), dtype=np.float64)

    for u, v, data in torus.graph.edges(data=True):
        mu    = data['mu']
        phase = data['phase']
        theta = float(phase[0]) if hasattr(phase, '__len__') else float(phase)
        cu    = torus.coord(u)
        if cu[mu] == (torus.coord(v)[mu] - 1) % L:
            phase_mat[u, mu] = theta
        else:
            phase_mat[v, mu] = theta

    return phase_mat


def _build_link_table_su2(torus) -> np.ndarray:
    """
    U_mat[node, mu] = SU(2) matrix (2×2) for link (node, μ)
    Returns (N, 4, 2, 2) complex array.
    SU(2) матрицы для всех линков тора.
    """
    L     = torus.L
    N     = torus.N
    U_mat = np.zeros((N, 4, 2, 2), dtype=complex)

    for u, v, data in torus.graph.edges(data=True):
        mu     = data['mu']
        params = data['phase']                  # (3,) Euler angles
        U      = su2_matrix(params)             # (2,2)
        cu     = torus.coord(u)
        if cu[mu] == (torus.coord(v)[mu] - 1) % L:
            U_mat[u, mu] = U
            U_mat[v, mu] = U.conj().T           # U†_{x+μ̂,μ} = U†
        else:
            U_mat[v, mu] = U
            U_mat[u, mu] = U.conj().T

    return U_mat


def _build_link_table_su3(torus) -> np.ndarray:
    """
    U_mat[node, mu] = SU(3) matrix (3×3) for link (node, μ)
    Returns (N, 4, 3, 3) complex array.
    SU(3) матрицы для всех линков тора.
    """
    from scipy.linalg import expm
    L     = torus.L
    N     = torus.N
    lam   = _gell_mann_matrices()
    U_mat = np.zeros((N, 4, 3, 3), dtype=complex)

    for u, v, data in torus.graph.edges(data=True):
        mu     = data['mu']
        params = data['phase']                  # (8,) Gell-Mann coefficients
        H      = 0.5j * np.einsum('k,kij->ij', params, lam)
        U      = expm(H)                        # (3,3)
        cu     = torus.coord(u)
        if cu[mu] == (torus.coord(v)[mu] - 1) % L:
            U_mat[u, mu] = U
            U_mat[v, mu] = U.conj().T
        else:
            U_mat[v, mu] = U
            U_mat[u, mu] = U.conj().T

    return U_mat


# ══════════════════════════════════════════════════════════════════════════════
# U(1) PLAQUETTES (fast vectorised)
# ══════════════════════════════════════════════════════════════════════════════

def compute_plaquettes_u1(torus, J: float = 1.0) -> np.ndarray:
    """
    U(1) plaquette phase W_{x,μν} = θ_{x,μ} + θ_{x+μ̂,ν} - θ_{x+μ̂+ν̂,μ} - θ_{x+ν̂,ν}

    Returns (N, 6) float array.
    """
    if torus.gauge_group != 'U1':
        raise ValueError("compute_plaquettes_u1 requires gauge_group='U1'")

    L      = torus.L
    coords = torus.coords
    phases = _build_link_table_u1(torus)

    def shift(c, mu):
        c_new = c.copy(); c_new[:, mu] = (c_new[:, mu] + 1) % L; return c_new

    def c2id(c):
        return c[:,0]*L**3 + c[:,1]*L**2 + c[:,2]*L + c[:,3]

    plaq = np.zeros((torus.N, 6), dtype=np.float64)
    for pidx, (mu, nu) in enumerate(MU_NU_PAIRS):
        id_x    = c2id(coords)
        id_xmu  = c2id(shift(coords, mu))
        id_xnu  = c2id(shift(coords, nu))
        id_xmnu = c2id(shift(shift(coords, mu), nu))
        plaq[:, pidx] = (phases[id_x, mu] + phases[id_xmu, nu]
                       - phases[id_xmnu, mu] - phases[id_xnu, nu])
    return plaq


# ══════════════════════════════════════════════════════════════════════════════
# SU(2) PLAQUETTES
# ══════════════════════════════════════════════════════════════════════════════

def compute_plaquettes_su2(torus) -> np.ndarray:
    """
    SU(2) Wilson plaquette:
      P_{x,μν} = U_{x,μ} · U_{x+μ̂,ν} · U†_{x+μ̂+ν̂,μ} · U†_{x+ν̂,ν}

    Returns (N, 6) float array of Re Tr[P_{x,μν}] / 2.
    The Wilson action is H = -J Σ Re Tr[P] / 2.

    Плакетный оператор SU(2): упорядоченное произведение матриц вдоль □.
    """
    if torus.gauge_group != 'SU2':
        raise ValueError("compute_plaquettes_su2 requires gauge_group='SU2'")

    L      = torus.L
    coords = torus.coords
    U      = _build_link_table_su2(torus)       # (N, 4, 2, 2)

    def shift(c, mu):
        c_new = c.copy(); c_new[:, mu] = (c_new[:, mu] + 1) % L; return c_new

    def c2id(c):
        return c[:,0]*L**3 + c[:,1]*L**2 + c[:,2]*L + c[:,3]

    plaq = np.zeros((torus.N, 6), dtype=np.float64)

    for pidx, (mu, nu) in enumerate(MU_NU_PAIRS):
        id_x    = c2id(coords)
        id_xmu  = c2id(shift(coords, mu))
        id_xnu  = c2id(shift(coords, nu))
        id_xmnu = c2id(shift(shift(coords, mu), nu))

        # P = U_{x,μ} @ U_{x+μ̂,ν} @ U†_{x+μ̂+ν̂,μ} @ U†_{x+ν̂,ν}
        # Vectorised over all N nodes simultaneously
        # Векторизованное вычисление для всех N узлов
        P = (U[id_x,   mu]                       # (N,2,2)
           @ U[id_xmu, nu]                        # (N,2,2)
           @ U[id_xmnu, mu].conj().swapaxes(-1,-2)  # U†
           @ U[id_xnu,  nu].conj().swapaxes(-1,-2))  # U†

        # Re Tr[P] / 2  (d_r = 2 for fundamental rep)
        plaq[:, pidx] = np.real(np.trace(P, axis1=-2, axis2=-1)) / 2.0

    return plaq


def wilson_loop_su2(torus, mu: int, nu: int,
                    R: int = 1, T: int = 1) -> float:
    """
    Rectangular R×T Wilson loop in the (μ,ν) plane for SU(2).
    W(R,T) = <Re Tr[ ∏ U ]> / 2

    Area law: W(R,T) ~ exp(-σ R T)  (confinement)
    Perimeter law: W(R,T) ~ exp(-μ (2R+2T))  (deconfinement)

    Прямоугольный вильсоновский контур R×T.
    """
    if torus.gauge_group != 'SU2':
        raise ValueError("wilson_loop_su2 requires gauge_group='SU2'")

    L      = torus.L
    coords = torus.coords
    U_mat  = _build_link_table_su2(torus)

    def shift_n(c, direction, n):
        c_new = c.copy()
        c_new[:, direction] = (c_new[:, direction] + n) % L
        return c_new

    def c2id(c):
        return c[:,0]*L**3 + c[:,1]*L**2 + c[:,2]*L + c[:,3]

    N   = torus.N
    W   = np.eye(2, dtype=complex)[np.newaxis].repeat(N, axis=0)  # (N,2,2)

    # Bottom edge: R steps in μ direction
    c = coords.copy()
    for _ in range(R):
        W = W @ U_mat[c2id(c), mu]
        c = shift_n(c, mu, 1)

    # Right edge: T steps in ν direction
    for _ in range(T):
        W = W @ U_mat[c2id(c), nu]
        c = shift_n(c, nu, 1)

    # Top edge: R steps back in μ direction (conjugate transpose)
    for _ in range(R):
        c = shift_n(c, mu, -1)
        W = W @ U_mat[c2id(c), mu].conj().swapaxes(-1,-2)

    # Left edge: T steps back in ν direction (conjugate transpose)
    for _ in range(T):
        c = shift_n(c, nu, -1)
        W = W @ U_mat[c2id(c), nu].conj().swapaxes(-1,-2)

    traces = np.real(np.trace(W, axis1=-2, axis2=-1)) / 2.0
    return float(np.mean(traces))


# ══════════════════════════════════════════════════════════════════════════════
# SU(3) PLAQUETTES
# ══════════════════════════════════════════════════════════════════════════════

def compute_plaquettes_su3(torus) -> np.ndarray:
    """
    SU(3) Wilson plaquette Re Tr[P_{x,μν}] / 3.

    Returns (N, 6) float array.
    Плакетный оператор SU(3).
    """
    if torus.gauge_group != 'SU3':
        raise ValueError("compute_plaquettes_su3 requires gauge_group='SU3'")

    L      = torus.L
    coords = torus.coords
    U      = _build_link_table_su3(torus)       # (N, 4, 3, 3)

    def shift(c, mu):
        c_new = c.copy(); c_new[:, mu] = (c_new[:, mu] + 1) % L; return c_new

    def c2id(c):
        return c[:,0]*L**3 + c[:,1]*L**2 + c[:,2]*L + c[:,3]

    plaq = np.zeros((torus.N, 6), dtype=np.float64)

    for pidx, (mu, nu) in enumerate(MU_NU_PAIRS):
        id_x    = c2id(coords)
        id_xmu  = c2id(shift(coords, mu))
        id_xnu  = c2id(shift(coords, nu))
        id_xmnu = c2id(shift(shift(coords, mu), nu))

        P = (U[id_x,    mu]
           @ U[id_xmu,  nu]
           @ U[id_xmnu, mu].conj().swapaxes(-1,-2)
           @ U[id_xnu,  nu].conj().swapaxes(-1,-2))

        # Re Tr[P] / 3  (d_r = 3 for fundamental representation)
        plaq[:, pidx] = np.real(np.trace(P, axis1=-2, axis2=-1)) / 3.0

    return plaq


def polyakov_loop_su3(torus, mu: int = 3) -> complex:
    """
    Polyakov loop in direction μ (default: temporal direction μ=3):
      P = (1/N_space) Σ_x Re Tr[ ∏_{t=0}^{L-1} U_{(x,t),μ} ] / 3

    Nonzero ⟺ deconfined phase (finite-temperature order parameter).
    Петля Полякова — порядковый параметр деконфайнмента.
    """
    if torus.gauge_group != 'SU3':
        raise ValueError("polyakov_loop_su3 requires gauge_group='SU3'")

    L      = torus.L
    coords = torus.coords
    U_mat  = _build_link_table_su3(torus)

    def c2id(c):
        return c[:,0]*L**3 + c[:,1]*L**2 + c[:,2]*L + c[:,3]

    # Spatial slice: fix μ coordinate to 0
    mask = coords[:, mu] == 0
    base = coords[mask]                         # (N_space, 4)
    N_sp = len(base)

    P_loop = np.eye(3, dtype=complex)[np.newaxis].repeat(N_sp, axis=0)

    for t in range(L):
        c    = base.copy()
        c[:, mu] = t
        ids  = c2id(c)
        P_loop = P_loop @ U_mat[ids, mu]

    traces = np.real(np.trace(P_loop, axis1=-2, axis2=-1)) / 3.0
    return complex(np.mean(traces))


# ══════════════════════════════════════════════════════════════════════════════
# UNIFIED HAMILTONIANS
# Унифицированные гамильтонианы
# ══════════════════════════════════════════════════════════════════════════════

def hamiltonian_u1(torus, J: float = 1.0) -> float:
    """H = -J Σ cos(W_{x,μν})"""
    return float(-J * np.sum(np.cos(compute_plaquettes_u1(torus))))


def hamiltonian_su2(torus, J: float = 1.0) -> float:
    """H = -J Σ Re Tr[P_{x,μν}] / 2"""
    return float(-J * np.sum(compute_plaquettes_su2(torus)))


def hamiltonian_su3(torus, J: float = 1.0) -> float:
    """H = -J Σ Re Tr[P_{x,μν}] / 3"""
    return float(-J * np.sum(compute_plaquettes_su3(torus)))


def hamiltonian(torus, J: float = 1.0) -> float:
    """Unified dispatcher: calls correct H based on torus.gauge_group."""
    g = torus.gauge_group
    if g == 'U1':  return hamiltonian_u1(torus, J)
    if g == 'SU2': return hamiltonian_su2(torus, J)
    if g == 'SU3': return hamiltonian_su3(torus, J)
    raise ValueError(f"Unknown gauge group: {g}")


def plaquette_energy_density(torus, J: float = 1.0) -> np.ndarray:
    """
    Local energy density ε(x) = -J Σ_{μ<ν} [plaquette observable at x].
    Works for U(1), SU(2), SU(3).
    Локальная плотность энергии — единый интерфейс для всех групп.
    """
    g = torus.gauge_group
    if g == 'U1':
        plaq = compute_plaquettes_u1(torus)
        return -J * np.sum(np.cos(plaq), axis=1)
    if g == 'SU2':
        plaq = compute_plaquettes_su2(torus)
        return -J * np.sum(plaq, axis=1)
    if g == 'SU3':
        plaq = compute_plaquettes_su3(torus)
        return -J * np.sum(plaq, axis=1)
    raise ValueError(f"Unknown gauge group: {g}")


# ══════════════════════════════════════════════════════════════════════════════
# TOPOLOGICAL CHARGE (U(1) only — well-defined for Abelian theory)
# Топологический заряд — только для U(1)
# ══════════════════════════════════════════════════════════════════════════════

def topological_charge_density(torus) -> np.ndarray:
    """
    q(x) = (1/16π²) Σ_{μνρσ} ε_{μνρσ} W_{x,μν} W_{x,ρσ}

    Discrete analogue of F∧F.  Q = Σ_x q(x) is integer for smooth configs.
    Дискретный аналог F∧F в U(1) теории.
    """
    if torus.gauge_group != 'U1':
        raise NotImplementedError(
            "Topological charge via F∧F is defined for U(1). "
            "For SU(2)/SU(3) use the clover definition (not yet implemented).")

    plaq     = compute_plaquettes_u1(torus)
    pair_idx = {(mu, nu): i for i, (mu, nu) in enumerate(MU_NU_PAIRS)}

    def levi_civita(p):
        inv = sum(1 for i in range(4) for j in range(i+1,4) if p[i] > p[j])
        return (-1)**inv

    q         = np.zeros(torus.N, dtype=np.float64)
    prefactor = 1.0 / (16.0 * np.pi**2)

    for perm in permutations(range(4)):
        mu, nu, rho, sigma = perm
        eps = levi_civita(perm)
        mn  = (min(mu,nu), max(mu,nu))
        rs  = (min(rho,sigma), max(rho,sigma))
        if mn not in pair_idx or rs not in pair_idx:
            continue
        sign_mn = 1 if mu < nu else -1
        sign_rs = 1 if rho < sigma else -1
        q += (eps * prefactor * sign_mn * sign_rs
              * plaq[:, pair_idx[mn]] * plaq[:, pair_idx[rs]])
    return q


def total_topological_charge(torus) -> float:
    """Q = Σ_x q(x)."""
    return float(np.sum(topological_charge_density(torus)))


# ══════════════════════════════════════════════════════════════════════════════
# UNIFIED ANALYSIS CLASS
# Унифицированный класс анализа
# ══════════════════════════════════════════════════════════════════════════════

class PlaquetteAnalysis:
    """
    Full plaquette analysis for a Z4Torus — supports U(1), SU(2), SU(3).

    Usage
    -----
    >>> torus_u1  = Z4Torus(L=6, gauge_group='U1')
    >>> torus_su2 = Z4Torus(L=4, gauge_group='SU2')
    >>> torus_su3 = Z4Torus(L=4, gauge_group='SU3')
    >>>
    >>> for t in [torus_u1, torus_su2, torus_su3]:
    ...     pa = PlaquetteAnalysis(t, J=1.0).run()
    ...     print(pa.summary())
    """

    def __init__(self, torus, J: float = 1.0):
        self.torus = torus
        self.J     = J
        self.energy_density: Optional[np.ndarray] = None
        self.q_density:      Optional[np.ndarray] = None
        self.H_total:        Optional[float]      = None
        self.Q_total:        Optional[float]      = None
        self.plaq_values:    Optional[np.ndarray] = None
        # SU(2) extras
        self.wilson_loops:   Optional[dict]       = None
        # SU(3) extras
        self.polyakov:       Optional[complex]    = None

    def run(self) -> 'PlaquetteAnalysis':
        g = self.torus.gauge_group

        self.energy_density = plaquette_energy_density(self.torus, self.J)
        self.H_total        = float(np.sum(self.energy_density))

        if g == 'U1':
            self.plaq_values = compute_plaquettes_u1(self.torus)
            self.q_density   = topological_charge_density(self.torus)
            self.Q_total     = float(np.sum(self.q_density))

        elif g == 'SU2':
            self.plaq_values = compute_plaquettes_su2(self.torus)
            # Wilson loops for R=T=1,2
            self.wilson_loops = {}
            for R in [1, 2]:
                for T in [1, 2]:
                    self.wilson_loops[(R,T)] = wilson_loop_su2(
                        self.torus, mu=0, nu=3, R=R, T=T)

        elif g == 'SU3':
            self.plaq_values = compute_plaquettes_su3(self.torus)
            self.polyakov    = polyakov_loop_su3(self.torus, mu=3)

        return self

    def summary(self) -> str:
        if self.H_total is None:
            return "PlaquetteAnalysis: not yet run. Call .run() first."
        g   = self.torus.gauge_group
        out = [
            f"PlaquetteAnalysis — Z4Torus(L={self.torus.L}, N={self.torus.N})",
            f"  Gauge group  : {g}",
            f"  J (coupling) : {self.J}",
            f"  H_total      : {self.H_total:.6f}",
            f"  H per site   : {self.H_total/self.torus.N:.6f}",
            f"  <P> (mean plaq): {self.plaq_values.mean():.6f}",
            f"  std(P)       : {self.plaq_values.std():.6f}",
        ]
        if g == 'U1' and self.Q_total is not None:
            out.append(f"  Q_total      : {self.Q_total:.4f}  (integer→smooth config)")
        if g == 'SU2' and self.wilson_loops:
            out.append("  Wilson loops <Re Tr W(R,T)>/2 :")
            for (R,T), val in sorted(self.wilson_loops.items()):
                out.append(f"    W({R},{T}) = {val:.6f}")
        if g == 'SU3' and self.polyakov is not None:
            out.append(f"  Polyakov loop: Re={self.polyakov.real:.6f}  "
                       f"(~0 confined, >0 deconfined)")
        return "\n".join(out)