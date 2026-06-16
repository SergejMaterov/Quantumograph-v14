"""
Quantumograph v14
=================
Finite quantum graph theory of spacetime on the Z⁴ torus.

Modules
-------
torus      : Z⁴ hypercubic torus with periodic boundary conditions
spectral   : Spectral dimension dₛ via heat kernel and random walks
plaquette  : Plaquette Hamiltonian and topological charge density
visualize  : Physically correct visualisation of all observables

Quick start
-----------
>>> from quantumograph import Z4Torus, SpectralAnalysis, PlaquetteAnalysis
>>> torus = Z4Torus(L=6, gauge_group='U1')
>>> sa = SpectralAnalysis(torus, mode='sparse', k=256).run()
>>> pa = PlaquetteAnalysis(torus, J=1.0).run()
>>> print(sa.summary())
>>> print(pa.summary())

Author  : © 2025–2026 Sergej Materov <sergejmaterov2@gmail.com>
License : CC BY-NC 4.0  (https://creativecommons.org/licenses/by-nc/4.0/)

© 2025–2026 Sergej Materov. Patent applications pending at Rospatent.
Commercial use prohibited without written permission.
"""

from .torus     import Z4Torus
from .spectral  import SpectralAnalysis, fss_extrapolate, check_4d_criterion
from .plaquette import (PlaquetteAnalysis,
                        hamiltonian, hamiltonian_u1, hamiltonian_su2, hamiltonian_su3,
                        compute_plaquettes_u1, compute_plaquettes_su2, compute_plaquettes_su3,
                        wilson_loop_su2, polyakov_loop_su3,
                        topological_charge_density, total_topological_charge,
                        su2_matrix, su3_matrix, su2_matrix_vec)
from .visualize import (plot_spectral_dimension, plot_full_analysis,
                        animate_torus, plot_fss)

__version__  = "0.1.1"
__author__   = "Sergej Materov"
__email__    = "sergejmaterov2@gmail.com"
__license__  = "CC BY-NC 4.0"

__all__ = [
    # Torus
    "Z4Torus",
    # Spectral
    "SpectralAnalysis", "fss_extrapolate", "check_4d_criterion",
    # Plaquette — analysis
    "PlaquetteAnalysis",
    # Plaquette — group elements
    "su2_matrix", "su3_matrix", "su2_matrix_vec",
    # Plaquette — observables
    "hamiltonian", "hamiltonian_u1", "hamiltonian_su2", "hamiltonian_su3",
    "compute_plaquettes_u1", "compute_plaquettes_su2", "compute_plaquettes_su3",
    "wilson_loop_su2", "polyakov_loop_su3",
    "topological_charge_density", "total_topological_charge",
    # Visualisation
    "plot_spectral_dimension", "plot_full_analysis",
    "animate_torus", "plot_fss",
]
