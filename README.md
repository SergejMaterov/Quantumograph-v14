# Quantumograph v14

**Finite quantum graph theory of spacetime**

> Supplementary code for the monograph *Quantumograph v14* by Sergej Materov.
> Patent applications pending at Rospatent (Russian Federation).

---

## English

### What is this?

This repository contains the reference implementation of **Quantumograph v14** —
a discrete quantum gravity model in which 4-dimensional spacetime emerges from
the spectral properties of a finite graph on the Z⁴ hypercubic torus.

The code is released as supplementary material to the monograph so that
results can be independently reproduced. The algorithms that are subject to
pending patent applications are **not included**.

### Theoretical foundations

| Concept | Implementation |
|---|---|
| Z⁴ torus T = Z_L⁴ with PBC | `quantumograph/torus.py` |
| Gauge link variables U(1)/SU(2)/SU(3) | `torus.py` → `_assign_link_phases()` |
| Laplacian spectrum {λᵢ} | `torus.laplacian_sparse()` |
| Heat kernel K(σ) = (1/N)Σ exp(-σλᵢ) | `spectral.heat_kernel_exact()` |
| Spectral dimension dₛ(σ) = -2 d ln K / d ln σ | `spectral.spectral_dimension()` |
| 4D criterion \|dₛ(∞) - 4.0\| < 0.1 | `spectral.check_4d_criterion()` |
| FSS extrapolation dₛ(L) → dₛ(∞) | `spectral.fss_extrapolate()` |
| **U(1)** plaquette W_{x,μν} = Σ θ on □ | `plaquette.compute_plaquettes_u1()` |
| **U(1)** Hamiltonian H = -J Σ cos(W) | `plaquette.hamiltonian_u1()` |
| **U(1)** topological charge q(x) = ε F∧F | `plaquette.topological_charge_density()` |
| **SU(2)** link U = exp(i α·σ/2), ZYZ Euler | `plaquette.su2_matrix()` |
| **SU(2)** plaquette Re Tr[P]/2 | `plaquette.compute_plaquettes_su2()` |
| **SU(2)** Wilson loop W(R,T) | `plaquette.wilson_loop_su2()` |
| **SU(3)** link U = exp(i α·λ/2), Gell-Mann | `plaquette.su3_matrix()` |
| **SU(3)** plaquette Re Tr[P]/3 | `plaquette.compute_plaquettes_su3()` |
| **SU(3)** Polyakov loop (deconfinement) | `plaquette.polyakov_loop_su3()` |

### Gauge groups

The torus supports three gauge groups selectable via `gauge_group` parameter:

| Group | Link variable | Parameters | d_r | Key observable |
|---|---|---|---|---|
| `U1` | U = exp(iθ) | 1 phase θ ∈ [0,2π) | 1 | Topological charge Q |
| `SU2` | U = exp(i α·σ/2) | 3 Euler angles | 2 | Wilson loop W(R,T) |
| `SU3` | U = exp(i α·λ/2) | 8 Gell-Mann params | 3 | Polyakov loop P |

Wilson action for all groups: H = -J Σ_{x,μ<ν} Re Tr[P_{x,μν}] / d_r

### Computation modes

| Mode | Method | Feasible N | Recommended for |
|---|---|---|---|
| `exact` | Full diagonalisation of L | N ≤ 3 000 | L ≤ 7, validation |
| `sparse` | ARPACK k leading eigenpairs | N ≤ 100 000 | L ≤ 10, desktop |
| `walk` | Monte Carlo random walks | unlimited | L ≥ 16, HPC |

### Quick start

```bash
pip install -r requirements.txt

# U(1), L=6 (~1296 nodes) — default
python run_demo.py

# SU(2), L=4 + Wilson loop table
python run_demo.py --gauge SU2 --L 4 --wilson

# SU(3), L=4 + Polyakov loop
python run_demo.py --gauge SU3 --L 4

# Compare all three gauge groups side by side
python run_demo.py --gauge all --L 4

# Finite-size scaling
python run_demo.py --fss --fss_L 3 4 5 6 7

# Animated t-slice scan
python run_demo.py --animate
```

### HPC usage

```bash
# Large torus, Monte Carlo walks (scales to any N)
python run_demo.py --L 20 --mode walk --n_walks 500000

# FSS on cluster: large L range, sparse spectrum
python run_demo.py --fss --fss_L 6 8 10 12 16 --mode sparse --k 1024

# SU(3) on large torus
python run_demo.py --gauge SU3 --L 8 --mode sparse --k 512
```

### Package API

```python
from quantumograph import (Z4Torus, SpectralAnalysis, PlaquetteAnalysis,
                            plot_full_analysis, wilson_loop_su2, polyakov_loop_su3)

# ── U(1) ──────────────────────────────────────────────────────────────
torus = Z4Torus(L=6, gauge_group='U1', seed=42)
sa = SpectralAnalysis(torus, mode='exact').run()
pa = PlaquetteAnalysis(torus, J=1.0).run()
print(sa.summary())   # dₛ(∞) = 4.000  ✓
print(pa.summary())   # H_total, Q_total

# ── SU(2) ─────────────────────────────────────────────────────────────
torus2 = Z4Torus(L=4, gauge_group='SU2', seed=42)
pa2 = PlaquetteAnalysis(torus2, J=1.0).run()
print(pa2.summary())  # H_total, Wilson loops W(1,1) … W(2,2)

# Direct Wilson loop
W = wilson_loop_su2(torus2, mu=0, nu=3, R=2, T=2)
print(f"W(2,2) = {W:.6f}")

# ── SU(3) ─────────────────────────────────────────────────────────────
torus3 = Z4Torus(L=4, gauge_group='SU3', seed=42)
pa3 = PlaquetteAnalysis(torus3, J=1.0).run()
print(pa3.summary())  # H_total, Polyakov loop

P = polyakov_loop_su3(torus3, mu=3)
print(f"Polyakov Re<P> = {P.real:.6f}")

# ── Visualise (works for all three groups) ────────────────────────────
sa3 = SpectralAnalysis(torus3, mode='exact').run()
plot_full_analysis(torus3, sa3, pa3)
```

### Repository structure

```
quantumograph/
├── quantumograph/
│   ├── torus.py        # Z⁴ torus, PBC, link variables U(1)/SU(2)/SU(3)
│   ├── spectral.py     # K(σ), dₛ(σ), FSS, Monte Carlo walks
│   ├── plaquette.py    # Plaquette Hamiltonians, Wilson/Polyakov loops
│   └── visualize.py    # 4-panel figures, animation, FSS plot
├── tests/
│   └── test_core.py    # 38 unit tests (pytest)
├── run_demo.py         # CLI: --gauge --L --mode --fss --wilson --animate
├── README.md
├── LICENSE             # CC BY-NC 4.0
├── setup.py
└── requirements.txt
```

### Running tests

```bash
pip install pytest
pytest tests/ -v
# Expected: 38 passed
```

### License
[![License: CC BY-NC-SA 4.0](https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png)](https://creativecommons.org/licenses/by-nc-sa/4.0/)

**CC BY-NC 4.0** — Creative Commons Attribution-NonCommercial 4.0 International.

You are free to use, share, and adapt this code for **non-commercial purposes**
provided you give appropriate credit. Commercial use requires written permission
from the author.

Patent applications for the core algorithms are pending at Rospatent (Russia).
After patent issuance the license will be updated to Apache 2.0.

© 2025–2026 Sergej Materov <sergejmaterov2@gmail.com>

---

## Русский

### Что это?

Репозиторий содержит эталонную реализацию **Quantumograph v14** —
модели дискретной квантовой гравитации, в которой 4-мерное пространство-время
возникает из спектральных свойств конечного графа на Z⁴-гиперкубическом торе.

Код опубликован как дополнительный материал к монографии для независимой
воспроизводимости результатов. Алгоритмы, являющиеся предметом патентных
заявок, **не включены**.

### Теоретические основы

| Концепция | Реализация |
|---|---|
| Z⁴-тор T = Z_L⁴ с периодическими граничными условиями | `quantumograph/torus.py` |
| Калибровочные линковые переменные U(1)/SU(2)/SU(3) | `torus.py` → `_assign_link_phases()` |
| Спектр лапласиана {λᵢ} | `torus.laplacian_sparse()` |
| Ядро теплопроводности K(σ) = (1/N)Σ exp(-σλᵢ) | `spectral.heat_kernel_exact()` |
| Спектральная размерность dₛ(σ) = -2 d ln K / d ln σ | `spectral.spectral_dimension()` |
| Критерий 4D: \|dₛ(∞) - 4.0\| < 0.1 | `spectral.check_4d_criterion()` |
| FSS-экстраполяция dₛ(L) → dₛ(∞) | `spectral.fss_extrapolate()` |
| **U(1)** плакета W_{x,μν} = Σ θ по □ | `plaquette.compute_plaquettes_u1()` |
| **U(1)** гамильтониан H = -J Σ cos(W) | `plaquette.hamiltonian_u1()` |
| **U(1)** топологический заряд q(x) = ε F∧F | `plaquette.topological_charge_density()` |
| **SU(2)** линк U = exp(i α·σ/2), углы Эйлера ZYZ | `plaquette.su2_matrix()` |
| **SU(2)** плакета Re Tr[P]/2 | `plaquette.compute_plaquettes_su2()` |
| **SU(2)** вильсоновский контур W(R,T) | `plaquette.wilson_loop_su2()` |
| **SU(3)** линк U = exp(i α·λ/2), матрицы Гелл-Манна | `plaquette.su3_matrix()` |
| **SU(3)** плакета Re Tr[P]/3 | `plaquette.compute_plaquettes_su3()` |
| **SU(3)** петля Полякова (деконфайнмент) | `plaquette.polyakov_loop_su3()` |

### Калибровочные группы

| Группа | Линк | Параметры | d_r | Ключевая наблюдаемая |
|---|---|---|---|---|
| `U1` | U = exp(iθ) | 1 фаза θ ∈ [0,2π) | 1 | Топологический заряд Q |
| `SU2` | U = exp(i α·σ/2) | 3 угла Эйлера | 2 | Вильсоновский контур W(R,T) |
| `SU3` | U = exp(i α·λ/2) | 8 параметров Гелл-Манна | 3 | Петля Полякова P |

Вильсоновское действие для всех групп: H = -J Σ_{x,μ<ν} Re Tr[P_{x,μν}] / d_r

### Режимы вычислений

| Режим | Метод | Максимальный N | Рекомендуется для |
|---|---|---|---|
| `exact` | Полная диагонализация L | N ≤ 3 000 | L ≤ 7, верификация |
| `sparse` | ARPACK k собственных значений | N ≤ 100 000 | L ≤ 10, десктоп |
| `walk` | Случайные блуждания Монте-Карло | без ограничений | L ≥ 16, HPC |

### Быстрый старт

```bash
pip install -r requirements.txt

# U(1), L=6 (~1296 узлов) — по умолчанию
python run_demo.py

# SU(2), L=4 + таблица вильсоновских контуров
python run_demo.py --gauge SU2 --L 4 --wilson

# SU(3), L=4 + петля Полякова
python run_demo.py --gauge SU3 --L 4

# Сравнение всех трёх групп
python run_demo.py --gauge all --L 4

# FSS-экстраполяция
python run_demo.py --fss --fss_L 3 4 5 6 7

# Анимация временных срезов
python run_demo.py --animate
```

### Структура репозитория

```
quantumograph/
├── quantumograph/
│   ├── torus.py        # Z⁴ тор, ПГУ, линки U(1)/SU(2)/SU(3)
│   ├── spectral.py     # K(σ), dₛ(σ), FSS, блуждания МК
│   ├── plaquette.py    # Гамильтонианы, контуры Вильсона, петля Полякова
│   └── visualize.py    # 4-панельные графики, анимация, FSS-plot
├── tests/
│   └── test_core.py    # 38 юнит-тестов (pytest)
├── run_demo.py         # CLI: --gauge --L --mode --fss --wilson --animate
├── README.md
├── LICENSE             # CC BY-NC 4.0
├── setup.py
└── requirements.txt
```

### Запуск тестов

```bash
pip install pytest
pytest tests/ -v
# Ожидается: 38 passed
```

### Лицензия
[![License: CC BY-NC-SA 4.0](https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png)](https://creativecommons.org/licenses/by-nc-sa/4.0/)

**CC BY-NC 4.0** — Creative Commons Attribution-NonCommercial 4.0.

Использование, распространение и адаптация разрешены только в **некоммерческих**
целях при указании авторства. Коммерческое использование требует письменного
разрешения автора.

Патентные заявки на основные алгоритмы находятся на рассмотрении в Роспатенте.
После выдачи патентов лицензия будет обновлена до Apache 2.0.

© 2025–2026 Сергей Матеров <sergejmaterov2@gmail.com>
