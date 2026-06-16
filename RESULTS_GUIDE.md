# Quantumograph v14 — Guide to Results / Руководство по интерпретации результатов

> This document explains what every number, plot, and observable means physically.
> Intended for readers of the monograph who want to understand the code output.
>
> Этот документ объясняет физический смысл каждого числа, графика и наблюдаемой.
> Предназначен для читателей монографии, желающих понять вывод кода.

---

## English

---

### 1. Torus construction output

```
[1/3] Torus built in 0.04s — 1296 nodes, 5184 edges
```

| Value | Meaning |
|---|---|
| `1296 nodes` | L⁴ = 6⁴ lattice sites — the discrete spacetime points |
| `5184 edges` | 4 × L⁴ = 4 × 1296 nearest-neighbour links (each node connects to 8 neighbours = 2 per direction × 4 directions, giving 4L⁴ undirected edges) |

**What this means physically:** You are looking at a finite piece of discrete spacetime. Each node is a "cubit" — an elementary quantum of spacetime. The edges are the quantum entanglement links between adjacent spacetime points. The torus topology (periodic boundary conditions) eliminates boundary effects, making the geometry translation-invariant — the same at every point, just like the continuum.

---

### 2. Spectral analysis output

```
SpectralAnalysis — Z4Torus(L=6, N=1296)
  Mode        : exact
  dₛ(∞)       : 4.0000
  Criterion   : |dₛ(∞) - 4.0| < 0.1  →  ✓ PASSES
  dₛ range    : [0.020, 4.845]
```

#### 2.1 dₛ(∞) — spectral dimension

This is **the most important number in the entire output**.

The spectral dimension dₛ measures how many effective dimensions a diffusing particle "sees" as it performs a random walk on the graph. It is computed from the heat kernel:

```
K(σ) = (1/N) Σᵢ exp(-σ λᵢ)    [return probability after diffusion time σ]
dₛ(σ) = -2 · d ln K(σ) / d ln σ
```

| dₛ value | Physical interpretation |
|---|---|
| dₛ → 0 | UV regime: diffusion has not yet explored the lattice structure |
| dₛ → 4 | **Emergent 4D spacetime** — the physical regime of the model |
| dₛ > 4 | IR artefact: finite torus size causes over-counting of return paths |
| dₛ(∞) ≈ 4.0 | The model correctly reproduces 4-dimensional spacetime at large scales |

**The criterion** `|dₛ(∞) - 4.0| < 0.1` is the quantitative test of whether this particular Z⁴ torus configuration exhibits 4-dimensional emergent geometry. A ✓ means the model passes.

#### 2.2 dₛ range

`dₛ range : [0.020, 4.845]` — the full range of dₛ(σ) as σ sweeps from UV to IR.

- The **minimum** (~0.02) is the UV value: at very short diffusion times, only a few nearest neighbours are explored, giving an apparent dimension close to 0.
- The **maximum** (~4.8) is the IR peak before finite-size effects dominate.
- The physical 4D value is found at the **crossing point** where dₛ(σ*) = 4.0, typically at σ* ≈ 0.39–0.40 in lattice units.

#### 2.3 Computation modes and their accuracy

| Mode | What it computes | When to trust it |
|---|---|---|
| `exact` | All N eigenvalues of L | Always accurate, use for L ≤ 7 |
| `sparse` | k smallest eigenvalues | Accurate if k ≥ N/4; increase k if dₛ looks wrong |
| `walk` | MC estimate of K(σ) | Accurate for n_walks ≥ 100 000; noisy at large σ |

---

### 3. Plaquette analysis output — U(1)

```
PlaquetteAnalysis — Z4Torus(L=6, N=1296)
  Gauge group  : U1
  J (coupling) : 1.0
  H_total      : 11.750
  H per site   : 0.009
  <P> (mean plaq): 0.000
  std(P)       : 3.609
  Q_total      : 29.703  (integer for smooth config)
```

#### 3.1 H_total and H per site

`H = -J Σ_{x,μ<ν} cos(W_{x,μν})`

| Value | Meaning |
|---|---|
| `H_total > 0` | The gauge field is in a **disordered (hot) phase** — random link phases give ⟨cos W⟩ ≈ 0, so H ≈ 0 per site |
| `H_total → -J·6·N` | The gauge field is in an **ordered (cold) phase** — all plaquettes are flat (W≈0), maximum coherence |
| `H per site ≈ 0.009` | Each lattice site contributes almost zero energy on average — consistent with a random (infinite temperature) initial configuration |

The random seed generates completely random link phases θ ∈ [0, 2π), giving cos(W) averaging to zero. This is the **quenched disordered** starting point, not the physical ground state. To find the ground state, one would run a Monte Carlo simulation (not included — subject to pending patents).

#### 3.2 ⟨P⟩ — mean plaquette

`<P> = 0.000` for U(1): the mean of W_{x,μν} over all plaquettes and nodes.

Since W is a sum of four random phases, its mean is zero by symmetry. The **std(P) ≈ 3.6** reflects the spread of plaquette values across the torus — comparable to 2π, as expected for uniformly random link phases.

#### 3.3 Q_total — topological charge

`Q_total : 29.703  (integer for smooth config)`

The topological charge is:
```
Q = Σ_x q(x) = (1/16π²) Σ_{x,μνρσ} ε_{μνρσ} W_{x,μν} W_{x,ρσ}
```

This is the discrete analogue of the instanton number ∫ F∧F / 16π² in continuum gauge theory.

| Q_total value | Meaning |
|---|---|
| Q ≈ integer | **Smooth (topologically non-trivial) configuration** — a lattice instanton |
| Q non-integer (e.g. 29.7) | **Rough (random) configuration** — the lattice discretisation is too coarse or the config is too disordered for topological quantisation to hold |
| Q = 0 | Topologically trivial vacuum |
| Q = 1, 2, … | Instanton, di-instanton, … sectors |

In the current output Q ≈ 29.7 is non-integer because the initial configuration is completely random (infinite temperature). Integer Q emerges after thermalization via Monte Carlo or gradient flow.

---

### 4. Plaquette analysis output — SU(2)

```
PlaquetteAnalysis — Z4Torus(L=4, N=256)
  Gauge group  : SU2
  H_total      : -6.306
  H per site   : -0.025
  <P> (mean plaq): 0.004
  Wilson loops <Re Tr W(R,T)>/2 :
    W(1,1) = -0.016
    W(1,2) = -0.051
    W(2,1) =  0.021
    W(2,2) =  0.086
```

#### 4.1 H_total for SU(2)

`H = -J Σ Re Tr[P_{x,μν}] / 2`

For SU(2), Re Tr[U]/2 ∈ [-1, +1], so H ∈ [-6JN, 0].

- `H = -6JN`: **ordered phase** (all links U = I, complete coherence)
- `H ≈ 0`: **disordered phase** (random SU(2) matrices, ⟨Re Tr U⟩ = 0)
- `H per site ≈ -0.025`: slightly negative — just below the disordered value, consistent with random SU(2) matrices where ⟨Re Tr U⟩ ≈ 0 (the Haar measure average over SU(2) gives exactly 0)

#### 4.2 Wilson loops W(R,T)

The Wilson loop measures the **phase accumulated by a quark-antiquark pair** separated by distance R and propagating for time T:

```
W(R,T) = ⟨Re Tr[ ∏_{path R×T} U ] / 2⟩
```

Two possible behaviours determine the phase of the gauge theory:

| Behaviour | Formula | Physical meaning |
|---|---|---|
| **Area law** | W(R,T) ~ exp(-σ · R · T) | **Confinement** — quarks are confined, string tension σ > 0 |
| **Perimeter law** | W(R,T) ~ exp(-μ · (2R+2T)) | **Deconfinement** — quarks are free, no string |

For a random SU(2) configuration (infinite temperature), W(R,T) oscillates near zero for all R,T — neither law applies cleanly. The ratio `ln W(2,2) / ln W(1,1)` distinguishes the two laws: area law predicts ratio → 4, perimeter law → 2.

---

### 5. Plaquette analysis output — SU(3)

```
PlaquetteAnalysis — Z4Torus(L=4, N=256)
  Gauge group  : SU3
  H_total      : -26.598
  H per site   : -0.104
  <P> (mean plaq): 0.017
  Polyakov loop: Re=0.224  (~0 confined, >0 deconfined)
```

#### 5.1 H_total for SU(3)

`H = -J Σ Re Tr[P_{x,μν}] / 3`

For SU(3), Re Tr[U]/3 ∈ [-1/2, +1] (the SU(3) constraint gives a tighter lower bound than SU(2)).

`H per site ≈ -0.104` — more negative than SU(2) because SU(3) has more internal degrees of freedom (8 generators vs 3), so random SU(3) matrices contribute a slightly non-zero mean plaquette.

#### 5.2 Polyakov loop — the deconfinement order parameter

```
P = (1/N_space) Σ_x Re Tr[ ∏_{t=0}^{L-1} U_{(x,t),3} ] / 3
```

The Polyakov loop winds around the entire temporal extent of the torus. It is the **order parameter for the confinement–deconfinement phase transition**:

| Re⟨P⟩ | Phase | Physical meaning |
|---|---|---|
| Re⟨P⟩ ≈ 0 | **Confined** | Quarks and gluons are bound inside hadrons; Z₃ centre symmetry is unbroken |
| Re⟨P⟩ > 0 | **Deconfined** (quark-gluon plasma) | Quarks move freely; Z₃ centre symmetry is spontaneously broken |

`Re=0.224` in the current output: a random SU(3) configuration gives a nonzero Polyakov loop because finite-L statistical fluctuations are large. True confinement (Re⟨P⟩ → 0) requires averaging over many thermalized configurations at low temperature (large L in the temporal direction).

---

### 6. Reading the four-panel visualisation

#### Panel [0,0] — Spectral dimension curve dₛ(σ)

- **Horizontal axis** (log scale): diffusion time σ, from UV (left) to IR (right)
- **Vertical axis**: dₛ(σ), the effective spacetime dimension seen at scale σ
- **Blue curve**: dₛ(σ) computed from the Laplacian spectrum
- **Dashed orange line**: dₛ = 4.0 target
- **Yellow band**: the 4D criterion region |dₛ - 4| < 0.1
- **Green/red dotted line**: dₛ(∞) estimate with ✓/✗ status
- **What to look for**: the curve should cross dₛ = 4 at some σ* ≈ 0.4 and be roughly flat there — this is the emergent 4D regime

#### Panel [0,1] — Energy density ε(x) on time slice t

- Each point is a lattice node at fixed time coordinate t, shown in the (x₀, x₁) plane
- Colour encodes the local plaquette energy ε(x) = -J Σ_{μ<ν} cos(W_{x,μν})
- **Dark (inferno colormap)**: high energy — frustrated, strongly fluctuating gauge field
- **Bright**: low energy — locally ordered, smooth gauge field
- **What to look for**: for a random config, colours are uniformly distributed; after thermalization, spatial correlations appear

#### Panel [1,0] — Group-specific observable

**U(1): Topological charge density q(x)**
- Colour encodes q(x) = (1/16π²) ε_{μνρσ} W_{μν} W_{ρσ}
- Positive (red): local instanton density
- Negative (blue): local anti-instanton density
- Regions where |q(x)| is large are topologically active

**SU(2): Wilson loop table W(R,T)**
- Grid: rows = R (spatial separation), columns = T (temporal extent)
- Green cell: W > 0 (constructive interference)
- Red cell: W < 0 (destructive, typical of random hot configs)
- Pattern of decrease with R·T reveals the string tension

**SU(3): Plaquette value distribution**
- Histogram of Re Tr[P_{x,μν}]/3 over all plaquettes
- Centred near 0 for random configs; shifts right (toward +1) as the system cools
- Polyakov loop value shown in title

#### Panel [1,1] — Torus graph coloured by ε(x)

- Nodes: lattice sites in the (x₀, x₁) slice at fixed (x₂, t)
- Edges: nearest-neighbour links within the slice
- Colour (plasma): local energy density ε(x)
- **What to look for**: spatial correlations in ε(x) indicate the formation of gauge field structures (flux tubes, vortices, instantons)

---

### 7. FSS extrapolation output

```
L=4  dₛ = 3.955
L=5  dₛ = 4.000
L=6  dₛ = 4.000
L=7  dₛ = 4.000
FSS extrapolation dₛ(∞) = 4.068  ✓
```

Finite-size scaling (FSS) extrapolates the spectral dimension measured on finite toruses L=3,4,…,7 to the infinite-volume limit L→∞ using the fit:

```
dₛ(L) = dₛ(∞) + a/L + b/L²
```

- **Why necessary**: any finite lattice has boundary effects (even with PBC, the finite spectral gap cuts off IR modes). The true emergent dimension is only visible at L→∞.
- **dₛ(∞) ≈ 4.07**: consistent with 4D spacetime within FSS uncertainties.
- **The FSS plot**: shows dₛ(L) vs 1/L with the fitted curve; extrapolation to 1/L=0 gives dₛ(∞). The green dot at 1/L=0 is the final result.

---

### 8. Summary: what a "passing" result looks like

A fully consistent Quantumograph v14 result satisfies all of the following:

| Observable | Expected value | Status |
|---|---|---|
| dₛ(∞) | 4.0 ± 0.1 | ✓ confirmed for L ≥ 5 |
| FSS dₛ(∞) | 4.0 ± 0.2 | ✓ confirmed for L=3..7 |
| σ* (crossing point) | ≈ 0.39–0.40 | ✓ stable across L |
| Node degree | exactly 8 | ✓ (PBC enforced) |
| Plaquette range U(1) | W ∈ (-8π, 8π) | ✓ |
| Plaquette range SU(2) | Re Tr/2 ∈ [-1, 1] | ✓ |
| Plaquette range SU(3) | Re Tr/3 ∈ [-0.5, 1] | ✓ |
| Q_total (smooth config) | near-integer | ✓ after gradient flow |

---

---

## Русский

---

### 1. Вывод при построении тора

```
[1/3] Torus built in 0.04s — 1296 nodes, 5184 edges
```

| Значение | Смысл |
|---|---|
| `1296 nodes` | L⁴ = 6⁴ узлов решётки — дискретные точки пространства-времени |
| `5184 edges` | 4 × L⁴ = 4 × 1296 рёбер ближайших соседей (каждый узел соединён с 8 соседями — по 2 в каждом из 4 направлений, итого 4L⁴ неориентированных рёбер) |

**Физический смысл:** перед вами конечный фрагмент дискретного пространства-времени. Каждый узел — «кубит»: элементарный квант пространства-времени. Рёбра — связи квантовой запутанности между соседними точками пространства-времени. Топология тора (периодические граничные условия) устраняет граничные эффекты, делая геометрию трансляционно-инвариантной и одинаковой во всех точках, как в континууме.

---

### 2. Вывод спектрального анализа

```
SpectralAnalysis — Z4Torus(L=6, N=1296)
  Mode        : exact
  dₛ(∞)       : 4.0000
  Criterion   : |dₛ(∞) - 4.0| < 0.1  →  ✓ PASSES
  dₛ range    : [0.020, 4.845]
```

#### 2.1 dₛ(∞) — спектральная размерность

Это **самое важное число во всём выводе**.

Спектральная размерность dₛ измеряет, сколько эффективных измерений «видит» диффундирующая частица при случайном блуждании по графу. Вычисляется через ядро теплопроводности:

```
K(σ) = (1/N) Σᵢ exp(-σ λᵢ)    [вероятность возврата после диффузии σ]
dₛ(σ) = -2 · d ln K(σ) / d ln σ
```

| Значение dₛ | Физическая интерпретация |
|---|---|
| dₛ → 0 | УФ-режим: диффузия ещё не исследовала структуру решётки |
| dₛ → 4 | **Эмерджентное 4D пространство-время** — физический режим модели |
| dₛ > 4 | ИК-артефакт: конечный размер тора вызывает пересчёт путей возврата |
| dₛ(∞) ≈ 4.0 | Модель корректно воспроизводит 4-мерное пространство-время на больших масштабах |

**Критерий** `|dₛ(∞) - 4.0| < 0.1` — количественная проверка того, что данная конфигурация Z⁴-тора демонстрирует  4-мерную эмерджентную геометрию. ✓ означает прохождение теста.

#### 2.2 Диапазон dₛ

`dₛ range : [0.020, 4.845]` — полный диапазон dₛ(σ) при сканировании от УФ до ИК.

- **Минимум** (~0.02): УФ-значение — при малых временах диффузии исследуются только ближайшие соседи.
- **Максимум** (~4.8): ИК-пик перед тем, как конечный размер тора начинает доминировать.
- Физическое значение 4D находится в **точке пересечения** dₛ(σ*) = 4.0, типично при σ* ≈ 0.39–0.40 в единицах решётки.

#### 2.3 Режимы вычисления и их точность

| Режим | Что вычисляется | Когда доверять |
|---|---|---|
| `exact` | Все N собственных значений L | Всегда точно, используйте при L ≤ 7 |
| `sparse` | k наименьших собственных значений | Точно при k ≥ N/4; увеличьте k если dₛ выглядит неверно |
| `walk` | МК-оценка K(σ) | Точно при n_walks ≥ 100 000; шумно при больших σ |

---

### 3. Вывод анализа плакет — U(1)

```
PlaquetteAnalysis — Z4Torus(L=6, N=1296)
  Gauge group  : U1
  J (coupling) : 1.0
  H_total      : 11.750
  H per site   : 0.009
  <P> (mean plaq): 0.000
  std(P)       : 3.609
  Q_total      : 29.703  (integer for smooth config)
```

#### 3.1 H_total и H на узел

`H = -J Σ_{x,μ<ν} cos(W_{x,μν})`

| Значение | Смысл |
|---|---|
| `H_total > 0` | Калибровочное поле в **разупорядоченной (горячей) фазе** — случайные фазы дают ⟨cos W⟩ ≈ 0 |
| `H_total → -J·6·N` | Калибровочное поле в **упорядоченной (холодной) фазе** — все плакеты плоские (W≈0) |
| `H per site ≈ 0.009` | Каждый узел вносит почти нулевую энергию, что полностью соответствует случайной (бесконечно горячей) начальной конфигурации |

Случайный начальный seed генерирует полностью случайные фазы θ ∈ [0, 2π), что даёт cos(W) со средним нулём. Это **закалённая разупорядоченная** начальная точка, не физическое основное состояние. Для нахождения основного состояния необходима симуляция методом Монте-Карло.

#### 3.2 ⟨P⟩ — средняя плакета

`<P> = 0.000` для U(1): среднее W_{x,μν} по всем плакетам и узлам равно нулю по симметрии.

`std(P) ≈ 3.6` отражает разброс значений плакет по тору — сравнимо с 2π, как и ожидается для равномерно случайных фаз.

#### 3.3 Q_total — топологический заряд

`Q_total : 29.703  (integer for smooth config)`

Топологический заряд:
```
Q = Σ_x q(x) = (1/16π²) Σ_{x,μνρσ} ε_{μνρσ} W_{x,μν} W_{x,ρσ}
```

Это дискретный аналог числа инстантонов ∫ F∧F / 16π² в континуальной теории.

| Значение Q_total | Смысл |
|---|---|
| Q ≈ целое число | **Гладкая (топологически нетривиальная) конфигурация** — решёточный инстантон |
| Q нецелое (напр. 29.7) | **Грубая (случайная) конфигурация** — дискретизация слишком груба для квантования топологии |
| Q = 0 | Топологически тривиальный вакуум |
| Q = 1, 2, … | Инстантон, ди-инстантон, … |

В текущем выводе Q ≈ 29.7 — нецелое, поскольку начальная конфигурация полностью случайна. Целочисленный Q возникает после термализации (Монте-Карло) или применения градиентного потока.

---

### 4. Вывод анализа плакет — SU(2)

```
PlaquetteAnalysis — Z4Torus(L=4, N=256)
  Gauge group  : SU2
  H_total      : -6.306
  H per site   : -0.025
  <P>          : 0.004
  Wilson loops <Re Tr W(R,T)>/2 :
    W(1,1) = -0.016
    W(1,2) = -0.051
    W(2,1) =  0.021
    W(2,2) =  0.086
```

#### 4.1 Вильсоновские контуры W(R,T)

Вильсоновский контур измеряет **фазу, накопленную парой кварк-антикварк**, разделённой на расстояние R и распространяющейся во времени T:

```
W(R,T) = ⟨Re Tr[ ∏_{путь R×T} U ] / 2⟩
```

| Поведение | Формула | Физический смысл |
|---|---|---|
| **Закон площади** | W(R,T) ~ exp(-σ · R · T) | **Конфайнмент** — кварки связаны, натяжение струны σ > 0 |
| **Закон периметра** | W(R,T) ~ exp(-μ · (2R+2T)) | **Деконфайнмент** — кварки свободны |

Отношение `ln W(2,2) / ln W(1,1)`: закон площади даёт → 4, закон периметра → 2.

Для случайной конфигурации W(R,T) флуктуирует около нуля — ни один закон не применим напрямую. Истинное поведение проявляется после термализации.

---

### 5. Вывод анализа плакет — SU(3)

```
PlaquetteAnalysis — Z4Torus(L=4, N=256)
  Gauge group  : SU3
  H_total      : -26.598
  H per site   : -0.104
  <P>          : 0.017
  Polyakov loop: Re=0.224  (~0 confined, >0 deconfined)
```

#### 5.1 Петля Полякова — параметр порядка деконфайнмента

```
P = (1/N_пр) Σ_x Re Tr[ ∏_{t=0}^{L-1} U_{(x,t),3} ] / 3
```

Петля Полякова обматывает тор по всей временно́й протяжённости. Она является **параметром порядка фазового перехода конфайнмент–деконфайнмент**:

| Re⟨P⟩ | Фаза | Физический смысл |
|---|---|---|
| Re⟨P⟩ ≈ 0 | **Конфайнмент** | Кварки и глюоны связаны в адроны; Z₃-центровая симметрия сохранена |
| Re⟨P⟩ > 0 | **Деконфайнмент** (кварк-глюонная плазма) | Кварки движутся свободно; Z₃-симметрия спонтанно нарушена |

`Re=0.224` в текущем выводе: случайная SU(3)-конфигурация даёт ненулевую петлю Полякова из-за больших статистических флуктуаций при конечном L. Истинный конфайнмент (Re⟨P⟩ → 0) требует усреднения по многим термализованным конфигурациям.

---

### 6. Чтение четырёхпанельной визуализации

#### Панель [0,0] — Кривая спектральной размерности dₛ(σ)

- **Горизонтальная ось** (логарифмическая): время диффузии σ, от УФ (слева) до ИК (справа)
- **Вертикальная ось**: dₛ(σ) — эффективная размерность пространства-времени на масштабе σ
- **Синяя кривая**: dₛ(σ), вычисленная из спектра лапласиана
- **Пунктирная оранжевая**: цель dₛ = 4.0
- **Жёлтая полоса**: область критерия 4D — |dₛ - 4| < 0.1
- **Зелёная/красная пунктирная**: оценка dₛ(∞) со статусом ✓/✗
- **Что искать**: кривая должна пересекать dₛ = 4 при σ* ≈ 0.4 и быть там приблизительно плоской — это и есть эмерджентный 4D режим

#### Панель [0,1] — Плотность энергии ε(x) на временно́м срезе t

- Каждая точка — узел решётки при фиксированной временно́й координате t, в плоскости (x₀, x₁)
- Цвет кодирует локальную плакетную энергию ε(x) = -J Σ_{μ<ν} cos(W_{x,μν})
- **Тёмный (colormap inferno)**: высокая энергия — расстроенное, сильно флуктуирующее поле
- **Яркий**: низкая энергия — локально упорядоченное, гладкое поле

#### Панель [1,0] — Группо-специфическая наблюдаемая

**U(1): Плотность топологического заряда q(x)**
- Красное: локальная плотность инстантонов
- Синее: локальная плотность антиинстантонов
- Области с большим |q(x)| — топологически активные

**SU(2): Таблица вильсоновских контуров W(R,T)**
- Строки = R (пространственное разделение), столбцы = T (временна́я протяжённость)
- Зелёная клетка: W > 0 (конструктивная интерференция)
- Красная клетка: W < 0 (деструктивная, типично для горячих конфигураций)

**SU(3): Распределение значений плакет**
- Гистограмма Re Tr[P_{x,μν}]/3 по всем плакетам
- Центрирована около 0 для случайных конфигураций; смещается вправо (к +1) при охлаждении системы
- Значение петли Полякова указано в заголовке

#### Панель [1,1] — Граф тора, окрашенный по ε(x)

- Узлы: узлы решётки в срезе (x₀, x₁) при фиксированных (x₂, t)
- Рёбра: связи ближайших соседей внутри среза
- Цвет (plasma): локальная плотность энергии ε(x)
- **Что искать**: пространственные корреляции в ε(x) указывают на формирование структур калибровочного поля (флюксовые трубки, вихри, инстантоны)

---

### 7. FSS-экстраполяция

```
L=4  dₛ = 3.955
L=5  dₛ = 4.000
L=6  dₛ = 4.000
L=7  dₛ = 4.000
FSS dₛ(∞) = 4.068  ✓
```

FSS экстраполирует спектральную размерность, измеренную на конечных торах L=3,4,...,7, в бесконечный объём L→∞ с помощью подгонки:

```
dₛ(L) = dₛ(∞) + a/L + b/L²
```

- **Зачем нужно**: любая конечная решётка имеет эффекты конечного размера — конечный спектральный зазор обрезает ИК-моды. Истинная эмерджентная размерность видна только при L→∞.
- **dₛ(∞) ≈ 4.07**: совместимо с 4D пространством-временем в пределах погрешностей FSS.
- **FSS-график**: показывает dₛ(L) как функцию 1/L с подогнанной кривой; экстраполяция к 1/L=0 даёт dₛ(∞). Зелёная точка при 1/L=0 — итоговый результат.

---

### 8. Итог: как выглядит успешный результат

Полностью согласованный результат Quantumograph v14 удовлетворяет всем следующим условиям:

| Наблюдаемая | Ожидаемое значение | Статус |
|---|---|---|
| dₛ(∞) | 4.0 ± 0.1 | ✓ подтверждено для L ≥ 5 |
| FSS dₛ(∞) | 4.0 ± 0.2 | ✓ подтверждено для L=3..7 |
| σ* (точка пересечения) | ≈ 0.39–0.40 | ✓ стабильно по L |
| Степень узла | ровно 8 | ✓ (ПГУ выполнены) |
| Диапазон плакет U(1) | W ∈ (-8π, 8π) | ✓ |
| Диапазон плакет SU(2) | Re Tr/2 ∈ [-1, 1] | ✓ |
| Диапазон плакет SU(3) | Re Tr/3 ∈ [-0.5, 1] | ✓ |
| Q_total (гладкая конф.) | близко к целому | ✓ после градиентного потока |

© 2025–2026 Sergej Materov <sergejmaterov2@gmail.com> — CC BY-NC 4.0
