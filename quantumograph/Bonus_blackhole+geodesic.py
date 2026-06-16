"""
Quantumograph v14 — Визуализация искривления метрики пространства-времени
Оптимизировано для 24-ядерного CPU:
  • Рёбра: LineCollection (1 вызов вместо 178k)
  • Проекция: полностью векторизована через numpy
  • Фон метрики: кешируется, не пересчитывается каждый кадр
  • Кеш кадров: активные узлы + позиции считаются заранее
Authors : Sergej Materov <sergejmaterov2@gmail.com>
License : CC BY-NC 4.0
"""

import os
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.widgets as widgets
from matplotlib.animation import FuncAnimation
from matplotlib.gridspec import GridSpec
from matplotlib.collections import LineCollection
from concurrent.futures import ProcessPoolExecutor, as_completed

# ─────────────────────────────────────────────
# 1. ПАРАМЕТРЫ
# ─────────────────────────────────────────────
NUM_CUBITS   = 2000
SEED         = 10
BH_MASS      = 1.8
EDGE_THRESH  = 0.6
TIME_WINDOW  = 1.8
ANIM_SPEED   = 60
RS           = 0.3
SHOW_GRID    = True
SHOW_ARROWS  = True
SHOW_CONTOUR = True

GRAPH_WORKERS = max(1, (os.cpu_count() or 2) - 1)

np.random.seed(SEED)
X = np.random.uniform(-1.5, 1.5, NUM_CUBITS)
Y = np.random.uniform(-1.5, 1.5, NUM_CUBITS)
Z = np.random.uniform(-1.5, 1.5, NUM_CUBITS)
T = np.random.uniform(0,   10,   NUM_CUBITS)
P = np.column_stack([X, Y, Z])
D = np.linalg.norm(P, axis=1)

timeline = np.linspace(0, 10, 120)

Quantumograph  = None
frame_cache    = None   # list[dict] — активные узлы + позиции + данные рёбер
bg_cache       = {}     # (bh_mass, rs) -> list of Artists для фона

# ─────────────────────────────────────────────
# 2. ВЕКТОРИЗОВАННАЯ ГРАВИТАЦИОННАЯ ПРОЕКЦИЯ
# ─────────────────────────────────────────────
def schwarzschild_project_vec(coords, bh_mass, rs):
    """
    coords: (N,3) float array
    Возвращает (N,2) float array — спроецированные координаты.
    Всё через numpy, без Python-цикла.
    """
    r      = np.linalg.norm(coords, axis=1)                    # (N,)
    r_safe = np.maximum(r, 1e-6)

    phi   = np.arctan2(coords[:, 1], coords[:, 0])             # (N,)
    theta = np.arccos(np.clip(coords[:, 2] / r_safe, -1, 1))   # (N,)

    rs_eff = rs * bh_mass
    shrink = np.sqrt(np.maximum(0.0, 1.0 - rs_eff / (r_safe + 1e-4)))
    fall   = 1.0 / (1.0 + bh_mass * 0.4 / (r_safe + 0.05)**2)
    r_vis  = r_safe * shrink * fall

    vx = r_vis * np.sin(theta) * np.cos(phi)
    vy = r_vis * np.sin(theta) * np.sin(phi)
    return np.column_stack([vx, vy])                            # (N,2)


# ─────────────────────────────────────────────
# 3. ПОСТРОЙКА ГРАФА (параллельная, без изменений)
# ─────────────────────────────────────────────
def _build_chunk(args):
    idx_chunk, bh_mass, edge_thresh, P_, T_, D_, N = args
    edges = []
    for i in idx_chunk:
        j_idx  = np.arange(i + 1, N)
        if j_idx.size == 0:
            continue
        diff_p = P_[j_idx] - P_[i]
        d4     = np.sqrt(np.sum(diff_p**2, axis=1) + (T_[j_idx] - T_[i])**2)
        grav   = 1.0 + (bh_mass / (D_[i] + 0.1)) * (bh_mass / (D_[j_idx] + 0.1))
        eff    = d4 / grav
        mask   = eff < edge_thresh
        if np.any(mask):
            jj = j_idx[mask]; ew = eff[mask]
            edges.extend((int(i), int(j), float(1.0 / (e + 0.01)))
                         for j, e in zip(jj, ew))
    return edges


def build_graph(bh_mass, edge_thresh):
    G = nx.Graph()
    for i in range(NUM_CUBITS):
        G.add_node(i, pos_4d=(X[i], Y[i], Z[i], T[i]))
    chunks = np.array_split(np.arange(NUM_CUBITS), max(1, GRAPH_WORKERS * 2))
    print(f"Строим граф в {GRAPH_WORKERS} процессах…")
    all_edges = []
    with ProcessPoolExecutor(max_workers=GRAPH_WORKERS) as ex:
        futs = [ex.submit(_build_chunk, (c, bh_mass, edge_thresh, P, T, D, NUM_CUBITS))
                for c in chunks if len(c) > 0]
        for f in as_completed(futs):
            all_edges.extend(f.result())
    G.add_weighted_edges_from(all_edges)
    print(f"Граф: {G.number_of_nodes()} узлов, {G.number_of_edges()} рёбер")
    return G


# ─────────────────────────────────────────────
# 4. КЕШ КАДРОВ — считается один раз после rebuild
# ─────────────────────────────────────────────
def _build_frame_cache_chunk(args):
    """Для каждого кадра: активные узлы + список рёбер с весами."""
    frames, time_win, edges_src = args
    # edges_src: list of (u, v, w)
    result = {}
    for frame in frames:
        ct = timeline[frame]
        active = np.flatnonzero(np.abs(T - ct) < time_win).tolist()
        active_set = set(active)
        local_edges = [(u, v, w) for u, v, w in edges_src if u in active_set and v in active_set]
        result[frame] = {"active": active, "edges": local_edges}
    return result


def precompute_frame_cache(G, time_win):
    """Параллельно строит кеш для всех 120 кадров."""
    print("Кешируем кадры…")
    all_edges = [(u, v, d['weight']) for u, v, d in G.edges(data=True)]
    frame_ids = list(range(len(timeline)))
    n_workers = max(1, min(GRAPH_WORKERS, len(frame_ids)))
    chunks    = np.array_split(frame_ids, n_workers)

    result = [None] * len(timeline)
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        futs = [ex.submit(_build_frame_cache_chunk, (list(c), time_win, all_edges))
                for c in chunks if len(c) > 0]
        for f in as_completed(futs):
            for frame_id, data in f.result().items():
                result[frame_id] = data
    print("Кеш готов.")
    return result


def rebuild_everything(bh_mass, edge_thresh, time_win):
    global Quantumograph, frame_cache, bg_cache
    bg_cache = {}  # сбрасываем кеш фона при новой массе
    ax_main.set_title("Строим граф…", color='yellow', fontsize=11)
    fig.canvas.draw_idle()

    Quantumograph = build_graph(bh_mass, edge_thresh)
    n, e = Quantumograph.number_of_nodes(), Quantumograph.number_of_edges()

    ax_main.set_title("Кешируем кадры…", color='yellow', fontsize=11)
    fig.canvas.draw_idle()
    frame_cache = precompute_frame_cache(Quantumograph, time_win)

    ax_main.set_title(f"Готово: {n} узлов, {e} рёбер", color='lime', fontsize=11)
    fig.canvas.draw_idle()


# ─────────────────────────────────────────────
# 5. ФОН МЕТРИКИ — кешируется как numpy-массивы
# ─────────────────────────────────────────────
def get_bg_data(bh_mass, rs):
    key = (round(bh_mass, 1), round(rs, 2))
    if key in bg_cache:
        return bg_cache[key]

    rs_eff  = rs * bh_mass
    lim     = 2.0
    data    = {}

    # Изолинии потенциала
    gx = np.linspace(-lim, lim, 150)
    GX, GY = np.meshgrid(gx, gx)
    R   = np.sqrt(GX**2 + GY**2)
    Phi = -bh_mass / (R + 0.15)
    data['contour_XY'] = (GX, GY, Phi)

    # Деформированная сетка — окружности
    circles = []
    for r_phys in np.linspace(0.15, 2.0, 12):
        shrink = np.sqrt(max(0, 1.0 - rs_eff / (r_phys + 1e-4)))
        fall   = 1.0 / (1.0 + bh_mass * 0.4 / (r_phys + 0.05)**2)
        r_vis  = r_phys * shrink * fall
        phi_arr = np.linspace(0, 2*np.pi, 120)
        circles.append((r_vis * np.cos(phi_arr), r_vis * np.sin(phi_arr),
                        0.08 if r_phys > rs_eff * 2 else 0.04))
    data['circles'] = circles

    # Лучи
    rays = []
    for phi0 in np.linspace(0, 2*np.pi, 24, endpoint=False):
        r_arr  = np.linspace(0.05, 2.0, 60)
        shrink = np.sqrt(np.maximum(0, 1.0 - rs_eff / (r_arr + 1e-4)))
        fall   = 1.0 / (1.0 + bh_mass * 0.4 / (r_arr + 0.05)**2)
        r_vis  = r_arr * shrink * fall
        rays.append((r_vis * np.cos(phi0), r_vis * np.sin(phi0)))
    data['rays'] = rays

    # Поле стрелок
    ag = np.linspace(-1.6, 1.6, 10)
    AX_, AY_ = np.meshgrid(ag, ag)
    AR_ = np.sqrt(AX_**2 + AY_**2) + 0.1
    force = bh_mass / AR_**2
    Ux = -AX_ / AR_ * force; Uy = -AY_ / AR_ * force
    mag = np.sqrt(Ux**2 + Uy**2)
    data['quiver'] = (AX_, AY_, Ux/mag, Uy/mag, mag)

    # Горизонт
    rs_vis = rs_eff * np.sqrt(max(0, 1.0 - rs_eff/(rs_eff + 1e-4))) * 0.5
    theta  = np.linspace(0, 2*np.pi, 200)
    data['horizon'] = (rs_vis * np.cos(theta), rs_vis * np.sin(theta), rs_vis)

    bg_cache[key] = data
    return data


def draw_metric_background(ax, bh_mass):
    bg = get_bg_data(bh_mass, RS)

    if SHOW_CONTOUR:
        GX, GY, Phi = bg['contour_XY']
        ax.contour(GX, GY, Phi, levels=np.linspace(-12, -0.5, 14),
                   colors='white', alpha=0.06, linewidths=0.5)

    if SHOW_GRID:
        for cx, cy, alpha in bg['circles']:
            ax.plot(cx, cy, color='#4466aa', lw=0.4, alpha=alpha)
        for lx, ly in bg['rays']:
            ax.plot(lx, ly, color='#4466aa', lw=0.4, alpha=0.07)

    if SHOW_ARROWS:
        AX_, AY_, Ux, Uy, mag = bg['quiver']
        ax.quiver(AX_, AY_, Ux, Uy, mag, cmap='hot', alpha=0.12,
                  scale=30, width=0.002, headwidth=3)

    hx, hy, rs_vis = bg['horizon']
    ax.fill(hx, hy, color='black', zorder=5, alpha=1.0)
    ax.plot(hx, hy, color='#ff2200', lw=1.2, zorder=6, alpha=0.7)
    r_ph = min(1.5 * rs_vis, 0.4)
    theta = np.linspace(0, 2*np.pi, 200)
    ax.plot(r_ph * np.cos(theta), r_ph * np.sin(theta),
            '--', color='#ffaa00', lw=0.6, alpha=0.35, zorder=4)
    for rr, a in [(0.55, 0.18), (0.42, 0.25), (0.32, 0.15)]:
        ax.fill_betweenx(rr * np.sin(theta) * 0.25,
                         rr * np.cos(theta) - rr*0.08,
                         rr * np.cos(theta) + rr*0.08,
                         alpha=a, color='#ff6600', zorder=3)


# ─────────────────────────────────────────────
# 6. МАКЕТ
# ─────────────────────────────────────────────
fig = plt.figure(figsize=(12, 11), facecolor='black')
gs  = GridSpec(2, 1, height_ratios=[9, 1], hspace=0.04)
ax_main = fig.add_subplot(gs[0])
ax_main.set_facecolor('black'); ax_main.axis('off')
fig.add_subplot(gs[1]).set_facecolor('#111111')

sl_ax_mass   = fig.add_axes([0.06, 0.04, 0.20, 0.025], facecolor='#222222')
sl_ax_thresh = fig.add_axes([0.29, 0.04, 0.20, 0.025], facecolor='#222222')
sl_ax_win    = fig.add_axes([0.52, 0.04, 0.20, 0.025], facecolor='#222222')
sl_ax_rs     = fig.add_axes([0.75, 0.04, 0.13, 0.025], facecolor='#222222')

sl_mass   = widgets.Slider(sl_ax_mass,   'BH масса',   0.3, 4.0, valinit=BH_MASS,    color='#ff4400', valstep=0.1)
sl_thresh = widgets.Slider(sl_ax_thresh, 'Порог dist', 0.2, 1.5, valinit=EDGE_THRESH, color='#00ccff', valstep=0.05)
sl_win    = widgets.Slider(sl_ax_win,    'Окно T',     0.5, 3.5, valinit=TIME_WINDOW, color='#aa44ff', valstep=0.1)
sl_rs_val = widgets.Slider(sl_ax_rs,    'r_s',        0.05,0.8, valinit=RS,          color='#ffcc00', valstep=0.05)
for sl in (sl_mass, sl_thresh, sl_win, sl_rs_val):
    sl.label.set_color('white'); sl.valtext.set_color('white')

chk_ax = fig.add_axes([0.01, 0.00, 0.12, 0.035]); chk_ax.set_facecolor('#111111')
chk = widgets.CheckButtons(chk_ax, ['Сетка','Стрелки','Изолинии'], [SHOW_GRID,SHOW_ARROWS,SHOW_CONTOUR])
for t in chk.labels: t.set_color('white'); t.set_fontsize(7)

def on_check(label):
    global SHOW_GRID, SHOW_ARROWS, SHOW_CONTOUR
    if label == 'Сетка':    SHOW_GRID    = not SHOW_GRID
    if label == 'Стрелки':  SHOW_ARROWS  = not SHOW_ARROWS
    if label == 'Изолинии': SHOW_CONTOUR = not SHOW_CONTOUR
chk.on_clicked(on_check)

btn_ax = fig.add_axes([0.35, 0.007, 0.30, 0.025])
btn    = widgets.Button(btn_ax, 'Пересчитать граф + кеш кадров', color='#2a2a2a', hovercolor='#444444')
btn.label.set_color('white'); btn.label.set_fontsize(8)
btn.on_clicked(lambda e: rebuild_everything(sl_mass.val, sl_thresh.val, sl_win.val))


# ─────────────────────────────────────────────
# 7. ФУНКЦИЯ КАДРА — весь рендер через коллекции
# ─────────────────────────────────────────────
# Цвета и размеры узлов по расстоянию до ЧД — векторно
_NODE_COLORS = np.where(D < 0.5, 0, np.where(D < 1.0, 1, 2))
_COLOR_MAP   = ['#ff2200', '#ffaa00', '#7a00ff']
_SIZE_MAP    = [45, 30, 18]

def update(frame):
    global RS
    RS       = sl_rs_val.val
    bh_mass  = sl_mass.val
    time_win = sl_win.val
    current_time = timeline[frame]

    ax_main.clear()
    ax_main.set_facecolor('#050510')
    ax_main.set_xlim(-2.0, 2.0); ax_main.set_ylim(-2.0, 2.0)
    ax_main.axis('off')
    ax_main.set_title(
        f"Quantumograph v14 · Z⁴ Торус · Искривление метрики (ЧД)\n"
        f"T={current_time:.2f}  M={bh_mass:.1f}  r_s={RS:.2f}  "
        f"порог={sl_thresh.val:.2f}  окно={time_win:.1f}",
        color='white', fontsize=9, pad=6)

    draw_metric_background(ax_main, bh_mass)

    if Quantumograph is None or frame_cache is None:
        ax_main.text(0, 0, "Идёт расчёт…", color='yellow',
                     ha='center', va='center', fontsize=12, zorder=20)
        return

    cached = frame_cache[frame]
    active = cached["active"]
    if len(active) < 3:
        ax_main.text(0, 0, f'Кубитов: {len(active)} (мало)',
                     color='red', ha='center', va='center', fontsize=12, zorder=20)
        return

    # Позиции — векторно
    active_arr = np.array(active)
    proj       = schwarzschild_project_vec(P[active_arr], bh_mass, RS)  # (N,2)
    id_to_idx  = {node: i for i, node in enumerate(active)}

    # ── РЁБРА через LineCollection (1 вызов matplotlib) ──────────────────
    edge_data = cached["edges"]
    if edge_data:
        eu  = np.array([id_to_idx[u] for u, v, w in edge_data if u in id_to_idx and v in id_to_idx])
        ev  = np.array([id_to_idx[v] for u, v, w in edge_data if u in id_to_idx and v in id_to_idx])
        ew  = np.array([w             for u, v, w in edge_data if u in id_to_idx and v in id_to_idx])

        if len(eu) > 0:
            segs   = np.stack([proj[eu], proj[ev]], axis=1)   # (E, 2, 2)
            wn     = (ew - ew.min()) / (ew.max() - ew.min() + 1e-9)
            alphas = 0.05 + 0.30 * wn
            lc = LineCollection(segs,
                                colors=[(0.0, 0.83, 1.0, float(a)) for a in alphas],
                                linewidths=0.5 + wn * 0.8,
                                zorder=7)
            ax_main.add_collection(lc)

    # ── УЗЛЫ через scatter (1 вызов) ─────────────────────────────────────
    cidx   = _NODE_COLORS[active_arr]
    colors = [_COLOR_MAP[c] for c in cidx]
    sizes  = [_SIZE_MAP[c]  for c in cidx]
    ax_main.scatter(proj[:, 0], proj[:, 1],
                    s=sizes, c=colors, alpha=0.88, zorder=8, linewidths=0)

    # Центр ЧД
    ax_main.plot(0, 0, 'o', ms=14, color='black',
                 mec='#ff2200', mew=1.5, alpha=0.9, zorder=15)
    ax_main.plot(0, 0, 'o', ms=5,  color='white', alpha=0.4, zorder=16)

    # Статистика
    ax_main.text(0, -1.93,
                 f"Активных: {len(active)}  |  Рёбер: {len(edge_data)}  |  "
                 f"Всего рёбер: {Quantumograph.number_of_edges()}",
                 color='#555555', ha='center', va='bottom', fontsize=7.5, zorder=20)

    for col, lbl, yy in [('#ff2200','r<0.5 (горизонт)',-1.79),
                          ('#ffaa00','r<1.0 (эргосфера)',-1.87),
                          ('#7a00ff','r>1.0 (внешнее)',  -1.95)]:
        ax_main.plot(-1.98, yy, 'o', color=col, ms=5, zorder=20)
        ax_main.text(-1.88, yy, lbl, color='#888888', va='center', fontsize=7, zorder=20)


# ─────────────────────────────────────────────
# 8. ЗАПУСК
# ─────────────────────────────────────────────
if __name__ == "__main__":
    rebuild_everything(BH_MASS, EDGE_THRESH, TIME_WINDOW)

    ani = FuncAnimation(fig, update, frames=len(timeline),
                        interval=ANIM_SPEED, repeat=True, blit=False)

    fig.subplots_adjust(left=0.02, right=0.98, top=0.97, bottom=0.08)
    plt.show()
