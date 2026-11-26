from math import comb
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import scipy.special as sp

# 1. Borwein加速法によるリーマンゼータ関数の高精度計算
N_TERMS = 45
n_arr = np.arange(1, N_TERMS + 1)
weights = np.array([
    sum(comb(N_TERMS - 1, j) for j in range(k, N_TERMS))
    / (2.0 ** (N_TERMS - 1))
    for k in range(N_TERMS)
])


def zeta_half_plane(s_arr):
    terms = (
        ((-1.0) ** (n_arr[None, :] - 1))
        * weights[None, :]
        / (n_arr[None, :] ** s_arr[:, None])
    )
    eta = np.sum(terms, axis=1)
    denom = 1.0 - 2.0 ** (1.0 - s_arr)
    denom = np.where(np.abs(denom) < 1e-12, 1e-12 + 0j, denom)
    return eta / denom


def lambda_adelic(s_grid):
    res = np.zeros_like(s_grid, dtype=np.complex128)
    mask_pos = s_grid.real >= 0.5
    s_pos = s_grid[mask_pos]
    res[mask_pos] = (
        (np.pi ** (-s_pos / 2.0))
        * sp.gamma(s_pos / 2.0)
        * zeta_half_plane(s_pos)
    )
    s_neg = s_grid[~mask_pos]
    s_ref = 1.0 - s_neg
    res[~mask_pos] = (
        (np.pi ** (-s_ref / 2.0))
        * sp.gamma(s_ref / 2.0)
        * zeta_half_plane(s_ref)
    )
    return res


# 2. グリッド生成（低次元〜高次元の分岐点シフトが見やすい t in [0, 18] にフォーカス）
re_vals = np.linspace(-1.5, 2.5, 131)
im_vals = np.linspace(0.0, 18.0, 181)
RE, IM = np.meshgrid(re_vals, im_vals)
S = RE + 1j * IM
LAM = lambda_adelic(S)

g = 3.0  # クラインの4次曲線固定

# 3. 描画設定（左：d=24 リーチ格子の3D地形、右：臨界線上の d=1, 8, 24 分岐点シフト比較）
fig = plt.figure(figsize=(15.5, 6.8), dpi=150)

# --- 左パネル：リーチ格子 (d=24, g=3) の 3D サーフェス ---
ax1 = fig.add_subplot(1, 2, 1, projection='3d')
d_leech = 24.0
sharp_24 = 2.0 * (2.0 * g + 1.0) * d_leech
q_24 = 1.0 + (1.0 / np.log2(3.0)) / d_leech
CAPA_24 = sharp_24 * np.log2(1.0 + q_24 / LAM)

# 1次元あたりの正規化容量 (CAPA / d) でプロットしてスケールを統一
Z_24_norm = np.clip(CAPA_24.real / d_leech, -35.0, 120.0)

surf1 = ax1.plot_surface(
    RE, IM, Z_24_norm, cmap='RdYlBu_r',
    edgecolor='none', alpha=0.92, rstride=1, cstride=1, antialiased=True
)
ax1.contour(RE, IM, Z_24_norm, zdir='z', offset=-45.0, levels=[-20, 0, 14.3, 40, 80], cmap='RdYlBu_r', linewidths=1.0)

idx_half = np.argmin(np.abs(re_vals - 0.5))
ax1.plot(
    RE[:, idx_half], IM[:, idx_half], Z_24_norm[:, idx_half] + 1.5,
    color='black', linewidth=1.8, label=r'Critical Line $\mathrm{Re}(s)=1/2$'
)

ax1.set_title(r'Leech Lattice ($d=24, g=3$): Normalized 3D Capacity / $d$', fontsize=11, fontweight='bold', pad=10)
ax1.set_xlabel(r'$\mathrm{Re}(s)$', fontsize=10, labelpad=6)
ax1.set_ylabel(r'$\mathrm{Im}(s)$', fontsize=10, labelpad=6)
ax1.set_zlabel(r'Capacity / $d$ (bits/dim)', fontsize=10, labelpad=6)
ax1.set_zlim(-45.0, 125.0)
ax1.view_init(elev=26, azim=-58)
fig.colorbar(surf1, ax=ax1, shrink=0.58, aspect=18, pad=0.08)
ax1.legend(loc='upper left', fontsize=8.5)

# --- 右パネル：臨界線 Re(s)=1/2 上での次元別 (d=1, 8, 24) 分岐点シフト拡大図 ---
ax2 = fig.add_subplot(1, 2, 2)
t_line = np.linspace(0.0, 16.5, 1000)
s_crit = 0.5 + 1j * t_line
lam_crit = lambda_adelic(s_crit).real

configs = [
    (1.0, 'd = 1 (Base)', '#2980b9', '-'),
    (8.0, 'd = 8 ($E_8$ Lattice)', '#e67e22', '--'),
    (24.0, 'd = 24 (Leech Lattice)', '#c0392b', '-.')
]

for d_val, lbl, col, ls in configs:
    q_d = 1.0 + (1.0 / np.log2(3.0)) / d_val
    sharp_d = 2.0 * (2.0 * g + 1.0) * d_val
    # 1次元あたりの正規化実部容量 Re[capa_eff] / d
    capa_norm = (sharp_d / d_val) * np.log2(np.abs(1.0 + q_d / lam_crit))
    ax2.plot(t_line, capa_norm, linestyle=ls, color=col, linewidth=2.2, label=f'{lbl}: $q={q_d:.4f}$')

    # 分岐点 (Lambda + q_d = 0) の正確な位置を検出してマーカー表示
    idx_branch = np.argmin(np.abs(lam_crit[:400] + q_d))
    t_b = t_line[idx_branch]
    ax2.axvline(t_b, color=col, linestyle=':', linewidth=1.2, alpha=0.8)
    ax2.annotate(
        rf'$t_b \approx {t_b:.2f}$', xy=(t_b, -22 + d_val * 0.4),
        xytext=(t_b + 0.35, -24 + d_val * 0.55),
        fontsize=8.5, fontweight='bold', color=col,
        arrowprops=dict(arrowstyle='->', color=col, lw=1.0)
    )

# 第1リーマン零点 gamma_1 の位置
ax2.axvline(14.1347, color='black', linestyle='--', linewidth=1.4, label=r'1st Riemann Zero $\gamma_1 \approx 14.135$')
ax2.axhline(0.0, color='gray', linewidth=0.9)

ax2.set_title(r'Branch Point Shift on Critical Line ($d=1 \to 8 \to 24$)', fontsize=11, fontweight='bold', pad=10)
ax2.set_xlabel(r'Imaginary Part $t$ along $\mathrm{Re}(s) = 1/2$', fontsize=10.5)
ax2.set_ylabel(r'Normalized Capacity $\mathrm{Re}[\mathrm{capa}_{\mathrm{eff}}] / d$ (bits/dim)', fontsize=10.5)
ax2.set_xlim(0.0, 16.5)
ax2.set_ylim(-30.0, 110.0)
ax2.grid(True, alpha=0.35)
ax2.legend(loc='upper left', fontsize=8.5, framealpha=0.92)

plt.tight_layout()
plt.savefig('zeta_capacity_high_dim_shift.png', dpi=200, bbox_inches='tight')
plt.show()
