from math import comb
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # 3D描画用
import numpy as np
import scipy.special as sp

# 1. Borwein加速法によるリーマンゼータ関数の高精度計算（Re(s) >= 0.5 用）
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


# 2. 関数等式 Lambda(s) = Lambda(1 - s) を組み込んだ完備ゼータ関数
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


# 3. 3D描画用に滑らかで見やすい解像度のグリッドを生成
re_vals = np.linspace(-2.0, 4.0, 141)
im_vals = np.linspace(0.0, 28.0, 181)
RE, IM = np.meshgrid(re_vals, im_vals)
S = RE + 1j * IM

# 4. 有効シャノン＝ハートレー情報容量の計算（d=1, g=0）
d = 1.0
g = 0.0
sharpness = 2.0 * (2.0 * g + 1.0) * d
q_shannon = 1.0 + (1.0 / np.log2(3.0)) / d

LAM = lambda_adelic(S)
CAPA = sharpness * np.log2(1.0 + q_shannon / LAM)

# 3Dの山の頂上が突き抜けすぎないよう適切な高さでクリップ
Z_RE = np.clip(CAPA.real, -6.0, 42.0)
Z_IM = CAPA.imag

# 5. 3Dサーフェスプロットの描画
fig = plt.figure(figsize=(15.0, 6.8), dpi=150)

# --- 左パネル：実部 Re[capa_eff(s)] の3D山脈 ---
ax1 = fig.add_subplot(1, 2, 1, projection='3d')
surf1 = ax1.plot_surface(
    RE, IM, Z_RE, cmap='RdYlBu_r',
    edgecolor='none', alpha=0.92, rstride=1, cstride=1, antialiased=True
)
# 底面に等高線を投影（山と谷の位置を把握しやすくする）
ax1.contour(RE, IM, Z_RE, zdir='z', offset=-8.0, levels=[-2, 2, 10, 20, 30], cmap='RdYlBu_r', linewidths=1.0)

# 臨界線 Re(s)=1/2 上の軌跡を黒い太線で強調
idx_half = np.argmin(np.abs(re_vals - 0.5))
ax1.plot(
    RE[:, idx_half], IM[:, idx_half], Z_RE[:, idx_half] + 0.3,
    color='black', linewidth=1.8, label=r'Critical Line $\mathrm{Re}(s)=1/2$'
)

ax1.set_title(r'3D Surface: Real Part $\mathrm{Re}[\mathrm{capa}_{\mathrm{eff}}(s)]$', fontsize=11.5, fontweight='bold', pad=10)
ax1.set_xlabel(r'$\mathrm{Re}(s)$', fontsize=10, labelpad=6)
ax1.set_ylabel(r'$\mathrm{Im}(s)$', fontsize=10, labelpad=6)
ax1.set_zlabel(r'Capacity (bits)', fontsize=10, labelpad=6)
ax1.set_zlim(-8.0, 43.0)
ax1.view_init(elev=28, azim=-55)  # 見やすい斜め上からの俯瞰アングル
fig.colorbar(surf1, ax=ax1, shrink=0.58, aspect=18, pad=0.08)
ax1.legend(loc='upper left', fontsize=8.5)

# --- 右パネル：虚部 Im[capa_eff(s)] の3Dらせん階段（位相ホロノミー） ---
ax2 = fig.add_subplot(1, 2, 2, projection='3d')
surf2 = ax2.plot_surface(
    RE, IM, Z_IM, cmap='twilight_shifted',
    edgecolor='none', alpha=0.92, rstride=1, cstride=1, antialiased=True
)
ax2.contour(RE, IM, Z_IM, zdir='z', offset=-11.0, levels=[-6, -3, 0, 3, 6], cmap='twilight_shifted', linewidths=1.0)

ax2.set_title(r'3D Surface: Imaginary Part $\mathrm{Im}[\mathrm{capa}_{\mathrm{eff}}(s)]$ (Phase)', fontsize=11.5, fontweight='bold', pad=10)
ax2.set_xlabel(r'$\mathrm{Re}(s)$', fontsize=10, labelpad=6)
ax2.set_ylabel(r'$\mathrm{Im}(s)$', fontsize=10, labelpad=6)
ax2.set_zlabel(r'Phase Shift', fontsize=10, labelpad=6)
ax2.set_zlim(-11.0, 11.0)
ax2.view_init(elev=28, azim=-55)
fig.colorbar(surf2, ax=ax2, shrink=0.58, aspect=18, pad=0.08)

plt.tight_layout()
plt.savefig('zeta_capacity_3d.png', dpi=200, bbox_inches='tight')
plt.show()
