from math import comb
import matplotlib.pyplot as plt
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

    # Re(s) >= 0.5 の領域
    mask_pos = s_grid.real >= 0.5
    s_pos = s_grid[mask_pos]
    res[mask_pos] = (
        (np.pi ** (-s_pos / 2.0))
        * sp.gamma(s_pos / 2.0)
        * zeta_half_plane(s_pos)
    )

    # Re(s) < 0.5 の領域（s <-> 1 - s で鏡像反転）
    s_neg = s_grid[~mask_pos]
    s_ref = 1.0 - s_neg
    res[~mask_pos] = (
        (np.pi ** (-s_ref / 2.0))
        * sp.gamma(s_ref / 2.0)
        * zeta_half_plane(s_ref)
    )

    return res


# 3. 複素平面グリッドの生成（Re(s) in [-2, 4], Im(s) in）
re_vals = np.linspace(-2.0, 4.0, 241)
im_vals = np.linspace(0.0, 30.0, 301)
RE, IM = np.meshgrid(re_vals, im_vals)
S = RE + 1j * IM

# 4. 有効シャノン＝ハートレー情報容量の計算（d=1, g=0）
d = 1.0
g = 0.0
sharpness = 2.0 * (2.0 * g + 1.0) * d
q_shannon = 1.0 + (1.0 / np.log2(3.0)) / d

LAM = lambda_adelic(S)
CAPA = sharpness * np.log2(1.0 + q_shannon / LAM)

# 特異点の発散を見やすくするため実部をクリップ
RE_CAPA = np.clip(CAPA.real, -10.0, 40.0)
IM_CAPA = CAPA.imag

# 5. ヒートマップ・等高線の描画
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 6.2), dpi=150, sharey=True)

# --- 左パネル：実部 Re[capa_eff(s)] ---
levels_re = np.linspace(-8, 36, 23)
cf1 = ax1.contourf(RE, IM, RE_CAPA, levels=levels_re, cmap='RdYlBu_r', extend='both')
cs1 = ax1.contour(
    RE, IM, RE_CAPA,
    levels=[-4, 0, 2.038, 4.082, 10, 20, 30],
    colors='black', linewidths=0.65, alpha=0.65
)
ax1.clabel(cs1, inline=True, fontsize=7.5, fmt='%.1f')

ax1.axvline(0.5, color='white', linestyle='--', linewidth=1.3, label=r'Critical Line $\mathrm{Re}(s)=1/2$')
ax1.axvline(0.0, color='gray', linestyle=':', linewidth=1.0)
ax1.axvline(1.0, color='gray', linestyle=':', linewidth=1.0)

# 非自明なリーマン零点（最初の3つ）のプロット
zeros_t = [14.1347, 21.0220, 25.0109]
for idx, zt in enumerate(zeros_t, 1):
    ax1.plot(0.5, zt, 'w*', markersize=9, markeredgecolor='black', markeredgewidth=0.8)
    ax1.text(
        0.68, zt - 0.35, rf'$\gamma_{idx}\approx{zt:.2f}$',
        color='white', fontsize=8.5, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.15', facecolor='black', alpha=0.55, lw=0)
    )

# 実軸上の評価点（log2(e), 2, 1-log2(e)）のプロット
eval_x = [np.log2(np.e), 2.0, 1.0 - np.log2(np.e)]
ax1.plot(eval_x, [0.0, 0.0, 0.0], 'ko', markersize=5.5, markerfacecolor='lime')
ax1.text(1.44, 0.8, r'$\log_2 e$', color='black', fontsize=8, fontweight='bold', ha='center')
ax1.text(2.05, 0.8, r'$s=2$', color='black', fontsize=8, fontweight='bold', ha='left')
ax1.text(-0.44, 0.8, r'$1-\log_2 e$', color='black', fontsize=8, fontweight='bold', ha='center')

cbar1 = fig.colorbar(cf1, ax=ax1, pad=0.02, aspect=28)
cbar1.set_label(r'$\mathrm{Re}[\mathrm{capa}_{\mathrm{eff}}(s)]$ (bits)', fontsize=9.5)
ax1.set_title(r'Real Part: $\mathrm{Re}[\mathrm{capa}_{\mathrm{eff}}(s)]$ ($d=1, g=0$)', fontsize=11, fontweight='bold')
ax1.set_xlabel(r'$\mathrm{Re}(s)$', fontsize=10.5)
ax1.set_ylabel(r'$\mathrm{Im}(s)$', fontsize=10.5)
ax1.legend(loc='upper right', fontsize=8.2, framealpha=0.88)

# --- 右パネル：虚部 Im[capa_eff(s)]（位相ホロノミー） ---
levels_im = np.linspace(-9.1, 9.1, 25)
cf2 = ax2.contourf(RE, IM, IM_CAPA, levels=levels_im, cmap='twilight_shifted')
cs2 = ax2.contour(RE, IM, IM_CAPA, levels=[-6, -3, 0, 3, 6], colors='black', linewidths=0.6, alpha=0.55)
ax2.clabel(cs2, inline=True, fontsize=7.5, fmt='%.1f')

ax2.axvline(0.5, color='black', linestyle='--', linewidth=1.3, label=r'Critical Line $\mathrm{Re}(s)=1/2$')
for zt in zeros_t:
    ax2.plot(0.5, zt, 'ko', markersize=6, markerfacecolor='yellow')

cbar2 = fig.colorbar(cf2, ax=ax2, pad=0.02, aspect=28)
cbar2.set_label(r'$\mathrm{Im}[\mathrm{capa}_{\mathrm{eff}}(s)]$ (Phase Shift)', fontsize=9.5)
ax2.set_title(r'Imaginary Part: $\mathrm{Im}[\mathrm{capa}_{\mathrm{eff}}(s)]$ (Phase Holonomy)', fontsize=11, fontweight='bold')
ax2.set_xlabel(r'$\mathrm{Re}(s)$', fontsize=10.5)
ax2.legend(loc='upper right', fontsize=8.2, framealpha=0.88)

plt.tight_layout()
plt.savefig('zeta_capacity_heatmap.png', dpi=200, bbox_inches='tight')
plt.show()
