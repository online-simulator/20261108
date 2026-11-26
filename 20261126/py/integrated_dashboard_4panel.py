import matplotlib.pyplot as plt
import numpy as np

# 1. ご提示CSVから導出された厳密定数
dS0 = -0.8366732334793517
base_decay = 1.0 + dS0  # 約 0.163327
capa0 = -dS0 / base_decay  # 約 5.122695
alpha_svp_asymp = 0.5464698129876524
kappa_kl = 0.401007  # カバチャンスキー＝レーベンシュタイン接吻数上限指数

fig, axes = plt.subplots(2, 2, figsize=(12.8, 9.4), dpi=150)
ax1, ax2, ax3, ax4 = axes.flatten()

# --- パネル1（左上）：各算術曲線のフルヴィッツ余裕度 (ln|Aut| / capa0) ---
curve_names = [
    'X0(37)\n(g=2)',
    'Fermat F4\n(g=3)',
    'Fermat F5\n(g=6)',
    'Klein\n(g=3)',
    'Fermat F6\n(g=10)',
    'Macbeath\n(g=7)',
    'Hurwitz 3\n(g=14)',
]
aut_orders = [4, 96, 150, 168, 216, 504, 1092]
margins = [np.log(a) / capa0 for a in aut_orders]
colors_p1 = [
    '#e74c3c' if m < 1.0 else ('#f39c12' if m < 1.01 else '#27ae60')
    for m in margins
]

bars = ax1.bar(
    curve_names,
    margins,
    color=colors_p1,
    edgecolor='black',
    linewidth=0.8,
    alpha=0.88,
)
ax1.axhline(
    1.0,
    color='black',
    linestyle='--',
    linewidth=1.4,
    label=r'Critical Threshold ($\ln|\mathrm{Aut}| = \mathrm{capa}_0$)',
)
for bar, m in zip(bars, margins):
    ax1.text(
        bar.get_x() + bar.get_width() / 2,
        m + 0.03,
        f'{m:.3f}',
        ha='center',
        va='bottom',
        fontsize=8,
        fontweight='bold',
    )

ax1.set_title(
    r'1. Hurwitz Symmetry Margin Index ($\ln|\mathrm{Aut}| / \mathrm{capa}_0$)',
    fontsize=10.5,
    fontweight='bold',
)
ax1.set_ylabel(r'Margin Ratio (Stable if $\geq 1.0$)', fontsize=9.5)
ax1.set_ylim(0, 1.55)
ax1.tick_params(axis='x', labelsize=8)
ax1.grid(axis='y', alpha=0.3)
ax1.legend(loc='upper left', fontsize=8.5)

# --- パネル2（右上）：実効Core-SVPふるい法指数とNIST基準の比較 ---
dims = np.array([1, 2, 3, 4, 5, 6, 8, 12, 24, 36])
alpha_csv = np.array([
    0.8261629401947822,
    0.6924278365287646,
    0.6452136247047902,
    0.6210782156897083,
    0.6064235810559829,
    0.5965805952646658,
    0.5841936653788219,
    0.5717134654391159,
    0.5591391284668298,
    0.5549266172070758,
])
alpha_sieve = 0.5 * alpha_csv

ax2.plot(
    dims,
    alpha_sieve,
    'o-',
    color='#2980b9',
    linewidth=2.2,
    markersize=6,
    label=r'Proposed Sieving Index: $0.5 \cdot \alpha(s, d, g)$',
)
ax2.axhline(
    0.292,
    color='#c0392b',
    linestyle='--',
    linewidth=1.4,
    label='NIST Classical Sieving (BDGL: 0.292)',
)
ax2.axhline(
    0.265,
    color='#8e44ad',
    linestyle='-.',
    linewidth=1.4,
    label='NIST Quantum Sieving (0.265)',
)
ax2.axhline(
    0.5 * alpha_svp_asymp,
    color='#2c3e50',
    linestyle=':',
    linewidth=1.4,
    label=f'Asymptotic Limit ({0.5 * alpha_svp_asymp:.4f})',
)

ax2.set_title(
    '2. Effective Core-SVP Sieving Exponent vs. NIST Standards',
    fontsize=10.5,
    fontweight='bold',
)
ax2.set_xlabel('Lattice Dimension ($d$)', fontsize=9.5)
ax2.set_ylabel('Exponent $c$ in $2^{c \cdot d}$', fontsize=9.5)
ax2.set_xticks([1, 4, 8, 12, 16, 20, 24, 28, 32, 36])
ax2.set_ylim(0.24, 0.43)
ax2.grid(True, alpha=0.3)
ax2.legend(loc='upper right', fontsize=8)

# --- パネル3（左下）：アデール空間圧縮利得 (1 - (0.401 + alpha)) ---
adelic_gain = 1.0 - (kappa_kl + alpha_csv)
ax3.plot(
    dims,
    adelic_gain,
    's-',
    color='#16a085',
    linewidth=2.2,
    markersize=6,
    label=r'Adelic Gain: $1 - (\kappa_{\mathrm{KL}} + \alpha)$',
)
ax3.axhline(0.0, color='black', linewidth=1.2, linestyle='-')
ax3.axhline(
    1.0 - (kappa_kl + alpha_svp_asymp),
    color='#d35400',
    linestyle='--',
    linewidth=1.3,
    label='Asymptotic Gain (+0.0525 bits/dim)',
)
ax3.fill_between(
    dims,
    0,
    adelic_gain,
    where=(adelic_gain >= 0),
    color='#16a085',
    alpha=0.18,
    interpolate=True,
)
ax3.fill_between(
    dims,
    0,
    adelic_gain,
    where=(adelic_gain < 0),
    color='#e74c3c',
    alpha=0.15,
    interpolate=True,
)

ax3.set_title(
    '3. Adelic Space Compression Gain (bits / dimension)',
    fontsize=10.5,
    fontweight='bold',
)
ax3.set_xlabel('Lattice Dimension ($d$)', fontsize=9.5)
ax3.set_ylabel('Compression Gain (bits / dim)', fontsize=9.5)
ax3.set_xticks([1, 4, 8, 12, 16, 20, 24, 28, 32, 36])
ax3.set_ylim(-0.25, 0.08)
ax3.grid(True, alpha=0.3)
ax3.legend(loc='lower right', fontsize=8)

# --- パネル4（右下）：量子並列繰り返し減衰指数 lambda_QPR(d, g) ---
g_vals = np.array([0, 1, 2, 3, 7, 14])
for d_val, col, ls in [(1, '#2980b9', '-'), (2, '#e67e22', '--'), (3, '#c0392b', '-.')]:
    lam_qpr = -(2 * g_vals + 1) * d_val * np.log(base_decay)
    ax4.plot(
        g_vals,
        lam_qpr,
        marker='o',
        linestyle=ls,
        color=col,
        linewidth=2.0,
        markersize=5.5,
        label=f'Dimension $d={d_val}$',
    )

ax4.axvline(3, color='#7f8c8d', linestyle=':', linewidth=1.2)
ax4.text(3.2, 130, 'Klein (g=3)', fontsize=8, color='#34495e', fontweight='bold')
ax4.axvline(7, color='#7f8c8d', linestyle=':', linewidth=1.2)
ax4.text(7.2, 130, 'Macbeath (g=7)', fontsize=8, color='#34495e', fontweight='bold')

ax4.set_title(
    r'4. Quantum Parallel Repetition Decay Rate $\lambda_{\mathrm{QPR}}(d, g)$',
    fontsize=10.5,
    fontweight='bold',
)
ax4.set_xlabel('Arithmetic Curve Genus ($g$)', fontsize=9.5)
ax4.set_ylabel(r'Decay Exponent per Round ($\lambda_{\mathrm{QPR}}$)', fontsize=9.5)
ax4.set_xticks(g_vals)
ax4.grid(True, alpha=0.3)
ax4.legend(loc='upper left', fontsize=8.5)

plt.tight_layout(pad=1.8)
plt.savefig('integrated_dashboard_4panel.png', dpi=200, bbox_inches='tight')
plt.show()
