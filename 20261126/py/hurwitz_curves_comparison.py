import matplotlib.pyplot as plt
import numpy as np

# 1. ご提示CSVから導出された基本定数
dS0 = -0.8366732334793517
capa0 = -dS0 / (1.0 + dS0)  # 約 5.122695

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.2), dpi=150)

# --- 左パネル：量子もつれ容量の飽和曲線と対称性天井 ---
n_vals = np.linspace(0, 8, 250)
hurwitz_curves = [
    ('Klein Quartic (g=3, |Aut|=168)', 3, 168, '#c0392b', '-'),
    ('Macbeath Curve (g=7, |Aut|=504)', 7, 504, '#2980b9', '--'),
    ('1st Hurwitz Triplet (g=14, |Aut|=1092)', 14, 1092, '#27ae60', '-.'),
]

for label, g, aut, col, ls in hurwitz_curves:
    sharp = (2 * g + 1) * 1.0
    # n=0 でのゼロ除算を回避しつつ容量を計算
    with np.errstate(divide='ignore', invalid='ignore'):
        capa = np.where(
            n_vals == 0,
            0.0,
            (-dS0)
            * n_vals
            * sharp
            * (1.0 - np.exp(-1.0 / (n_vals * sharp * (1.0 + dS0)))),
        )
    ax1.plot(n_vals, capa, linestyle=ls, color=col, linewidth=2.3, label=label)
    ax1.axhline(np.log(aut), color=col, linestyle=':', linewidth=1.3, alpha=0.85)

ax1.axhline(
    capa0,
    color='black',
    linewidth=1.6,
    linestyle='-',
    label=r'Zeta Asymptote $\mathrm{capa}_0 \approx 5.123$',
)

# 対称性天井のテキストラベル
ax1.text(
    4.2,
    np.log(168) + 0.14,
    r'Klein $\ln(168) \approx 5.124$ (Critical)',
    fontsize=8.5,
    color='#c0392b',
    fontweight='bold',
)
ax1.text(
    4.2,
    np.log(504) + 0.14,
    r'Macbeath $\ln(504) \approx 6.223$',
    fontsize=8.5,
    color='#2980b9',
    fontweight='bold',
)
ax1.text(
    4.2,
    np.log(1092) + 0.14,
    r'Hurwitz Triplet $\ln(1092) \approx 6.996$',
    fontsize=8.5,
    color='#27ae60',
    fontweight='bold',
)

ax1.set_title(
    'Hurwitz Curves: Capacity Saturation vs. Symmetry Ceilings',
    fontsize=11,
    fontweight='bold',
)
ax1.set_xlabel('Parallel Repetitions ($n$)', fontsize=10.5)
ax1.set_ylabel('Capacity of Entanglement ($d=1$)', fontsize=10.5)
ax1.set_ylim(0, 7.6)
ax1.grid(True, alpha=0.3)
ax1.legend(loc='lower right', fontsize=8.2, framealpha=0.92)

# --- 右パネル：スピロ充填可能セル数 log10(Nsol) の推移 ---
x_pts = np.array([
    1.4426950408889634,
    1.498425536800067,
    1.5541560327111708,
    1.6098865286222743,
    1.6656170245333781,
    1.7213475204444817,
    1.7770780163555853,
    1.832808512266689,
    1.8885390081777926,
    1.9442695040888962,
    2.0,
])
y_g1_d1 = np.array([
    2.662906276629065,
    2.3585866617637925,
    2.0630001212809894,
    1.776161045596143,
    1.4979877599986904,
    1.2283302805968233,
    0.966991321745623,
    0.7137420985925562,
    0.46833414007718555,
    0.23050805777504957,
    0.0,
])
y_unit = y_g1_d1 / 3.0  # (2g+1) あたりの単位スケーリング

for label, g, aut, col, ls in hurwitz_curves:
    y_vals = y_unit * (2 * g + 1)
    ax2.plot(
        x_pts,
        y_vals,
        linestyle=ls,
        marker='o',
        markersize=4.5,
        color=col,
        linewidth=2.2,
        label=label,
    )

ax2.set_title(
    r'Szpiro Packable Cells $\log_{10}(N_{\mathrm{sol}})$ ($d=1$)',
    fontsize=11,
    fontweight='bold',
)
ax2.set_xlabel('Parameter $x$ ($s = x$)', fontsize=10.5)
ax2.set_ylabel(r'$\log_{10}(N_{\mathrm{sol}})$', fontsize=10.5)
ax2.set_xlim(1.42, 2.02)
ax2.set_ylim(0, 27.5)
ax2.grid(True, alpha=0.3)
ax2.legend(loc='upper right', fontsize=8.2, framealpha=0.92)

plt.tight_layout()
plt.savefig('hurwitz_curves_comparison.png', dpi=200, bbox_inches='tight')
plt.show()
