import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

# ==============================================================================
# 0. 共通の基本定数・マスター方程式の定義
# ==============================================================================
plt.style.use('seaborn-v0_8-whitegrid')

gamma_1 = 14.134725141734693790
def integrand(t):
    return np.log(t / (2 * np.pi)) / (t * np.log(t))**2

integral_val, _ = quad(integrand, gamma_1, np.inf)
L_inf = np.sqrt(2 * np.pi * integral_val)  # ~0.2391632423
d_crit = 1.0 / np.sqrt(L_inf)              # ~2.04481
d_inf_GL1 = np.sqrt(1.5 / L_inf)           # ~2.50437
d_Gauss = np.sqrt(2 * np.pi)               # ~2.50663

def calc_d_eff(C_val, I_val, n=1, r=0):
    """GL_n 拡張有効次元 d_eff の計算"""
    val = n * (1.0 + r) * (1.0 + np.log(I_val) / np.log(C_val)) / L_inf
    return np.sqrt(max(val, 0.0))


# ==============================================================================
# Figure 1: 超大規模統一数論相図 (siegel_zero_phase_diagram_ultrascale.png)
# ==============================================================================
def generate_figure_1():
    imag_quad_data = [
        (3, 1, "D=-3", (-8, 8), None),
        (11, 1, "D=-11", None, (-15, 10)),
        (43, 1, "D=-43", None, (-12, 8)),
        (163, 1, "D=-163 (Max Heegner)", None, (-45, -14)),
        (427, 2, "D=-427 (h=2)", None, (-30, -14)),
        (907, 3, "D=-907 (h=3)", None, (6, -6)),
        (239, 15, "D=-239 (h=15)", None, (6, 4)),
        (959, 36, "D=-959 (h=36)", None, (-45, 7)),
        (18483, 8, "", None, None),
        (999983, 552, "D=-999983 (h=552)", (-55, 8), None)
    ]
    real_quad_data = [
        (5, 0.4812, "D=5", (6, -4), None),
        (8, 0.8814, "D=8", (-25, -8), (6, -8)),
        (12, 1.3170, "D=12", None, (-26, 7)),
        (73, 7.6667, "D=73 (R=7.67)", None, (-55, 7)),
        (401, 25.12, "", None, None)
    ]
    cubic_data = [
        (49, 0.5255, "Cubic D=49", None, (6, -5)),
        (81, 0.8493, "Cubic D=81", None, (-48, -12)),
        (148, 1.6624, "", None, None),
        (1957, 9.3125, "Cubic D=1957", (6, 4), None)
    ]
    elliptic_data = [
        (11, 1.0000, 0, "11a1 (r=0)", None, (6, -10)),
        (681, 9.0000, 0, "681b1 (Sha=9)", None, (6, -6)),
        (37, 0.05111, 1, "37a1 (r=1)", (7, -4), None),
        (57, 0.21463, 1, "57a1 (1st Cross)", None, (-55, 8)),
        (91, 0.09427, 1, "91a1 (r=1)", None, (-18, -14)),
        (197, 0.24498, 1, "197a1 (r=1)", None, (-58, 4)),
        (389, 0.1524, 2, "389a1 (r=2)", (-55, 4), (7, 4)),
        (5077, 0.4171, 3, "5077a1 (r=3, GGZ)", (-85, 4), None),
        (234446, 1.5043, 4, "234446a1 (r=4)", (-85, 4), None),
        (19047851, 30.289, 5, "19047851a1 (r=5)", (-85, 4), None),
        (5187563742, 150.42, 6, "Min r=6 Curve", (-80, 4), None),
        (1e116, 1e5, 28, "Elkies r=28 Curve (N~10^116, d~11.2)", (-165, -14), None)
    ]
    gln_data = [
        (1e30, 1.0, 10, 0, "GL_10 Cusp Form (n=10, d~6.47)", (-110, 8)),
        (1e80, 1.0, 100, 0, "Rankin-Selberg GL_10 x GL_10 (n=100, d~20.45)", (-180, -14))
    ]

    N_grid = np.logspace(np.log10(2.5), 125, 1000)
    d_GL1_typ = np.sqrt(np.maximum(0, (1.0 + np.log(np.sqrt(N_grid) / np.pi) / np.log(N_grid)) / L_inf))
    d_GL2_r1_typ = np.sqrt(np.maximum(0, 2.0 * (1.0 + np.log(0.15 * np.log(N_grid)) / np.log(N_grid)) / L_inf))
    d_GL2_r2_typ = np.sqrt(np.maximum(0, 3.0 * (1.0 + np.log(0.20 * np.log(N_grid)) / np.log(N_grid)) / L_inf))
    d_GL3_base = np.full_like(N_grid, np.sqrt(3.0 / L_inf))
    d_GL10_base = np.full_like(N_grid, np.sqrt(10.0 / L_inf))

    fig, ax = plt.subplots(figsize=(16, 9.5), dpi=150)

    def plot_layers(target_ax, is_inset=False):
        target_ax.axhspan(0, d_crit, color='#ffcccc', alpha=0.35, label='Localized Phase ($d < d_{crit} \\approx 2.045$)' if not is_inset else "")
        target_ax.axhspan(d_crit, 3.0, color='#d9f2e6', alpha=0.35, label='Diffusive Phase ($d > d_{crit}$)' if not is_inset else "")
        target_ax.axhspan(3.0, 23.0, color='#cce6ff', alpha=0.28, label='Topological / Higher-Rank & $GL_n$ Barrier ($d > 3.0$)' if not is_inset else "")

        target_ax.axhline(d_crit, color='#d62728', linestyle='--', linewidth=2.0, zorder=4)
        target_ax.axhline(d_inf_GL1, color='#2ca02c', linestyle=':', linewidth=1.6, zorder=3)
        target_ax.axhline(d_Gauss, color='#e377c2', linestyle='--', linewidth=1.5, alpha=0.9, zorder=4,
                           label='Gaussian Limit $d_{\\mathrm{Gauss}} = \\sqrt{2\\pi} \\approx 2.5066$' if not is_inset else "")

        target_ax.plot(N_grid, d_GL1_typ, color='#2ca02c', linewidth=1.8, alpha=0.85, label='$GL_1$ Typical Scaling ($d_\\infty \\approx 2.5044$)' if not is_inset else "")
        target_ax.plot(N_grid, d_GL2_r1_typ, color='#1f77b4', linewidth=1.8, linestyle='-.', alpha=0.85, label='$GL_2$ Rank $r=1$ Scaling' if not is_inset else "")
        target_ax.plot(N_grid, d_GL2_r2_typ, color='#9467bd', linewidth=1.8, linestyle='-.', alpha=0.85, label='$GL_2$ Rank $r=2$ Scaling' if not is_inset else "")
        if not is_inset:
            target_ax.plot(N_grid, d_GL3_base, color='#ff7f0e', linewidth=1.8, linestyle='--', alpha=0.9, label='$GL_3$ Baseline ($n=3 \\Rightarrow d \\approx 3.54$)')
            target_ax.plot(N_grid, d_GL10_base, color='#8c564b', linewidth=1.8, linestyle='--', alpha=0.85, label='$GL_{10}$ Baseline ($n=10 \\Rightarrow d \\approx 6.47$)')

        sz = 55 if not is_inset else 65
        for data, col, mk, lbl_leg, txt_col, fw in [
            (imag_quad_data, '#006600', 'o', '$GL_1$ Imaginary Quad ($D<0$)', '#004d00', 'bold'),
            (real_quad_data, '#009999', 's', '$GL_1$ Real Quad ($D>0$)', '#006666', 'normal'),
            (cubic_data, '#e68a00', '^', '$GL_1$ Real Cubic ($n=3$)', '#995c00', 'normal')
        ]:
            xs = [N for N, _, _, _, _ in data]
            ys = [calc_d_eff(N, I_val, n=1, r=0) for N, I_val, _, _, _ in data]
            target_ax.scatter(xs, ys, color=col, marker=mk, s=sz, zorder=6, edgecolors='white', linewidth=0.8, label=lbl_leg if not is_inset else "")
            for (N, _, lbl, off_m, off_i), y in zip(data, ys):
                off = off_i if is_inset else off_m
                if lbl and off:
                    target_ax.annotate(lbl, (N, y), xytext=off, textcoords='offset points', fontsize=7.5, color=txt_col, fontweight=fw)

        x_ec = [N for N, _, _, _, _, _ in elliptic_data]
        y_ec = [calc_d_eff(N, I_val, n=1, r=r) for N, I_val, r, _, _, _ in elliptic_data]
        target_ax.scatter(x_ec, y_ec, color='#0033cc', marker='D', s=sz+10, zorder=7, edgecolors='white', linewidth=0.9, label='$GL_2$ Elliptic Curves ($r=0..28$)' if not is_inset else "")
        for (N, _, r, lbl, off_m, off_i), y in zip(elliptic_data, y_ec):
            off = off_i if is_inset else off_m
            if lbl and off:
                target_ax.annotate(lbl, (N, y), xytext=off, textcoords='offset points', fontsize=7.8, color='#002288', fontweight='bold')

        if not is_inset:
            x_gln = [C for C, _, _, _, _, _ in gln_data]
            y_gln = [calc_d_eff(C, I_val, n=n, r=r) for C, I_val, n, r, _, _ in gln_data]
            target_ax.scatter(x_gln, y_gln, color='#cc0066', marker='*', s=160, zorder=8, edgecolors='white', linewidth=0.9, label='Higher $GL_n$ & Rankin-Selberg')
            for (C, _, _, _, lbl, off_m), y in zip(gln_data, y_gln):
                target_ax.annotate(lbl, (C, y), xytext=off_m, textcoords='offset points', fontsize=8.2, color='#99004d', fontweight='bold')

    plot_layers(ax, is_inset=False)
    ax.set_xscale('log')
    ax.set_xlim(1.5, 1e125)
    ax.set_ylim(0.8, 22.5)
    ax.set_xlabel('Analytic Conductor $C(\\pi)$ / Discriminant $|D|$ (Log Scale up to $10^{125}$)', fontsize=12, fontweight='bold', labelpad=8)
    ax.set_ylabel('Effective Dimension $d_{\\mathrm{eff}}$ / $d_{\\mathrm{rank}}$', fontsize=12, fontweight='bold', labelpad=8)
    ax.set_title('Ultra-Scale Unified Phase Diagram ($GL_1 \\to GL_2(r \\leq 28) \\to GL_n$ & Gaussian Limit $\\sqrt{2\\pi}$)', fontsize=14, fontweight='bold', pad=14)
    ax.legend(loc='upper left', fontsize=7.8, frameon=True, facecolor='white', framealpha=0.95, ncol=2)

    axins = ax.inset_axes([0.38, 0.34, 0.44, 0.40])
    plot_layers(axins, is_inset=True)
    axins.text(7.8, d_Gauss + 0.025, '$d_{\\mathrm{Gauss}} = \\sqrt{2\\pi} \\approx 2.5066$', color='#b83b92', fontweight='bold', fontsize=7.8)
    axins.set_xscale('log')
    axins.set_xlim(7, 1300)
    axins.set_ylim(1.75, 2.68)
    axins.set_title('Zoom-in: Critical Transition Zone ($7 \\leq N \\leq 1300,\\ 1.75 \\leq d \\leq 2.68$)', fontsize=8.8, fontweight='bold', pad=5,
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.95))
    axins.tick_params(labelsize=8.0)
    axins.grid(True, which='both', linestyle=':', alpha=0.6)
    mark_inset(ax, axins, loc1=2, loc2=3, fc="none", ec="0.45", linestyle="--", linewidth=1.1)

    plt.savefig('siegel_zero_phase_diagram_ultrascale.png', dpi=300, bbox_inches='tight', pad_inches=0.2)
    plt.close(fig)
    print("Saved: siegel_zero_phase_diagram_ultrascale.png")

# ==============================================================================
# Figure 2: GL_2 楕円曲線シフト検証 (elliptic_curves_real_data_phase_shift.png)
# ==============================================================================
def generate_figure_2():
    real_curves = [
        ("37a1 (Min r=1)",      37,     1, 0.05111, "#ff0055", (-45, -36)),
        ("43a1 (r=1)",          43,     1, 0.06282, "#ff6600", (-22, -48)),
        ("53a1 (r=1)",          53,     1, 0.08857, "#d98a00", (-25, 16)),
        ("389a1 (Min r=2)",     389,    2, 0.15246, "#9467bd", (-18, 16)),
        ("5077a1 (Min r=3)",    5077,   3, 0.41714, "#0088cc", (-95, 14)),
        ("681b1 (r=0, Sha=9)",  681,    0, 9.00000, "#2ca02c", (16, 12)),
    ]

    C_ref = 1000.0
    I_small = 1e-2
    factor_small = 1.0 + np.log(I_small) / np.log(C_ref)

    n_vals = np.linspace(1, 10, 300)
    r_vals = np.linspace(0, 8, 300)
    N_grid, R_grid = np.meshgrid(n_vals, r_vals)
    W_grid = N_grid * (1.0 + R_grid)

    d_eff_base = np.sqrt(W_grid / L_inf)
    d_eff_small = np.sqrt(W_grid * factor_small / L_inf)
    log10_P_small = -(W_grid**2) * np.log10(C_ref * I_small)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16.5, 7.5), dpi=150)

    levels_P = np.linspace(-100, 0, 26)
    cf = ax1.contourf(N_grid, R_grid, log10_P_small, levels=levels_P, cmap='magma_r', extend='min')
    cbar = fig.colorbar(cf, ax=ax1, pad=0.02)
    cbar.set_label('$\\log_{10} P_{\\mathrm{tunnel}}^{(GL_n)}$  Background ($C=10^3, I=10^{-2}$)', fontsize=11, fontweight='bold')

    ax1.contour(N_grid, R_grid, d_eff_base, levels=[d_crit], colors=['#00ffcc'], linestyles=['--'], linewidths=[1.6], alpha=0.7)
    ax1.contour(N_grid, R_grid, d_eff_base, levels=[3.0], colors=['#ffffff'], linestyles=['--'], linewidths=[1.6], alpha=0.7)
    ax1.contour(N_grid, R_grid, d_eff_small, levels=[d_crit], colors=['#00ffcc'], linestyles=['-'], linewidths=[2.6])
    ax1.contour(N_grid, R_grid, d_eff_small, levels=[3.0], colors=['#ffffff'], linestyles=['-'], linewidths=[2.6])

    r_annotations = {
        0: ("681b1 ($\\mathrm{Sha}=9$): $d_{\\mathrm{eff}}=3.35$ (>3D Barrier)", (18, 2), "#2ca02c"),
        1: ("37a1 ($d=1.74 < d_{crit}$)\n43a1 ($d=2.16$), 53a1 ($d=2.55$)", (18, -6), "#ff3366"),
        2: ("389a1 (Min $r=2$): $d_{\\mathrm{eff}}=4.13$ (from 5.01)", (18, -4), "#cc99ff"),
        3: ("5077a1 (Min $r=3$): $d_{\\mathrm{eff}}=5.48$ (from 5.78)", (18, -4), "#66ccff"),
    }
    for r_val, (txt, off, col) in r_annotations.items():
        ax1.scatter([2.0], [r_val], color=col, edgecolor='black', s=110, marker='*', zorder=8)
        ax1.annotate(txt, (2.0, r_val), xytext=off, textcoords='offset points',
                     fontsize=8.0, fontweight='bold', color='white',
                     bbox=dict(boxstyle='round,pad=0.25', facecolor='black', alpha=0.78, edgecolor=col),
                     arrowprops=dict(arrowstyle='-', color=col, lw=1.0))

    ax1.plot([], [], color='#00ffcc', ls='--', lw=1.6, label='$d_{crit}$ ($I=1.0$)')
    ax1.plot([], [], color='#00ffcc', ls='-', lw=2.6, label='$d_{crit}$ Shifted ($I=10^{-2}$)')
    ax1.plot([], [], color='#ffffff', ls='-', lw=2.6, label='3D Barrier ($I=10^{-2}$)')
    ax1.scatter([], [], color='#ff3366', edgecolor='black', s=110, marker='*', label='Real Elliptic Curves ($GL_2$)')
    ax1.legend(loc='upper right', fontsize=8.5, frameon=True, facecolor='black', framealpha=0.8, labelcolor='white')

    ax1.set_title('Real Elliptic Curves on Shifted Phase Diagram ($GL_2, n=2$)', fontsize=12.0, fontweight='bold', pad=10)
    ax1.set_xlabel('Degree $n$ of Automorphic Representation ($GL_n$)', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Zero Multiplicity $r_{\\mathrm{zero}}$ (Algebraic Rank $r$)', fontsize=11, fontweight='bold')
    ax1.set_xlim(1, 10)
    ax1.set_ylim(-0.3, 8)

    d_axis = np.linspace(0.5, 6.5, 400)
    for W_c, lbl_w, ls in [(2, '$W=2$ ($r=0$)', '--'), (4, '$W=4$ ($r=1$)', '-'), (6, '$W=6$ ($r=2$)', ':'), (8, '$W=8$ ($r=3$)', '-.')]:
        ax2.plot(d_axis, -W_c * L_inf * (d_axis**2) * np.log10(C_ref), color='gray', alpha=0.45, linestyle=ls, linewidth=1.4, label=lbl_w)

    for lbl, N_val, r_val, I_val, col, off in real_curves:
        W_val = 2 * (1 + r_val)
        d_ideal = np.sqrt(W_val / L_inf)
        log10_P_ideal = -(W_val**2) * np.log10(N_val)

        factor_real = 1.0 + np.log(I_val) / np.log(N_val)
        d_real = np.sqrt(max(0.0, W_val * factor_real / L_inf))
        log10_P_real = -(W_val**2) * np.log10(N_val * I_val)

        ax2.scatter([d_ideal], [log10_P_ideal], color=col, alpha=0.45, s=55, edgecolor='black', zorder=5)
        ax2.scatter([d_real], [log10_P_real], color=col, s=160, marker='*', edgecolor='black', linewidth=0.9, zorder=8)
        ax2.annotate('', xy=(d_real, log10_P_real), xytext=(d_ideal, log10_P_ideal),
                     arrowprops=dict(arrowstyle='->', color=col, lw=1.8, ls='--'))
        ax2.annotate(f"{lbl}\n($d={d_real:.2f}$)", (d_real, log10_P_real), xytext=off, textcoords='offset points',
                     fontsize=8.0, fontweight='bold', color=col,
                     bbox=dict(boxstyle='round,pad=0.22', facecolor='white', alpha=0.92, edgecolor=col),
                     arrowprops=dict(arrowstyle='->', color=col, lw=1.1))

    ax2.axvline(d_crit, color='#d62728', linestyle='--', linewidth=1.8, alpha=0.85)
    ax2.text(d_crit + 0.05, -75, '$d_{crit} \\approx 2.045$', color='#d62728', fontweight='bold', fontsize=9)
    ax2.axvline(3.0, color='#3366cc', linestyle='--', linewidth=1.8, alpha=0.85)
    ax2.text(3.05, -75, '3D / Kakeya Barrier ($d=3.0$)', color='#3366cc', fontweight='bold', fontsize=9)

    ax2.axvspan(0.5, d_crit, color='#ffcccc', alpha=0.25, label='Localized Phase ($d < d_{crit}$)')
    ax2.axvspan(d_crit, 3.0, color='#d9f2e6', alpha=0.25, label='Diffusive Phase')
    ax2.axvspan(3.0, 6.5, color='#cce6ff', alpha=0.20, label='Higher-Rank / $GL_n$ Barrier')

    ax2.scatter([], [], color='gray', alpha=0.5, s=55, edgecolor='black', label='Ideal Base ($I=1.0$)')
    ax2.scatter([], [], color='black', s=150, marker='*', label='Actual Curve ($I = \\mathrm{Reg}\\cdot|\\mathrm{Sha}|$)')

    ax2.set_title('Real Elliptic Curves: Shift from Ideal ($I=1$) to Actual ($I=\\mathrm{Reg}\\cdot|\\mathrm{Sha}|$)', fontsize=11.5, fontweight='bold', pad=10)
    ax2.set_xlabel('Effective Dimension $d_{\\mathrm{eff}}^{(GL_2)}$', fontsize=11, fontweight='bold')
    ax2.set_ylabel('$\\log_{10} P_{\\mathrm{tunnel}}^{(GL_2)} = -W^2 \\log_{10}(N \\cdot I)$', fontsize=11, fontweight='bold')
    ax2.set_xlim(1.0, 6.3)
    ax2.set_ylim(-235, 22)
    ax2.legend(loc='lower left', fontsize=8.0, frameon=True, facecolor='white', framealpha=0.92, ncol=2)

    plt.tight_layout()
    plt.savefig('elliptic_curves_real_data_phase_shift.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print("Saved: elliptic_curves_real_data_phase_shift.png")


# ==============================================================================
# Figure 3: GL_3 尖点形式・総実3次体検証 (gl3_cubic_verification_phase_diagram.png)
# ==============================================================================
def generate_figure_3():
    higher_dim_data = [
        ("Cubic D=49 (Min)",      3, 49,    0, 0.5255, "^", "#e6550d", (-10, -34)),
        ("Cubic D=81",            3, 81,    0, 0.8493, "^", "#fd8d3c", (-68, -22)),
        ("Cubic D=148",           3, 148,   0, 1.6624, "^", "#d9781e", (18, -6)),
        ("Cubic D=1957",          3, 1957,  0, 9.3125, "^", "#d94801", (16, -6)),
        ("Sym^2(11a1) [GL_3]",    3, 121,   0, 1.4820, "P", "#006d2c", (14, 12)),
        ("GL_3 Maass (Lvl 1)",    3, 2560,  0, 1.0000, "P", "#238b45", (16, -12)),
        ("37a1 (GL_2, r=1)",      2, 37,    1, 0.05111, "*", "#ff0055", (-10, 14)),
        ("389a1 (GL_2, r=2)",     2, 389,   2, 0.15246, "*", "#9467bd", (-85, 16)),
    ]

    C_ref = 1000.0
    I_small = 1e-2
    factor_small = 1.0 + np.log(I_small) / np.log(C_ref)

    n_vals = np.linspace(1, 10, 300)
    r_vals = np.linspace(0, 8, 300)
    N_grid, R_grid = np.meshgrid(n_vals, r_vals)
    W_grid = N_grid * (1.0 + R_grid)

    d_eff_base = np.sqrt(W_grid / L_inf)
    d_eff_small = np.sqrt(W_grid * factor_small / L_inf)
    log10_P_small = -(W_grid**2) * np.log10(C_ref * I_small)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16.5, 7.5), dpi=150)

    levels_P = np.linspace(-100, 0, 26)
    cf = ax1.contourf(N_grid, R_grid, log10_P_small, levels=levels_P, cmap='magma_r', extend='min')
    cbar = fig.colorbar(cf, ax=ax1, pad=0.02)
    cbar.set_label('$\\log_{10} P_{\\mathrm{tunnel}}^{(GL_n)}$  Background ($C=10^3, I=10^{-2}$)', fontsize=11, fontweight='bold')

    ax1.contour(N_grid, R_grid, d_eff_base, levels=[d_crit], colors=['#00ffcc'], linestyles=['--'], linewidths=[1.6], alpha=0.7)
    ax1.contour(N_grid, R_grid, d_eff_base, levels=[3.0], colors=['#ffffff'], linestyles=['--'], linewidths=[1.6], alpha=0.7)
    ax1.contour(N_grid, R_grid, d_eff_small, levels=[d_crit], colors=['#00ffcc'], linestyles=['-'], linewidths=[2.6])
    ax1.contour(N_grid, R_grid, d_eff_small, levels=[3.0], colors=['#ffffff'], linestyles=['-'], linewidths=[2.6])

    ax1.scatter([3.0], [0.0], color='#fd8d3c', edgecolor='black', s=130, marker='^', zorder=8)
    ax1.annotate("Cubic Fields ($n=3, r=0$):\nD=49 ($d=1.87 < d_{crit}$)\nD=81 ($d=2.00$), D=148 ($d=2.15$)",
                 (3.0, 0.0), xytext=(16, 16), textcoords='offset points',
                 fontsize=8.0, fontweight='bold', color='white',
                 bbox=dict(boxstyle='round,pad=0.25', facecolor='black', alpha=0.78, edgecolor='#fd8d3c'))

    ax1.annotate("$GL_3$ Cusp / $\\mathrm{Sym}^2(11a1)$:\n$d_{\\mathrm{eff}} \\approx 3.54 \\sim 3.68$ (>3D Barrier)",
                 (3.0, 0.0), xytext=(16, 72), textcoords='offset points',
                 fontsize=8.0, fontweight='bold', color='white',
                 bbox=dict(boxstyle='round,pad=0.25', facecolor='black', alpha=0.78, edgecolor='#238b45'),
                 arrowprops=dict(arrowstyle='->', color='#238b45', lw=1.5))

    ax1.scatter([2.0], [1.0], color='#ff0055', edgecolor='black', s=110, marker='*', zorder=8)
    ax1.annotate("37a1 ($GL_2, r=1$): $d=1.74$", (2.0, 1.0), xytext=(14, -4), textcoords='offset points',
                 fontsize=8.0, fontweight='bold', color='white',
                 bbox=dict(boxstyle='round,pad=0.2', facecolor='black', alpha=0.75, edgecolor='#ff0055'))

    ax1.plot([], [], color='#00ffcc', ls='--', lw=1.6, label='$d_{crit}$ ($I=1.0$)')
    ax1.plot([], [], color='#00ffcc', ls='-', lw=2.6, label='$d_{crit}$ Shifted ($I=10^{-2}$)')
    ax1.plot([], [], color='#ffffff', ls='-', lw=2.6, label='3D Barrier ($I=10^{-2}$)')
    ax1.scatter([], [], color='#fd8d3c', edgecolor='black', s=110, marker='^', label='Totally Real Cubic ($n=3$)')
    ax1.legend(loc='upper right', fontsize=8.5, frameon=True, facecolor='black', framealpha=0.8, labelcolor='white')

    ax1.set_title('Higher-Degree $GL_3$ & Totally Real Cubic Fields on Phase Diagram', fontsize=12.0, fontweight='bold', pad=10)
    ax1.set_xlabel('Degree $n$ of Automorphic Representation ($GL_n$)', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Zero Multiplicity $r_{\\mathrm{zero}}$ at Central Point $s=1$', fontsize=11, fontweight='bold')
    ax1.set_xlim(1, 10)
    ax1.set_ylim(-0.3, 8)

    d_axis = np.linspace(0.5, 6.5, 400)
    for W_c, lbl_w, ls in [(3, '$W=3$ ($GL_3, r=0$)', '-.'), (4, '$W=4$ ($GL_2, r=1$)', '-'), (6, '$W=6$ ($GL_2, r=2$)', ':')]:
        ax2.plot(d_axis, -W_c * L_inf * (d_axis**2) * np.log10(C_ref), color='gray', alpha=0.45, linestyle=ls, linewidth=1.4, label=lbl_w)

    for lbl, n_deg, N_val, r_val, I_val, mk, col, off in higher_dim_data:
        W_val = n_deg * (1 + r_val)
        n_scale = 1.0 if "Cubic" in lbl else float(n_deg)
        W_eff = n_scale * (1 + r_val)

        d_ideal = np.sqrt(W_eff / L_inf)
        log10_P_ideal = -(W_val**2) * np.log10(N_val)

        factor_real = 1.0 + np.log(I_val) / np.log(N_val)
        d_real = np.sqrt(max(0.0, W_eff * factor_real / L_inf))
        log10_P_real = -(W_val**2) * np.log10(N_val * I_val)

        ax2.scatter([d_ideal], [log10_P_ideal], color=col, alpha=0.4, s=50, edgecolor='black', zorder=5)
        ax2.scatter([d_real], [log10_P_real], color=col, s=145, marker=mk, edgecolor='black', linewidth=0.9, zorder=8)
        ax2.annotate('', xy=(d_real, log10_P_real), xytext=(d_ideal, log10_P_ideal),
                     arrowprops=dict(arrowstyle='->', color=col, lw=1.6, ls='--'))
        ax2.annotate(f"{lbl}\n($d={d_real:.2f}$)", (d_real, log10_P_real), xytext=off, textcoords='offset points',
                     fontsize=8.0, fontweight='bold', color=col,
                     bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.92, edgecolor=col),
                     arrowprops=dict(arrowstyle='->', color=col, lw=1.1))

    ax2.axvline(d_crit, color='#d62728', linestyle='--', linewidth=1.8, alpha=0.85)
    ax2.text(d_crit + 0.05, -75, '$d_{crit} \\approx 2.045$', color='#d62728', fontweight='bold', fontsize=9)
    ax2.axvline(3.0, color='#3366cc', linestyle='--', linewidth=1.8, alpha=0.85)
    ax2.text(3.05, -75, '3D / Kakeya Barrier ($d=3.0$)', color='#3366cc', fontweight='bold', fontsize=9)

    ax2.axvspan(0.5, d_crit, color='#ffcccc', alpha=0.25, label='Localized Phase ($d < d_{crit}$)')
    ax2.axvspan(d_crit, 3.0, color='#d9f2e6', alpha=0.25, label='Diffusive Phase')
    ax2.axvspan(3.0, 6.5, color='#cce6ff', alpha=0.20, label='Higher-Rank / $GL_n$ Barrier')

    ax2.scatter([], [], color='gray', alpha=0.45, s=50, edgecolor='black', label='Ideal Base ($I=1.0$)')
    ax2.scatter([], [], color='#fd8d3c', s=130, marker='^', edgecolor='black', label='Totally Real Cubic ($n=3$)')
    ax2.scatter([], [], color='#006d2c', s=130, marker='P', edgecolor='black', label='$GL_3$ Automorphic Form')

    ax2.set_title('Verification of $GL_3$ & Cubic Fields: Localized Cubic vs. Barrier-Protected $GL_3$ Cusp Forms', fontsize=11.0, fontweight='bold', pad=10)
    ax2.set_xlabel('Effective Dimension $d_{\\mathrm{eff}}^{(GL_n)}$', fontsize=11, fontweight='bold')
    ax2.set_ylabel('$\\log_{10} P_{\\mathrm{tunnel}}^{(GL_n)} = -W^2 \\log_{10}(C \\cdot I)$', fontsize=11, fontweight='bold')
    ax2.set_xlim(1.3, 4.35)
    ax2.set_ylim(-115, 12)
    ax2.legend(loc='lower left', fontsize=7.8, frameon=True, facecolor='white', framealpha=0.92, ncol=2)

    plt.tight_layout()
    plt.savefig('gl3_cubic_verification_phase_diagram.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print("Saved: gl3_cubic_verification_phase_diagram.png")


if __name__ == "__main__":
    generate_figure_1()
    generate_figure_2()
    generate_figure_3()
    print("All 3 figures generated successfully.")
