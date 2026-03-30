"""
pk_engine.py — Pharmacokinetic Analysis Engine
All original drawing logic, color constants, and metric calculations preserved
exactly from Bioavailability.py. Only adaptation: inputs come from a dict
instead of input() calls, and the figure is returned as JPEG bytes.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
import numpy as np
from scipy.interpolate import make_interp_spline
from scipy.stats import linregress
from scipy.integrate import trapezoid
import io

# ── ORIGINAL CONSTANT (unchanged) ────────────────────────────────────────────
CRAYON_BG = '#FFFEF5'   # warm off-white paper


# ── ORIGINAL HELPERS (unchanged) ─────────────────────────────────────────────

def parse_input_string(input_str):
    clean = input_str.replace(',', ' ')
    return [float(x) for x in clean.split() if x.strip()]


def get_action_times(t_smooth, c_smooth, mec):
    if len(t_smooth) == 0 or mec <= 0:
        return 0, 0, 0
    above = np.where(c_smooth >= mec)[0]
    if len(above) < 2:
        return 0, 0, 0
    return t_smooth[above[0]], t_smooth[above[-1]], t_smooth[above[-1]] - t_smooth[above[0]]


def analyze_pk_data(time_raw, conc_raw, dose, mec):
    """
    AUC  = Σ [(C₁+C₂)/2 · (t₂−t₁)]   trapezoidal rule
    Ke   = |slope|  ln(C) vs t  terminal log-linear regression
    t½   = 0.693 / Ke
    Cl   = Dose / AUC
    Vd   = Cl  / Ke
    F    = (AUC_test · Dose_ref) / (AUC_ref · Dose_test) × 100
    """
    if len(time_raw) == 0:
        return {}
    data  = sorted(zip(time_raw, conc_raw))
    t_arr = np.array([x[0] for x in data])
    c_arr = np.array([x[1] for x in data])

    if len(t_arr) > 3:
        try:
            spl      = make_interp_spline(t_arr, c_arr, k=3)
            t_smooth = np.linspace(t_arr.min(), t_arr.max(), 1200)
            c_smooth = np.clip(spl(t_smooth), 0, None)
        except Exception:
            t_smooth = np.linspace(t_arr.min(), t_arr.max(), 1200)
            c_smooth = np.interp(t_smooth, t_arr, c_arr)
    else:
        t_smooth = np.linspace(t_arr.min(), t_arr.max(), 1200)
        c_smooth = np.interp(t_smooth, t_arr, c_arr)

    auc  = float(trapezoid(c_arr, t_arr))
    aumc = float(trapezoid(c_arr * t_arr, t_arr))
    mrt  = aumc / auc if auc > 0 else 0.0
    c_max = float(np.max(c_smooth))
    t_max = float(t_smooth[np.argmax(c_smooth)])

    idx_max = int(np.argmax(c_arr))
    ke = 0.0
    if idx_max < len(c_arr) - 1:
        t_term = t_arr[idx_max:]
        c_safe = np.where(c_arr[idx_max:] > 0, c_arr[idx_max:], 1e-9)
        slope, *_ = linregress(t_term, np.log(c_safe))
        ke = abs(slope)

    t_half = 0.693 / ke if ke > 0 else 0.0
    cl     = dose  / auc if auc > 0 else 0.0
    vd     = cl    / ke  if ke  > 0 else 0.0

    return dict(t_arr=t_arr, c_arr=c_arr, t_smooth=t_smooth, c_smooth=c_smooth,
                auc=auc, aumc=aumc, mrt=mrt, c_max=c_max, t_max=t_max,
                ke=ke, t_half=t_half, cl=cl, vd=vd)


def draw_gradient_bars(ax, bar_labels, bar_values, cmap_name, title_text, edge):
    """Horizontal gradient bars — dark→light per bar, no overlap."""
    n      = len(bar_labels)
    maxval = max(bar_values) if bar_values else 1
    ax.set_xlim(0, maxval * 1.40)
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_facecolor(CRAYON_BG)
    ax.set_yticks(range(n))
    ax.set_yticklabels(bar_labels, fontsize=8.5, weight='bold')
    ax.tick_params(axis='x', labelsize=7.5)
    ax.set_title(title_text, fontsize=9, weight='bold', color='#1A237E', pad=5)
    for sp in ax.spines.values():
        sp.set_linewidth(1.3)
        sp.set_color('#5D4037')

    cmap   = matplotlib.colormaps.get_cmap(cmap_name)
    bar_h  = 0.52
    shades = [0.82, 0.56, 0.32]
    n_strips = 120

    for i, (lbl, val) in enumerate(zip(bar_labels, bar_values)):
        if val <= 0:
            continue
        s = shades[i % len(shades)]
        for j in range(n_strips):
            x0 = val * j / n_strips
            x1 = val * (j + 1) / n_strips
            frac = j / (n_strips - 1)
            color_val = s * (1.0 - frac * 0.72)
            ax.fill_betweenx([i - bar_h/2, i + bar_h/2], x0, x1,
                             color=cmap(color_val), zorder=3, linewidth=0)
        rect = mpatches.FancyBboxPatch(
                   (0, i - bar_h/2), val, bar_h,
                   boxstyle='square,pad=0',
                   linewidth=1.3, edgecolor=edge, facecolor='none', zorder=4)
        ax.add_patch(rect)
        ax.text(val + maxval * 0.018, i, f'{val:.2f}',
                va='center', fontsize=8.0, weight='bold', color='#212121')


# ── MAIN GENERATION FUNCTION ─────────────────────────────────────────────────

def generate_pk_image(params: dict) -> bytes:
    """
    Generate the PK analysis chart and return it as JPEG bytes.

    params keys:
        drug_name       : str
        mec             : float  (0 = skip)
        msc             : float  (0 = skip)
        route1_name     : str
        dose1           : float
        conc1_str       : str   (space/comma separated)
        time1_str       : str   (space/comma separated)
        has_second_route: bool
        route2_name     : str   (optional)
        dose2           : float (optional)
        conc2_str       : str   (optional)
        time2_str       : str   (optional)
    """
    drug_name = params['drug_name']
    mec       = float(params['mec'])
    msc       = float(params['msc'])

    route1_name = params['route1_name']
    dose1       = float(params['dose1'])
    conc1_raw   = parse_input_string(params['conc1_str'])
    time1_raw   = parse_input_string(params['time1_str'])

    show_mec_msc = (mec > 0) and (msc > 0)

    pk1 = analyze_pk_data(time1_raw, conc1_raw, dose1, mec)
    onset1, _, duration1 = get_action_times(pk1['t_smooth'], pk1['c_smooth'], mec)
    loading_dose1     = pk1['vd'] * mec if mec > 0 else 0
    maintenance_rate1 = pk1['cl'] * mec if mec > 0 else 0

    has_second_route = params.get('has_second_route', False)
    pk2 = None
    route2_name = "—"
    dose2 = 0
    onset2 = duration2 = 0
    bioavailability = None
    loading_dose2 = maintenance_rate2 = 0

    if has_second_route:
        route2_name = params['route2_name']
        dose2       = float(params['dose2'])
        conc2_raw   = parse_input_string(params['conc2_str'])
        time2_raw   = parse_input_string(params['time2_str'])

        pk2 = analyze_pk_data(time2_raw, conc2_raw, dose2, mec)
        onset2, _, duration2 = get_action_times(pk2['t_smooth'], pk2['c_smooth'], mec)
        loading_dose2     = pk2['vd'] * mec if mec > 0 else 0
        maintenance_rate2 = pk2['cl'] * mec if mec > 0 else 0

        f_ratio = 0.0
        if pk1['auc'] > 0 and pk2['auc'] > 0:
            if "iv" in route2_name.lower():
                f_ratio = (pk1['auc'] * dose2) / (pk2['auc'] * dose1)
            elif "iv" in route1_name.lower():
                f_ratio = (pk2['auc'] * dose1) / (pk1['auc'] * dose2)
            else:
                f_ratio = (pk2['auc'] * dose1) / (pk1['auc'] * dose2)
        bioavailability = f_ratio * 100

    # Toxicity
    if show_mec_msc:
        is_toxic_1 = pk1['c_max'] > msc
        status_1   = "Toxic  ☠" if is_toxic_1 else "Safe  ♥"
        is_toxic_2 = (pk2['c_max'] > msc) if (has_second_route and pk2) else False
        status_2   = ("Toxic  ☠" if is_toxic_2 else "Safe  ♥") if has_second_route else "—"
    else:
        is_toxic_1 = is_toxic_2 = False
        status_1 = "Safe  ♥"
        status_2 = "Safe  ♥" if has_second_route else "—"

    # Fixed colors (unchanged)
    COLOR_R1_LINE   = '#C2185B'
    COLOR_R1_FILL   = '#FFB3C6'
    COLOR_R2_LINE   = '#0277BD'
    COLOR_R2_FILL   = '#B3E5FC'
    color_drug_name = '#D4AC0D'
    t_color_r1      = '#D32F2F' if is_toxic_1 else '#1565C0'
    t_color_r2      = '#D32F2F' if is_toxic_2 else '#1565C0'

    # Figure (unchanged)
    fig = plt.figure(figsize=(8.27, 11.69), facecolor=CRAYON_BG)
    gs  = fig.add_gridspec(3, 2,
                           height_ratios=[3.8, 4.6, 1.55],
                           hspace=0.28, wspace=0.22,
                           left=0.08, right=0.97,
                           top=0.97, bottom=0.07)

    ax1      = fig.add_subplot(gs[0, :])
    ax_table = fig.add_subplot(gs[1, :])
    ax2      = fig.add_subplot(gs[2, 0])
    ax3      = fig.add_subplot(gs[2, 1])

    # ── GRAPH (unchanged) ─────────────────────────────────────────────────────
    ax1.set_facecolor(CRAYON_BG)
    for sp in ax1.spines.values():
        sp.set_linewidth(2.2)
        sp.set_color('#3E2723')

    t1s, c1s = pk1['t_smooth'], pk1['c_smooth']
    t1_max, c1_max = pk1['t_max'], pk1['c_max']
    t_span = t1s.max() - t1s.min()

    ax1.axvspan(t1s.min(), t1_max,  color='#E8F5E9', alpha=0.60, zorder=0)
    ax1.axvspan(t1_max,   t1s.max(), color='#FFF9C4', alpha=0.60, zorder=0)

    ax1.text((t1s.min() + t1_max)/2, 0.89, "◀  Absorption",
             transform=ax1.get_xaxis_transform(),
             ha='center', va='bottom', fontsize=9, weight='bold',
             color='#2E7D32', style='italic',
             path_effects=[pe.withStroke(linewidth=2, foreground=CRAYON_BG)])
    ax1.text((t1_max + t1s.max())/2, 0.89, "Elimination  ▶",
             transform=ax1.get_xaxis_transform(),
             ha='center', va='bottom', fontsize=9, weight='bold',
             color='#E65100', style='italic',
             path_effects=[pe.withStroke(linewidth=2, foreground=CRAYON_BG)])

    ax1.fill_between(t1s, c1s, color=COLOR_R1_FILL,
                     alpha=0.70, zorder=2, label=f'AUC₁  ({route1_name})')

    ax1.plot(t1s, c1s, color=COLOR_R1_LINE, lw=7.5, alpha=0.15, zorder=3)
    ax1.plot(t1s, c1s, color='white',       lw=4.8, alpha=0.20, zorder=3)
    ax1.plot(t1s, c1s, color='#5D1A2E',     lw=4.0, zorder=4)
    ax1.plot(t1s, c1s, color=COLOR_R1_LINE, lw=2.4, zorder=5,
             label=f'{route1_name}  Curve')

    ax1.vlines(t1_max, 0, c1_max, colors='#880E4F',
               linestyles='--', lw=1.8, zorder=6,
               label=f'Tmax₁  =  {t1_max:.2f} h')
    ax1.hlines(c1_max, t1s.min(), t1_max, colors='#880E4F',
               linestyles='--', lw=1.8, zorder=6)

    ax1.text(t1_max + t_span*0.012, c1_max*0.04,
             f'Tmax₁\n{t1_max:.2f} h',
             ha='left', va='bottom', fontsize=7.2, color='#880E4F',
             weight='bold', zorder=10,
             bbox=dict(facecolor='#FCE4EC', edgecolor='#880E4F',
                       boxstyle='round,pad=0.3', linewidth=1.3, alpha=0.92))

    ax1.plot(t1_max, c1_max, 'x', ms=11, color='#FF6D00',
             markeredgewidth=3.0, zorder=8,
             label=f'Cmax₁  =  {c1_max:.2f} mg/L')

    ax1.annotate(f'Cmax₁  {c1_max:.2f} mg/L',
                 xy=(t1_max, c1_max),
                 xytext=(t1_max + t_span*0.07, c1_max*1.08),
                 fontsize=7.8, color=COLOR_R1_LINE, weight='bold', zorder=10,
                 bbox=dict(facecolor='#FCE4EC', edgecolor=COLOR_R1_LINE,
                           boxstyle='round,pad=0.35', linewidth=1.3, alpha=0.95),
                 arrowprops=dict(arrowstyle='->', color=COLOR_R1_LINE, lw=1.4,
                                 connectionstyle='arc3,rad=-0.15'))

    if pk1['t_half'] > 0:
        c_half    = c1_max / 2.0
        t_at_half = t1_max + pk1['t_half']
        if t_at_half <= t1s.max():
            ax1.annotate('', xy=(t_at_half, c_half), xytext=(t1_max, c_half),
                         arrowprops=dict(arrowstyle='<->', color='#6A1B9A', lw=1.4))
            ax1.text((t1_max + t_at_half)/2, c_half*1.10,
                     f't½₁ = {pk1["t_half"]:.2f} h',
                     ha='center', va='bottom', fontsize=7.5,
                     color='#6A1B9A', weight='bold', zorder=10,
                     bbox=dict(facecolor='#EDE7F6', edgecolor='#6A1B9A',
                               boxstyle='round,pad=0.32', linewidth=1.3, alpha=0.95))

    if has_second_route and pk2:
        t2s, c2s  = pk2['t_smooth'], pk2['c_smooth']
        t2_max    = pk2['t_max']
        c2_max    = pk2['c_max']
        t2_span   = t2s.max() - t2s.min()

        ax1.fill_between(t2s, c2s, color=COLOR_R2_FILL,
                         alpha=0.55, zorder=2, label=f'AUC₂  ({route2_name})')

        ax1.plot(t2s, c2s, color=COLOR_R2_LINE, lw=7.5, alpha=0.15,
                 linestyle='--', zorder=3)
        ax1.plot(t2s, c2s, color='white', lw=4.8, alpha=0.20,
                 linestyle='--', zorder=3)
        ax1.plot(t2s, c2s, color='#0D2B45', lw=4.0,
                 linestyle='--', zorder=4)
        ax1.plot(t2s, c2s, color=COLOR_R2_LINE, lw=2.4,
                 linestyle='--', zorder=5, label=f'{route2_name}  Curve')

        ax1.vlines(t2_max, 0, c2_max, colors='#01579B',
                   linestyles='--', lw=1.8, zorder=6,
                   label=f'Tmax₂  =  {t2_max:.2f} h')
        ax1.hlines(c2_max, t2s.min(), t2_max, colors='#01579B',
                   linestyles='--', lw=1.8, zorder=6)

        _tmax2_ha  = 'right' if abs(t2_max - t1_max) < t_span * 0.12 else 'left'
        _tmax2_off = -t_span*0.012 if _tmax2_ha == 'right' else t_span*0.012
        _tmax2_yoff = c2_max * 0.04 if abs(c2_max - c1_max) > c1_max * 0.05 else c2_max * 0.18
        ax1.text(t2_max + _tmax2_off, _tmax2_yoff,
                 f'Tmax₂\n{t2_max:.2f} h',
                 ha=_tmax2_ha, va='bottom', fontsize=7.2, color='#01579B',
                 weight='bold', zorder=10,
                 bbox=dict(facecolor='#E3F2FD', edgecolor='#01579B',
                           boxstyle='round,pad=0.3', linewidth=1.3, alpha=0.92))

        ax1.plot(t2_max, c2_max, 'x', ms=11, color='#00BCD4',
                 markeredgewidth=3.0, zorder=8,
                 label=f'Cmax₂  =  {c2_max:.2f} mg/L')

        _c2_ydelta = -c2_max * 0.18 if c2_max > c1_max * 0.80 else c2_max * 0.10
        ax1.annotate(f'Cmax₂  {c2_max:.2f} mg/L',
                     xy=(t2_max, c2_max),
                     xytext=(t2_max - t2_span*0.10, c2_max + _c2_ydelta),
                     fontsize=7.8, color=COLOR_R2_LINE, weight='bold', zorder=10,
                     bbox=dict(facecolor='#E3F2FD', edgecolor=COLOR_R2_LINE,
                               boxstyle='round,pad=0.35', linewidth=1.3, alpha=0.95),
                     arrowprops=dict(arrowstyle='->', color=COLOR_R2_LINE, lw=1.4,
                                     connectionstyle='arc3,rad=0.15'))

        if pk2['t_half'] > 0:
            c_half2    = c2_max / 2.0
            t_at_half2 = t2_max + pk2['t_half']
            _thalf2_y  = c_half2 * 0.72
            if t_at_half2 <= t2s.max():
                ax1.annotate('', xy=(t_at_half2, _thalf2_y),
                             xytext=(t2_max, _thalf2_y),
                             arrowprops=dict(arrowstyle='<->', color='#0277BD', lw=1.4))
                ax1.text((t2_max + t_at_half2)/2, _thalf2_y * 1.10,
                         f't½₂ = {pk2["t_half"]:.2f} h',
                         ha='center', va='bottom', fontsize=7.5,
                         color='#0277BD', weight='bold', zorder=10,
                         bbox=dict(facecolor='#E3F2FD', edgecolor='#0277BD',
                                   boxstyle='round,pad=0.32', linewidth=1.3, alpha=0.95))

    if show_mec_msc:
        ax1.axhspan(mec, msc, color='#A5D6A7', alpha=0.18, zorder=0)

    all_c_max = max(c1_max, pk2['c_max'] if (has_second_route and pk2) else c1_max)
    ax1.set_ylim(bottom=0, top=all_c_max * 1.40)

    ax1.set_title(f"Plasma Drug Concentration  vs  Time  —  {drug_name}",
                  fontsize=12.5, weight='bold', pad=8, color='#1A237E',
                  path_effects=[pe.withStroke(linewidth=2, foreground=CRAYON_BG)])
    ax1.set_ylabel("Concentration (mg/L)", fontsize=10, weight='bold', color='#3E2723')
    ax1.set_xlabel("Time (hours)",          fontsize=10, weight='bold', color='#3E2723', labelpad=8)
    ax1.tick_params(axis='x', labelsize=9, pad=4)
    ax1.grid(True, alpha=0.20, linestyle='--', color='#BCAAA4', linewidth=0.8)

    if show_mec_msc:
        ax1.axhline(mec, color='#388E3C', lw=1.5, linestyle='-.', zorder=7,
                    label=f'MEC  =  {mec:.2f} mg/L')
        ax1.axhline(msc, color='#D32F2F', lw=1.5, linestyle='-.', zorder=7,
                    label=f'MSC  =  {msc:.2f} mg/L')
        ax1.text(t1s.max(), mec, f' MEC {mec:.2f}',
                 va='bottom', ha='right', fontsize=7.5, color='#388E3C',
                 weight='bold', zorder=8)
        ax1.text(t1s.max(), msc, f' MSC {msc:.2f}',
                 va='bottom', ha='right', fontsize=7.5, color='#D32F2F',
                 weight='bold', zorder=8)

    legend = ax1.legend(loc="upper right", fontsize=7.5, fancybox=True,
                        framealpha=0.96, edgecolor='#5D4037', shadow=True,
                        ncol=1, borderpad=0.7)
    legend.get_frame().set_linewidth(1.5)
    legend.get_frame().set_facecolor(CRAYON_BG)

    ax1.text(0.988, 0.012, "  @NautamKerai  ",
             transform=ax1.transAxes, va='bottom', ha='right',
             fontsize=8.5, weight='bold', color='white', zorder=15,
             clip_on=False,
             bbox=dict(facecolor='#1565C0', edgecolor='#FFFFFF',
                       boxstyle='round,pad=0.50', linewidth=2.2, alpha=0.97))

    # ── TABLE (unchanged) ─────────────────────────────────────────────────────
    ax_table.set_xlim(0, 1)
    ax_table.set_ylim(0, 1)
    ax_table.set_xticks([])
    ax_table.set_yticks([])
    ax_table.set_facecolor(CRAYON_BG)
    for sp in ax_table.spines.values():
        sp.set_visible(True)
        sp.set_color('#3E2723')
        sp.set_linewidth(1.8)

    r2_label = route2_name.upper() if has_second_route else "—"

    ax_table.text(0.5, 0.972,
                  f"COMPARATIVE PHARMACOKINETIC ANALYSIS  :  [ {drug_name.upper()} ]",
                  ha='center', va='center', fontsize=10.5, weight='bold',
                  color=color_drug_name, transform=ax_table.transAxes)

    h_y = 0.918
    ax_table.text(0.04, h_y, "PARAMETER",
                  weight='bold', fontsize=9, color='#1A237E',
                  transform=ax_table.transAxes)
    ax_table.text(0.52, h_y, route1_name.upper(),
                  weight='bold', fontsize=9, ha='center', color=COLOR_R1_LINE,
                  transform=ax_table.transAxes)
    ax_table.text(0.80, h_y, r2_label,
                  weight='bold', fontsize=9, ha='center', color=COLOR_R2_LINE,
                  transform=ax_table.transAxes)
    ax_table.plot([0.03, 0.97], [h_y-0.020, h_y-0.020],
                  color='#3E2723', lw=1.4, transform=ax_table.transAxes)

    rows = []
    rows.append(("Dose (mg)",
                 f"{dose1:.1f} mg", f"{dose2:.1f} mg" if has_second_route else "—", False))

    rows.append(("PHARMACODYNAMIC CONSTANTS", "", "", True))
    rows.append(("Onset of Action",
                 f"{onset1:.3f} h", f"{onset2:.3f} h" if has_second_route else "—", False))
    rows.append(("Duration of Action",
                 f"{duration1:.3f} h", f"{duration2:.3f} h" if has_second_route else "—", False))

    rows.append(("PHARMACOKINETIC CONSTANTS", "", "", True))
    rows.append(("AUC (mg·h/L)",
                 f"{pk1['auc']:.4f}", f"{pk2['auc']:.4f}" if has_second_route else "—", False))
    rows.append(("AUMC (mg·h²/L)",
                 f"{pk1['aumc']:.4f}", f"{pk2['aumc']:.4f}" if has_second_route else "—", False))
    rows.append(("MRT — Mean Residence Time (h)",
                 f"{pk1['mrt']:.4f}", f"{pk2['mrt']:.4f}" if has_second_route else "—", False))
    rows.append(("Cmax — Peak Concentration (mg/L)",
                 f"{pk1['c_max']:.4f}", f"{pk2['c_max']:.4f}" if has_second_route else "—", False))
    rows.append(("Tmax — Time to Peak (h)",
                 f"{pk1['t_max']:.4f}", f"{pk2['t_max']:.4f}" if has_second_route else "—", False))
    rows.append(("Ke — Elimination Rate (h⁻¹)",
                 f"{pk1['ke']:.5f}", f"{pk2['ke']:.5f}" if has_second_route else "—", False))
    rows.append(("t½ — Half-Life (h)",
                 f"{pk1['t_half']:.4f}", f"{pk2['t_half']:.4f}" if has_second_route else "—", False))
    rows.append(("Cl — Clearance (L/h)",
                 f"{pk1['cl']:.5f}", f"{pk2['cl']:.5f}" if has_second_route else "—", False))
    rows.append(("Vd — Volume of Distribution (L)",
                 f"{pk1['vd']:.5f}", f"{pk2['vd']:.5f}" if has_second_route else "—", False))

    rows.append(("DOSING STRATEGY", "", "", True))
    rows.append(("Loading Dose (mg)",
                 f"{loading_dose1:.3f}", f"{loading_dose2:.3f}" if has_second_route else "—", False))
    rows.append(("Maintenance Rate (mg/h)",
                 f"{maintenance_rate1:.3f}", f"{maintenance_rate2:.3f}" if has_second_route else "—", False))

    rows.append(("STATUS", status_1, status_2, False))

    row_h  = min(0.042, (h_y - 0.092) / len(rows))
    curr_y = h_y - 0.052
    alt_bg = False

    for lbl, v1, v2, is_header in rows:
        if is_header:
            curr_y -= 0.004
            ax_table.fill_betweenx([curr_y-0.008, curr_y+0.028], 0.03, 0.97,
                                    color='#E3F2FD', alpha=0.85,
                                    transform=ax_table.transAxes)
            ax_table.text(0.042, curr_y+0.007, f"▸  {lbl}",
                          fontsize=8.5, weight='bold', color='#0D47A1',
                          transform=ax_table.transAxes)
            curr_y -= row_h
            alt_bg = False
            continue

        if alt_bg:
            ax_table.fill_betweenx([curr_y-0.010, curr_y+0.022], 0.03, 0.97,
                                    color='#FFF8E1', alpha=0.50,
                                    transform=ax_table.transAxes)
        alt_bg = not alt_bg

        ax_table.text(0.042, curr_y, lbl,
                      fontsize=8.0, fontfamily='monospace', color='#212121',
                      transform=ax_table.transAxes, clip_on=True)

        cv1 = t_color_r1
        cv2 = t_color_r2
        fw  = 'normal'
        if lbl == "STATUS":
            fw  = 'bold'
            cv1 = '#B71C1C' if "Toxic" in v1 else '#1B5E20'
            cv2 = '#B71C1C' if "Toxic" in v2 else ('#1B5E20' if has_second_route else '#888')

        ax_table.text(0.52, curr_y, v1,
                      fontsize=8.0, fontfamily='monospace',
                      ha='center', color=cv1, weight=fw,
                      transform=ax_table.transAxes, clip_on=True)
        ax_table.text(0.80, curr_y, v2,
                      fontsize=8.0, fontfamily='monospace',
                      ha='center', color=cv2, weight=fw,
                      transform=ax_table.transAxes, clip_on=True)

        ax_table.plot([0.03, 0.97], [curr_y-0.010, curr_y-0.010],
                      color='#CFD8DC', lw=0.5, alpha=0.6,
                      transform=ax_table.transAxes)
        curr_y -= row_h

    if has_second_route and bioavailability is not None:
        bio_y    = max(curr_y - 0.020, 0.018)
        bio_text = f"BIOAVAILABILITY  [ F ]  =  {bioavailability:.2f} %"
        ax_table.text(0.5, bio_y, bio_text,
                      ha='center', va='center',
                      fontsize=9.5, weight='bold', color='#4A148C',
                      bbox=dict(facecolor='#F3E5F5', edgecolor='#7B1FA2',
                                boxstyle='round,pad=0.6', linewidth=2.0),
                      transform=ax_table.transAxes, clip_on=True)

    # ── BENCHMARK BARS (unchanged) ────────────────────────────────────────────
    if show_mec_msc:
        b_labels = ['MEC', 'Cmax', 'MSC']
        b_vals_1 = [mec, pk1['c_max'], msc]
    else:
        b_labels = ['Cmax']
        b_vals_1 = [pk1['c_max']]

    draw_gradient_bars(ax2, b_labels, b_vals_1,
                       cmap_name='Blues',
                       title_text=f"Safety Benchmark  ({route1_name})",
                       edge='#0277BD')

    if has_second_route and pk2:
        b_vals_2 = [mec, pk2['c_max'], msc] if show_mec_msc else [pk2['c_max']]
        draw_gradient_bars(ax3,
                           b_labels if show_mec_msc else ['Cmax'],
                           b_vals_2,
                           cmap_name='Greens',
                           title_text=f"Safety Benchmark  ({route2_name})",
                           edge='#2E7D32')
    else:
        ax3.set_visible(False)

    # ── RETURN AS BYTES (replaces plt.savefig to file) ────────────────────────
    buf = io.BytesIO()
    plt.savefig(buf, format='jpeg', dpi=300, bbox_inches='tight', facecolor=CRAYON_BG)
    plt.close(fig)
    buf.seek(0)
    return buf.read()
