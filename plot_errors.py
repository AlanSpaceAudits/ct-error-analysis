"""
Generate per-aircraft error analysis plots matching the CelTheo scatter style.

One plot per aircraft showing:
  - Green dot: FE predicted elevation ± measurement error
  - Magenta dot: GE predicted elevation ± measurement error
  - Orange line: inscribed angle (curvature drop) as reference

Legend shows aircraft details, predictions, and error values.
Style matches 00_scatter_graph.py from cel_theo_graphs.

Figures saved to error_graphs/.

Usage:
    pip install matplotlib geographiclib
    python plot_errors.py
"""

import math
import os
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from geographiclib.geodesic import Geodesic


# ============================================================
# Output directory
# ============================================================
OUT_DIR = 'error_graphs'
os.makedirs(OUT_DIR, exist_ok=True)


# ============================================================
# Plot style variables — matching 00_scatter_graph.py
# ============================================================
FIG_WIDTH = 9.0
FIG_HEIGHT = 6.0
RIGHT_MARGIN = 0.68

MARKER_SIZE = 8
ERR_LINEWIDTH = 1.5
ERR_CAPSIZE = 4
TARGET_LINEWIDTH = 2.0
TARGET_BAND_LINEWIDTH = 1.0

LEGEND_FONT_SIZE = 8.5
LEGEND_LABELSPACING = 0.4
LEGEND_HANDLELEN = 1.5
LEGEND_HANDLETEXTPAD = 0.6

Y_PAD = 0.15
X_RANGE = 0.15

TARGET_COLOR = "orange"
FE_COLOR = "green"
GE_COLOR = "magenta"


# ============================================================
# Constants — same as error_analysis.py
# ============================================================
BALLOON_LAT, BALLOON_LON = 25.756, -80.384
GEOID_UNDULATION = -24.022
EARTH_RADIUS = 6371000.0

ALTIMETER_TABLE = [
    (0, 20), (1000, 20), (2000, 30), (4000, 35), (6000, 40),
    (8000, 60), (10000, 80), (12000, 90), (14000, 100),
    (16000, 110), (18000, 120), (20000, 130), (25000, 155),
    (30000, 180), (35000, 205), (40000, 230), (45000, 255),
]


def altimeter_tolerance_ft(alt_ft):
    for i in range(len(ALTIMETER_TABLE) - 1):
        lo_alt, lo_tol = ALTIMETER_TABLE[i]
        hi_alt, hi_tol = ALTIMETER_TABLE[i + 1]
        if lo_alt <= alt_ft <= hi_alt:
            frac = (alt_ft - lo_alt) / (hi_alt - lo_alt)
            return lo_tol + frac * (hi_tol - lo_tol)
    return ALTIMETER_TABLE[-1][1]


def fmt_arcmin(deg):
    """Format degrees to M'SS.S\" string."""
    total_arcmin = abs(deg) * 60
    arcmin = int(total_arcmin)
    arcsec = (total_arcmin - arcmin) * 60
    return f"{arcmin}'{arcsec:04.1f}\""


# ============================================================
# Aircraft data
# ============================================================
AIRCRAFT = [
    {
        'name': 'Spirit Airlines',
        'flight_level': 'FL340',
        'star': 'HIP 44471',
        'date': '1/6/2026',
        'time_local': '8:11:30 PM',
        'utc_offset': '-5',
        'obs': (27.06811, -82.22231, -1.07),
        'tgt': (27.2931, -81.9698, 10375.54),
        'alt_ft': 34040,
        'gs_ms': None,
        'vr_ms': 4.88,
        'delta_s': 2,
        'ge_pred': 16.2199,
        'fe_pred': 16.3714,
    },
    {
        'name': 'Allegiant 453',
        'flight_level': 'FL310',
        'star': 'HIP 47693',
        'date': '1/16/2026',
        'time_local': '8:59:08 PM',
        'utc_offset': '-5',
        'obs': (27.06894, -82.22377, 4.42),
        'tgt': (27.0581, -81.5134, 9466.84),
        'alt_ft': 31050,
        'gs_ms': 215.0,
        'vr_ms': 5.20,
        'delta_s': 4,
        'ge_pred': 7.3687,
        'fe_pred': 7.6468,
    },
    {
        'name': 'Jet Speed 691',
        'flight_level': 'FL390',
        'star': 'HIP 23078',
        'date': '1/29/2026',
        'time_local': '9:44:35 PM',
        'utc_offset': '-5',
        'obs': (27.06765, -82.22612, 6.2),
        'tgt': (26.5664, -82.3093, 11898.24),
        'alt_ft': 39030,
        'gs_ms': 287.6,
        'vr_ms': -14.31,
        'delta_s': 3,
        'ge_pred': 11.7265,
        'fe_pred': 11.9578,
    },
    {
        'name': 'AAL1517',
        'flight_level': 'FL308',
        'star': 'HIP 39850',
        'date': '1/29/2026',
        'time_local': '9:56:07 PM',
        'utc_offset': '-5',
        'obs': (27.06765, -82.22612, 6.2),
        'tgt': (26.3003, -81.9388, 9392.17),
        'alt_ft': 30810,
        'gs_ms': 275.2,
        'vr_ms': -4.88,
        'delta_s': 5,
        'ge_pred': 5.6207,
        'fe_pred': 5.9734,
    },
]


def compute_errors(ac):
    """Compute the full error budget for one aircraft."""
    obs_lat, obs_lon, obs_alt = ac['obs']
    tgt_lat, tgt_lon, tgt_alt = ac['tgt']

    inv = Geodesic.WGS84.Inverse(obs_lat, obs_lon, tgt_lat, tgt_lon)
    nom_dist = inv['s12']
    nom_alt_diff = tgt_alt - obs_alt

    balloon_to_ac = Geodesic.WGS84.Inverse(
        BALLOON_LAT, BALLOON_LON, tgt_lat, tgt_lon)['s12'] / 1000

    nacp9 = 30.0
    cpr = 5.0
    latency_err = ac['gs_ms'] * 0.6 if ac['gs_ms'] else 0
    temporal_horiz = ac['gs_ms'] * ac['delta_s'] if ac['gs_ms'] else 0
    rss_horiz = math.sqrt(nacp9**2 + cpr**2 + latency_err**2)
    total_horiz = rss_horiz + temporal_horiz

    alt_tol_m = altimeter_tolerance_ft(ac['alt_ft']) * 0.3048
    adsb_alt_m = 12.5 * 0.3048
    sse_m = 90 * 0.3048
    geom_alt = (EARTH_RADIUS * tgt_alt) / (EARTH_RADIUS - tgt_alt)
    geop_err = geom_alt - tgt_alt
    temp_err_K = 0.42 * (balloon_to_ac / 100)
    alt_repr_m = temp_err_K * 4.1
    temporal_vert = abs(ac['vr_ms'] * ac['delta_s'])

    rss_alt = math.sqrt(alt_tol_m**2 + adsb_alt_m**2 + sse_m**2
                        + alt_repr_m**2)
    systematic_alt = geop_err + abs(GEOID_UNDULATION)
    total_alt_min = -rss_alt - temporal_vert
    total_alt_max = rss_alt + systematic_alt + temporal_vert

    close_dist = max(nom_dist - total_horiz, 1000)
    elev_close = math.degrees(math.atan2(
        nom_alt_diff + total_alt_max, close_dist))
    far_dist = nom_dist + total_horiz
    elev_far = math.degrees(math.atan2(
        nom_alt_diff + total_alt_min, far_dist))
    uncertainty = (elev_close - elev_far) / 2

    central_angle_deg = math.degrees(nom_dist / EARTH_RADIUS)
    inscribed_angle_deg = central_angle_deg / 2
    drop_m = nom_dist**2 / (2 * EARTH_RADIUS)
    error_drop_m = drop_m * (uncertainty / inscribed_angle_deg) if inscribed_angle_deg > 0 else 0

    return {
        'dist_km': nom_dist / 1000,
        'uncertainty': uncertainty,
        'inscribed_angle': inscribed_angle_deg,
        'drop_m': drop_m,
        'error_drop_m': error_drop_m,
        'ratio': uncertainty / inscribed_angle_deg if inscribed_angle_deg > 0 else 0,
    }


# ============================================================
# Generate one plot per aircraft
# ============================================================

for ac in AIRCRAFT:
    r = compute_errors(ac)

    fe_val = ac['fe_pred']
    ge_val = ac['ge_pred']
    err = r['uncertainty']
    inscribed = r['inscribed_angle']
    drop_m = r['drop_m']
    err_drop_m = r['error_drop_m']
    ratio = r['ratio']

    # X positions — FE left, GE right (no x-axis, just spacing)
    x_fe = -0.05
    x_ge = 0.05

    fig, ax = plt.subplots(figsize=(FIG_WIDTH, FIG_HEIGHT))
    fig.subplots_adjust(right=RIGHT_MARGIN)

    # ── Inscribed angle as reference band ──
    # Show the GE-FE gap: horizontal lines at FE and GE predictions
    # with the inscribed angle (drop) as the space between them
    mid = (fe_val + ge_val) / 2
    inscribed_line = ax.axhline(
        mid, color=TARGET_COLOR, linewidth=TARGET_LINEWIDTH
    )
    ax.axhline(
        mid + inscribed / 2,
        color=TARGET_COLOR, linestyle="--",
        linewidth=TARGET_BAND_LINEWIDTH, alpha=0.7,
    )
    ax.axhline(
        mid - inscribed / 2,
        color=TARGET_COLOR, linestyle="--",
        linewidth=TARGET_BAND_LINEWIDTH, alpha=0.7,
    )

    # ── FE dot + error bar ──
    fe_errbar = ax.errorbar(
        x_fe, fe_val, yerr=err,
        fmt="o", markersize=MARKER_SIZE,
        color=FE_COLOR, ecolor=FE_COLOR,
        elinewidth=ERR_LINEWIDTH, capsize=ERR_CAPSIZE,
    )

    # ── GE dot + error bar ──
    ge_errbar = ax.errorbar(
        x_ge, ge_val, yerr=err,
        fmt="o", markersize=MARKER_SIZE,
        color=GE_COLOR, ecolor=GE_COLOR,
        elinewidth=ERR_LINEWIDTH, capsize=ERR_CAPSIZE,
    )

    # ── Y limits with padding ──
    y_values = [
        fe_val - err, fe_val + err,
        ge_val - err, ge_val + err,
    ]
    y_min = min(y_values)
    y_max = max(y_values)
    span = (y_max - y_min) if y_max != y_min else 0.1
    pad = Y_PAD * span
    ax.set_ylim(y_min - pad, y_max + pad)

    ax.set_xlim(-X_RANGE, X_RANGE)
    ax.set_xticks([])
    ax.set_xlabel("")
    ax.set_ylabel("Elevation Angle (°)")
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    ax.set_title(ac['name'])

    # ── Legend ──
    dummy = Line2D([0], [0], color="none")

    ge_fe_diff = fe_val - ge_val

    legend_handles = [
        dummy,               # Aircraft name
        dummy,               # Star
        dummy,               # Time
        dummy,               # Date
        dummy,               # UTC offset
        dummy,               # Blank spacer
        inscribed_line,      # Inscribed angle line
        dummy,               # Inscribed angle error band
        fe_errbar[0],        # FE prediction
        ge_errbar[0],        # GE prediction
        dummy,               # Blank spacer
        dummy,               # Ratio
    ]

    legend_labels = [
        f"{ac['name']} ({ac['flight_level']})",
        f"Occl. Star: {ac['star']}",
        f"Occl. Time: {ac['time_local']}",
        f"Occl. Date: {ac['date']}",
        f"UTC Offset: {ac['utc_offset']}",
        "",
        (f"GE–FE diff (inscribed ∠): {ge_fe_diff:.3f}° "
         f"({fmt_arcmin(ge_fe_diff)}), {drop_m:.0f} m"),
        f"Measurement error: ±{err:.3f}° ({fmt_arcmin(err)}), ±{err_drop_m:.0f} m",
        f"θFE elev: {fe_val:.4f}° ±{err:.3f}°",
        f"θGE elev: {ge_val:.4f}° ±{err:.3f}°",
        "",
        f"Error ÷ Drop: {ratio:.1f}×  (~{r['dist_km']:.0f} km)",
    ]

    ax.legend(
        legend_handles,
        legend_labels,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        borderaxespad=0.0,
        frameon=True,
        fontsize=LEGEND_FONT_SIZE,
        labelspacing=LEGEND_LABELSPACING,
        handlelength=LEGEND_HANDLELEN,
        handletextpad=LEGEND_HANDLETEXTPAD,
    )

    plt.tight_layout()

    safe_name = ac['name'].replace(" ", "_").replace("'", "")
    save_path = os.path.join(OUT_DIR, f'{safe_name}_error.png')
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f'Saved: {save_path}')
    plt.close(fig)

print('\nDone. All per-aircraft plots saved to error_graphs/')
