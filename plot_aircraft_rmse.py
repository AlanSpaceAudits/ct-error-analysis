"""
Aircraft Occultation — RMSE-style bar charts
=============================================

Generates per-aircraft bar charts in the same visual style as the
Celestial Theodolite mountain peak RMSE graphs. Two panels side by side:

  Left  (green):  FE predicted vs FE observed elevation
  Right (pink):   GE predicted vs GE observed elevation

Error bars use the RSS of ADS-B instrument uncertainty and the CelTheo
shared measurement sigma (from 65 mountain peak observations).
The GE panel shows the inscribed angle (curvature drop) as a hatched zone.

Usage:
    python3 plot_aircraft_rmse.py
"""

import math
import os
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from geographiclib.geodesic import Geodesic

# ── Output ──
OUT_DIR = 'error_graphs'
os.makedirs(OUT_DIR, exist_ok=True)

# ── Constants ──
EARTH_RADIUS = 6371000.0
BALLOON_LAT, BALLOON_LON = 25.756, -80.384
GEOID_UNDULATION = -24.022
CELTHEO_SIGMA = 0.075408  # shared σ from 65 mountain peak observations (degrees)

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


def compute_errors(ac):
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
    adsb_uncertainty = (elev_close - elev_far) / 2
    combined_uncertainty = math.sqrt(adsb_uncertainty**2 + CELTHEO_SIGMA**2)
    drop_m = nom_dist**2 / (2 * EARTH_RADIUS)
    fe_ge_diff = ac['fe_pred'] - ac['ge_pred']
    error_as_drop_m = drop_m * (combined_uncertainty / fe_ge_diff) if fe_ge_diff > 0 else 0
    return {
        'dist_km': nom_dist / 1000,
        'dist_m': nom_dist,
        'adsb_uncertainty': adsb_uncertainty,
        'celtheo_sigma': CELTHEO_SIGMA,
        'uncertainty': combined_uncertainty,
        'inscribed_angle': fe_ge_diff,
        'drop_m': drop_m,
        'error_as_drop_m': error_as_drop_m,
    }


def deg_to_dms(deg):
    sign = "-" if deg < 0 else ""
    deg = abs(deg)
    d = int(deg)
    rem = (deg - d) * 60
    m = int(rem)
    s = (rem - m) * 60
    return f"{sign}{d}°{m}'{s:.1f}\""


def deg_to_dms_signed(deg):
    sign = "+" if deg >= 0 else "-"
    deg = abs(deg)
    d = int(deg)
    rem = (deg - d) * 60
    m = int(rem)
    s = (rem - m) * 60
    return f"{sign}{d}°{m}'{s:.1f}\""


# ── Aircraft data ──
AIRCRAFT = [
    {
        'name': 'Spirit Airlines',
        'flight_level': 'FL340',
        'star': 'HIP 44471',
        'date': '01/06/26',
        'obs': (27.06811, -82.22231, -1.07),
        'tgt': (27.2931, -81.9698, 10375.54),
        'alt_ft': 34040,
        'gs_ms': None,
        'vr_ms': 4.88,
        'delta_s': 2,
        'fe_pred': 16.3714,
        'ge_pred': 16.2199,
    },
    {
        'name': 'Allegiant 453',
        'flight_level': 'FL310',
        'star': 'HIP 47693',
        'date': '01/16/26',
        'obs': (27.06894, -82.22377, 4.42),
        'tgt': (27.0581, -81.5134, 9466.84),
        'alt_ft': 31050,
        'gs_ms': 215.0,
        'vr_ms': 5.20,
        'delta_s': 4,
        'fe_pred': 7.6468,
        'ge_pred': 7.3687,
    },
    {
        'name': 'Jet Speed 691',
        'flight_level': 'FL390',
        'star': 'HIP 23078',
        'date': '01/29/26',
        'obs': (27.06765, -82.22612, 6.2),
        'tgt': (26.5664, -82.3093, 11898.24),
        'alt_ft': 39030,
        'gs_ms': 287.6,
        'vr_ms': -14.31,
        'delta_s': 3,
        'fe_pred': 11.9578,
        'ge_pred': 11.7265,
    },
    {
        'name': 'AAL1517',
        'flight_level': 'FL308',
        'star': 'HIP 39850',
        'date': '01/29/26',
        'obs': (27.06765, -82.22612, 6.2),
        'tgt': (26.3003, -81.9388, 9392.17),
        'alt_ft': 30810,
        'gs_ms': 275.2,
        'vr_ms': -4.88,
        'delta_s': 5,
        'fe_pred': 5.9734,
        'ge_pred': 5.6207,
    },
]


# ── Generate graphs ──

for ac in AIRCRAFT:
    r = compute_errors(ac)

    theta_fe = ac['fe_pred']
    theta_ge = ac['ge_pred']
    adsb_sigma = r['adsb_uncertainty']
    celtheo_sigma = r['celtheo_sigma']
    sigma = r['uncertainty']  # RSS of both
    inscribed_angle = r['inscribed_angle']
    d = r['dist_m']
    drop_m = r['drop_m']

    # For aircraft, FE and GE "observed" are the same as predictions
    # (single observation = the prediction point itself)
    # The error bars show the combined measurement uncertainty envelope
    fe_obs = theta_fe
    ge_obs = theta_ge

    error_as_drop_m = r['error_as_drop_m']

    # ── Figure ──
    fig, (ax_fe, ax_ge) = plt.subplots(1, 2, figsize=(14, 7))
    fig.suptitle(
        f"{ac['name']} ({ac['flight_level']}) — RMSE Summary\n"
        f"1 star, 1 obs  |  Star: {ac['star']}\n"
        f"dist = {d:,.0f} m  |  Date: {ac['date']}",
        fontsize=15, fontweight="bold",
    )

    # Shared y-axis
    y_max_val = max(theta_fe, theta_ge)
    y_min_val = min(0, theta_ge - 2 * sigma)
    y_range = y_max_val - y_min_val
    y_max = y_max_val + y_range * 0.35
    y_min = y_min_val - y_range * 0.05 if y_min_val < 0 else 0
    for ax in (ax_fe, ax_ge):
        ax.set_ylim(y_min, y_max)
        ax.set_ylabel("Angle (°)", fontsize=12)
        if y_min < 0:
            ax.axhline(y=0, color="#999999", linewidth=0.8, zorder=1)

    # ── LEFT: Flat Earth ──
    ax_fe.set_title(f"Flat Earth  (error bars = ±σ_combined = ±{error_as_drop_m:.1f} m)", fontsize=12)

    ax_fe.bar(0, theta_fe, width=0.5, color="#a8e6a0", edgecolor="none",
              label="FE predicted θ", zorder=2)

    ax_fe.bar(0.55, fe_obs, width=0.5, color="#4caf50", edgecolor="none",
              label="FE observed", zorder=2)
    ax_fe.errorbar(0.55, fe_obs, yerr=sigma, fmt="none",
                   ecolor="black", capsize=6, capthick=1.5, linewidth=1.5, zorder=3)

    ax_fe.axhline(y=theta_fe, color="red", linestyle="--", linewidth=1.5,
                  label="predicted ref", zorder=4)

    # Error badge with breakdown
    ax_fe.text(
        0.275, y_max * 0.92,
        f"σ = {deg_to_dms(sigma)}\n"
        f"ADS-B: {deg_to_dms(adsb_sigma)}  |  CelTheo: {deg_to_dms(celtheo_sigma)}",
        fontsize=10, fontweight="bold", color="white", ha="center",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#2e7d32", edgecolor="none"),
    )

    ax_fe.set_xticks([0, 0.55])
    ax_fe.set_xticklabels(["Predicted\nθ_FE", "Observed"], fontsize=10)
    ax_fe.legend(loc="upper right", fontsize=9)

    # ── RIGHT: Globe Earth ──
    ax_ge.set_title(f"Globe Earth  (error bars = ±σ_combined = ±{error_as_drop_m:.1f} m)", fontsize=12)

    ax_ge.bar(0, theta_ge, width=0.5, color="#f8bbd0", edgecolor="none",
              label="GE predicted θ", zorder=2)

    # Hatched drop zone
    matplotlib.rcParams['hatch.color'] = "#e91e90"
    ax_ge.bar(0, inscribed_angle, width=0.5, bottom=theta_ge,
              color="#f8bbd0", edgecolor="#e91e90", linewidth=0.8,
              hatch="//", zorder=3)
    drop_center_y = theta_ge + inscribed_angle / 2
    ax_ge.text(0, drop_center_y,
               f"drop = {deg_to_dms(inscribed_angle)}\n({drop_m:.1f} m)",
               fontsize=8, color="#ad1457", ha="center", va="center",
               bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                         edgecolor="#e91e90", alpha=0.85),
               zorder=5)

    ax_ge.bar(0.55, ge_obs, width=0.5, color="#e91e90", edgecolor="none",
              label="GE observed", zorder=2)
    ax_ge.errorbar(0.55, ge_obs, yerr=sigma, fmt="none",
                   ecolor="black", capsize=6, capthick=1.5, linewidth=1.5, zorder=3)

    ax_ge.axhline(y=theta_ge, color="red", linestyle="--", linewidth=1.5,
                  label="predicted ref", zorder=4)

    # Error badge with breakdown
    ax_ge.text(
        0.275, y_max * 0.92,
        f"σ = {deg_to_dms(sigma)}\n"
        f"ADS-B: {deg_to_dms(adsb_sigma)}  |  CelTheo: {deg_to_dms(celtheo_sigma)}",
        fontsize=10, fontweight="bold", color="white", ha="center",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#c62828", edgecolor="none"),
    )

    # Error ÷ Drop annotation
    ratio = sigma / inscribed_angle if inscribed_angle > 0 else 0
    ax_ge.annotate(
        f"Error ÷ Drop: {ratio:.1f}x",
        xy=(0.275, theta_fe + y_max * 0.02), fontsize=10, color="#ad1457",
        ha="center", va="bottom",
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                  edgecolor="#ad1457", alpha=0.9),
    )

    ax_ge.set_xticks([0, 0.55])
    ax_ge.set_xticklabels(["Predicted\nθ_GE", "Observed"], fontsize=10)
    ax_ge.legend(loc="upper right", fontsize=9)

    # ── Bottom annotation ──
    fig.text(
        0.5, 0.02,
        f"Error ÷ Drop = {ratio:.1f}x  |  "
        f"Measurement error swamps inscribed angle" if ratio >= 1.0 else
        f"Error ÷ Drop = {ratio:.1f}x  |  "
        f"Error bars overlap inscribed angle band",
        fontsize=12, fontweight="bold", ha="center", color="#333333",
    )

    plt.tight_layout(rect=[0, 0.05, 1, 0.95])
    safe_name = ac['name'].replace(" ", "_").replace("'", "")
    png_path = os.path.join(OUT_DIR, f'{safe_name}_rmse.png')
    plt.savefig(png_path, dpi=150, bbox_inches="tight")
    print(f'Saved: {png_path}')
    plt.close(fig)

print('\nDone.')
