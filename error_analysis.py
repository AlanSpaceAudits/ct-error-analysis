"""
Full-spectrum error analysis for AstroLive aircraft occultation data.

AstroLive recorded four aircraft stellar occultations using a telescope with
timestamped video and an RTL-SDR ADS-B receiver. The aircraft positions from
those ADS-B broadcasts were fed into a raytracer to predict star positions and
compare globe vs flat earth models. This script quantifies every error source
in the aircraft position chain to determine whether the input data is precise
enough to support the claimed results.

The aircraft positions come from ADS-B (1090 MHz) decoded by SDRAngel. Each
position passes through: GPS fix -> transponder -> ADS-B broadcast -> RTL-SDR
reception -> SDRAngel decode -> CSV -> raytracer. Every link adds uncertainty.

The altitude chain is even worse: barometric pressure on the aircraft ->
ISA pressure-altitude -> interpolation from a weather balloon 234 km away ->
geopotential height -> [geometric conversion SKIPPED in the raytracer].

Primary sources for every number:
  - RTCA DO-260B: NACp accuracy categories, CPR encoding, latency limits
  - 14 CFR Part 43, Appendix E, Table I: Altimeter scale error tolerances
  - 14 CFR Part 91, Appendix G: RVSM altimetry system error (ASE) limits
  - 14 CFR 25.1325(e): Static source error (SSE) allowance
  - Sun et al. (2010), J. Geophys. Res. 115, D23104: Radiosonde spatial
    representativeness — 0.42 K per 100 km at individual pressure levels
  - NOAA GEOID18: Geoid undulation at observer location

Usage:
    pip install geographiclib
    python error_analysis.py
"""

import math
from geographiclib.geodesic import Geodesic


# ============================================================
# FAA Part 43 Appendix E Table I — altimeter scale error
# ============================================================
# Each tuple is (altitude_ft, allowed_error_ft).
# These are the maximum permissible scale errors when an
# altimeter is bench-tested per 14 CFR Part 43 Appendix E.
# We interpolate between table entries for the actual flight level.
# Source: https://www.law.cornell.edu/cfr/text/14/appendix-E_to_part_43
ALTIMETER_TABLE = [
    (0, 20), (1000, 20), (2000, 30), (4000, 35), (6000, 40),
    (8000, 60), (10000, 80), (12000, 90), (14000, 100),
    (16000, 110), (18000, 120), (20000, 130), (25000, 155),
    (30000, 180), (35000, 205), (40000, 230), (45000, 255),
]


def altimeter_tolerance_ft(alt_ft):
    """
    Interpolate 14 CFR Part 43 Appendix E Table I to get the
    allowed altimeter scale error (±ft) at a given altitude.
    Linear interpolation between the table entries above.
    """
    for i in range(len(ALTIMETER_TABLE) - 1):
        lo_alt, lo_tol = ALTIMETER_TABLE[i]
        hi_alt, hi_tol = ALTIMETER_TABLE[i + 1]
        if lo_alt <= alt_ft <= hi_alt:
            frac = (alt_ft - lo_alt) / (hi_alt - lo_alt)
            return lo_tol + frac * (hi_tol - lo_tol)
    # Above the table — return the highest entry
    return ALTIMETER_TABLE[-1][1]


# ============================================================
# Event data — from ADS-B screenshots and raytracer CSV
# ============================================================

# Weather balloon launch site (WMO station 72202, Miami FL).
# This is where the radiosonde data comes from. All four events
# are 167-234 km away from this station.
BALLOON_LAT, BALLOON_LON = 25.756, -80.384

# Geoid undulation at the observer location (27.068°N, 82.222°W).
# Measured using NOAA GEOID18 interactive calculator.
# N = -24.022 m means MSL is 24 m below the WGS84 ellipsoid here.
# The raytracer uses WGS84 geometry but the observer altitude is MSL,
# and it never accounts for this 24 m offset.
# Source: https://geodesy.noaa.gov/GEOID/GEOID18/computation.html
GEOID_UNDULATION = -24.022  # meters

# Mean Earth radius for geopotential-to-geometric height conversion.
# h_geometric = (R * h_geopotential) / (R - h_geopotential)
EARTH_RADIUS = 6371000.0  # meters

# Each aircraft dict contains:
#   name:         Flight identifier
#   flight_level: Reported flight level from ADS-B
#   obs:          Observer position (lat, lon, altitude_m) — from GPS at observation site
#   tgt:          Aircraft position (lat, lon, altitude_m) — lat/lon from ADS-B screenshot,
#                 altitude from raytracer CSV (geopotential height via balloon interpolation)
#   alt_ft:       Aircraft altitude in feet from ADS-B broadcast (barometric, Mode S encoded)
#   gs_ms:        Ground speed in m/s from ADS-B (None if not available)
#   vr_ms:        Vertical rate in m/s from ADS-B (positive = climbing)
#   trk:          Track heading in degrees from ADS-B
#   delta_s:      Time gap in seconds between occultation frame and ADS-B "Updated" timestamp.
#                 The ADS-B position is AFTER the occultation — the plane has already moved.
AIRCRAFT = [
    {
        'name': 'Spirit Airlines',
        'flight_level': 'FL340',
        'obs': (27.06811, -82.22231, -1.07),  # Observer slightly below MSL
        'tgt': (27.2931, -81.9698, 10375.54),  # 34,040 ft geopotential
        'alt_ft': 34040,
        'gs_ms': None,   # Ground speed not available in ADS-B screenshot
        'vr_ms': 4.88,   # Climbing at +960 ft/min = +4.88 m/s
        'trk': 290,      # Heading west-northwest
        'delta_s': 2,    # ADS-B timestamp is 2s after occultation
    },
    {
        'name': 'Allegiant 453',
        'flight_level': 'FL310',
        'obs': (27.06894, -82.22377, 4.42),
        'tgt': (27.0581, -81.5134, 9466.84),  # 31,050 ft geopotential
        'alt_ft': 31050,
        'gs_ms': 215.0,  # 418 knots = 215.0 m/s
        'vr_ms': 5.20,   # Climbing at +1024 ft/min = +5.20 m/s
        'trk': 322,      # Heading northwest
        'delta_s': 4,    # ADS-B timestamp is 4s after occultation
    },
    {
        'name': 'Jet Speed 691',
        'flight_level': 'FL390',
        'obs': (27.06765, -82.22612, 6.2),
        'tgt': (26.5664, -82.3093, 11898.24),  # 39,030 ft geopotential
        'alt_ft': 39030,
        'gs_ms': 287.6,  # 559 knots = 287.6 m/s
        'vr_ms': -14.31, # Descending at -2816 ft/min = -14.31 m/s
        'trk': 125,      # Heading southeast
        'delta_s': 3,    # ADS-B timestamp is 3s after occultation
    },
    {
        'name': 'American Airlines 1517',
        'flight_level': 'FL308',
        'obs': (27.06765, -82.22612, 6.2),
        'tgt': (26.3003, -81.9388, 9392.17),  # 30,810 ft geopotential
        'alt_ft': 30810,
        'gs_ms': 275.2,  # 535 knots = 275.2 m/s
        'vr_ms': -4.88,  # Descending at -960 ft/min = -4.88 m/s
        'trk': 114,      # Heading east-southeast
        'delta_s': 5,    # ADS-B timestamp is 5s after occultation
    },
]


def analyze(ac):
    """
    Run the full error budget for a single aircraft event.

    Computes the nominal geometry (distance, elevation angle), then
    quantifies every error source in horizontal position and altitude,
    and propagates them into an elevation angle uncertainty range.
    """
    obs_lat, obs_lon, obs_alt = ac['obs']
    tgt_lat, tgt_lon, tgt_alt = ac['tgt']

    # ── Nominal geometry ──
    # WGS84 geodesic inverse: exact surface distance between observer and aircraft
    inv = Geodesic.WGS84.Inverse(obs_lat, obs_lon, tgt_lat, tgt_lon)
    nom_dist = inv['s12']  # surface distance in meters
    nom_alt_diff = tgt_alt - obs_alt  # height of aircraft above observer
    # Geometric elevation angle: arctan(height / distance), no refraction
    nom_elev = math.degrees(math.atan2(nom_alt_diff, nom_dist))

    # ── How far is the balloon from the aircraft and observer? ──
    # The weather balloon (Miami) provides the temperature/pressure profiles
    # used to convert barometric altitude to geopotential height. The further
    # away the balloon is, the less representative its data is.
    balloon_to_ac = Geodesic.WGS84.Inverse(
        BALLOON_LAT, BALLOON_LON, tgt_lat, tgt_lon)['s12'] / 1000  # km
    balloon_to_obs = Geodesic.WGS84.Inverse(
        BALLOON_LAT, BALLOON_LON, obs_lat, obs_lon)['s12'] / 1000  # km

    print(f'\n{"=" * 90}')
    print(f'  {ac["name"]} — {ac["flight_level"]}')
    print(f'{"=" * 90}')
    print(f'  Nominal: dist={nom_dist/1000:.2f} km, '
          f'alt_diff={nom_alt_diff:.1f} m, elev={nom_elev:.4f}°')
    print(f'  Balloon→aircraft: {balloon_to_ac:.1f} km | '
          f'Balloon→observer: {balloon_to_obs:.1f} km')

    # ══════════════════════════════════════════════════════════
    # HORIZONTAL POSITION ERRORS
    # ══════════════════════════════════════════════════════════
    print(f'\n  ── HORIZONTAL POSITION ERRORS ──')

    # 1. GPS accuracy per DO-260B NACp categories (Table 2-73)
    # NACp=9: EPU < 30m at 95% confidence — typical for commercial jets
    # NACp=8: EPU < 92.6m (0.05 NM) — FAA minimum for ADS-B Out (14 CFR 91.227)
    # We use NACp=9 as the baseline (benefit of the doubt)
    nacp9 = 30.0   # meters, 95% confidence radius
    nacp8 = 92.6   # meters, shown for comparison
    print(f'  GPS/NACp=9 (95%):          ±{nacp9:.0f} m')
    print(f'  GPS/NACp=8 (95%):          ±{nacp8:.1f} m')

    # 2. CPR (Compact Position Reporting) quantization
    # ADS-B encodes lat/lon with 17 bits per coordinate (Nb=17, airborne).
    # Each ~6° zone is divided into 2^17 = 131,072 bins.
    # Resolution: 6° / 131072 × 111,320 m/° ≈ 5.1 m per bin
    # This is derived from the bit width, not directly stated in DO-260B.
    cpr = 5.0  # meters
    print(f'  CPR quantization:          ±{cpr:.0f} m')

    # 3. Uncompensated latency in the ADS-B broadcast itself
    # DO-260B §2.2.3.2.7: total latency GPS fix → broadcast ≤ 2.0 seconds,
    # of which up to 0.6 seconds may be UNCOMPENSATED — the transponder
    # doesn't extrapolate the position forward to cover this delay.
    # So the broadcast position is where the plane was 0.6s ago.
    if ac['gs_ms']:
        latency_err = ac['gs_ms'] * 0.6  # worst case uncompensated
        print(f'  Uncompensated latency:     ±{latency_err:.0f} m '
              f'(0.6s × {ac["gs_ms"]:.0f} m/s)')

        # 4. Temporal mismatch: the big one
        # The ADS-B "Updated" timestamp in SDRAngel is 2-5 seconds AFTER the
        # occultation frame. The raytracer uses these delayed positions with
        # NO extrapolation — no correction for ground speed or heading.
        # The position used is literally where the plane was seconds later.
        temporal_horiz = ac['gs_ms'] * ac['delta_s']
        print(f'  TEMPORAL MISMATCH:         {temporal_horiz:.0f} m '
              f'({ac["delta_s"]}s × {ac["gs_ms"]:.0f} m/s)')
    else:
        latency_err = 0
        temporal_horiz = 0
        print(f'  Uncompensated latency:     N/A (no GS)')
        print(f'  TEMPORAL MISMATCH:         N/A (no GS, '
              f'delta={ac["delta_s"]}s)')

    # RSS the random errors (GPS + CPR + latency), then add temporal
    # mismatch on top as a systematic offset
    rss_horiz = math.sqrt(nacp9**2 + cpr**2 + latency_err**2)
    total_horiz = rss_horiz + temporal_horiz
    print(f'  RSS (random, NACp=9):      ±{rss_horiz:.0f} m')
    print(f'  TOTAL (random + temporal): {total_horiz:.0f} m')

    # ══════════════════════════════════════════════════════════
    # ALTITUDE ERRORS
    # ══════════════════════════════════════════════════════════
    print(f'\n  ── ALTITUDE ERRORS ──')

    # 1. Altimeter instrument scale error
    # 14 CFR Part 43, Appendix E, Table I — bench test tolerance
    # Interpolated for the actual flight level of this aircraft
    alt_tol = altimeter_tolerance_ft(ac['alt_ft'])
    alt_tol_m = alt_tol * 0.3048
    print(f'  Altimeter scale error:     ±{alt_tol:.0f} ft '
          f'(±{alt_tol_m:.1f} m) [14 CFR Part 43 App E]')

    # 2. ADS-B altitude quantization (Mode S encoding)
    # The 12-bit altitude field uses a Q-bit to set resolution:
    # Q=1 (below 50,175 ft): 25 ft steps, altitude = 25N - 1000
    # Q=0 (above 50,175 ft): 100 ft Gray code (legacy Gillham)
    # All four aircraft are well below 50,175 ft → Q=1 → 25 ft steps
    # Rounding to nearest 25 ft step gives max error of ±12.5 ft
    adsb_alt_ft = 12.5
    adsb_alt_m = adsb_alt_ft * 0.3048
    print(f'  ADS-B alt quantization:    ±{adsb_alt_ft} ft '
          f'(±{adsb_alt_m:.1f} m) [25 ft steps, Mode S]')

    # 3. Static Source Error (SSE)
    # 14 CFR 25.1325(e): airflow over the fuselage at cruise speed creates
    # pressure distortions at the static ports, causing the altimeter to
    # read wrong. Allowed error: ±30 ft per 100 knots of airspeed.
    # At typical cruise IAS of ~300 knots: ±90 ft
    # This doesn't affect a weather balloon — it floats vertically through
    # still air with no aerodynamic disturbance around its sensor.
    sse_ft = 90  # ±30 ft/100 kn × ~300 kn cruise IAS
    sse_m = sse_ft * 0.3048
    print(f'  Static source error:       ±{sse_ft} ft '
          f'(±{sse_m:.1f} m) [14 CFR 25.1325(e), ~300 kn]')

    # 4. RVSM total altimetry system error
    # 14 CFR Part 91, Appendix G: for RVSM airspace (FL290-FL410),
    # mean ASE + 3σ must be ≤ ±200 ft. This is the total system budget
    # including instrument, encoder, SSE, and avionics.
    # Shown for reference — the components above are already counted separately.
    rvsm_ft = 200
    rvsm_m = rvsm_ft * 0.3048
    print(f'  RVSM total ASE (3σ):       ±{rvsm_ft} ft '
          f'(±{rvsm_m:.1f} m) [14 CFR 91 App G]')

    # 5. Geopotential-to-geometric height conversion (SKIPPED)
    # The raytracer converts geopotential → geometric ONLY for weather
    # balloon targets (line 116 of CelestialTheodoliteDataProcessor.py).
    # Aircraft targets skip this entirely. The formula:
    # h_geometric = (R × h_geopotential) / (R - h_geopotential)
    # At 10-12 km altitude this adds 14-22 m. The aircraft is always
    # higher than assumed — systematic positive bias.
    geom_alt = (EARTH_RADIUS * tgt_alt) / (EARTH_RADIUS - tgt_alt)
    geop_err = geom_alt - tgt_alt
    print(f'  Geopotential→geometric:    +{geop_err:.1f} m '
          f'[SKIPPED in processing]')

    # 6. Geoid undulation (not accounted for anywhere)
    # NOAA GEOID18 at observer location: N = -24.022 m
    # The observer altitude is MSL (the geoid). The raytracer does WGS84
    # ellipsoid geometry. These are different surfaces — 24 m apart here.
    # The code never reconciles this. Three reference frames are mixed:
    # observer on the geoid, aircraft on geopotential, geometry on the ellipsoid.
    print(f'  Geoid undulation:          {GEOID_UNDULATION} m '
          f'[NOAA GEOID18, not accounted for]')

    # 7. Balloon spatial representativeness
    # Sun et al. (2010): temperature SD at individual pressure levels
    # increases by 0.42 K per 100 km of separation distance.
    # Conservative estimate: assume the error only affects a single ~1 km
    # layer at aircraft altitude. Via the hypsometric equation:
    # Δz = (R_d / g₀) × ΔT × ln(P_bottom / P_top)
    # For a single 1 km layer: ~4.1 m per 1 K of temperature error.
    # If the error extends over a broader column (likely at 200+ km),
    # it scales proportionally — up to ~39 m for a full-column error.
    temp_err_K = 0.42 * (balloon_to_ac / 100)  # K, scaled by distance
    single_level_factor = 4.1  # meters per K, conservative single-layer
    alt_repr_m = temp_err_K * single_level_factor
    print(f'  Balloon representativeness: ±{alt_repr_m:.1f} m '
          f'({temp_err_K:.2f}K × {single_level_factor} m/K, conservative) '
          f'[Sun et al. 2010]')

    # 8. Temporal mismatch (vertical component)
    # Same 2-5 second delay as horizontal, but applied to vertical rate.
    # Aircraft climbing or descending means the altitude is also wrong.
    temporal_vert = ac['vr_ms'] * ac['delta_s']
    print(f'  TEMPORAL MISMATCH (vert):  {temporal_vert:+.1f} m '
          f'({ac["delta_s"]}s × {ac["vr_ms"]:+.1f} m/s)')

    # ── Combine altitude errors ──
    # RSS the random components (altimeter + ADS-B quant + SSE + balloon)
    # Add systematic offsets separately (geopotential skip + geoid)
    # Temporal mismatch is directional — adds to one end of the range
    rss_alt = math.sqrt(alt_tol_m**2 + adsb_alt_m**2 + sse_m**2
                        + alt_repr_m**2)
    systematic_alt = geop_err + abs(GEOID_UNDULATION)
    total_alt_min = -rss_alt - abs(temporal_vert)
    total_alt_max = rss_alt + systematic_alt + abs(temporal_vert)

    print(f'  RSS (random):              ±{rss_alt:.1f} m')
    print(f'  Systematic (geop+geoid):   {geop_err:.1f} + '
          f'{abs(GEOID_UNDULATION):.1f} = {systematic_alt:.1f} m')
    print(f'  TOTAL ALT RANGE:           {total_alt_min:+.1f} to '
          f'{total_alt_max:+.1f} m')

    # ══════════════════════════════════════════════════════════
    # IMPACT ON ELEVATION ANGLE
    # ══════════════════════════════════════════════════════════
    # The elevation angle is what the raytracer actually uses.
    # We compute the worst-case bounds:
    # - Maximum elevation: aircraft as close as possible + as high as possible
    # - Minimum elevation: aircraft as far as possible + as low as possible
    # The half-width of that range is the elevation angle uncertainty.
    print(f'\n  ── IMPACT ON ELEVATION ANGLE ──')

    # Closest possible position: subtract horizontal error from distance
    # (floor at 1 km to avoid divide-by-zero edge cases)
    close_dist = max(nom_dist - total_horiz, 1000)
    close_alt = nom_alt_diff + total_alt_max  # highest possible altitude
    elev_close = math.degrees(math.atan2(close_alt, close_dist))

    # Farthest possible position: add horizontal error to distance
    far_dist = nom_dist + total_horiz
    far_alt = nom_alt_diff + total_alt_min  # lowest possible altitude
    elev_far = math.degrees(math.atan2(far_alt, far_dist))

    # Uncertainty is half the total range
    uncertainty = (elev_close - elev_far) / 2
    arcmin = (elev_close - elev_far) * 60  # convert to arcminutes

    print(f'  Nominal elevation:         {nom_elev:.4f}°')
    print(f'  Max elevation (close+high):{elev_close:.4f}°')
    print(f'  Min elevation (far+low):   {elev_far:.4f}°')
    print(f'  TOTAL RANGE:               {elev_far:.4f}° to '
          f'{elev_close:.4f}°')
    print(f'  UNCERTAINTY:               ±{uncertainty:.4f}° '
          f'({arcmin:.2f} arcmin)')
    print(f'  vs claimed globe error:    ~0.01-0.06°')

    return {
        'name': ac['name'],
        'fl': ac['flight_level'],
        'horiz': total_horiz,
        'alt_min': total_alt_min,
        'alt_max': total_alt_max,
        'elev_uncertainty': uncertainty,
        'elev_range': (elev_far, elev_close),
        'nom_elev': nom_elev,
    }


def main():
    print('=' * 90)
    print('  FULL-SPECTRUM ERROR ANALYSIS — AIRCRAFT POSITION UNCERTAINTY')
    print('=' * 90)
    print('  Sources: DO-260B, 14 CFR 43/91/25, Sun et al. 2010, NOAA GEOID18')

    results = [analyze(ac) for ac in AIRCRAFT]

    # ── Summary table ──
    print(f'\n\n{"=" * 90}')
    print('  SUMMARY')
    print(f'{"=" * 90}')
    print(f'  {"Aircraft":<28s} {"Horiz (m)":>10s} {"Vert Range (m)":>16s} '
          f'{"Elev ±°":>10s} {"Claimed":>10s}')
    print(f'  {"─" * 80}')
    for r in results:
        print(f'  {r["name"] + " " + r["fl"]:<28s} '
              f'{r["horiz"]:>10.0f} '
              f'{r["alt_min"]:>+8.0f}/{r["alt_max"]:>+.0f} '
              f'{r["elev_uncertainty"]:>10.4f} '
              f'{"0.01-0.06":>10s}')

    print(f'\n  In every case, measurement uncertainty exceeds the claimed ')
    print(f'  result by 3-6×. The error bars swallow the difference.')


if __name__ == '__main__':
    main()
