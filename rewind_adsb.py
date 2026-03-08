"""
Rewind ADS-B positions to the exact occultation time and recompute geometry.

The ADS-B "Updated" timestamp is always AFTER the occultation frame.
This script extrapolates the aircraft position backwards using
ground speed, track heading, and vertical rate to estimate where
the plane actually was at the moment of occultation.
"""

import math
from geographiclib.geodesic import Geodesic

def rewind_position(lat, lon, alt_m, gs_ms, trk_deg, vr_ms, delta_s):
    """
    Given an ADS-B position at time T+delta, rewind to time T.
    Move backwards along track heading by gs * delta.
    Adjust altitude by vr * delta.
    """
    # Distance to rewind (meters)
    dist_back = gs_ms * delta_s

    # Bearing is opposite of track heading
    bearing_back = (trk_deg + 180.0) % 360.0

    # Use geodesic forward calculation to get rewound lat/lon
    result = Geodesic.WGS84.Direct(lat, lon, bearing_back, dist_back)
    new_lat = result['lat2']
    new_lon = result['lon2']

    # Rewind altitude (subtract the climb/descent that happened in delta)
    new_alt = alt_m - (vr_ms * delta_s)

    return new_lat, new_lon, new_alt


def compute_geometry(obs_lat, obs_lon, obs_alt, tgt_lat, tgt_lon, tgt_alt):
    """
    Compute surface distance and geometric elevation angle from observer to target.
    """
    inv = Geodesic.WGS84.Inverse(obs_lat, obs_lon, tgt_lat, tgt_lon)
    surface_dist = inv['s12']  # meters

    # Geometric elevation angle (no refraction)
    alt_diff = tgt_alt - obs_alt
    elev_rad = math.atan2(alt_diff, surface_dist)
    elev_deg = math.degrees(elev_rad)

    return surface_dist, elev_deg, alt_diff


def analyze_event(name, obs_lat, obs_lon, obs_alt,
                  adsb_lat, adsb_lon, adsb_alt,
                  gs_ms, trk_deg, vr_ms, delta_s):

    print(f"\n{'='*70}")
    print(f"  {name}")
    print(f"{'='*70}")

    # Original (uncorrected) geometry
    orig_dist, orig_elev, orig_alt_diff = compute_geometry(
        obs_lat, obs_lon, obs_alt,
        adsb_lat, adsb_lon, adsb_alt
    )

    # Rewind position
    rw_lat, rw_lon, rw_alt = rewind_position(
        adsb_lat, adsb_lon, adsb_alt,
        gs_ms, trk_deg, vr_ms, delta_s
    )

    # Corrected geometry
    corr_dist, corr_elev, corr_alt_diff = compute_geometry(
        obs_lat, obs_lon, obs_alt,
        rw_lat, rw_lon, rw_alt
    )

    print(f"\n  ADS-B Position (T+{delta_s}s):  {adsb_lat:.6f}, {adsb_lon:.6f}, {adsb_alt:.1f}m")
    print(f"  Rewound Position (T):       {rw_lat:.6f}, {rw_lon:.6f}, {rw_alt:.1f}m")
    print(f"  GS: {gs_ms:.1f} m/s | VR: {vr_ms:+.2f} m/s | Trk: {trk_deg}° | Delta: {delta_s}s")

    print(f"\n  {'':30s} {'Original':>12s}  {'Corrected':>12s}  {'Difference':>12s}")
    print(f"  {'─'*70}")

    print(f"  {'Lat (°)':30s} {adsb_lat:12.6f}  {rw_lat:12.6f}  {rw_lat - adsb_lat:+12.6f}")
    print(f"  {'Lon (°)':30s} {adsb_lon:12.6f}  {rw_lon:12.6f}  {rw_lon - adsb_lon:+12.6f}")
    print(f"  {'Alt (m)':30s} {adsb_alt:12.1f}  {rw_alt:12.1f}  {rw_alt - adsb_alt:+12.1f}")
    print(f"  {'Surface Distance (m)':30s} {orig_dist:12.1f}  {corr_dist:12.1f}  {corr_dist - orig_dist:+12.1f}")
    print(f"  {'Surface Distance (km)':30s} {orig_dist/1000:12.3f}  {corr_dist/1000:12.3f}  {(corr_dist - orig_dist)/1000:+12.3f}")
    print(f"  {'Alt above observer (m)':30s} {orig_alt_diff:12.1f}  {corr_alt_diff:12.1f}  {corr_alt_diff - orig_alt_diff:+12.1f}")
    print(f"  {'Geometric Elevation (°)':30s} {orig_elev:12.6f}  {corr_elev:12.6f}  {corr_elev - orig_elev:+12.6f}")

    return {
        'name': name,
        'orig_dist': orig_dist,
        'corr_dist': corr_dist,
        'orig_elev': orig_elev,
        'corr_elev': corr_elev,
        'elev_diff': corr_elev - orig_elev,
        'dist_diff': corr_dist - orig_dist,
        'rw_lat': rw_lat,
        'rw_lon': rw_lon,
        'rw_alt': rw_alt,
    }


# ============================================================
# Event data from ADS-B screenshots and CSV
# ============================================================

results = []

# Allegiant 453
# ADS-B: 27.0581, -81.5134, CSV alt: 9466.84m
# GS: 418 kn = 215.0 m/s, VR: +1024 ft/m = +5.20 m/s, Trk: 322, Delta: 4s
results.append(analyze_event(
    name="Allegiant 453 (AAY453)",
    obs_lat=27.06894, obs_lon=-82.22377, obs_alt=4.42,
    adsb_lat=27.0581, adsb_lon=-81.5134, adsb_alt=9466.84,
    gs_ms=215.0, trk_deg=322, vr_ms=5.20, delta_s=4
))

# EJM691 (Jet Speed 691)
# ADS-B: 26.5664, -82.3093, CSV alt: 11898.24m
# GS: 559 kn = 287.6 m/s, VR: -2816 ft/m = -14.31 m/s, Trk: 125, Delta: 3s
results.append(analyze_event(
    name="EJM691 (Jet Speed 691)",
    obs_lat=27.06765, obs_lon=-82.22612, obs_alt=6.2,
    adsb_lat=26.5664, adsb_lon=-82.3093, adsb_alt=11898.24,
    gs_ms=287.6, trk_deg=125, vr_ms=-14.31, delta_s=3
))

# AAL1517 (American Airlines 1517)
# ADS-B: 26.3003, -81.9388, CSV alt: 9392.17m
# GS: 535 kn = 275.2 m/s, VR: -960 ft/m = -4.88 m/s, Trk: 114, Delta: 5s
results.append(analyze_event(
    name="AAL1517 (American Airlines 1517)",
    obs_lat=27.06765, obs_lon=-82.22612, obs_alt=6.2,
    adsb_lat=26.3003, adsb_lon=-81.9388, adsb_alt=9392.17,
    gs_ms=275.2, trk_deg=114, vr_ms=-4.88, delta_s=5
))

# Spirit Airlines — no GS available, skip horizontal correction
# But we can still correct altitude
print(f"\n{'='*70}")
print(f"  Spirit Airlines (NKS395) — altitude only, no GS")
print(f"{'='*70}")
print(f"  VR: +4.88 m/s, Delta: 2s")
print(f"  Alt correction: -{4.88*2:.1f}m")
print(f"  CSV alt: 10375.54m -> corrected: {10375.54 - 4.88*2:.1f}m")
obs_lat, obs_lon, obs_alt = 27.06811, -82.22231, -1.07
tgt_lat, tgt_lon = 27.2931, -81.9698
orig_dist, orig_elev, _ = compute_geometry(obs_lat, obs_lon, obs_alt, tgt_lat, tgt_lon, 10375.54)
corr_alt = 10375.54 - 4.88*2
corr_dist, corr_elev, _ = compute_geometry(obs_lat, obs_lon, obs_alt, tgt_lat, tgt_lon, corr_alt)
print(f"  Geometric Elevation: {orig_elev:.6f}° -> {corr_elev:.6f}° (diff: {corr_elev - orig_elev:+.6f}°)")
print(f"  (Horizontal position unchanged — GS not available)")
