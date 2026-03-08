# ADS-B Error Analysis — AstroLive Aircraft Occultations

AstroLive recorded four aircraft stellar occultations using a telescope with timestamped video and an RTL-SDR ADS-B receiver (SDRAngel). The aircraft positions from those ADS-B broadcasts were fed into a raytracer to compare globe vs flat earth model predictions. The experiment claims globe model errors of 0.01-0.06° in elevation angle.

This script quantifies every error source in the aircraft position chain — from GPS fix to broadcast encoding to the timestamp mismatch — using primary source specifications. The result: the measurement uncertainty is 0.16-0.35°, which is 3-6x larger than the claimed precision. The error bars on the input data swallow the difference the experiment is trying to measure.

## What it calculates

### Horizontal position errors
- **GPS accuracy (NACp)** — DO-260B Table 2-73 defines accuracy categories. Most commercial jets are NACp=9 (±30m at 95% confidence). The FAA minimum for ADS-B Out is NACp=8 (±92.6m).
- **CPR quantization** — ADS-B encodes positions with 17 bits per coordinate. Each ~6° zone gets 131,072 bins, giving ~5m resolution per bin. Derived from the bit width, not directly stated in DO-260B.
- **Uncompensated latency** — DO-260B §2.2.3.2.7 allows up to 0.6 seconds of uncompensated delay between the GPS fix and the broadcast. At cruise speed that's 129-173m of position error baked into the broadcast itself.
- **Temporal mismatch** — The ADS-B "Updated" timestamp in every screenshot is 2-5 seconds after the occultation frame. The raytracer uses these delayed positions with no correction for ground speed or heading. At cruise speed: 860-1,376m of drift. This is the largest error source.

### Altitude errors
- **Altimeter scale error** — 14 CFR Part 43 Appendix E Table I. At FL308-FL390: ±55 to ±70m allowed.
- **ADS-B altitude quantization** — Mode S encodes altitude in 25 ft steps (Q-bit = 1 below 50,175 ft). Rounding to nearest step gives ±12.5 ft (±3.8m).
- **Static Source Error (SSE)** — 14 CFR 25.1325(e): ±30 ft per 100 knots. At ~300 kn cruise IAS: ±90 ft (±27m). This is airflow distortion at the static ports — doesn't affect a weather balloon floating in still air.
- **RVSM total system error** — 14 CFR Part 91 Appendix G: mean + 3σ ≤ ±200 ft (±61m) for the entire altimetry system.
- **Geopotential-to-geometric conversion (skipped)** — The raytracer only does this conversion for weather balloon targets, not aircraft. At 10-12 km altitude the difference is 14-22m. The aircraft is always higher than assumed.
- **Geoid undulation** — NOAA GEOID18 gives N = -24.022m at the observer location. MSL is 24m below the WGS84 ellipsoid. The raytracer does ellipsoid geometry with an MSL altitude and never accounts for the offset.
- **Balloon representativeness** — The radiosonde is from Miami (WMO 72202), 167-234 km away. Per Sun et al. (2010), temperature SD increases by 0.42 K per 100 km. Conservative single-level estimate: ±3-4m altitude error.
- **Temporal mismatch (vertical)** — Same 2-5s delay applied to vertical rate. Aircraft climbing/descending at 5-14 m/s means 10-43m of altitude drift.

### Combined result
The script RSS-combines random errors, adds systematic offsets, and propagates everything into an elevation angle uncertainty range. The half-width of that range is the measurement precision.

## Results

| Aircraft | Horiz. Error | Vert. Range | Elevation Uncertainty | Claimed |
|----------|-------------|-------------|----------------------|---------|
| Spirit Airlines FL340 | ±30m | -77 to +118m | **±0.16°** | 0.01-0.06° |
| Allegiant 453 FL310 | 993m | -84 to +122m | **±0.19°** | 0.01-0.06° |
| Jet Speed 691 FL390 | 1,038m | -117 to +163m | **±0.35°** | 0.01-0.06° |
| AAL1517 FL308 | 1,544m | -87 to +125m | **±0.17°** | 0.01-0.06° |

## Requirements

```
pip install geographiclib
```

Python 3.6+ (f-strings).

## Usage

```
python error_analysis.py
```

No external data files needed — all event data is hardcoded from the original ADS-B screenshots and raytracer CSV.

## Primary sources

| Document | What it defines |
|----------|----------------|
| [RTCA DO-260B](https://mode-s.org/1090mhz/content/ads-b/7-uncertainty.html) | NACp accuracy categories, CPR encoding, latency limits |
| [14 CFR Part 43, App E](https://www.law.cornell.edu/cfr/text/14/appendix-E_to_part_43) | Altimeter scale error tolerances (Table I) |
| [14 CFR Part 91, App G](https://www.law.cornell.edu/cfr/text/14/appendix-G_to_part_91) | RVSM altimetry system error limits |
| [14 CFR 25.1325](https://www.ecfr.gov/current/title-14/section-25.1325) | Static source error allowance |
| [14 CFR 91.227](https://www.ecfr.gov/current/title-14/chapter-I/subchapter-F/part-91/subpart-E/section-91.227) | ADS-B Out equipment requirements |
| [Sun et al. (2010)](https://www.arl.noaa.gov/documents/JournalPDFs/SunEtAl.JGR2010.JDO14457.pdf) | Radiosonde spatial representativeness (0.42 K/100 km) |
| [NOAA GEOID18](https://geodesy.noaa.gov/GEOID/GEOID18/computation.html) | Geoid undulation at observer location |
| [ADS-B altitude encoding](https://mode-s.org/1090mhz/content/ads-b/3-airborne-position.html) | Mode S Q-bit, 25 ft step encoding |
