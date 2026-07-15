# ZwickerLoudness.jl

[![CI](https://github.com/slink/ZwickerLoudness.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/slink/ZwickerLoudness.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/slink/ZwickerLoudness.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/slink/ZwickerLoudness.jl)
[![DOI](https://zenodo.org/badge/1152458596.svg)](https://zenodo.org/badge/latestdoi/1152458596)

Zwicker loudness calculation per **ISO 532-1:2017**: the clause 5 method for
stationary sounds, and (since v0.3.0) the clause 6 method for time-varying
sounds.

Converts one-third-octave band SPL values into perceptual loudness in **sones** and **phons**, with 240-bin specific loudness at 0.1 Bark resolution.

## Installation

```julia
using Pkg
Pkg.add("ZwickerLoudness")
```

## Quick Start

```julia
using ZwickerLoudness

# 28 one-third-octave band SPL values [dB], 25 Hz to 12.5 kHz
spl = [60, 62, 65, 68, 70, 72, 74, 75, 73, 71,
       69, 67, 65, 63, 61, 59, 57, 55, 53, 50,
       47, 44, 41, 38, 35, 32, 29, 26.0]

result = zwicker_loudness(spl)

result.loudness           # total loudness [sone]
result.loudness_level     # loudness level [phon]
result.specific_loudness  # N'(z) at 0.1 Bark resolution [sone/Bark] (240 values)

# Diffuse field listening condition
result_diffuse = zwicker_loudness(spl; field_type=:diffuse)
```

Also accepts 31-band input (standard 20 Hz -- 20 kHz); bands 2--29 are extracted automatically.

## How It Works

The ISO 532-1:2017 pipeline:

```
28 one-third-octave bands (25 Hz -- 12.5 kHz)
    |
    v
Low-frequency correction (DLL table, equal-loudness contours)
    |
    v
20 critical band excitation levels (energy summation + ear transmission)
    |
    v
Core specific loudness (Zwicker power-law, s=0.25, korry correction)
    |
    v
240-bin spreading (USL/RNS slope tables, 0.1 Bark resolution)
    |
    v
Integration --> total loudness [sone] --> loudness level [phon]
```

**Key psychoacoustic effects modeled:**

- **Low-frequency correction** -- equal-loudness contour weighting for bands below 315 Hz
- **Threshold in quiet** -- frequency-dependent minimum audibility per critical band
- **Ear transmission** -- outer/middle ear transfer function (A0 correction)
- **Loudness growth** -- power-law compression with critical band adaptation (DCB)
- **Spreading function** -- upward spread of masking via tabulated upper slopes (USL)
- **Free/diffuse field** -- selectable listening condition (DDF correction)

## API

### `zwicker_loudness(spl; field_type=:free) -> ZwickerResult`

Compute Zwicker loudness from one-third-octave band SPL values.

**Input:** Vector of SPL values [dB] for 28 bands (25 Hz -- 12.5 kHz). Also accepts 31 bands (20 Hz -- 20 kHz), shorter vectors (padded with silence), or longer vectors (truncated).

**Keyword arguments:**

| Argument | Default | Description |
|----------|---------|-------------|
| `field_type` | `:free` | Listening condition: `:free` or `:diffuse` |

**Returns** a `ZwickerResult`:

| Field | Type | Description |
|-------|------|-------------|
| `loudness` | `Float64` | Total loudness [sone] |
| `loudness_level` | `Float64` | Loudness level [phon] |
| `specific_loudness` | `Vector{Float64}` | N'(z) at 0.1 Bark resolution [sone/Bark], 240 values |

## Conformance

Validated against ISO 532-1:2017 Annex B **Signal 1** (the only 28-band reference in Annex B) and cross-checked per-bin against the [MoSQITo](https://github.com/Eomys/MoSQITo) Python reference implementation.

| Test | Expected | Result |
|------|----------|--------|
| Signal 1 total `N` | 83.296 sone (+/-5%) | 83.296 sone |
| Signal 1 `N'(z)` (240 bins) | MoSQITo CSV reference | max abs diff < 0.001 sone/Bark |

> **Scope note.** Annex B Signals 2-5 are time-domain `.wav` files that require a third-octave filter bank (per IEC 61260) to convert to 28-band SPL before the Method 1 kernel implemented here can run. That preprocessing is out of scope for this package (no external dependencies beyond the `Statistics` stdlib); pair it with a third-octave filter bank to validate against Signals 2-5.

## Time-Varying Loudness (ISO 532-1:2017 clause 6)

Since v0.3.0, this package also implements **ISO 532-1:2017 clause 6, "method
for time-varying sounds"** -- colloquially "Method 2" in this project's own
shorthand, but that is *not* the standard's own terminology (the standard's
clause heading is literally "Method for time-varying sounds"; ISO never
labels it "Method 2").

> **Not to be confused with ISO 532-2.** `ISO 532-2:2017` is a completely
> different, unrelated document ("Acoustics -- Methods for calculating
> loudness -- Part 2: Moore-Glasberg method"). This package implements
> **ISO 532-1** only (both its clause 5 stationary method, above, and its
> clause 6 time-varying method, below) -- never ISO 532-2.

Clause 6 extends Method 1 to arbitrary non-stationary sounds (stationary
sounds being a special case) by adding nonlinear temporal decay, temporal
weighting, and time-domain percentile statistics (`N5`/`N10`) around the
same per-block core-loudness/spreading kernel Method 1 already implements.

### Kernel input contract

`zwicker_loudness_time_varying` is the **kernel** stage only: it takes a
`28 x nblocks` matrix of one-third-octave band SPL values, one column per
block at the standard's fixed **2000 Hz (0.5 ms)** block rate -- it does
**not** accept a raw audio signal. Producing that band-level time series
from a waveform (the signal-domain "front end": time-domain third-octave
filtering + square-and-smooth) is a separate concern, out of scope for this
zero-external-dependency package; it is planned for
[ZwickerLoudnessAudio.jl](https://github.com/slink/ZwickerLoudnessAudio.jl)'s
`loudness_zwtv(signal, fs)`.

```julia
using ZwickerLoudness

band_levels = fill(60.0, 28, 400)  # 400 blocks @ 2000 Hz = 0.2 s
result = zwicker_loudness_time_varying(band_levels)

result.loudness_over_time  # N(t) [sone], 2 ms axis (100 values for 400 input blocks)
result.N5                  # loudness exceeded 5% of the time [sone]
result.N10                 # loudness exceeded 10% of the time [sone]
result.specific_loudness   # N'(z, t) [sone/Bark], 240 x 100
result.time_axis           # [s], 2 ms steps
```

Real output for the snippet above:

```julia-repl
julia> length(result.loudness_over_time)
100

julia> result.loudness_over_time[1:5]
5-element Vector{Float64}:
  0.09412794403247204
  7.066389778630747
 11.197814560228862
 13.719501002870494
 15.32682682674735

julia> result.loudness_over_time[end]
30.93988817101539

julia> result.N5
30.787930815252135

julia> result.N10
30.612887643196228

julia> size(result.specific_loudness)
(240, 100)
```

(The rising values in `loudness_over_time[1:5]` reflect the nonlinear decay
network's attack behavior settling from silence toward the steady-state
60 dB loudness; `N5`/`N10` sit close to the settled value since most of the
0.2 s block is at steady state.)

### `zwicker_loudness_time_varying(band_levels; field_type=:free) -> ZwickerTimeVaryingResult`

**Input:** `28 x nblocks` matrix of one-third-octave band SPL values [dB],
one column per 0.5 ms (2000 Hz) block. The 2000 Hz block rate is fixed by
the standard, not a caller-supplied parameter.

**Keyword arguments:** same `field_type` as `zwicker_loudness`.

**Returns** a `ZwickerTimeVaryingResult`:

| Field | Type | Description |
|-------|------|--------------|
| `loudness_over_time` | `Vector{Float64}` | Time-varying total loudness `N(t)` [sone], 2 ms axis |
| `N5` | `Float64` | Loudness reached or exceeded 5% of the time [sone] |
| `N10` | `Float64` | Loudness reached or exceeded 10% of the time [sone] |
| `specific_loudness` | `Matrix{Float64}` | `N'(z, t)` [sone/Bark], 240 x `length(loudness_over_time)`; decimated but never temporally weighted |
| `time_axis` | `Vector{Float64}` | Time [s] for each output column, nominal 2 ms steps |

### Percentile definition

`N5`/`N10` are the ISO 532-1 / DIN 45631/A1 definition: the loudness value
reached or exceeded 5% (respectively 10%) of the signal's duration -- the
95th and 90th percentiles of `loudness_over_time`. MoSQITo computes no
percentiles anywhere in its loudness code, so there is no upstream
implementation to transcribe; this package computes them via Julia's
`Statistics.quantile`, which (like numpy's `percentile` default) uses the
**Hyndman-Fan type 7** ("linear") interpolation scheme. This matters because
MATLAB's `prctile` (used by SQAT, a common crosscheck target) instead
follows a Hazen-style (type 5) scheme that can diverge from type 7 by a few
percent on the same data -- a difference in *percentile definition*, not a
transcription error. If you are crosschecking `N5`/`N10` against a
MATLAB-based tool, expect this gap and attribute it accordingly.

### Conformance

The clause 6 kernel is unit-tested stage-by-stage (nonlinear decay, temporal
weighting, decimation, percentiles) against vendored MoSQITo `loudness_zwtv`
@ `d990c33f94f1` intermediates on synthesized signals (transcription-class
agreement, ~1e-9 to 1e-15), and additionally against **ISO 532-1:2017 Annex
B.4 "Test signal 6"**: since the fixture deliberately does not vendor the
Annex B band-level matrix (the underlying `.wav` is ISO-copyrighted material
MoSQITo merely redistributes -- see `test/fixtures/NOTICE.md`), this
end-to-end check regenerates `band_levels` locally from a MoSQITo checkout at
test time and self-skips (`@test_skip`, CI stays green) when that local
reference material isn't present. See
`test/test_method2_annexb_conformance.jl`.

## References

- ISO 532-1:2017 -- Acoustics: Methods for calculating loudness, Part 1: Zwicker method (clause 5: stationary sounds; clause 6: time-varying sounds)
- Zwicker, E. (1991). "Program for calculating loudness according to DIN 45631 (ISO 532B)". *J. Acoust. Soc. Jpn. (E)* 12, 1.
- Zwicker, E. & Fastl, H. (2007). *Psychoacoustics: Facts and Models*. Springer.
- MoSQITo project, `loudness_zwtv` @ commit `d990c33f94f1` (<https://github.com/Eomys/MoSQITo>, Apache License 2.0) -- transcription reference and crosscheck oracle for the clause 6 kernel.

## License

MIT
