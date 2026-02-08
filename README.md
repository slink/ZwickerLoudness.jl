# ZwickerLoudness.jl

[![CI](https://github.com/slink/ZwickerLoudness.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/slink/ZwickerLoudness.jl/actions/workflows/CI.yml)

Zwicker loudness calculation for stationary sounds per **ISO 532-1:2017** (Method 1).

Converts 28 one-third-octave band SPL values into perceptual loudness in **sones** and **phons**, with 240-bin specific loudness at 0.1 Bark resolution.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/slink/ZwickerLoudness.jl")
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

Validated against ISO 532-1:2017 Annex B test signals and cross-checked with the [MoSQITo](https://github.com/Eomys/MoSQITo) Python reference implementation.

| Test Signal | Expected [sone] | Result [sone] | Tolerance |
|-------------|-----------------|---------------|-----------|
| Annex B Signal 1 (broadband) | 83.296 | 83.296 | +/-5% |

## References

- ISO 532-1:2017 -- Acoustics: Methods for calculating loudness, Part 1: Zwicker method
- Zwicker, E. (1991). "Program for calculating loudness according to DIN 45631 (ISO 532B)". *J. Acoust. Soc. Jpn. (E)* 12, 1.
- Zwicker, E. & Fastl, H. (2007). *Psychoacoustics: Facts and Models*. Springer.

## License

MIT
