# ZwickerLoudness.jl

[![CI](https://github.com/slink/ZwickerLoudness.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/slink/ZwickerLoudness.jl/actions/workflows/CI.yml)

Zwicker loudness calculation for stationary sounds per **ISO 532B** (Method B, simplified).

Converts 28 one-third-octave band SPL values into perceptual loudness in **sones** and **phons**, accounting for the frequency-dependent sensitivity of human hearing.

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
result.specific_loudness  # N'(z) per critical band [sone/Bark] (24 values)
```

Also accepts 31-band input (standard 20 Hz -- 20 kHz); bands 2--29 are extracted automatically.

## How It Works

The ISO 532B pipeline:

```
28 one-third-octave bands (25 Hz -- 12.5 kHz)
    |
    v
24 Bark critical bands (energy summation)
    |
    v
Ear transfer function (outer/middle ear correction)
    |
    v
Core specific loudness (Zwicker power-law model)
    |
    v
Upward spread of masking
    |
    v
Integration --> total loudness [sone] --> loudness level [phon]
```

**Key psychoacoustic effects modeled:**

- **Threshold in quiet** -- frequency-dependent minimum audibility (ISO 226)
- **Ear transfer function** -- resonance gain of the ear canal, peaking at 2--4 kHz
- **Loudness growth** -- power-law compression (Stevens' law), calibrated so 1 kHz at 40 dB SPL = 1 sone = 40 phon
- **Upward spread of masking** -- low-frequency sounds reduce the perceived loudness of higher-frequency sounds

## API

### `zwicker_loudness(spl::Vector{Float64}) -> ZwickerResult`

Compute Zwicker loudness from one-third-octave band SPL values.

**Input:** Vector of SPL values [dB] for 28 bands (25 Hz -- 12.5 kHz). Also accepts 31 bands (20 Hz -- 20 kHz), shorter vectors (padded with silence), or longer vectors (truncated).

**Returns** a `ZwickerResult`:

| Field | Type | Description |
|-------|------|-------------|
| `loudness` | `Float64` | Total loudness [sone] |
| `loudness_level` | `Float64` | Loudness level [phon] |
| `specific_loudness` | `Vector{Float64}` | N'(z) per critical band [sone/Bark], 24 values |

## References

- ISO 532:1975 (Method B) -- Acoustics: Method for calculating loudness level
- DIN 45631 -- Berechnung des Lautstarkepegels und der Lautheit
- Zwicker, E. & Fastl, H. (2007). *Psychoacoustics: Facts and Models*. Springer.

## License

MIT
