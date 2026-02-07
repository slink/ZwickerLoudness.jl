module ZwickerLoudness

export ZwickerResult, zwicker_loudness

"""
    ZwickerResult

Result of Zwicker loudness calculation (ISO 532B).
"""
struct ZwickerResult
    loudness::Float64                  # total loudness [sone]
    loudness_level::Float64            # loudness level [phon]
    specific_loudness::Vector{Float64} # N'(z) per critical band [sone/Bark]
end

# =========================================================================== #
#  ISO 532B Tabulated Constants
# =========================================================================== #

# Standard 1/3-octave center frequencies (bands 1-28: 25 Hz to 12.5 kHz)
const THIRD_OCTAVE_CENTERS_28 = Float64[
    25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200,
    250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000,
    2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500
]

# Mapping: which 1/3-octave bands (1-indexed into the 28-band array) contribute
# to each of the 24 critical bands (Bark scale). Some critical bands span
# multiple 1/3-octave bands; energy is summed.
#
# Based on Zwicker & Fastl "Psychoacoustics" Table 6.1 and ISO 532B Table 1.
# Each entry is a vector of 1/3-octave band indices that map to that Bark band.
const BARK_BAND_MAPPING = [
    [1, 2],       # Bark 1:  0-100 Hz   (25, 31.5 Hz)
    [3, 4],       # Bark 2:  100-200 Hz  (40, 50 Hz)
    [5, 6],       # Bark 3:  200-300 Hz  (63, 80 Hz)
    [7],          # Bark 4:  300-400 Hz  (100 Hz)
    [8],          # Bark 5:  400-510 Hz  (125 Hz)
    [9],          # Bark 6:  510-630 Hz  (160 Hz)
    [10],         # Bark 7:  630-770 Hz  (200 Hz)
    [11],         # Bark 8:  770-920 Hz  (250 Hz)
    [12],         # Bark 9:  920-1080 Hz (315 Hz)
    [13],         # Bark 10: 1080-1270 Hz (400 Hz)
    [14],         # Bark 11: 1270-1480 Hz (500 Hz)
    [15],         # Bark 12: 1480-1720 Hz (630 Hz)
    [16],         # Bark 13: 1720-2000 Hz (800 Hz)
    [17],         # Bark 14: 2000-2320 Hz (1000 Hz)
    [18],         # Bark 15: 2320-2700 Hz (1250 Hz)
    [19],         # Bark 16: 2700-3150 Hz (1600 Hz)
    [20],         # Bark 17: 3150-3700 Hz (2000 Hz)
    [21],         # Bark 18: 3700-4400 Hz (2500 Hz)
    [22],         # Bark 19: 4400-5300 Hz (3150 Hz)
    [23],         # Bark 20: 5300-6400 Hz (4000 Hz)
    [24],         # Bark 21: 6400-7700 Hz (5000 Hz)
    [25],         # Bark 22: 7700-9500 Hz (6300 Hz)
    [26],         # Bark 23: 9500-12000 Hz (8000 Hz)
    [27, 28],     # Bark 24: 12000-15500 Hz (10000, 12500 Hz)
]

# Threshold in quiet per critical band [dB SPL]
# From ISO 226 / Zwicker & Fastl Table 8.1
const THRESHOLD_IN_QUIET = Float64[
    60.0, 38.0, 26.0, 18.0, 14.0, 11.0, 8.5, 6.5, 5.5, 4.5,
    4.0,  3.5,  3.0,  2.5,  2.5,  3.0,  4.0, 5.5, 7.0, 9.5,
    13.0, 18.5, 26.0, 38.0
]

# Ear transfer function correction [dB] per critical band
# Free-field to eardrum transfer function (frontal incidence).
# Peaks around 2-4 kHz (Bark bands 17-20) due to ear canal resonance.
# At 1 kHz (Bark 14) correction is ~0 dB (reference frequency).
# Based on ISO 532B / DIN 45631 diffuse-field listening condition.
const EAR_TRANSFER = Float64[
    -1.0, -0.5,  0.0,  0.5,  0.5,  0.5,  0.0,  0.0, -0.5,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.5,  1.5,  3.0,  5.0,  5.5,  4.5,
     2.0, -1.0, -4.5, -8.0
]

# Critical band level factor 's' per Bark band
# Controls the steepness of the loudness growth function
# From Zwicker & Fastl / ISO 532B
const LEVEL_FACTOR_S = Float64[
    0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18,
    0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.20,
    0.21, 0.22, 0.24, 0.28
]

# Upper slope of masking pattern [dB/Bark] as function of excitation level
# For the upward spread of masking calculation.
# slope = USL_BASE + USL_SLOPE * (excitation_level - 40)
# Capped between USL_MIN and USL_MAX
const USL_BASE = 27.0   # slope at 40 dB
const USL_SLOPE = -0.2  # slope change per dB above 40
const USL_MIN = 10.0    # minimum slope (loud sounds)
const USL_MAX = 37.0    # maximum slope (quiet sounds)

# =========================================================================== #
#  Core Algorithm
# =========================================================================== #

"""
    map_to_critical_bands(spl_third::Vector{Float64}) -> Vector{Float64}

Map 28 1/3-octave band levels to 24 critical band levels by energy summation.
"""
function map_to_critical_bands(spl_third::Vector{Float64})
    cb_levels = zeros(Float64, 24)
    for (z, bands) in enumerate(BARK_BAND_MAPPING)
        # Energy sum of contributing 1/3-octave bands
        p2_sum = 0.0
        for b in bands
            p2_sum += 10.0^(spl_third[b] / 10.0)
        end
        cb_levels[z] = 10.0 * log10(max(p2_sum, 1e-30))
    end
    return cb_levels
end

"""
    apply_ear_transfer(cb_levels::Vector{Float64}) -> Vector{Float64}

Apply outer/middle ear transfer function correction to critical band levels.
Returns excitation levels.
"""
function apply_ear_transfer(cb_levels::Vector{Float64})
    return cb_levels .+ EAR_TRANSFER
end

"""
    core_loudness(excitation::Float64, z::Int) -> Float64

Compute specific loudness N' [sone/Bark] for critical band z given
excitation level [dB SPL]. Uses the Zwicker power-law model per
DIN 45631 / ISO 532B Method B.

The formula (working in dB domain):
    N' = 0.08 × 10^(0.023 × LTQ) × [(0.5 + 0.5 × 10^((LE-LTQ)/(10/s)))^0.23 - 1]

Calibrated so that 1 kHz at 40 dB SPL produces 1 sone total loudness.
"""
function core_loudness(excitation::Float64, z::Int)
    LTQ = THRESHOLD_IN_QUIET[z]
    s = LEVEL_FACTOR_S[z]

    # Below threshold → zero loudness
    excitation <= LTQ && return 0.0

    # Normalization factor: converts from excitation to specific loudness
    # ISO 532B: N' = C × (E_TQ/s)^0.23 × [(0.5 + 0.5×E/E_TQ)^0.23 - 1]
    # In dB: (E_TQ/s)^0.23 = 10^(0.023×LTQ) / s^0.23
    # C = 0.1133 calibrated so that 1 kHz pure tone at 40 dB SPL = 1.0 sone
    norm = 0.1133 * 10.0^(0.023 * LTQ) / s^0.23

    # Power-law loudness growth
    # E/E_TQ = 10^((LE-LTQ)/10) in linear power domain
    inner = 0.5 + 0.5 * 10.0^((excitation - LTQ) / 10.0)
    N_prime = norm * (inner^0.23 - 1.0)

    return max(N_prime, 0.0)
end

"""
    upper_slope(excitation::Float64) -> Float64

Compute the upper slope [dB/Bark] of the masking pattern for a given
excitation level. The slope decreases with increasing level (louder
sounds mask more broadly).
"""
function upper_slope(excitation::Float64)
    slope = USL_BASE + USL_SLOPE * (excitation - 40.0)
    return clamp(slope, USL_MIN, USL_MAX)
end

"""
    apply_masking!(specific_loudness::Vector{Float64}, excitation_levels::Vector{Float64})

Apply upward spread of masking: each critical band's excitation spreads
into higher bands, reducing the specific loudness contribution of higher
bands that are masked by lower-frequency sounds.

Modifies `specific_loudness` in place.
"""
function apply_masking!(specific_loudness::Vector{Float64}, excitation_levels::Vector{Float64})
    n_bands = length(specific_loudness)

    for z in 1:(n_bands - 1)
        excitation_levels[z] <= THRESHOLD_IN_QUIET[z] && continue

        slope = upper_slope(excitation_levels[z])

        # Spread into higher bands
        for zz in (z + 1):n_bands
            delta_z = zz - z  # Bark distance
            masking_reduction = slope * delta_z  # dB reduction per Bark

            # Masking threshold from band z at position zz
            masked_level = excitation_levels[z] - masking_reduction

            # If the masked level exceeds the actual excitation at zz,
            # the signal at zz is (partially) masked. We reduce its
            # specific loudness contribution proportionally.
            if masked_level > excitation_levels[zz]
                # Signal is fully masked at this distance
                # Reduce specific loudness at zz
                mask_excess = masked_level - excitation_levels[zz]
                # Exponential reduction factor
                reduction = 10.0^(-mask_excess / 20.0)
                specific_loudness[zz] *= reduction
            end
        end
    end
end

"""
    sones_to_phons(loudness::Float64) -> Float64

Convert loudness in sones to loudness level in phons.
N = 2^((LN - 40)/10)  →  LN = 40 + 10 × log₂(N)
"""
function sones_to_phons(loudness::Float64)
    loudness <= 0.0 && return 0.0
    return 40.0 + 10.0 * log2(loudness)
end

"""
    phons_to_sones(phons::Float64) -> Float64

Convert loudness level in phons to loudness in sones.
LN = 40 + 10 × log₂(N)  →  N = 2^((LN - 40)/10)
"""
function phons_to_sones(phons::Float64)
    return 2.0^((phons - 40.0) / 10.0)
end

# =========================================================================== #
#  Main Entry Point
# =========================================================================== #

"""
    zwicker_loudness(spl_third_octave::Vector{Float64}) -> ZwickerResult

Compute Zwicker loudness (ISO 532B, stationary sound, free-field) from
1/3-octave band SPL values.

# Input
- `spl_third_octave`: SPL values [dB] for 28 1/3-octave bands
  (center frequencies 25 Hz to 12.5 kHz). Accepts vectors of length 28
  directly, or length 31 (the standard 20 Hz–20 kHz set, from which
  bands 2-29 are extracted).

# Returns
A `ZwickerResult` with:
- `loudness`: total loudness [sone]
- `loudness_level`: loudness level [phon]
- `specific_loudness`: N'(z) per critical band [sone/Bark] (24 values)
"""
function zwicker_loudness(spl_third_octave::Vector{Float64})
    n = length(spl_third_octave)

    # Handle different input lengths
    if n == 28
        spl_28 = spl_third_octave
    elseif n == 31
        # Standard 31-band (20 Hz–20 kHz): extract bands 2-29 (25 Hz–12.5 kHz)
        spl_28 = spl_third_octave[2:29]
    elseif n > 28
        # Take the first 28 bands
        spl_28 = spl_third_octave[1:28]
    else
        # Pad with silence (below threshold)
        spl_28 = vcat(spl_third_octave, fill(-Inf, 28 - n))
    end

    # Step 1: Map 1/3-octave bands to 24 critical bands
    cb_levels = map_to_critical_bands(spl_28)

    # Step 2: Apply ear transfer function
    excitation_levels = apply_ear_transfer(cb_levels)

    # Step 3: Compute core specific loudness per band
    specific_loudness = [core_loudness(excitation_levels[z], z) for z in 1:24]

    # Step 4: Apply upward spread of masking
    apply_masking!(specific_loudness, excitation_levels)

    # Step 5: Integrate specific loudness → total loudness [sone]
    # Each critical band is 1 Bark wide, so N = Σ N'(z) × Δz = Σ N'(z)
    loudness = sum(specific_loudness)

    # Step 6: Convert to phons
    loudness_level = sones_to_phons(loudness)

    return ZwickerResult(loudness, loudness_level, specific_loudness)
end

end # module
