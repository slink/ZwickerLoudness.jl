module ZwickerLoudness

using Statistics: quantile

export ZwickerResult, zwicker_loudness

"""
    ZwickerResult

Result of an ISO 532-1:2017 Method 1 loudness calculation.

# Fields
- `loudness::Float64`: total loudness `N` in sone.
- `loudness_level::Float64`: loudness level `L_N` in phon.
- `specific_loudness::Vector{Float64}`: specific loudness `N'(z)` in sone/Bark,
  240 values at 0.1 Bark resolution (z = 0.1, 0.2, …, 24.0 Bark).
"""
struct ZwickerResult
    loudness::Float64
    loudness_level::Float64
    specific_loudness::Vector{Float64}
end

# =========================================================================== #
#  ISO 532-1:2017 Tables (Zwicker:1991 / ISO 532-1:2017 Annex A)
#
#  Provenance: these tables (RAP/DLL/LTQ/A0/DDF/DCB/ZUP/RNS/USL) and the
#  Method 1 pipeline they parameterize (correct_low_frequencies ->
#  compute_excitation_levels -> compute_core_loudness -> compute_spreading)
#  trace to Zwicker, E. (1991), "Program for calculating loudness according
#  to DIN 45631 (ISO 532B)", J. Acoust. Soc. Jpn. (E) 12, 1 -- the original
#  BASIC listing published alongside ISO 532B/DIN 45631, whose constants and
#  step structure ISO 532-1:2017 Annex A (the normative "computer program")
#  carries forward essentially unchanged. This package's implementation was
#  additionally cross-checked, table-for-table and step-for-step, against
#  the MoSQITo project's vectorized transcription of the same algorithm
#  (`_main_loudness.py` / `_calc_slopes.py`, Apache-2.0, see
#  `src/method2.jl`'s module docstring and `test/fixtures/NOTICE.md` for the
#  pinned commit) during the time-varying (clause 6) work, which confirmed
#  table-for-table and formula-for-formula equivalence with the code below.
# =========================================================================== #

const RAP = Float64[45, 55, 65, 71, 80, 90, 100, 120]

const DLL = Float64[
    -32 -24 -16 -10 -5  0 -7 -3  0 -2  0;
    -29 -22 -15 -10 -4  0 -7 -2  0 -2  0;
    -27 -19 -14  -9 -4  0 -6 -2  0 -2  0;
    -25 -17 -12  -9 -3  0 -5 -2  0 -2  0;
    -23 -16 -11  -7 -3  0 -4 -1  0 -1  0;
    -20 -14 -10  -6 -3  0 -4 -1  0 -1  0;
    -18 -12  -9  -6 -2  0 -3 -1  0 -1  0;
    -15 -10  -8  -4 -2  0 -3 -1  0 -1  0
]

const LTQ = Float64[30, 18, 12, 8, 7, 6, 5, 4, 3, 3,
                      3,  3,  3, 3, 3, 3, 3, 3, 3, 3]

const A0 = Float64[0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   -0.5, -1.6, -3.2, -5.4, -5.6, -4, -1.5, 2, 5, 12]

const DDF = Float64[0, 0, 0.5, 0.9, 1.2, 1.6, 2.3, 2.8, 3, 2,
                    0, -1.4, -2, -1.9, -1, 0.5, 3, 4, 4.3, 4]

const DCB = Float64[-0.25, -0.6, -0.8, -0.8, -0.5, 0, 0.5, 1.1, 1.5, 1.7,
                     1.8, 1.8, 1.7, 1.6, 1.4, 1.2, 0.8, 0.5, 0, -0.5]

const ZUP = Float64[0.9, 1.8, 2.8, 3.5, 4.4, 5.4, 6.6, 7.9, 9.2, 10.6,
                    12.3, 13.8, 15.2, 16.7, 18.1, 19.3, 20.6, 21.8, 22.7, 23.6, 24.0]

const RNS = Float64[21.5, 18, 15.1, 11.5, 9, 6.1, 4.4, 3.1,
                    2.13, 1.36, 0.82, 0.42, 0.30, 0.22, 0.15, 0.10, 0.035, 0]

const USL = Float64[
    13.0  8.2  6.3  5.5  5.5  5.5  5.5  5.5;
     9.0  7.5  6.0  5.1  4.5  4.5  4.5  4.5;
     7.8  6.7  5.6  4.9  4.4  3.9  3.9  3.9;
     6.2  5.4  4.6  4.0  3.5  3.2  3.2  3.2;
     4.5  3.8  3.6  3.2  2.9  2.7  2.7  2.7;
     3.7  3.0  2.8  2.35 2.2  2.2  2.2  2.2;
     2.9  2.3  2.1  1.9  1.8  1.7  1.7  1.7;
     2.4  1.7  1.5  1.35 1.3  1.3  1.3  1.3;
     1.95 1.45 1.3  1.15 1.1  1.1  1.1  1.1;
     1.5  1.2  0.94 0.86 0.82 0.82 0.82 0.82;
     0.72 0.67 0.64 0.63 0.62 0.62 0.62 0.62;
     0.59 0.53 0.51 0.50 0.42 0.42 0.42 0.42;
     0.40 0.33 0.26 0.24 0.24 0.22 0.22 0.22;
     0.27 0.21 0.20 0.18 0.17 0.17 0.17 0.17;
     0.16 0.15 0.14 0.12 0.11 0.11 0.11 0.11;
     0.12 0.11 0.10 0.08 0.08 0.08 0.08 0.08;
     0.09 0.08 0.07 0.06 0.06 0.06 0.06 0.05;
     0.06 0.05 0.03 0.02 0.02 0.02 0.02 0.02
]

# =========================================================================== #
#  Helpers
# =========================================================================== #

# Find RNS index: position in decreasing RNS array where n fits.
# equal_too=false: find first j where RNS[j] <= n
# equal_too=true:  find first j where RNS[j] < n (strictly)
function _rns_index(n::Float64; equal_too::Bool=false)
    for j in 1:17
        if equal_too
            RNS[j] < n && return j
        else
            RNS[j] <= n && return j
        end
    end
    return 18
end

function _usl_value(rns_idx::Int, ig::Int)
    return USL[rns_idx, min(ig, 8)]
end

# =========================================================================== #
#  Step 1: Low-Frequency Correction (bands 1-11, 25 Hz to 250 Hz)
# =========================================================================== #

function correct_low_frequencies(spl_11::AbstractVector{Float64})
    n = length(spl_11)
    corrected = copy(spl_11)
    ti = zeros(Float64, n)

    for i in 1:n
        k = 1
        for k_test in 1:7
            if spl_11[i] > (RAP[k_test] - DLL[k_test, i])
                k = k_test + 1
            else
                break
            end
        end
        corrected[i] = spl_11[i] + DLL[k, i]
        ti[i] = 10.0^(corrected[i] / 10.0)
    end

    return corrected, ti
end

# =========================================================================== #
#  Step 2: Compute Excitation Levels (20 critical bands)
# =========================================================================== #

function compute_excitation_levels(spl_28::AbstractVector{Float64}, field_type::Symbol)
    _, ti = correct_low_frequencies(spl_28[1:11])

    # Sum intensities into three critical band levels
    gi = zeros(Float64, 3)
    gi[1] = sum(ti[1:6])
    gi[2] = sum(ti[7:9])
    gi[3] = sum(ti[10:11])

    lcb = zeros(Float64, 3)
    for i in 1:3
        gi[i] > 0.0 && (lcb[i] = 10.0 * log10(gi[i]))
    end

    # Build 20 excitation levels: bands 9-28, replacing first 3 with LCB
    le = zeros(Float64, 20)
    for i in 1:20
        le[i] = spl_28[i + 8]
    end
    le[1] = lcb[1]
    le[2] = lcb[2]
    le[3] = lcb[3]

    # Subtract ear transmission correction
    for i in 1:20
        le[i] -= A0[i]
    end

    # Add diffuse field correction if needed
    if field_type == :diffuse
        for i in 1:20
            le[i] += DDF[i]
        end
    end

    return le
end

# =========================================================================== #
#  Step 3: Core Loudness (20 bands + 1 zero terminator = 21 values)
# =========================================================================== #

function compute_core_loudness(le::Vector{Float64})
    s = 0.25
    nm = zeros(Float64, 21)

    for i in 1:20
        if le[i] > LTQ[i]
            le_corr = le[i] - DCB[i]
            mp1 = 0.0635 * 10.0^(0.025 * LTQ[i])
            mp2 = (1.0 - s + s * 10.0^(0.1 * (le_corr - LTQ[i])))^0.25 - 1.0
            nm[i] = max(mp1 * mp2, 0.0)
        end
    end

    # Korry correction for lowest critical band
    korry = 0.4 + 0.32 * nm[1]^0.2
    if korry <= 1.0
        nm[1] *= korry
    end

    return nm
end

# =========================================================================== #
#  Step 4: Spreading Function (240-bin specific loudness + total loudness)
# =========================================================================== #

function compute_spreading(nm::Vector{Float64})
    N = 0.0
    ns = zeros(Float64, 240)
    dec = 8

    zup_ea = [round(Int, ZUP[i] * 10) for i in 1:21]

    n1 = 0.0
    z1 = 0.0

    for i in 1:21
        ig = clamp(i - 1, 1, 8)

        lo_bin = (i == 1) ? 1 : zup_ea[i-1] + 1
        hi_bin = zup_ea[i]

        if round(n1, digits=dec) <= round(nm[i], digits=dec)
            # Rising/flat: rectangular integration
            N += nm[i] * (ZUP[i] - z1)
            for j in lo_bin:hi_bin
                ns[j] = nm[i]
            end
            n1 = nm[i]
            z1 = ZUP[i]
        else
            # Falling: slope decay from n1 toward nm[i]
            j_rns = _rns_index(n1)
            usl_val = _usl_value(j_rns, ig)

            n2 = max(RNS[j_rns], nm[i])
            dz = (n1 - n2) / usl_val
            z2 = z1 + dz
            if z2 > ZUP[i]
                z2 = ZUP[i]
                dz = z2 - z1
                n2 = n1 - dz * usl_val
            end
            N += dz * (n1 + n2) / 2.0

            # Fill 0.1-Bark bins
            z = lo_bin * 0.1
            done = false
            for j in lo_bin:hi_bin
                if round(z2, digits=dec) > round(z, digits=dec)
                    ns[j] = max(n1 - (z - z1) * usl_val, 0.0)
                else
                    # End of slope segment; start new one
                    n1 = n2
                    z1 = z2

                    if round(n1, digits=dec) <= round(nm[i], digits=dec)
                        # Transitioned to flat
                        N += nm[i] * (ZUP[i] - z1)
                        for k in j:hi_bin
                            ns[k] = nm[i]
                        end
                        n1 = nm[i]
                        z1 = ZUP[i]
                        done = true
                        break
                    end

                    # New decay segment
                    j_rns = _rns_index(n1; equal_too=true)
                    usl_val = _usl_value(j_rns, ig)
                    n2 = max(RNS[j_rns], nm[i])
                    dz = (n1 - n2) / usl_val
                    z2 = z1 + dz
                    if z2 > ZUP[i]
                        z2 = ZUP[i]
                        dz = z2 - z1
                        n2 = n1 - dz * usl_val
                    end
                    N += dz * (n1 + n2) / 2.0

                    ns[j] = max(n1 - (z - z1) * usl_val, 0.0)
                end
                z += 0.1
            end

            if !done
                n1 = n2
                z1 = z2
            end
        end
    end

    N = max(N, 0.0)
    if N <= 16.0
        N = floor(N * 1000.0 + 0.5) / 1000.0
    else
        N = floor(N * 100.0 + 0.5) / 100.0
    end

    return N, ns
end

# =========================================================================== #
#  Sone-to-Phon Conversion (ISO 532-1:2017 two-branch formula)
# =========================================================================== #

"""
    sones_to_phons(N) -> Float64

Convert total loudness `N` in sone to loudness level `L_N` in phon using the
two-branch ISO 532-1:2017 formula:

- `N ≥ 1`: `L_N = 40 + 10·log₂(N)` (Stevens' power law).
- `0 < N < 1`: `L_N = 40·(N + 0.0005)^0.35`, floored at 3 phon.
- `N ≤ 0`: returns `0.0`.
"""
function sones_to_phons(N::Float64)
    N <= 0.0 && return 0.0
    if N >= 1.0
        return 40.0 + 10.0 * log2(N)
    else
        LN = 40.0 * (N + 0.0005)^0.35
        return max(LN, 3.0)
    end
end

# =========================================================================== #
#  Main Entry Point
# =========================================================================== #

"""
    zwicker_loudness(spl; field_type=:free) -> ZwickerResult

Compute Zwicker loudness from one-third-octave band SPL values per
ISO 532-1:2017 Method 1.

`spl` is a vector of band SPL values in dB (re. 20 µPa). Accepts:
- 28 values for the 25 Hz – 12.5 kHz range (preferred),
- 31 values for the 20 Hz – 20 kHz range; bands 2–29 are used,
- shorter inputs (padded with -60 dB silence) or longer (truncated to 28).

# Keyword arguments
- `field_type::Symbol = :free`: listening condition, either `:free` or `:diffuse`.

# Returns
A `ZwickerResult` with total loudness in sone, loudness level in phon,
and 240-bin specific loudness `N'(z)` in sone/Bark at 0.1 Bark resolution.

# Throws
- `ArgumentError` if `field_type` is not `:free` or `:diffuse`.
- `ArgumentError` if any of the first 11 band levels exceeds 120 dB (outside
  the Zwicker method's valid range).

# Example
```julia
spl = Float64[60, 62, 65, 68, 70, 72, 74, 75, 73, 71,
              69, 67, 65, 63, 61, 59, 57, 55, 53, 50,
              47, 44, 41, 38, 35, 32, 29, 26]
result = zwicker_loudness(spl)
result.loudness         # total loudness [sone]
result.loudness_level   # loudness level [phon]
result.specific_loudness  # 240-bin N'(z) [sone/Bark]
```
"""
function zwicker_loudness(spl_third_octave::AbstractVector{<:Real}; field_type::Symbol=:free)
    if field_type !== :free && field_type !== :diffuse
        throw(ArgumentError("field_type must be :free or :diffuse, got :$field_type"))
    end

    spl = Float64.(spl_third_octave)
    n = length(spl)

    if n == 28
        spl_28 = spl
    elseif n == 31
        spl_28 = spl[2:29]
    elseif n > 28
        spl_28 = spl[1:28]
    else
        spl_28 = vcat(spl, fill(-60.0, 28 - n))
    end

    if maximum(spl_28[1:11]) > 120.0
        throw(ArgumentError("1/3 octave band levels exceed 120 dB in bands 1-11; Zwicker method not valid."))
    end

    le = compute_excitation_levels(spl_28, field_type)
    nm = compute_core_loudness(le)
    N, ns = compute_spreading(nm)
    LN = sones_to_phons(N)

    return ZwickerResult(N, LN, ns)
end

include("method2.jl")

end # module
