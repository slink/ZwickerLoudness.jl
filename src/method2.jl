# =========================================================================== #
#  ISO 532-1:2017 clause 6, "Method for time-varying sounds"
#  (colloquially "Method 2" -- NOT to be confused with ISO 532-2, the
#  unrelated Moore-Glasberg standard).
#
#  Kernel stage: takes a 28-band, 2000 Hz (0.5 ms) third-octave level time
#  series (produced by a front end such as ZwickerLoudnessAudio.jl's
#  `loudness_zwtv`) and computes time-varying loudness N(t) plus N5/N10
#  percentiles.
#
#  Pipeline, per `.superpowers/sdd/zwtv-pins.md` §1 (verified against
#  MoSQITo `loudness_zwtv.py:111-136` @ d990c33f94f1):
#    1. Per-block core loudness `nm` (21 x nblocks): REUSES this package's
#       existing Method 1 stages `compute_excitation_levels` +
#       `compute_core_loudness` (src/ZwickerLoudness.jl) -- MoSQITo's own
#       vectorized `_main_loudness.py` was verified to implement the
#       identical scalar per-block algorithm, just vectorized over time;
#       re-transcribing it would duplicate already-tested code.
#    2. Nonlinear decay (`_nonlinear_decay`, this file) on the PRE-spreading
#       `nm` (21 x nblocks) -- transcribed from
#       `_nonlinear_decay.py` @ d990c33f94f1.
#    3. Spreading (REUSES `compute_spreading`) AFTER decay, producing total
#       loudness N(t) and specific loudness (240 x nblocks) from the
#       decayed `nm`.
#    4. Temporal weighting (`_temporal_weighting` / `_lowpass_intp`, this
#       file) on TOTAL loudness ONLY -- transcribed from
#       `_temporal_weighting.py` / `_lowpass_intp.py` @ d990c33f94f1.
#       Specific loudness is never temporally weighted
#       (`loudness_zwtv.py:122-130`: it is only ever decimated).
#    5. Final /4 decimation from the 0.5 ms block rate to the 2 ms N(t)
#       output axis (`loudness_zwtv.py:128-131`) -- a KERNEL-side property,
#       applied here to total loudness, specific loudness, and time axis
#       alike.
#    6. N5/N10 percentiles: not computed anywhere in MoSQITo -- this is the
#       ISO 532-1 / DIN 45631/A1 definition (N5 = loudness reached or
#       exceeded 5 % of the time = 95th percentile of N(t), N10 = 90th),
#       via `Statistics.quantile` (Hyndman-Fan type 7, matching both
#       Julia's and numpy's `percentile` defaults -- see the percentile pin
#       test for the measured gap against MATLAB `prctile`'s Hazen-style
#       scheme).
# =========================================================================== #

export ZwickerTimeVaryingResult, zwicker_loudness_time_varying

"""
    ZwickerTimeVaryingResult

Result of an ISO 532-1:2017 clause 6 ("method for time-varying sounds")
loudness calculation.

# Fields
- `loudness_over_time::Vector{Float64}`: time-varying total loudness `N(t)`
  in sone, on the 2 ms output axis.
- `N5::Float64`: loudness in sone reached or exceeded 5 % of the time
  (95th percentile of `loudness_over_time`, ISO 532-1 / DIN 45631/A1
  definition).
- `N10::Float64`: loudness in sone reached or exceeded 10 % of the time
  (90th percentile of `loudness_over_time`).
- `specific_loudness::Matrix{Float64}`: specific loudness `N'(z, t)` in
  sone/Bark, 240 rows (0.1 Bark resolution, z = 0.1, 0.2, …, 24.0 Bark) by
  `length(loudness_over_time)` columns (2 ms axis). Never temporally
  weighted (only decimated), per the reference algorithm.
- `time_axis::Vector{Float64}`: time in seconds for each `loudness_over_time`
  / `specific_loudness` column, at the nominal 2 ms output rate.
"""
struct ZwickerTimeVaryingResult
    loudness_over_time::Vector{Float64}
    N5::Float64
    N10::Float64
    specific_loudness::Matrix{Float64}
    time_axis::Vector{Float64}
end

# =========================================================================== #
#  Nonlinear decay (attack/release network on pre-spreading core loudness)
#  Transcribed from `_nonlinear_decay.py` @ d990c33f94f1.
# =========================================================================== #

# Upstream quirk (preserved for fidelity, see `_nl_loudness` in the pinned
# source): the recursion's "previous column" state for the very first
# upsampled sample wraps around to the LAST column via Python's negative-index
# semantics (`uo_mat[:, col - 1]` at `col == 0` reads `uo_mat[:, -1]`), which
# at that point still holds its raw (unprocessed) upsampled-input value. This
# is very likely an unintended artifact of the vectorized rewrite rather than
# a deliberate boundary condition, but the vendored fixtures were generated
# against it, so it is reproduced here rather than "fixed".
function _nonlinear_decay(nm::AbstractMatrix{Float64})
    nbands, nblocks = size(nm)
    nl_iter = 24
    sample_rate = 2000.0
    t_short = 0.005
    t_long = 0.015
    t_var = 0.075

    delta_t = 1.0 / (sample_rate * nl_iter)
    P = (t_var + t_long) / (t_var * t_short)
    Q = 1.0 / (t_short * t_var)
    lambda1 = -P / 2 + sqrt(P * P / 4 - Q)
    lambda2 = -P / 2 - sqrt(P * P / 4 - Q)
    den = t_var * (lambda1 - lambda2)
    e1 = exp(lambda1 * delta_t)
    e2 = exp(lambda2 * delta_t)
    # B[1..6] correspond to MoSQITo's B[0..5] respectively.
    B = (
        (e1 - e2) / den,
        ((t_var * lambda2 + 1) * e1 - (t_var * lambda1 + 1) * e2) / den,
        ((t_var * lambda1 + 1) * e1 - (t_var * lambda2 + 1) * e2) / den,
        (t_var * lambda1 + 1) * (t_var * lambda2 + 1) * (e1 - e2) / den,
        exp(-delta_t / t_long),
        exp(-delta_t / t_var),
    )

    # Virtual upsampling by nl_iter=24 via linear interpolation toward the
    # NEXT block's value (0.0 "next" for the final block), flattened
    # time-major/substep-minor (matches numpy's C-order reshape of the
    # (nbands, nblocks, nl_iter) array MoSQITo builds).
    ncols = nblocks * nl_iter
    ui = Matrix{Float64}(undef, nbands, ncols)
    @inbounds for t in 1:nblocks
        base = (t - 1) * nl_iter
        nextval_row = t < nblocks ? nm[:, t + 1] : nothing
        for b in 1:nbands
            x = nm[b, t]
            nxt = t < nblocks ? nextval_row[b] : 0.0
            d = (nxt - x) / nl_iter
            ui[b, base + 1] = x
            for k in 2:nl_iter
                ui[b, base + k] = ui[b, base + k - 1] + d
            end
        end
    end

    uo = copy(ui)
    u2 = zeros(Float64, nbands, ncols)

    @inbounds for col in 1:ncols
        prev = col == 1 ? ncols : col - 1
        for b in 1:nbands
            uo_prev = uo[b, prev]
            u2_prev = u2[b, prev]
            ui_col = ui[b, col]

            cur_uo = ui_col  # uo initialized as a copy of ui; default unless a branch below fires
            uo2 = uo_prev * B[3] - u2_prev * B[4]
            if uo_prev > u2_prev && uo2 >= ui_col
                cur_uo = uo2
            end
            uo2b = uo_prev * B[5]
            if uo_prev <= u2_prev && uo2b >= ui_col
                cur_uo = uo2b
            end
            uo[b, col] = cur_uo

            cur_u2 = cur_uo  # unconditional u2[:, col] = uo[:, col] in the reference
            u22 = uo_prev * B[1] - u2_prev * B[2]
            if ui_col < uo_prev && uo_prev > u2_prev && u22 <= cur_uo
                cur_u2 = u22
            end
            u2_2 = (u2_prev - ui_col) * B[6] + ui_col
            if ui_col >= uo_prev &&
               !(abs(ui_col - uo_prev) < 1e-5 && cur_uo <= u2_prev)
                cur_u2 = u2_2
            end
            u2[b, col] = cur_u2
        end
    end

    out = Matrix{Float64}(undef, nbands, nblocks)
    @inbounds for t in 1:nblocks
        base = (t - 1) * nl_iter + 1
        for b in 1:nbands
            out[b, t] = uo[b, base]
        end
    end
    return out
end

# =========================================================================== #
#  Temporal weighting of total loudness (two first-order low-pass filters)
#  Transcribed from `_temporal_weighting.py` / `_lowpass_intp.py` @ d990c33f94f1.
# =========================================================================== #

# 1st-order low-pass with linear-interpolation upsampling for precision, per
# `_lowpass_intp.py`. `loudness` is the scalar total-loudness time series
# (never applied to specific loudness -- see module docstring).
function _lowpass_intp(loudness::AbstractVector{Float64}, tau::Float64)
    sample_rate = 2000.0
    lp_iter = 24
    n = length(loudness)
    a1 = exp(-1.0 / (sample_rate * lp_iter * tau))
    b0 = 1.0 - a1

    ncols = n * lp_iter
    ui = Vector{Float64}(undef, ncols)
    @inbounds for t in 1:n
        base = (t - 1) * lp_iter
        x = loudness[t]
        nxt = t < n ? loudness[t + 1] : 0.0
        d = (nxt - x) / lp_iter
        ui[base + 1] = x
        for k in 2:lp_iter
            ui[base + k] = ui[base + k - 1] + d
        end
    end

    # Causal 1st-order IIR: y[n] = b0*x[n] + a1*y[n-1], zero initial state
    # (matches `scipy.signal.lfilter(..., zi=None)`), run continuously across
    # the WHOLE flattened series (state carries across block boundaries) --
    # the reference keeps only the first-substep output of each block.
    out = Vector{Float64}(undef, n)
    yprev = 0.0
    @inbounds for i in 1:ncols
        yi = b0 * ui[i] + a1 * yprev
        yprev = yi
        t, k = fldmod1(i, lp_iter)
        if k == 1
            out[t] = yi
        end
    end
    return out
end

function _temporal_weighting(loudness::AbstractVector{Float64})
    filt1 = _lowpass_intp(loudness, 3.5e-3)
    filt2 = _lowpass_intp(loudness, 70e-3)
    return 0.47 .* filt1 .+ 0.53 .* filt2
end

# =========================================================================== #
#  Final decimation (0.5 ms block rate -> 2 ms output rate)
#  Transcribed from `loudness_zwtv.py:128-131` @ d990c33f94f1
#  (`dec_factor = 4`, a KERNEL-side constant distinct from the front end's
#  `_third_octave_levels.py:36` `dec_factor`).
# =========================================================================== #

function _decimate_loudness(loudness::AbstractVector{Float64}, specific_loudness::AbstractMatrix{Float64})
    nblocks = length(loudness)
    idx = 1:4:nblocks
    time_axis = Float64[(k - 1) / 2000.0 for k in idx]
    return loudness[idx], specific_loudness[:, idx], time_axis
end

# =========================================================================== #
#  Main entry point
# =========================================================================== #

"""
    zwicker_loudness_time_varying(band_levels; field_type=:free) -> ZwickerTimeVaryingResult

Compute time-varying Zwicker loudness per ISO 532-1:2017 clause 6, "method
for time-varying sounds" (colloquially "Method 2" -- a distinct document
from ISO 532-2, the unrelated Moore-Glasberg method).

`band_levels` is a `28 × nblocks` matrix of one-third-octave band SPL values
in dB (re. 20 µPa), one column per 0.5 ms (2000 Hz) time block -- the exact
output contract of a front end such as ZwickerLoudnessAudio.jl's
`loudness_zwtv` third-octave stage. The 2000 Hz block rate is fixed by the
reference algorithm, not a caller-supplied parameter.

# Keyword arguments
- `field_type::Symbol = :free`: listening condition, either `:free` or `:diffuse`.

# Returns
A [`ZwickerTimeVaryingResult`](@ref) with time-varying total loudness `N(t)`
in sone (2 ms axis), `N5`/`N10` percentile loudness in sone, time-varying
240-bin specific loudness `N'(z, t)` in sone/Bark, and the corresponding
time axis in seconds.

# Throws
- `ArgumentError` if `field_type` is not `:free` or `:diffuse`.
- `ArgumentError` if `band_levels` does not have exactly 28 rows, or has
  zero columns.
- `ArgumentError` if any of the first 11 bands exceeds 120 dB in any block
  (outside the Zwicker method's valid range).

# Example
```julia
band_levels = fill(60.0, 28, 400)  # 400 blocks @ 2000 Hz = 0.2 s
result = zwicker_loudness_time_varying(band_levels)
result.loudness_over_time  # N(t) [sone], 2 ms axis
result.N5                  # loudness exceeded 5% of the time [sone]
```
"""
function zwicker_loudness_time_varying(band_levels::AbstractMatrix{<:Real}; field_type::Symbol=:free)
    if field_type !== :free && field_type !== :diffuse
        throw(ArgumentError("field_type must be :free or :diffuse, got :$field_type"))
    end

    nbands, nblocks = size(band_levels)
    if nbands != 28
        throw(ArgumentError("band_levels must have 28 rows (one per one-third-octave band), got $nbands"))
    end
    if nblocks == 0
        throw(ArgumentError("band_levels must have at least one time block"))
    end

    bl = Float64.(band_levels)

    if maximum(view(bl, 1:11, :)) > 120.0
        throw(ArgumentError("1/3 octave band levels exceed 120 dB in bands 1-11; Zwicker method not valid."))
    end

    # Stage 1: per-block core loudness (REUSES Method 1 stages, not re-transcribed).
    nm_pre_decay = Matrix{Float64}(undef, 21, nblocks)
    for j in 1:nblocks
        le = compute_excitation_levels(bl[:, j], field_type)
        nm_pre_decay[:, j] = compute_core_loudness(le)
    end

    # Stage 2: nonlinear decay, PRE-spreading.
    nm_post_decay = _nonlinear_decay(nm_pre_decay)

    # Stage 3: spreading AFTER decay (REUSES Method 1's compute_spreading).
    total_pre_weight = Vector{Float64}(undef, nblocks)
    specific_pre_weight = Matrix{Float64}(undef, 240, nblocks)
    for j in 1:nblocks
        N, ns = compute_spreading(nm_post_decay[:, j])
        total_pre_weight[j] = N
        specific_pre_weight[:, j] = ns
    end

    # Stage 4: temporal weighting on TOTAL loudness only.
    total_weighted = _temporal_weighting(total_pre_weight)

    # Stage 5: final /4 decimation to the 2 ms output axis.
    loudness_over_time, specific_loudness, time_axis = _decimate_loudness(total_weighted, specific_pre_weight)

    # Stage 6: N5/N10 percentiles (ISO/DIN definition, type-7 quantile).
    N5 = quantile(loudness_over_time, 0.95)
    N10 = quantile(loudness_over_time, 0.90)

    return ZwickerTimeVaryingResult(loudness_over_time, N5, N10, specific_loudness, time_axis)
end
