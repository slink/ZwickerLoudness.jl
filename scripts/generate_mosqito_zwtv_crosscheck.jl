# One-time generator for test/fixtures/zwtv_kernel_fixtures.jl.
# Synthesizes deterministic time-varying test signals, hands them to
# MoSQITo's loudness_zwtv (via uv, private per-stage functions + the
# public API for a consistency check), vendors the kernel-level fixtures:
# per-case (band_levels, N_t, N5/N10) plus full stage intermediates for
# one case, plus the percentile-definition pin. The kernel/front-end
# boundary is MoSQITo's private `_third_octave_levels` output (see
# src/method2.jl's module docstring); the licensing decisions this rig
# encodes are recorded in test/fixtures/NOTICE.md ("Licensing posture").
# Requires uv and network on first run.
# Usage: julia scripts/generate_mosqito_zwtv_crosscheck.jl
using Statistics: std

# Calibrated pure tone, `duration` seconds at `fc` Hz, `spl` dB re 20 uPa.
# Amplitude calibration transcribed from MoSQITo utils/am_sine_generator.py
# @ d990c33f94f1 (std-normalization, not RMS -- ddof=0), same convention
# PsychoacousticMetrics.jl's am_sine uses (test/support/am_generator.jl).
function tone(fs::Real, fc::Real, duration::Real, spl::Real)
    n = round(Int, duration * fs)
    t = (0:n-1) ./ fs
    y = sin.(2π * fc .* t)
    A = 2e-5 * 10.0^(spl / 20)
    return y .* (A / std(y; corrected=false))
end

# Tone burst: `burst_duration` s of a `spl`-calibrated tone (level defined
# over the burst itself, not the padded silence) starting at `burst_start`
# s within a `total_duration` s silent buffer.
function tone_burst(fs::Real, fc::Real, total_duration::Real, burst_duration::Real,
                     burst_start::Real, spl::Real)
    n_total = round(Int, total_duration * fs)
    n_burst = round(Int, burst_duration * fs)
    i0 = round(Int, burst_start * fs)
    burst = tone(fs, fc, burst_duration, spl)
    sig = zeros(n_total)
    sig[i0+1:i0+n_burst] .= burst
    return sig
end

# Amplitude-modulated tone: (1 + sin(2*pi*fm*t)) * sin(2*pi*fc*t), calibrated
# to `spl` dB by std-normalization (same convention as `tone`/am_sine).
function am_tone(fs::Real, fc::Real, fm::Real, duration::Real, spl::Real)
    n = round(Int, duration * fs)
    t = (0:n-1) ./ fs
    xmod = sin.(2π * fm .* t)
    xc = sin.(2π * fc .* t)
    y = (1 .+ xmod) .* xc
    A = 2e-5 * 10.0^(spl / 20)
    return y .* (A / std(y; corrected=false))
end

# (name, fs, field_type, signal)
# "steady_tone_1k_60db": baseline steady-state case (loudness should settle
# near the Method 1 single-band value once transients decay).
# "tone_burst_1k_70db_30ms": exercises the nonlinear decay (attack/release)
# and temporal weighting stages -- the reason for the full stage dump.
# "am_tone_1k_10hz_65db": non-constant envelope for temporal-weighting
# behavior beyond a single transient.
cases = [
    ("steady_tone_1k_60db", 48000.0, "free",
        tone(48000, 1000.0, 0.10, 60.0)),
    ("tone_burst_1k_70db_30ms", 48000.0, "free",
        tone_burst(48000, 1000.0, 0.15, 0.03, 0.06, 70.0)),
    ("am_tone_1k_10hz_65db", 48000.0, "free",
        am_tone(48000, 1000.0, 10.0, 0.20, 65.0)),
]
const FULL_DUMP_CASE = "tone_burst_1k_70db_30ms"

jsonvec(v) = string("[", join(string.(v), ","), "]")
tmp = tempname() * ".json"
open(tmp, "w") do io
    entries = String[]
    for (name, fs, field_type, sig) in cases
        push!(entries,
            """{"name": "$name", "fs": $fs, "field_type": "$field_type", "signal": $(jsonvec(sig))}""")
    end
    print(io, "[", join(entries, ","), "]")
end

out = joinpath(@__DIR__, "..", "test", "fixtures", "zwtv_kernel_fixtures.jl")

# Optional: one ISO 532-1 Annex B time-varying signal, end-to-end, DERIVED
# NUMBERS ONLY (see crosscheck_zwtv.py's Annex B block and "Licensing
# posture" in test/fixtures/NOTICE.md for the reasoning). Points at a
# *local* MoSQITo checkout at the pin -- never committed, and the path is
# only used transiently at generation time. Set to `nothing` to skip.
annexb_wav = joinpath(
    "/tmp", "mosqito-pinned", "validations", "sq_metrics", "loudness_zwtv",
    "input", "ISO_532-1", "Annex B.4", "Test signal 6 (tone 250 Hz 30 dB - 80 dB).wav",
)
annexb_name = "annexb_test_signal_6_tone_250hz"
annexb_provenance = "MoSQITo validations/sq_metrics/loudness_zwtv/input/ISO_532-1/Annex B.4/" *
    "Test signal 6 (tone 250 Hz 30 dB - 80 dB).wav @ d990c33f94f1 (Apache-2.0); " *
    "underlying signal is ISO 532-1:2017 Annex B.4 Test signal 6, ISO-copyrighted."

if isfile(annexb_wav)
    run(`uv run --with mosqito==1.2.1 --with numpy --with matplotlib python $(joinpath(@__DIR__, "crosscheck_zwtv.py")) $tmp $out $FULL_DUMP_CASE $annexb_wav $annexb_name $annexb_provenance`)
else
    @warn "Annex B wav not found locally at $annexb_wav; skipping Annex B derived-numbers fixture (clone MoSQITo @ d990c33f94f1 to /tmp/mosqito-pinned to include it)."
    run(`uv run --with mosqito==1.2.1 --with numpy --with matplotlib python $(joinpath(@__DIR__, "crosscheck_zwtv.py")) $tmp $out $FULL_DUMP_CASE`)
end
