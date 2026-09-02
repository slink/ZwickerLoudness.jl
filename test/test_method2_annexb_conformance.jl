using Test
using ZwickerLoudness

# =========================================================================== #
#  Kernel-level ISO 532-1:2017 Annex B.4 conformance (Task 3).
#
#  Mechanism: test/fixtures/zwtv_kernel_fixtures.jl's `ZWTV_ANNEXB_DERIVED`
#  vendors only MoSQITo's *derived numbers* (coarsened N(t) + exact N5/N10)
#  for Annex B.4 "Test signal 6" -- deliberately NOT the `band_levels` matrix
#  that produced them, since that .wav is ISO-copyrighted material MoSQITo
#  merely redistributes (see .superpowers/sdd/zwtv-pins.md §2 and
#  test/fixtures/NOTICE.md). So a kernel-level conformance check ("run OUR
#  kernel on the SAME band levels MoSQITo used, compare against MoSQITo's own
#  numbers") cannot use a vendored fixture for its input -- it must
#  regenerate `band_levels` locally, from a local copy of the Annex B .wav,
#  every time it runs.
#
#  This testset does exactly that: if a local MoSQITo checkout at the pinned
#  commit (with the Annex B .wav) is present at /tmp/mosqito-pinned, it
#  shells out to `scripts/dump_annexb_band_levels.py` (via `uv run --with
#  mosqito==1.2.1`: the PyPI release whose package code is identical to the
#  pinned commit, so the oracle stays reproducible whatever PyPI publishes
#  later) to regenerate `band_levels` from the .wav into a throwaway file
#  (never committed, deleted immediately after use), feeds that matrix
#  through THIS package's `zwicker_loudness_time_varying`, and compares the
#  result against the vendored `ZWTV_ANNEXB_DERIVED` numbers.
#
#  CI has no access to /tmp/mosqito-pinned (it is a local, gitignored,
#  throwaway clone -- see zwtv-pins.md §0), so this testset self-skips via
#  `@test_skip` with an explanatory message whenever that material is
#  absent, keeping CI green while local development (with the reference
#  material cloned) gets the full end-to-end gate. This is the deliberate
#  choice recorded in the Task 3 report over a Pkg.test-internal-only
#  design: regenerating band_levels via a Python subprocess is exactly what
#  Task 1's rig already does at fixture-generation time, so reusing that
#  mechanism here (rather than reimplementing a front end inside this
#  zero-external-dependency kernel package) is both faithful and minimal.
# =========================================================================== #

const _ANNEXB_MOSQITO_CHECKOUT = joinpath("/tmp", "mosqito-pinned")
const _ANNEXB_WAV = joinpath(
    _ANNEXB_MOSQITO_CHECKOUT, "validations", "sq_metrics", "loudness_zwtv",
    "input", "ISO_532-1", "Annex B.4", "Test signal 6 (tone 250 Hz 30 dB - 80 dB).wav",
)
const _ANNEXB_DUMP_SCRIPT = joinpath(@__DIR__, "..", "scripts", "dump_annexb_band_levels.py")

_annexb_skip_reason() =
    if Sys.which("uv") === nothing
        "Annex B kernel conformance skipped: `uv` not found on PATH (required to run " *
        "the local MoSQITo reference driver). Install uv to enable this gate locally; " *
        "CI intentionally lacks it and stays green via this skip."
    elseif !isfile(_ANNEXB_WAV)
        "Annex B kernel conformance skipped: local MoSQITo reference material not " *
        "found at \"$_ANNEXB_WAV\". Clone MoSQITo @ d990c33f94f1 to " *
        "$_ANNEXB_MOSQITO_CHECKOUT to enable (see .superpowers/sdd/zwtv-pins.md §2 for " *
        "why the .wav itself is never vendored in this repository). CI intentionally " *
        "lacks this material and stays green via this skip."
    else
        nothing
    end

@testset "ZwickerLoudness.jl - ISO 532-1:2017 clause 6 (Annex B.4 kernel conformance)" begin
    reason = _annexb_skip_reason()
    if reason !== nothing
        @test_skip reason
    else
        if !isdefined(Main, :ZWTV_ANNEXB_DERIVED)
            include(joinpath(@__DIR__, "fixtures", "zwtv_kernel_fixtures.jl"))
        end

        tmp_jl = tempname() * ".jl"
        try
            run(`uv run --with mosqito==1.2.1 --with numpy --with matplotlib python $_ANNEXB_DUMP_SCRIPT $_ANNEXB_WAV $tmp_jl`)
            include(tmp_jl)
            band_levels = ANNEXB_LOCAL_BAND_LEVELS  # 28 x nblocks @ 2000 Hz, regenerated locally

            @test size(band_levels, 1) == 28
            # 4x the kernel's own final decimation factor: the regenerated
            # band_levels are at the 0.5 ms input rate, N(t) is at 2 ms.
            @test size(band_levels, 2) == ZWTV_ANNEXB_DERIVED.n_full * 4

            result = zwicker_loudness_time_varying(band_levels; field_type=:free)
            @test length(result.loudness_over_time) == ZWTV_ANNEXB_DERIVED.n_full

            coarse = ZWTV_ANNEXB_DERIVED.coarse_factor
            N_t_coarse_local = result.loudness_over_time[1:coarse:end]
            # Measured (Task 3): max abs diff ~1.1e-14, max rel diff ~1.5e-14
            # against MoSQITo's own end-to-end N(t) on the identical .wav --
            # the same transcription-level (not ISO-tolerance-class) agreement
            # measured for the synthetic fixtures in Task 2, since this is
            # literally the same kernel algorithm run on the same band-level
            # input MoSQITo itself computed. rtol=1e-9 kept consistent with
            # the rest of the kernel test suite (not loosened to match the
            # measured precision).
            @test isapprox(N_t_coarse_local, ZWTV_ANNEXB_DERIVED.N_t_coarse; rtol=1e-9, atol=1e-9)
            @test isapprox(result.N5, ZWTV_ANNEXB_DERIVED.N5_iso; rtol=1e-9, atol=1e-9)
            @test isapprox(result.N10, ZWTV_ANNEXB_DERIVED.N10_iso; rtol=1e-9, atol=1e-9)
        finally
            isfile(tmp_jl) && rm(tmp_jl)
        end
    end
end
