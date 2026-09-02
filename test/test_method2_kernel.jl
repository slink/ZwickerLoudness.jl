using Test
using ZwickerLoudness
using Statistics
using SHA

# Guarded include -- Julia 1.10's const-reinclude landmine means re-including
# this file within the same session throws; `isdefined` makes `runtests.jl`
# (which may already have loaded it, or be re-run interactively) safe either way.
if !isdefined(Main, :ZWTV_KERNEL_FIXTURE_CASES)
    include(joinpath(@__DIR__, "fixtures", "zwtv_kernel_fixtures.jl"))
end

# SHA-256 of a matrix's bytes in ROW-MAJOR order (matches numpy's default
# `.tobytes()` on the (28, nblocks) array used to generate `band_levels_sha256`
# in the fixture file -- Julia matrices are column-major internally, so this
# must iterate row-then-column explicitly, not rely on `vec`/`reshape`).
function _sha256_rowmajor(m::AbstractMatrix{Float64})
    io = IOBuffer()
    nrows, ncols = size(m)
    for i in 1:nrows, j in 1:ncols
        write(io, m[i, j])
    end
    return bytes2hex(sha256(take!(io)))
end

@testset "ZwickerLoudness.jl - ISO 532-1:2017 clause 6 (time-varying kernel)" begin

    # ============================================================ #
    #  Fixture integrity (band_levels_sha256 cross-check)
    # ============================================================ #
    @testset "fixture band_levels integrity (sha256)" begin
        for case in ZWTV_KERNEL_FIXTURE_CASES
            @test _sha256_rowmajor(case.band_levels) == case.band_levels_sha256
        end
    end

    # ============================================================ #
    #  Per-stage: nonlinear decay vs MoSQITo `_nl_loudness` intermediate
    #  (tone_burst_1k_70db_30ms case; pre-decay `nm` regenerated from
    #  already-tested Method 1 stages, per zwtv-pins.md §5 -- not a
    #  separately vendored fixture).
    # ============================================================ #
    @testset "_nonlinear_decay vs fixture intermediate" begin
        case = only(filter(c -> c.name == ZWTV_STAGE_INTERMEDIATES.name, ZWTV_KERNEL_FIXTURE_CASES))
        bl = case.band_levels
        nblocks = size(bl, 2)

        nm_pre_decay = Matrix{Float64}(undef, 21, nblocks)
        for j in 1:nblocks
            le = ZwickerLoudness.compute_excitation_levels(bl[:, j], case.field_type)
            nm_pre_decay[:, j] = ZwickerLoudness.compute_core_loudness(le)
        end

        nm_post_decay = ZwickerLoudness._nonlinear_decay(nm_pre_decay)
        ref = ZWTV_STAGE_INTERMEDIATES.nm_post_decay

        @test size(nm_post_decay) == size(ref)
        # Measured deviation ~1e-15/1e-16 (transcription class); rtol=1e-9
        # per the plan's "start rtol 1e-9, measure" instruction -- NOT
        # loosened, since the measured agreement is far tighter.
        @test isapprox(nm_post_decay, ref; rtol=1e-9, atol=1e-12)
    end

    # ============================================================ #
    #  Per-stage: temporal weighting vs MoSQITo `_temporal_weighting`
    #  intermediate (post-weighting, pre-decimation N(t)).
    # ============================================================ #
    @testset "_temporal_weighting vs fixture intermediate" begin
        case = only(filter(c -> c.name == ZWTV_STAGE_INTERMEDIATES.name, ZWTV_KERNEL_FIXTURE_CASES))
        bl = case.band_levels
        nblocks = size(bl, 2)

        nm_pre_decay = Matrix{Float64}(undef, 21, nblocks)
        for j in 1:nblocks
            le = ZwickerLoudness.compute_excitation_levels(bl[:, j], case.field_type)
            nm_pre_decay[:, j] = ZwickerLoudness.compute_core_loudness(le)
        end
        nm_post_decay = ZwickerLoudness._nonlinear_decay(nm_pre_decay)

        total_pre_weight = Vector{Float64}(undef, nblocks)
        for j in 1:nblocks
            N, _ns = ZwickerLoudness.compute_spreading(nm_post_decay[:, j])
            total_pre_weight[j] = N
        end

        N_pre_decimation = ZwickerLoudness._temporal_weighting(total_pre_weight)
        ref = ZWTV_STAGE_INTERMEDIATES.N_pre_decimation

        @test length(N_pre_decimation) == length(ref)
        @test isapprox(N_pre_decimation, ref; rtol=1e-9, atol=1e-12)
    end

    # ============================================================ #
    #  End-to-end kernel vs fixture N(t), for all synthetic cases.
    # ============================================================ #
    @testset "zwicker_loudness_time_varying: end-to-end N(t)" begin
        for case in ZWTV_KERNEL_FIXTURE_CASES
            result = zwicker_loudness_time_varying(case.band_levels; field_type=case.field_type)

            @test length(result.loudness_over_time) == length(case.N_t)
            @test isapprox(result.loudness_over_time, case.N_t; rtol=1e-9, atol=1e-12)

            @test size(result.specific_loudness, 1) == 240
            @test size(result.specific_loudness, 2) == length(result.loudness_over_time)
            @test length(result.time_axis) == length(result.loudness_over_time)
        end
    end

    # ============================================================ #
    #  Kernel time_axis: the vendored fixture's time_axis is a FRONT-END
    #  artifact (numpy `linspace(0, len(sig)/fs, n_time)` over the ORIGINAL
    #  signal's sample count/duration, per `_third_octave_levels.py:259`)
    #  which the kernel cannot reproduce bit-exactly since it never sees the
    #  original signal or fs -- only the 28 x nblocks band-level matrix at
    #  the fixed, standard-mandated 2000 Hz block rate. The kernel instead
    #  synthesizes its own IDEAL, exact-multiples-of-2ms axis from the block
    #  count alone. Measured: for `steady_tone_1k_60db`, the fixture's step
    #  is 0.00201005... s vs the kernel's exact 0.002 s (~0.5% high because
    #  the front end's `linspace` divides the signal's total duration over
    #  (n_time-1) intervals rather than stepping at a fixed 1/2000 s rate).
    #  This is a front-end/API-contract difference, not a transcription bug,
    #  so it is checked for SHAPE/monotonicity/nominal rate here, not
    #  bit-exact equality against the fixture.
    # ============================================================ #
    @testset "kernel time_axis: shape and nominal 2 ms rate" begin
        for case in ZWTV_KERNEL_FIXTURE_CASES
            result = zwicker_loudness_time_varying(case.band_levels; field_type=case.field_type)
            @test result.time_axis[1] == 0.0
            @test issorted(result.time_axis)
            @test length(result.time_axis) == cld(size(case.band_levels, 2), 4)
            if length(result.time_axis) > 1
                steps = diff(result.time_axis)
                @test all(s -> isapprox(s, 0.002; rtol=1e-9), steps)
            end
        end
    end

    # ============================================================ #
    #  Percentiles: N5/N10 (ISO/DIN definition, type-7 quantile) vs fixture.
    #  NOT bit-exact -- fixture's N5_iso/N10_iso were computed by numpy
    #  `percentile(..., method="linear")` on the numpy-side N(t); Julia's
    #  `Statistics.quantile` on the same vendored N_t differs by up to 1 ULP
    #  (documented in the fixture header). rtol ~4*eps() per the plan.
    # ============================================================ #
    @testset "N5/N10 vs fixture (ISO/DIN type-7, ~1 ULP)" begin
        for case in ZWTV_KERNEL_FIXTURE_CASES
            result = zwicker_loudness_time_varying(case.band_levels; field_type=case.field_type)
            @test isapprox(result.N5, case.N5_iso; rtol=4 * eps())
            @test isapprox(result.N10, case.N10_iso; rtol=4 * eps())

            # Internally self-consistent: the kernel's own N5/N10 are exactly
            # `quantile` of its own returned `loudness_over_time` (not
            # re-derived some other way).
            @test result.N5 == quantile(result.loudness_over_time, 0.95)
            @test result.N10 == quantile(result.loudness_over_time, 0.90)
        end
    end

    # ============================================================ #
    #  Percentile-definition pin: bit-exact both ways (numpy "linear" ==
    #  Julia default `quantile`, Hyndman-Fan type 7). The type-5/Hazen
    #  companion fields (`n5_type5`/`n10_type5`) are recorded in the fixture
    #  for the future MATLAB-`prctile` gap attribution only -- not asserted
    #  here (no Julia-side implementation of that scheme in this package).
    # ============================================================ #
    @testset "percentile pin: type-7 bit-exact vs numpy" begin
        @test quantile(ZWTV_PERCENTILE_PIN_VECTOR, 0.95) === ZWTV_PERCENTILE_PIN.n5_type7
        @test quantile(ZWTV_PERCENTILE_PIN_VECTOR, 0.90) === ZWTV_PERCENTILE_PIN.n10_type7
    end

    # ============================================================ #
    #  Specific loudness is never temporally weighted (only decimated) --
    #  verify that decimating a freshly-computed pre-weighting specific
    #  loudness matches the kernel's returned `specific_loudness` exactly
    #  for one case.
    # ============================================================ #
    @testset "specific loudness: decimated but never temporally weighted" begin
        case = only(filter(c -> c.name == "steady_tone_1k_60db", ZWTV_KERNEL_FIXTURE_CASES))
        bl = case.band_levels
        nblocks = size(bl, 2)

        nm_pre_decay = Matrix{Float64}(undef, 21, nblocks)
        for j in 1:nblocks
            le = ZwickerLoudness.compute_excitation_levels(bl[:, j], case.field_type)
            nm_pre_decay[:, j] = ZwickerLoudness.compute_core_loudness(le)
        end
        nm_post_decay = ZwickerLoudness._nonlinear_decay(nm_pre_decay)

        specific_pre_weight = Matrix{Float64}(undef, 240, nblocks)
        for j in 1:nblocks
            _N, ns = ZwickerLoudness.compute_spreading(nm_post_decay[:, j])
            specific_pre_weight[:, j] = ns
        end
        expected_specific = specific_pre_weight[:, 1:4:nblocks]

        result = zwicker_loudness_time_varying(bl; field_type=case.field_type)
        @test result.specific_loudness == expected_specific
    end

    # ============================================================ #
    #  Input validation
    # ============================================================ #
    # ============================================================ #
    #  Diffuse-field branch: "steady_tone_1k_60db_diffuse" is the SAME
    #  signal as "steady_tone_1k_60db" (identical band_levels -- the front
    #  end is field-independent), so the two fixtures differ ONLY in the
    #  kernel's field_type path, and the end-to-end loop above already
    #  pinned N(t) for both against MoSQITo.
    # ============================================================ #
    @testset "diffuse-field case: same band levels, different N(t)" begin
        free = only(filter(c -> c.name == "steady_tone_1k_60db", ZWTV_KERNEL_FIXTURE_CASES))
        diff = only(filter(c -> c.name == "steady_tone_1k_60db_diffuse", ZWTV_KERNEL_FIXTURE_CASES))
        @test diff.field_type == :diffuse
        @test diff.band_levels_sha256 == free.band_levels_sha256
        @test diff.band_levels == free.band_levels
        @test diff.N5_iso != free.N5_iso
        r_free = zwicker_loudness_time_varying(free.band_levels; field_type=:free)
        r_diff = zwicker_loudness_time_varying(free.band_levels; field_type=:diffuse)
        @test r_free.loudness_over_time != r_diff.loudness_over_time
        @test isapprox(r_diff.loudness_over_time, diff.N_t; rtol=1e-9, atol=1e-12)
    end

    @testset "invalid field_type throws" begin
        @test_throws ArgumentError zwicker_loudness_time_varying(fill(60.0, 28, 10); field_type=:invalid)
    end

    @testset "wrong band count throws" begin
        @test_throws ArgumentError zwicker_loudness_time_varying(fill(60.0, 27, 10))
        @test_throws ArgumentError zwicker_loudness_time_varying(fill(60.0, 31, 10))
    end

    @testset "zero-column input throws" begin
        @test_throws ArgumentError zwicker_loudness_time_varying(zeros(28, 0))
    end

    @testset "excessive low-frequency level throws" begin
        bl = fill(-60.0, 28, 5)
        bl[3, 2] = 125.0
        @test_throws ArgumentError zwicker_loudness_time_varying(bl)
    end

    # ============================================================ #
    #  Sanity
    # ============================================================ #
    @testset "silence: zero loudness throughout" begin
        result = zwicker_loudness_time_varying(fill(-60.0, 28, 20))
        @test all(result.loudness_over_time .== 0.0)
        @test result.N5 == 0.0
        @test result.N10 == 0.0
    end

    @testset "monotonicity: louder band levels give higher N5" begin
        r_quiet = zwicker_loudness_time_varying(fill(40.0, 28, 40))
        r_loud = zwicker_loudness_time_varying(fill(80.0, 28, 40))
        @test r_loud.N5 > r_quiet.N5
    end
end
