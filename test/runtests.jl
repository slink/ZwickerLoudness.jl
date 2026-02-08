using Test
using ZwickerLoudness

@testset "ZwickerLoudness.jl - ISO 532-1:2017" begin

    # ============================================================ #
    #  Annex B Conformance: Signal 1 (broadband spectrum)
    # ============================================================ #
    @testset "Annex B Signal 1: broadband spectrum" begin
        spl = Float64[-60, -60, 78, 79, 89, 72, 80, 89, 75, 87,
                       85, 79, 86, 80, 71, 70, 72, 71, 72, 74,
                       69, 65, 67, 77, 68, 58, 45, 30]
        result = zwicker_loudness(spl)
        @test result.loudness ≈ 83.296 rtol=0.05
        @test result.loudness_level ≈ 103.802 rtol=0.05
        @test length(result.specific_loudness) == 240
    end

    # ============================================================ #
    #  Annex B Conformance: 1 kHz at 60 dB (~Signal 3)
    # ============================================================ #
    @testset "1 kHz at 60 dB -> ~4.019 sone" begin
        spl = fill(-60.0, 28)
        spl[17] = 60.0
        result = zwicker_loudness(spl)
        @test result.loudness ≈ 4.019 rtol=0.10
    end

    # ============================================================ #
    #  Sone-to-Phon: Two Branches
    # ============================================================ #
    @testset "sone-to-phon: N >= 1" begin
        @test ZwickerLoudness.sones_to_phons(1.0) ≈ 40.0
        @test ZwickerLoudness.sones_to_phons(2.0) ≈ 50.0
        @test ZwickerLoudness.sones_to_phons(4.0) ≈ 60.0
    end

    @testset "sone-to-phon: N < 1 (power-law branch)" begin
        @test ZwickerLoudness.sones_to_phons(0.5) ≈ 40.0 * 0.5005^0.35 atol=0.5
    end

    @testset "sone-to-phon: minimum 3 phon" begin
        @test ZwickerLoudness.sones_to_phons(0.001) >= 3.0
    end

    @testset "sone-to-phon: zero" begin
        @test ZwickerLoudness.sones_to_phons(0.0) ≈ 0.0
    end

    # ============================================================ #
    #  Low-Frequency Correction
    # ============================================================ #
    @testset "low-frequency correction" begin
        spl_11 = fill(-60.0, 11)
        spl_11[6] = 45.0
        corrected, _ = ZwickerLoudness.correct_low_frequencies(spl_11)
        @test corrected[6] ≈ 45.0 atol=0.1

        spl_11_2 = fill(-60.0, 11)
        spl_11_2[1] = 50.0
        corrected2, _ = ZwickerLoudness.correct_low_frequencies(spl_11_2)
        @test corrected2[1] < 50.0
    end

    # ============================================================ #
    #  Free vs Diffuse Field
    # ============================================================ #
    @testset "free vs diffuse field differ" begin
        spl = fill(70.0, 28)
        r_free = zwicker_loudness(spl; field_type=:free)
        r_diffuse = zwicker_loudness(spl; field_type=:diffuse)
        @test r_free.loudness != r_diffuse.loudness
    end

    @testset "default is free field" begin
        spl = fill(70.0, 28)
        @test zwicker_loudness(spl).loudness ≈ zwicker_loudness(spl; field_type=:free).loudness
    end

    # ============================================================ #
    #  Specific Loudness Shape
    # ============================================================ #
    @testset "specific loudness: 240 bins, non-negative" begin
        spl = fill(70.0, 28)
        result = zwicker_loudness(spl)
        @test length(result.specific_loudness) == 240
        @test all(result.specific_loudness .>= 0.0)
    end

    # ============================================================ #
    #  Input Length Handling
    # ============================================================ #
    @testset "28-band input" begin
        @test zwicker_loudness(fill(60.0, 28)).loudness > 0.0
    end

    @testset "31-band input matches 28-band" begin
        r_31 = zwicker_loudness(fill(60.0, 31))
        r_28 = zwicker_loudness(fill(60.0, 28))
        @test r_31.loudness ≈ r_28.loudness atol=1e-6
    end

    @testset "short input padded" begin
        @test zwicker_loudness(fill(60.0, 10)).loudness > 0.0
    end

    # ============================================================ #
    #  Sanity
    # ============================================================ #
    @testset "silence" begin
        result = zwicker_loudness(fill(-60.0, 28))
        @test result.loudness ≈ 0.0 atol=1e-10
    end

    @testset "monotonicity" begin
        @test zwicker_loudness(fill(80.0, 28)).loudness > zwicker_loudness(fill(40.0, 28)).loudness
    end

    @testset "invalid field_type throws" begin
        @test_throws ArgumentError zwicker_loudness(fill(60.0, 28); field_type=:invalid)
    end
end
