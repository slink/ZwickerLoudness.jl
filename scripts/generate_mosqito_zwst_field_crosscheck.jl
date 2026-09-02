# One-time generator for test/fixtures/zwst_field_fixtures.jl.
# Builds a few synthetic 28-band one-third-octave SPL vectors (no signal
# synthesis, no filter bank -- Method 1's input contract is band levels)
# and hands them to MoSQITo's stationary private stages via
# scripts/crosscheck_zwst_field.py, once per field type (free/diffuse).
# Requires uv and network on first run.
# Usage: julia scripts/generate_mosqito_zwst_field_crosscheck.jl

# (name, spl[28]) -- 25 Hz .. 12.5 kHz band centers, dB re 20 uPa.
# "flat_70db": every band equal -- diffuse-field correction visible in
#   every band at once.
# "tone_band17_60db": single 1 kHz band on a -60 dB floor -- isolates the
#   correction near the field-type curves' crossover.
# "slope_80_to_40": linear tilt -- exercises the low-frequency correction
#   (bands 1-11) at high level together with the high-band correction.
# "hump_40_plus_30": Gaussian hump centred on band 14 (500 Hz) -- a
#   speech-like shape with moderate levels everywhere.
cases = [
    ("flat_70db", fill(70.0, 28)),
    ("tone_band17_60db", [i == 17 ? 60.0 : -60.0 for i in 1:28]),
    ("slope_80_to_40", collect(range(80.0, 40.0; length=28))),
    ("hump_40_plus_30", [40.0 + 30.0 * exp(-((i - 14) / 5.0)^2) for i in 1:28]),
]

jsonvec(v) = string("[", join(string.(v), ","), "]")
tmp = tempname() * ".json"
open(tmp, "w") do io
    entries = [string("""{"name": "$name", "spl": """, jsonvec(spl), "}") for (name, spl) in cases]
    print(io, "[", join(entries, ","), "]")
end

out = joinpath(@__DIR__, "..", "test", "fixtures", "zwst_field_fixtures.jl")
run(`uv run --with mosqito==1.2.1 --with numpy --with matplotlib python $(joinpath(@__DIR__, "crosscheck_zwst_field.py")) $tmp $out`)
