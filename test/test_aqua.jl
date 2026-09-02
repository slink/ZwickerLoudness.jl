using Aqua
using ZwickerLoudness

# Standard Aqua.jl quality battery: undocumented exports, ambiguities, stale
# deps, unbound type parameters, undefined exports, piracy, etc. See
# https://github.com/JuliaTesting/Aqua.jl for what each check covers.
@testset "Aqua.jl quality checks" begin
    Aqua.test_all(ZwickerLoudness)
end
