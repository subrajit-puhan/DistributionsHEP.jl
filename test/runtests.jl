using Test

@testset "DistributionsHEP tests" verbose=true begin 
    include("test-chebyshev.jl")
    include("test-argusBG.jl")
    include("test-crystalball.jl")
    include("test-double-sided-crystal-ball.jl")
<<<<<<< HEAD
    include("test-RelativisticBreitWigner.jl")
=======
    include("test-voigt.jl")
>>>>>>> da142b2 (First commit)
end
