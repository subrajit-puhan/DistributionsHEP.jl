module DistributionsHEP

using Random
using Distributions
import Distributions: @check_args

include("chebychev.jl")
include("argusBG.jl")

end
