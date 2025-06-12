module GDIFF

using LinearAlgebra, SpecialFunctions,Interpolations
using SparseArrays

include("diffusion_grt.jl")
export diffusion_grt

include("linear_interpolation.jl")
export interp1_linear,linear_interpolation_1D

end
