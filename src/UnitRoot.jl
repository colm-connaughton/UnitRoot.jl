__precompile__()
module UnitRoot

using KernelDensity
using Distributions
using LinearAlgebra

#ar_lagmatrix
include("adf.jl")
include("kpss.jl")
include("dfgls.jl")

export df_test,
adf_test,
kpss_test,
dfgls_test,
simul_arma

end # module
