#__precompile__()
module UnitRoot

using KernelDensity
using Distributions

#ar_lagmatrix
include("adf.jl")
include("kpss.jl")
include("dfgls.jl")

export df_test,
adf_test,
kpss_test,
dfgls_test

end # module
