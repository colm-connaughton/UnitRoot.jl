<a id='UnitRoot.simul_arma' href='#UnitRoot.simul_arma'>#</a>
**`UnitRoot.simul_arma`** &mdash; *Function*.



Simulates a stochastic process {yₜ} following an ARMA(1,1) with parameters φ and Θ for t = 1, ..., T periods.

$$
  y_t = \mu_0 + \mu_1 + \phi*y_{t-1} + \epsilon_t + \theta*\epsilon_{t-1}, \quad \epsilon_t \sim \mathbf{N}(0, \sigma^2)
$$

**Required Input Parameters**

none

**Optional Input Parameters**

`mu::Float64`: AR(1) process constant term. (default=0.0)

`lin_trend::Float64`: AR(1) process linear time trend. (default=0.0)

`ar::Float64`: AR(1) autorregressive parameter. (default=0.5)

`ma::Float64`: MA(1) autorregressive parameter. (default=0.5)

`T::Int64`: Length of simulated time series. (default=100)

`y_0::Float64`: Initial value of the series. (default=0.0)

`stderr::Float64` : Normally distributed innovation standard error. (default=1.0)

**Output**

`y::Vector{Float64}`: Simulated ARMA(1,1) time series of length T.


<a target='_blank' href='https://github.com/p-chaim/UnitRoot.jl/tree/61100876c17fa088d9cdb045021baf18a0be3314/src/adf.jl#L136-L167' class='documenter-source'>source</a><br>

