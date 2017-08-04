@doc """
Calculates classic Dickey-Fuller (1979) statistic.

```math
  \\Delta y_t = \\alpha*y_{t-1} + \\sum^{pmax}_{j=1} + \\epsilon_t, \\quad \\epsilon_t\\sim\\mathcal{N}(0,\\sigma^2)
```

Required Input Parameters
----------------------------------
`y::Vector{Float64}`: Observed time series in vector form.

Optional Input Parameters
----------------------------------
`trend::String`: Deterministic trend specification. It determines the asymptotic distribution of the t-statistic. Possible values are `"constant"`, which considers a constant deterministic term; `"linear"`, which also considers a linear time trend; and `"none"`, which considers no deterministic detrending.

Output
----------------------------------
`df_stat::Float64`: Dickey-Fuller (1979) test statistic.

`df_pval::Float64`: p-value obtained from 10000 simulations.
"""->
function df_test(y::Vector{Float64}; trend::String="constant", pval::Bool=true)
  #constante = false; trend = false;
  #Test regression: Δy_t = μ + αy_t-1 + ɛ_t, test if α = 0
  #                 Δy_t = μ + βt + αy_t-1 + ɛ_t
  if trend=="none"
    x = ar_lagmatrix(y, 1; trend="none")
    (pval==true) && (println("no deteriministic detrending."))
  elseif trend=="constant"
    (pval==true) && (println("constant deteriministic detrending."))
    x = ar_lagmatrix(y, 1; trend="constant")
  elseif trend=="linear"
    (pval==true) && (println("constant plus linear deteriministic detrending."))
    x = ar_lagmatrix(y, 1; trend="linear")
  else
    error(" trend must be either none, constant, or linear")
  end
  T = size(x, 1)
  Δy = diff(y)
  β, t_stat, bic = df_ols(Δy, x); # Test regression
  df_stat = t_stat[end]
  (pval==false) && return df_stat
  df_pval = mean(df_stat.>simulate_dfdist(T, trend=trend))
  return df_stat, df_pval
end
#
#
function simulate_dfdist(T; trend::String="constant")
 nsim = 10000
 t_dist = Vector{Float64}(nsim)
 for sim ∈ 1:nsim
   t_h0 = simul_arma(mu = 0.0, lin_trend = 0.0, ar = 1.0, ma = 0.0, T = T, y_0 = 0.0, stderr = 1.0) # simulate distribution of t under I(1) null hypothesis
   t_dist[sim] = df_test(t_h0[:]; pval=false) #here pval=true will break your computer
 end
 return t_dist
end
#
function df_critval(T; confidence_level::Float64=0.05, trend::String="constant")
  quantile(simulate_dfdist(T, trend=trend), confidence_level; sorted=false)
end
#
#
@doc """
Calculates the Augmented Dickey-Fuller (1979) ADF statistic.

```math
  \\Delta y_t = \\alpha*y_{t-1} + \\sum^{pmax}_{j=1} \\beta_j*\\Delta y_{t-j}
```

Required Input Parameters
----------------------------------
`y::Vector{Float64}`: Observed time series in vector form.

Optional Input Parameters
----------------------------------
`pmax::Int`: Maximum number of lags to be considered in the correction for serial correlation.

`trend::String`: Deterministic trend specification. It determines the asymptotic distribution of the t-statistic. Possible values are `"constant"`, which considers a constant deterministic term; `"linear"`, which also considers a linear time trend; and `"none"`, which considers no deterministic detrending.

Output
----------------------------------
`df_stat::Float64`: ADF test statistic.

`df_pval::Float64`: p-value obtained from 10000 simulations.
""" ->
function adf_test(y::Vector{Float64}; pmax::Integer=15, trend::String="constant")
  bics = Vector{Float64}(pmax)
  for p ∈ 1:pmax
    β_adf, adf_stat, bics[p,:] = adfreg(y, p; trend=trend)
  end
  bic_star, bic_index = findmin(bics)
  β_adf, adf_stat, bics = adfreg(y, bic_index; trend=trend)
  adf_pval = mean(adf_stat.>simulate_dfdist(length(y), trend=trend))
  return adf_stat, adf_pval
end
# http://www.eviews.com/help/helpintro.html#page/content/advtimeser-Unit_Root_Testing.html
function adfreg(y::Vector{Float64}, p;  trend::String="constant")
  x =  [ar_lagmatrix(y, 1, trend="none")[p+1:end] ar_lagmatrix(diff(y), p, trend=trend)] #ar_lagmatrix(y, 1, trend="none")
  dy = diff(y)[p+1:end]
  β, t_stat, bic = df_ols(dy, x)
  return β[1], t_stat[1], bic
end
#
#
function ar_lagmatrix(x::Vector{Float64}, p::Integer; trend::String="none")
  T = length(x)
  xlag = Array{Float64}(T-p, p)
  for ip ∈ 1:p
    xlag[:, ip] = x[(p-ip+1):(end-ip)]
  end
  (trend=="none") && (return xlag)
  if trend == "constant"
    xlag = hcat( ones(T-p), xlag)
    return xlag
  end
  if trend == "linear"
    xlag = hcat( ones(T-p), collect(linspace(1,T-p,T-p)), xlag )
    return xlag
  end
end

function df_ols(y::Vector{Float64}, x::Array{Float64,})
  T, k = size(x)
  β = inv(x'x)*(x'y)
  ϵ = y - x*β
  yhat = x*β
  sse = ϵ'ϵ
  σ2_y = sse/(T-k)
  Σₓ = inv(x'x)
  Σᵦ = σ2_y.*Σₓ
  t_stat = β./sqrt(diag(Σᵦ))
  bic = (T-k).*log(sse./(T-k)) .+ k.*log(T-k)
  return β, t_stat, bic
end
#
@doc """
Simulates a stochastic process {yₜ} following an ARMA(1,1) with parameters φ and Θ for t = 1, ..., T periods.

```math
  y_t = \\mu_0 + \\mu_1 + \\phi*y_{t-1} + \\epsilon_t + \\theta*\\epsilon_{t-1}, \\quad \\epsilon_t \\sim \\mathbf{N}(0, \\sigma^2)
```

Required Input Parameters
----------------------------------
none

Optional Input Parameters
----------------------------------
`mu::Float64`: AR(1) process constant term. (default=0.0)

`lin_trend::Float64`: AR(1) process linear time trend. (default=0.0)

`ar::Float64`: AR(1) autorregressive parameter. (default=0.5)

`ma::Float64`: MA(1) autorregressive parameter. (default=0.5)

`T::Int64`: Length of simulated time series. (default=100)

`y_0::Float64`: Initial value of the series. (default=0.0)

`stderr::Float64` : Normally distributed innovation standard error. (default=1.0)

Output
----------------------------------
`y::Vector{Float64}`: Simulated ARMA(1,1) time series of length T.

""" ->
function simul_arma(;mu::Float64 = 0., lin_trend::Float64 = 0., ar::Float64 = .5, ma::Float64 = .5, T::Integer = 100, y_0::Float64 = 0.0, stderr::Float64 = 1.0)
  #y_t = μ + ϕy_t-1 + u_t/ u_t = ɛ_t + Θɛ_t-1
  y = zeros(T); ɛ = rand(Normal(0.0,stderr),T); u = zeros(T);
  for t in 2:T
    u[t, :] = ɛ[t, :] + ma*ɛ[t-1, :] #MA
    y[t, :] = mu + lin_trend*t + ar*y[t-1, :] + u[t, :] #AR
  end
  return y
end
