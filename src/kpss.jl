@doc """
Calculates the stationarity KPSS test statistic of Nyblom (1989) and Kwiatkowski et all (1992).

```math
 y_t = \\xi_t + r_t + \\epsilon_t
```
```math
 r_t = r_{t-1} + u_t, \\quad u_t\\sim\\mathcal{N}(0, \\sigma^2_u)
```

The KPSS test has a LM-type statistic given by

```math
  LM = \\frac{\\sum^T_t S_t}{\\epsilon'\\epsilon},
```

where Sₜ is the partial sum of resids up to time t

```math
  S_t = \\sum^t_{j=1} \\epsilon_j
```

Required Input Parameters
----------------------------------
`y::Vector{Float64}`: Observed time series in vector form.

Optional Input Parameters
----------------------------------

`trend::String`: Deterministic trend specification. It determines the asymptotic distribution of the t-statistic. Possible values are `"constant"`, which considers a constant deterministic term; amd `"linear"`, which also considers a linear time trend.

Output
----------------------------------
`kpss_stat::Float64`: KPSS test statistic.

`kpss_pval::Float64`: p-value obtained from 10000 simulations.

```@example 1
  series = simul_arma()
  kpss_test(series)
```

"""->
function kpss_test(y::Vector{Float64}; trend::String="constant")
  T = size(y,1)
  if trend=="constant"
    x = ones(T)
    println("constant deteriministic detrending.")
  elseif trend=="linear"
    x = [ones(T) collect(1:T)]
    println("constant plus linear deteriministic detrending.")
  else
    error(" trend must be either, constant, or linear")
  end
  β, t_stat, resid = kpss_ols(y, x[:,:])
  S = cumsum(resid, 1)
  nominator = sum(S.^2)/(T^2)
  denominator = resid'resid/T
  kpss_stat = nominator./denominator
  kpss_pval = mean(kpss_stat.<kpss_dist(T, trend=trend))
  return kpss_stat[1], kpss_pval
end
#
#
function wiener2(nobs)
    e0 = rand(Normal(), nobs);
    x = ones(nobs,1);
    e1 = e0 - mean(e0);
    b = inv(x'x)*(x'e0);
    e2 = e0-b.*e0;
    y1 = cumsum(e1);
    y2 = cumsum(e2);
    intW2_1 = nobs^(-2.)*sum(y1.^2, 1);
    intW2_2 = nobs^(-2.)*sum(y2.^2, 1);
    return intW2_1, intW2_2
end

function kpss_dist(T; trend="constant")
  nsim = 10000
  kpss = Vector{Float64}(nsim)
  for isim ∈ 1:nsim
    c_tmp, trnd_tmp = wiener2(T)
    (trend=="constant") && (kpss[isim, :] = c_tmp)
    (trend=="linear") && (kpss[isim, :] = trnd_tmp)
  end
  return kpss
end

function kpss_ols(y::Array{Float64,}, x::Array{Float64,})
  T, k = size(x)
  β = inv(x'x)*(x'y)
  ϵ = y - x*β
  yhat = x*β
  sse = ϵ'ϵ
  σ2_y = sse/(T-k)
  Σₓ = inv(x'x)
  Σᵦ = σ2_y.*Σₓ
  t_stat = β./sqrt.(diag(Σᵦ))
  bic = (T-k).*log(sse./(T-k)) .+ k.*log(T-k)
  resid = y - x*β
  return β, t_stat, resid
end


#kpss_test(randn(200))
