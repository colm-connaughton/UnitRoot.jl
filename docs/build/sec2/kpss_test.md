<a id='UnitRoot.kpss_test' href='#UnitRoot.kpss_test'>#</a>
**`UnitRoot.kpss_test`** &mdash; *Function*.



Calculates the stationarity KPSS test statistic of Nyblom (1989) and Kwiatkowski et all (1992).

$$
 y_t = \xi_t + r_t + \epsilon_t
$$

$$
 r_t = r_{t-1} + u_t, \quad u_t\sim\mathcal{N}(0, \sigma^2_u)
$$

The KPSS test has a LM-type statistic given by

$$
  LM = \frac{\sum^T_t S_t}{\epsilon'\epsilon},
$$

where Sₜ is the partial sum of resids up to time t

$$
  S_t = \sum^t_{j=1} \epsilon_j
$$

**Required Input Parameters**

`y::Vector{Float64}`: Observed time series in vector form.

**Optional Input Parameters**

`trend::String`: Deterministic trend specification. It determines the asymptotic distribution of the t-statistic. Possible values are `"constant"`, which considers a constant deterministic term; amd `"linear"`, which also considers a linear time trend.

**Output**

`kpss_stat::Float64`: KPSS test statistic.

`kpss_pval::Float64`: p-value obtained from 10000 simulations.

```@example 1
  series = simul_arma()
  kpss_test(series)
```


<a target='_blank' href='https://github.com/p-chaim/UnitRoot.jl/tree/61100876c17fa088d9cdb045021baf18a0be3314/src/kpss.jl#L1-L43' class='documenter-source'>source</a><br>

