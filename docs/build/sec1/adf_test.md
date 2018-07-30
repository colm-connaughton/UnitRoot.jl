<a id='UnitRoot.adf_test' href='#UnitRoot.adf_test'>#</a>
**`UnitRoot.adf_test`** &mdash; *Function*.



Calculates the Augmented Dickey-Fuller (1979) ADF statistic.

$$
  \Delta y_t = \alpha*y_{t-1} + \sum^{pmax}_{j=1} \beta_j*\Delta y_{t-j}
$$

**Required Input Parameters**

`y::Vector{Float64}`: Observed time series in vector form.

**Optional Input Parameters**

`pmax::Int`: Maximum number of lags to be considered in the correction for serial correlation.

`trend::String`: Deterministic trend specification. It determines the asymptotic distribution of the t-statistic. Possible values are `"constant"`, which considers a constant deterministic term; `"linear"`, which also considers a linear time trend; and `"none"`, which considers no deterministic detrending.

**Output**

`df_stat::Float64`: ADF test statistic.

`df_pval::Float64`: p-value obtained from 10000 simulations.


<a target='_blank' href='https://github.com/p-chaim/UnitRoot.jl/tree/61100876c17fa088d9cdb045021baf18a0be3314/src/adf.jl#L63-L85' class='documenter-source'>source</a><br>

