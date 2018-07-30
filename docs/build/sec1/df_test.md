<a id='UnitRoot.df_test' href='#UnitRoot.df_test'>#</a>
**`UnitRoot.df_test`** &mdash; *Function*.



Calculates classic Dickey-Fuller (1979) statistic.

$$
  \Delta y_t = \alpha*y_{t-1} + \epsilon_t, \quad \epsilon_t\sim\mathcal{N}(0,\sigma^2)
$$

**Required Input Parameters**

`y::Vector{Float64}`: Observed time series in vector form.

**Optional Input Parameters**

`trend::String`: Deterministic trend specification. It determines the asymptotic distribution of the t-statistic. Possible values are `"constant"`, which considers a constant deterministic term; `"linear"`, which also considers a linear time trend; and `"none"`, which considers no deterministic detrending.

**Output**

`df_stat::Float64`: Dickey-Fuller (1979) test statistic.

`df_pval::Float64`: p-value obtained from 10000 simulations.


<a target='_blank' href='https://github.com/p-chaim/UnitRoot.jl/tree/61100876c17fa088d9cdb045021baf18a0be3314/src/adf.jl#L1-L21' class='documenter-source'>source</a><br>

