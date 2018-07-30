@doc """

Calculates dfgls point optimal unit-root test statistic of Elliott, Rothenberg, and Stock (1996) value of the largest autorregressive coefficient under the alternative hypothesis given by ρ = 1+c/T, where c=-7 if trend="constant" and c=-13.5 if trend="linear")

"""->
function dfgls_test(y::Vector{Float64}; nar::Integer=5, trend::String="constant")
  #rewrite of mw's dfgls.prc
  t = length(y)
  dftau = NaN*zeros(1, 2)

  # @ -- GLS estimate of constant term -- @
  c = -7
  rho = (1+c/t)
  yrho = [y[1]; y[2:t]-rho*y[1:t-1]] #@ Quasi Differenced Y series @
  xcrho = [1.0; (1-rho)*ones(t-1)] #@ Quasi Difference Constant @
  mu = ((xcrho'*xcrho).^(-1)).*(xcrho'yrho)

  c = -13.5
  rho = (1+c/t)
  yrho = [y[1]; y[2:t]-rho*y[1:t-1]] #@ Quasi Differenced Y series @
  xcrho = [1.0; (1-rho)*ones(t-1)] #@ Quasi Difference Constant @
  trnd = collect(1:t)
  xtrho = [trnd[1]; trnd[2:t]-rho*trnd[1:t-1]]
  x = [xcrho xtrho]
  mut = inv(x'x)*(x'yrho)
  #  @ -- detrend y -- @
  y1 = y - mu.*ones(t)                #@ GLS Demeaned @
  y2 = y - [ones(t,1) trnd]*mut        #@ GLS Detrended @
  dy1 = diff(y1)
  dy2 = diff(y2)
  #compute adf
  x1 = ar_lagmatrix(dy1, nar)
  x2 = ar_lagmatrix(dy2, nar)

  z = [y1[1:size(x1,1)] x1]
  y = z[:, 1]
  x = z[:, 2:end]
  xxi = inv(x'x)
  xy = x'y
  bet = xxi*xy
  ssr = (y'y)-(xy)'bet
  s2 = ssr/(size(y,1)-size(x,2))
  b = bet[1]
  vb = s2.*xxi[1,1]
  dfmu = (b-1)./sqrt(vb)

  z = [y1[1:size(x1,1)] x1]
  y = z[:, 1]
  x = z[:, 2:end]
  xxi = inv(x'x)
  xy = x'y
  bet = xxi*xy
  ssr = (y'y)-(xy)'bet
  s2 = ssr/(size(y,1)-size(x,2))
  b = bet[1]
  vb = s2.*xxi[1,1]
  dftau = (b-1)./sqrt(vb)
  (trend=="constant") && return dfmu[1]
  (trend=="linear") && return dftau[1]
  (trend=="both") && return dfmu[1], dftau[1]
  end
