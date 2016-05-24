##Prior Models
"prior.quant" <-
  function(prob = 10^-(1:3), shape, scale, trendsd = 0)
{
  if(length(prob) != 3 || mode(prob) != "numeric")
    stop("`prob' must be a numeric vector of length three")
  if(min(prob) <= 0 || max(prob) >= 1)
    stop("`prob' must contain values in (0,1)")
  if(any(diff(prob) >= 0))
    stop("`prob' must be a vector of decreasing values")
  if(length(shape) != 3 || mode(shape) != "numeric")
    stop("`shape' must be a numeric vector of length three")
  if(length(scale) != 3 || mode(scale) != "numeric")
    stop("`scale' must be a numeric vector of length three")
  if(any(c(shape, scale) <= 0))
    stop("`shape'  and `scale' must contain positive values")
  if(length(trendsd) != 1 || mode(trendsd) != "numeric")
    stop("`trendsd' must be a numeric vector of length one")
  if(trendsd < 0)
    stop("`trendsd' must be positive")
  structure(list(prior = "dprior.quant", prob = prob, shape = shape,
                 scale = scale, trendsd = trendsd), class = "evprior") 
}

"prior.prob" <-
  function(quant, alpha, trendsd = 0)
{
  if(length(quant) != 3 || mode(quant) != "numeric")
    stop("`quant' must be a numeric vector of length three")
  if(any(diff(quant) <= 0))
    stop("`quant' must be a vector of increasing values")
  if(length(alpha) != 4 || mode(alpha) != "numeric")
    stop("`alpha' must be a numeric vector of length four")
  if(any(alpha <= 0))
    stop("`alpha' must contain positive values")
  if(length(trendsd) != 1 || mode(trendsd) != "numeric")
    stop("`trendsd' must be a numeric vector of length one")
  if(trendsd < 0)
    stop("`trendsd' must be positive")
  structure(list(prior = "dprior.prob", quant = quant, alpha = alpha,
                 trendsd = trendsd), class = "evprior")
}

"prior.norm" <-
  function(mean, cov, trendsd = 0)
{
  if(length(mean) != 3 || mode(mean) != "numeric")
    stop("`mean' must be a numeric vector of length three")
  if(!is.matrix(cov) || any(dim(cov) != 3) || mode(cov) != "numeric")
    stop("`cov' must be a symmetric three by three matrix")
  if(length(trendsd) != 1 || mode(trendsd) != "numeric")
    stop("`trendsd' must be a numeric vector of length one")
  if(trendsd < 0)
    stop("`trendsd' must be positive")
  if(any(abs(cov - t(cov)) > .Machine$double.eps^0.5))
    warning("`cov' may not be symmetric")
  eg <- eigen(cov, symmetric = TRUE, only.values = TRUE)$values
  if(any(eg <= 0))
    warning("`cov' may not be positive definite")
  
  icov <- solve(cov)
  icov <- icov[row(icov) >= col(icov)] 
  structure(list(prior = "dprior.norm", mean = mean, icov = icov,
                 trendsd = trendsd), class = "evprior")
}

"prior.loglognorm" <-
  function(mean, cov, trendsd = 0)
{
  if(length(mean) != 3 || mode(mean) != "numeric")
    stop("`mean' must be a numeric vector of length three")
  if(!is.matrix(cov) || any(dim(cov) != 3) || mode(cov) != "numeric")
    stop("`cov' must be a symmetric three by three matrix")
  if(length(trendsd) != 1 || mode(trendsd) != "numeric")
    stop("`trendsd' must be a numeric vector of length one")
  if(trendsd < 0)
    stop("`trendsd' must be positive")
  if(any(abs(cov - t(cov)) > .Machine$double.eps^0.5))
    warning("`cov' may not be symmetric")
  eg <- eigen(cov, symmetric = TRUE, only.values = TRUE)$values
  if(any(eg <= 0))
    warning("`cov' may not be positive definite")
  
  icov <- solve(cov)
  icov <- icov[row(icov) >= col(icov)] 
  structure(list(prior = "dprior.loglognorm", mean = mean, icov = icov,
                 trendsd = trendsd), class = "evprior")
}

## PRIOR DENSITIES

"dprior.quant" <-
  ## Computes log prior density for (mu,sigma,xi) based on gamma distributions for quantile differences corresponding to three specified probabilities. Includes optional normal trend for location.
  function(par, prob = 10^-(1:3), shape, scale, trendsd)
{
  .C("dprior_quant",
     par, prob, shape, scale, trendsd, dns = double(1),
     PACKAGE = "evdbayes")$dns
}

"dprior.prob" <-
  ## Computes log prior density for (mu,sigma,xi) based on beta distributions for probability ratios corresponding to three specified quantiles. Includes optional normal trend for location.
  function(par, quant, alpha, trendsd)
{
  .C("dprior_prob",
     par, quant, alpha, trendsd, dns = double(1),
     PACKAGE = "evdbayes")$dns
}

"dprior.norm" <-
                                        # Computes log prior density for (mu,sigma,xi) based on a trivariate normal distribution for (mu,log(sigma),xi). Includes optional normal trend for location.
  function(par, mean, icov, trendsd)
{
  .C("dprior_norm",
     par, mean, icov, trendsd, dns = double(1),
     PACKAGE = "evdbayes")$dns
}

"dprior.loglognorm" <-
  ## Computes log prior density for (loc,scale,shape) based on a trivariate
  ## normal distribution for (log(loc),log(scale),shape). Includes optional
  ## normal trend for location.
  function(par, mean, icov, trendsd)
{
  .C("dprior_loglognorm",
     par, mean, icov, trendsd, dns = double(1),
     PACKAGE = "evdbayes")$dns
}
