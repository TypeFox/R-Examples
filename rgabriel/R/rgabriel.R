###rgabriel
###Function for performing the procedure described in Gabriel's 1978 paper
###28 Dec 2013

rgabriel <- function (x, f, a = "alpha level") 
{
  if (a == "alpha level") {
    a  <-  0.05
  }
  if (class(f) != "factor") {
    f  <-  factor(f)
  }
  Level.Name  <-  levels(f)
  k  <-  length(Level.Name)
  input  <-  cbind(f[!sapply(is.na(x), all)], x[!sapply(is.na(x), all)])
  f  <-  factor(input[, 1])
  x  <-  input[, 2]
  n  <-  numeric(length  <-  nlevels(f))
  for (i in 1:length(n)) {
    n[i]  <-  length(f[f == i])
  }
  for (i in 1:length(n)) {
    if (i == 1) {
      df  <-  numeric(length = length(n))
    }
    df[i]  <-  n[i] - 1
  }
  s <- tapply(x,f,sd)
  dfstar <- sum(df)
  # the following code is written by Yihui according to Stoline's 1979 paper(equation 2.2). 
  psmm_x = function(x, c, r, nu) {
    snu = sqrt(nu)
    sx  = snu * x  # for the scaled Chi distribution
    lgx = log(snu) - lgamma(nu/2) + (1 - nu/2) * log(2) + (nu - 1) * log(sx) + (-sx^2/2)
    exp(r * log(2 * pnorm(c * x) - 1) + lgx)
  }
  psmm = function(x, r, nu) {
    res = integrate(psmm_x, 0, Inf, c = x, r = r, nu = nu)
    res$value
  }
  qsmm = function(q, r, nu) {
    r = (r * (r - 1)/2)
    if (!is.finite(nu)) return(qnorm(1 - .5 *(1 - q^(1/r))))
    res = uniroot(function(c, r, nu, q) {
      psmm(c, r = r, nu = nu) - q
    }, c(0, 100), r = r, nu = nu, q = q)
    res$root
  }
  SR <- qsmm(1-a, k, dfstar)
  vstar <- SR*s/sqrt(2*n)
  return(vstar)
}
