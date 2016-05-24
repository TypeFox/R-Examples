rtnorm.slice = function(lo, hi, mu, sd, niter=20) {
  lo=lo-mu; hi=hi-mu
  if(is.finite(lo) && is.finite(hi)) {
    x = runif(1, lo, hi)
  } else if (is.finite(lo)) {
    x = lo + 0.25*sd
  } else {
    x = hi-0.25*sd
  }
  for (i in 1:niter) {
    u = runif(1, 0, exp(-0.5*(x/sd)^2))
    bd = sqrt(-2*sd*sd*log(u))
    x = runif(1, max(lo, -bd), min(hi, bd))
  }
  return(x+mu)
}

rtnorm = function(lo, hi, mu, sd, maxn=100) {
  if (max(((lo-mu)/sd), ((hi-mu)/sd))<4) {
    out = qnorm(runif(1, pnorm(lo, mu, sd), pnorm(hi, mu, sd)), mu, sd)
    if(!is.finite(out)) out = rtnorm.slice(lo, hi, mu, sd)
  } else {
    out = rtnorm.slice(lo, hi, mu, sd)
  }
  return(out)
}