chi2ncp <- function(alpha, beta, df = 1) {
  cv <- qchisq(alpha, df = df, lower.tail = FALSE) # critical value for chisq statistic
  ncpmax <- 1 # find upper limit by doubling until bracket achieved
  while (pchisq(cv, df = df, lower.tail = FALSE, ncp = ncpmax) < beta) ncpmax <- ncpmax*2  
  return(uniroot(function(ncp)
                 return(pchisq(cv, df = df, lower.tail = FALSE, ncp = ncp) - beta),
                 c(0, ncpmax))$root) # line search
}
