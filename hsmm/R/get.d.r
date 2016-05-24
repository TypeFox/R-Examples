# Calculate d
get.d <- function(rd, J, M, param){
  d <- c()
  for (j in 1:J){
    # rd = "non.parametric"
    if (rd == "nonp")
      d <- c(d, param$np[,j])
    # rd = "geometric"
    if (rd == "geom")
      d <- c(d, dgeom(c(0:(M - 1)), param$p[j]))
    # rd = "negative.binomial"
    if (rd == "nbinom")
      d <- c(d, dnbinom(c(0:(M - 1)), size = param$r[j], prob = param$pi[j]))
    # rd = "logarithmic"                                                                                                                     
    if (rd == "log")
      d <- c(d, dlog(c(0:(M - 1)), param$p[j]))
    # rd = "logarithmic.geometric"
#    if (rd == "logarithmic.geometric")
#      d <- c(d, dloggeom(c(0:(M - 1)), param$p.loggeom[j], param$theta[j]))
    # rd = "Poisson"
    if (rd == "pois")
      d <- c(d, dpois(c(0:(M - 1)), param$lambda[j]))
    }
  return(d)
  }
