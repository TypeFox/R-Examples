"pdfpe3" <-
function(x,para) {
  if(! are.parpe3.valid(para)) return()
  # SMALL IS USED TO TEST WHETHER SKEWNESS IS EFFECTIVELY ZERO
  SMALL <- sqrt(.Machine$double.eps)
  names(para$para) <- NULL
  MU    <- para$para[1] # location
  SIGMA <- para$para[2] # scale
  GAMMA <- para$para[3] # shape

  ALPHA <- 4/GAMMA^2
  ops <- options(warn = -1)
  ARG <- gamma(ALPHA)
  options(ops)
  
  # distribution is normal
  if(abs(GAMMA) <= SMALL | ARG == Inf) return(dnorm((x - MU)/SIGMA))

  # GAMMA != 0, distribution is nonnormal
  BETA  <-      0.5*SIGMA * abs(GAMMA)
  XI    <- MU -   2*SIGMA /     GAMMA

  Y <- sign(GAMMA)*(x - XI)
  f <- (Y)^(ALPHA - 1) * exp(-Y/BETA) / (BETA^ALPHA * ARG)

  names(f) <- NULL
  f[! is.finite(f)] <- NA
  f[is.na(f)] <- 0 # decision Dec. 2015
  return(f)
}
