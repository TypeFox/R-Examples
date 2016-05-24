pred.logser=function(x, alpha, J, S){
  if(missing(alpha)&missing(J)&missing(S)){
    stop("Please provide at least two of these: alpha, J, S")
  }
  if(missing(alpha)){
      f1 <- function(a){ S + a * log((a/(a + J))) }
      sol <- uniroot(f1, interval = c(1/J, J))
      alpha <- sol$root
  }
  if(missing(J)){
    J <- alpha*exp(S/alpha) - alpha
  }
  X <- J/(alpha + J)
  alpha*(X^x)/x
}
