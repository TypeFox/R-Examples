###
### qimom.R
###

qimom <- function(p, V1=1, tau=1, nu=1) {

  ans <- double(length(p))
  ans[p<=.5] <- -sqrt(tau*V1/qgamma(2*p[p<=.5], nu/2, 1))
  ans[p>.5] <- sqrt(tau*V1/qgamma(2*(1-p[p>.5]), nu/2, 1))
  ans
}

