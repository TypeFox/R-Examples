XI.gen <- function(k, kd) {

  tempA <- - (1.0 + sqrt(kd)) / ((kd - 1.0)^(1.5))
  tempB <- tempA + sqrt(kd / (kd - 1.0))

  XI <- matrix(data = tempA, nrow = k-1L, ncol = k)

  XI[,1L] <- 1.0/sqrt(kd - 1.0)

  for( ii in 2L:k ) XI[ii-1L,ii] <- tempB

  return(XI)
}
