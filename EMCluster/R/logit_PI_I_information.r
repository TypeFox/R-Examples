### Get the Cov(logit pi_k) and 95% CE.

### Paritial all logit(pi_k). Return a K * (K - 1)
partial.logit.p <- function(PI){
  K <- length(PI)
  d.dPI <- 1 / PI + 1 / (1 - PI)
  ret <- diag(d.dPI[-K])
  ret <- rbind(ret, rep(-d.dPI[K], K - 1))
  ret
} # End of partial.logit.p().

get.cov.logit.PI <- function(x, emobj, cov.param = NULL){
  K <- emobj$nclass
  PI <- emobj$pi

  nabla.logit.PI <- partial.logit.p(PI)
  if(is.null(cov.param)){
    cov.param <- get.cov.param(x, emobj)$cov
  }

  cov.logit.PI <- nabla.logit.PI %*% cov.param[1:(K-1), 1:(K-1)] %*%
                   t(nabla.logit.PI)
  cov.logit.PI
} # End of get.cov.logit.PI().

