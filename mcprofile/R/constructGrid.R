constructGrid <- function(object, CM, control){
  est <- coefficients(object)
  Kest <- CM  %*% est
  sdKest <- sqrt(diag(CM %*% vcov(object) %*% t(CM)))
  n <- NROW(object$y)
  p <- length(est)
  
  if (n > p){
    zmax <- sqrt(qf(1 - control$alpha/nrow(CM), 1, n - p))
  } else {
    zmax <- sqrt(qchisq(1 - control$alpha/nrow(CM), 1))
  }
  del <- control$del(zmax)
  sst <- c(-1*(control$maxsteps:1), 0:control$maxsteps)

  grid <- sapply(1:nrow(CM), function(k){
    sapply(1:length(sst), function(i){
      Kest[k] + sst[i] * del * sdKest[k]
    })
  })  
  return(grid)
}
