FindBestTransform <-
function(x){
  alphas <- seq(.01, 1, len=50)
  gof <- rep(NA, length(alphas))
  for(alpha in alphas){
    gof[alphas==alpha] <- GoodnessOfFit(x^alpha,type="mle")
  }
  return(alphas[which.min(abs(gof-(nrow(x)-1)*(ncol(x)-1)))])
}

