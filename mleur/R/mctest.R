mctest <-
function(y, type=c("p", "n"), NumRep=1000, bootQ = FALSE){
  n <- length(y)
  tobs <- testStatUM(y, type=type)
  tstat <- numeric(NumRep)
  if (bootQ) {
      yc <- y-mean(y)
      phiHat <- ar1est(yc)
      res <- yc[-1]-phiHat*yc[-n]
      for (i in 1:NumRep) {
        z <- cumsum(sample(res, size=n, replace=TRUE))
        tstat[i] <- testStatUM(w, type=type)
    }
  }
  else {
    for (i in 1:NumRep) {
      z <- cumsum(rnorm(n))
      w <- z - mean(z)
      tstat[i] <- testStatUM(w, type=type)
    }
  }
  pval <- sum(tstat<tobs)/(NumRep+1)
  pval
  }

