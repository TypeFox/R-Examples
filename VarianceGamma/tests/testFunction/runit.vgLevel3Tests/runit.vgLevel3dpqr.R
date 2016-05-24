
## Inversion test: check qdist(pdist(x))=x
test.vgL3dpqrInversionqp <- function () {
  for (i in 1:nrow(testParam)) {
    param <- testParam[i,]
    x <- rvg(nqp, param = param)
    for (j in 1:length(x)) {
      xFunc <- qvg(pvg(x[j], param = param), param = param)
      checkTrue(abs(x[j] - xFunc) < errorThresholddpqrI,
        msg = paste(param[1], param[2], param[3], param[4], x[j]))
    }
  }
}

test.vgL3dpqrInversionpq <- function () {
  for (i in 1:nrow(testParam)) {
    qs <- c(0.001,0.01,0.025,0.05,0.1,0.2,0.4,0.5,0.6,0.8,0.9,0.95,0.975,0.99,
    0.999)
    param <- testParam[i,]
    for (j in 1:length(qs)){                                
      qFunc <- pvg(q = qvg(p = qs[j],param = param), param = param)
      checkTrue(abs(qs[j] - qFunc) < errorThresholddpqrI,
        msg = paste(param[1], param[2], param[3], param[4], qs[j]))
    }
  }
}


## dpqr random test based on the test in base R
## RNG tests using DKW inequality for rate of convergence
##
## P(sup | F_n - F | > t) < 2 exp(-2nt^2)
##
## The 2 in front of exp() was derived by Massart. It is the best possible
## constant valid uniformly in t,n,F. For large n*t^2 this agrees with the
## large-sample approximation to the Kolmogorov-Smirnov statistic.
## 
test.vgL3dpqrRandom <- function() {
  for (i in 1:nrow(testParam)) {
    param <- testParam[i,]
    x <- rvg(n = n, param = param)
    tx <- table(x)
    xi <- as.numeric(names(tx))
    f <- pvg(xi, param = param)
    fhat <- cumsum(tx)/n
    sup <- max(abs(fhat-f))
    superror <- signif(sup,2)
    pvalue <- min(1,round(2*exp(-2*n*sup*sup),4))
    checkTrue((sup < sqrt(log(thresholddpqrR/2)/(-2*n))), 
      msg = paste(param[1], param[2], param[3], param[4], superror, pvalue))
  }
}    


