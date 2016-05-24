
logLik.sma <- function(object, ...){
  
  n.groups <- length(object$groups)
  ns <- as.numeric(object$n)
  if (n.groups==1)
  {
    s <- cov(object$data, use="complete.obs")
    z <- c(s[1,1],s[2,2],s[1,2]) * (ns-1) / ns
    logL <- get.lr(ns,z,object$coef[[1]][2,1],object$method)
  }  
  else
  {
    logLs <- rep(NA,n.groups)
    for (i.group in 1:n.groups)
    {
      s.i <- cov(object$data[object$data[,3]==object$groups[i.group],1:2], use="complete.obs")
      s.i <- s.i * (ns[i.group]-1) / ns[i.group]
      z.i <- c(s.i[1,1],s.i[2,2],s.i[1,2])
      logLs[i.group] <- get.lr(ns[i.group],z.i,object$coef[[i.group]][2,1],object$method)
    }
    logL <- sum(logLs)
  }
  if(object$gt=="elevcom"||object$gt=="shiftcom")
    attributes(logL)$df <- 4*n.groups + 1
  else
    attributes(logL)$df <- 5*n.groups
  attributes(logL)$nobs <- sum(ns)
  class(logL) <- "logLik"
  return( logL )
}

get.lr = function(n, z, b, method, lambda=1)
{
  if (method == 1 | method == "SMA") {
    if (b == 0) {
      b <- 10^-6
    }
    l1b <- (b^2 * z[2] + 2 * b * z[3] + z[1])/2/abs(b)
    l2b <- (b^2 * z[2] - 2 * b * z[3] + z[1])/2/abs(b)
  }
  else if (method == 2 | method == "MA" | method == 3 | method == "lamest") {
    l1b <- (lambda^2 * z[2] + 2 * lambda * b * z[3] + 
              b^2 * z[1])/(lambda + b^2)
    l2b <- (b^2 * z[2] - 2 * b * z[3] + z[1])/(lambda + b^2)
  }
  logL <- - n * ( log(l1b * l2b)/2 + 1/2 + log(2*pi) )
  return(logL)
}  


### test code ####  
# library(smatr)
# data(leaflife)
# ft=sma(longev ~ lma, log='xy', data=leaflife)
# print( logLik(ft) )
# 'log Lik.' 80.14958 (df=5)


# leaf.low.soilp <- subset(leaflife, soilp == 'low')

# com.test <- sma(longev~lma*rain, log="xy", data=leaf.low.soilp)
# print( logLik(com.test) )
# 'log Lik.' 45.13477 (df=10)

# elev.res <- sma(longev~lma+rain, log="xy", data=leaf.low.soilp)
# print( logLik(elev.res) )
# 'log Lik.' 43.68789 (df=9)

# shift.res <- sma(longev~lma+rain, type="shift", log="xy", data=leaf.low.soilp)
# print( logLik(shift.res) )
# 'log Lik.' 43.68789 (df=9)
# print(AIC(shift.res))
# [1] -69.37578
# print(BIC(shift.res))
# [1] -57.71325