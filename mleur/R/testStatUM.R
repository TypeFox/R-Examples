testStatUM <-
function(y, type=c("p", "n")){
#computes our test statistics based on exact MLE
#p - pivotal, n - normalized
  type <- match.arg(type, type)
  n <- length(y)
  yc <- y-mean(y)
  phiHat <- ar1est(yc)
  if (type=="n")
    ans <- n*(phiHat-1)
  else {
    numer <- (phiHat-1)*sqrt(sum(yc[-n]^2))
    denom <- sqrt(sum((yc[-1]-phiHat*yc[-n])^2)/(n-2))
    ans <- numer/denom
    }
ans
}

