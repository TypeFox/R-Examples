normalApproximation <-
function(data, center, deviation, cutoff = 2.5) {
  ndx <- which(abs(data - center)/deviation < cutoff)
  center <- mean(data[ndx])
  devOld <- sd(data[ndx])
  devNew <- devOld + 1
  ndx <- which(abs(data - center)/devOld < cutoff)
  eps <- 1e-6
  repeat {
    if(length(ndx) < 2) {
      devNew <- .Machine$double.eps
      cutoff <- -cutoff
      break
    }
    center <- mean(data[ndx])
    devNew <- sd(data[ndx])
    if(abs(devNew - devOld) < eps) break
    ndx <- which(abs(data - center)/devNew < cutoff)
    devOld <- devNew
  }
  ans <- list()
  ans$center <- center
  ans$deviation <- devNew
  ans$cutoff <- cutoff
  return(ans)
}
