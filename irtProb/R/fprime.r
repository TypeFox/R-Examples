`fprime` <-
function(x, FUN="FP", h=0.001, names=paste("x",c(1:length(x)),sep="")) {
 # For match.fun function see Venables and Ripley (2000, p. 68)
 FUN           <- match.fun(FUN)
 dim           <- length(x)
 fprime        <- numeric(dim)
 for (i in 1:dim) {
  dh        <- numeric(dim)
  dh[i]     <- dh[i] + h
  # Yakowitz et Szidarovszky (1986, p. 84)
  fprime[i] <- (FUN(x + dh) - FUN(x - dh))/(2*h)
  }
 names(fprime) <-  names
 (fprime)
 }

