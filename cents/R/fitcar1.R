#Source: fitcar1.R
#
fitcar1 <- function(z, cL=-Inf, cU=Inf, verboseQ=FALSE){
  w <- z
  mw <- mean(w, na.rm=TRUE)
  w[is.na(w)] <- mw
  phi <- ar1est0(w-mw)
  err <- iter <- LL <- 1
  while (err > 1e-6 && iter < 50) {
    w <- z
    LL0 <- LL
    w[1] <- ifelse(is.na(w[1]), mw, w[1])
    for (i in 2:length(w))
      w[i] <- ifelse(is.na(w[i]), mw+phi*(w[i-1]-mw), w[i])
    ans <- fitar1(w)
    mw <- ans[1]
    phi <- ans[2]
    LL <- ans[3]
    err <- abs(LL-LL0)
    iter <- iter+1
    if (verboseQ){
      cat("iteration:", iter, "\n")
      cat("LL = ", LL, "\n")
      cat("phi =", phi, "\n")
      cat("mw =", mw, "\n")
    }
  }
  c(mw,phi,LL)
}