#Source: fitar1.R
#exact MLE for mean and phi in AR1
#
fitar1 <- function(z, meanZeroQ=FALSE){
  n <- length(z)
  stopifnot(n > 3)
  if (meanZeroQ) {
    phi <- ar1est0(z)
    S <- z[1]*(1-phi^2)+sum((z[2:n]-phi*z[1:(n-1)])^2)
    LL <- (-0.5*n)*(log(S/n)-log(1-phi^2))
    ans <- c(0, phi, LL)
  }
  else {
    yc <- y <- z-mean(z)
    ym <- 0
    A <- sum(y[-1])
    B <- sum(y[-n])
    LL <- (-n)*log(sum(yc^2)/n)
    iter <- err <- 1
    while (err > 1e-6 && iter < 20) {
      LL0 <- LL
      iter <- iter+1
      phi <- ar1est0(yc)
      A1 <- y[1]*(1-phi^2)+(1-phi)*(A-phi*B)
      A2 <- 1-phi^2+(n-1)*(1-phi)^2
      ym <- A1/A2
      yc <- y-ym
      S <- yc[1]*(1-phi^2)+sum((yc[2:n]-phi*yc[1:(n-1)])^2)
      LL <- (-0.5*n)*(log(S/n)-log(1-phi^2))
      err <- abs(LL-LL0)
    }
    ans <- c(mean(z)+ym, phi, LL)
  }
  ans
}