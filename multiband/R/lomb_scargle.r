#' Lomb Scargle method
#' 
#' \code{lomb_scargle} applies the Lomb Scargle method for performing weighted nonlinear least squares in a single band using magnitude measurements recorded at irregular intervals.
#' 
#' @param t times
#' @param m magnitudes
#' @param sigma errors
#' @param omega frequency
#' @param weights_flag boolean (TRUE to use inverse errors as weights; FALSE uses uniform weights)
#' @export
#' @examples
#' n <- 1e3
#' set.seed(12345)
#' sigma <- runif(n)
#' t <- cumsum(runif(n))
#' A <- 2.3
#' rho <- 0.1
#' omega <- 0.2
#' beta0 <- 0.3
#' m <- beta0 + A*sin(omega*t + rho) + sigma*rnorm(n)
#' plot(t,m,xlab='time',ylab='magnitude',main='Simulated Light Curve',pch=16)
#' 
#' ## Try several omega
#' nOmega <- 500
#' omega_seq <- seq(0.1,0.3,length=nOmega)
#' sol_ls <- vector(mode="list",length=nOmega)
#' RSS_seq <- double(nOmega)
#' for (i in 1:nOmega) {
#'   sol_ls[[i]] <- lomb_scargle(t,m,sigma,omega_seq[i])
#'   RSS_seq[i] <- sol_ls[[i]]$RSS
#' }
#' 
#' plot(omega_seq,RSS_seq,xlab=expression(omega),ylab='RSS',pch=16)
#' ix_min <- which(RSS_seq==min(RSS_seq))
#' sol_final <- sol_ls[[ix_min]]
lomb_scargle <- function(t,m,sigma,omega,weights_flag=TRUE) {
  n <- length(t)
  w <- 1/(sigma**2)

  ## Build 3-by-3 linear system of equations Bz = d
  X <- matrix(0,n,3)
  X[,1] <- 1
  X[,2] <- sin(omega*t)
  X[,3] <- cos(omega*t)

  if (weights_flag) {
    Wd <- matrix(0,n,3)
    Wd[,1] <- w
    Wd[,2] <- w
    Wd[,3] <- w
    B <- t(X)%*%(Wd*X)
    d <- t(X)%*%(m*w)
   } else {
    B <- t(X)%*%X
    d <- t(X)%*%m
  }

  z <- solve(B,d)
  
  rho <- atan2(z[3],z[2])
  A <- sqrt(z[2]**2 + z[3]**2)
  r <- m - c(X%*%z)

 if (weights_flag) {
   RSS <- sum(w*(r**2))
 } else {
   RSS <- sum(r**2)
 }
  return(list(beta0=z[1],A=A,rho=rho,RSS=RSS))
}
