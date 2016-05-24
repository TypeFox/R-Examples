## Draw points from the smoothed log-concave
'rslcd' <- function (n = 1, lcd, A=hatA(lcd), method = c("Independent","MH")) {
  method <- match.arg(method)
  if(class(lcd) != "LogConcDEAD") {
    stop("error: lcd must be of class LogConcDEAD")
  }
  d <- ncol(lcd$x)
  if(!(is.matrix(A) && (all(eigen(A)$values > .Machine$double.eps^0.5*4)) && (dim(A)[1]==d) && (dim(A)[2]==d))) {
    stop("error: Hat matrix A must be ",d," by ",d," positive definite")
  }
  samples <- rlcd(n, lcd, method) + rmvnorm(n, rep(0,d), A)
  return(samples)
}














