acf.custom <- function(..., main, ylim) {
  obj <- acf(..., plot=FALSE)
  obj$acf <- obj$acf[-1,,,drop=FALSE]
  obj$lag <- obj$lag[-1,,,drop=FALSE]
  plot(obj, main=main, ylim=ylim)
}

pacf.custom <- function(..., main, ylim) {
  obj <- pacf(..., plot=FALSE)
  obj$acf <- obj$acf[-1,,,drop=FALSE]
  obj$lag <- obj$lag[-1,,,drop=FALSE]
  plot(obj, main=main, ylim=ylim)
}

mutual.custom <- function(..., main, ylim) {
  mi <- mutual(..., plot=FALSE)[-1]
  plot(1:length(mi), mi, type = "h", xlab = "lag", ylab = "AMI",
         main = main, ylim=ylim)
  abline(h = 0)
}