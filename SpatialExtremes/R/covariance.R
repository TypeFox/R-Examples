covariance <- function(fitted, nugget, sill, range, smooth, smooth2 = NULL,
                       cov.mod = "whitmat", plot = TRUE, dist, xlab,
                       ylab, col = 1, ...){

  if (!missing(fitted)){
    cov.mod <- fitted$cov.mod
    smooth <- fitted$param["smooth"]
    range <- fitted$param["range"]
    nugget <- fitted$param["nugget"]
    sill <- 1 - fitted$param["nugget"]

    if (cov.mod == "caugen")
      smooth2 <- fitted$param["smooth2"]
  }

  if (cov.mod == "gauss")
    stop("''covariance'' is not implemented for the Smith's model")

  if (!(cov.mod %in% c("whitmat", "cauchy", "powexp", "bessel", "caugen")))
    stop("Invalid covariance model. ''cov.mod'' must be one of 'whitmat', 'cauchy', 'powexp', 'bessel', 'caugen'")
  
  if (cov.mod == "whitmat"){
    if ((smooth <= 0) || (range <= 0) || (smooth > 150) || (sill <= 0) || (nugget < 0))
      stop("invalid parameter for the whittle-matern covariance function")
    
    cov.fun <- function(dist) {
      idx <- dist == 0
      ans <- rep(nugget + sill, length(dist))
      ans[!idx] <- sill * 2^(1-smooth) / gamma(smooth) * (dist[!idx] / range)^smooth *
        besselK(dist[!idx] / range, smooth)
      dim(ans) <- dim(dist)
      return(ans)
    }
  }

  if (cov.mod == "cauchy"){
    if ((smooth <= 0) || (range <= 0) || (sill <= 0) || (nugget < 0))
      stop("invalid parameter for the cauchy covariance function")
    
    cov.fun <- function(dist){
      idx <- dist == 0
      ans <- rep(nugget + sill, length(dist))
      ans[!idx] <- sill * (1 + (dist[!idx] / range)^2)^-smooth
      dim(ans) <- dim(dist)
      return(ans)
    }
  }

  if (cov.mod == "caugen"){
    if ((smooth <= 0) || (range <= 0) || (sill <= 0) || (smooth2 <= 0) ||
        (smooth2 > 2) || (nugget < 0))
      stop("invalid parameter for the generalized cauchy covariance function")

    cov.fun <- function(dist){
      idx <- dist == 0
      ans <- rep(nugget + sill, length(dist))
      ans[!idx] <- sill * (1 + (dist[!idx] / range)^smooth2)^(-smooth/smooth2)
      dim(ans) <- dim(dist)
      return(ans)
    }
  }

  if (cov.mod == "powexp"){
    if ((smooth < 0) || (smooth > 2) || (range <= 0) || (sill <= 0) || (nugget < 0))
      stop("invalid parameter for the powered exponential covariance function")

    cov.fun <- function(dist){
      idx <- dist == 0
      ans <- rep(nugget + sill, length(dist))
      ans[!idx] <- sill * exp(-(dist[!idx] / range)^smooth)
      dim(ans) <- dim(dist)
      return(ans)
    }
  }

  if (cov.mod == "bessel"){
    if ((range <= 0) || (sill <= 0) || (nugget < 0))
      stop("invalid parameter for the Bessel covariance function")

    cov.fun <- function(dist){
      idx <- dist == 0
      ans <- rep(nugget + sill, length(dist))
      ans[!idx] <- sill * (2 * range / dist[!idx])^smooth * gamma(smooth + 1) *
        besselJ(dist[!idx] / range, smooth)
      dim(ans) <- dim(dist)
      return(ans)
    }
  }

  if (plot){

    if (missing(xlab))
      xlab <- "h"

    if (missing(ylab))
      ylab <- expression(gamma(h))

    if (is.null(list(...)$xlim)){
      tmp.fun <- function(dist) (cov.fun(dist) - 0.05)^2
      xlimsup <- optimize(tmp.fun, c(1e-6, 10 * sqrt(sill) * range))$minimum
    }

    else
      xlimsup <- list(...)$xlim[2]
    
    plot(cov.fun, from = 1e-3, to = xlimsup, xlab = xlab, ylab = ylab,
         col = col, ...)

    if (nugget > 0)
      points(0, nugget + sill, col = col)
  }

  if (!missing(dist)){
    cov.val <- cov.fun(dist)
    return(list(cov.fun = cov.fun, cov.val = cov.val))
  }

  else
    invisible(cov.fun)  
}
