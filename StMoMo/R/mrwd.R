
#' Fit a Multivariate Random Walk with Drift
#' 
#' Fits a Multivariate Random Walk with Drift to \code{x}, a
#' multivariate time series.
#' 
#' For further information on the Multivariate Random Walk with
#' drift see Appendix B in Haberman and Renshaw (2011).
#' 
#' @param x numeric matrix with a multivariate time series.
#' Series are arranged in rows with columns representing time.
#' 
#' @return an object of class \code{"mrwd"} with components:
#' \item{drift}{ a vector with the estimated drift.}
#' \item{sigma}{ a matrix with the estimated variance covariance 
#'  matrix.}
#' \item{fitted}{ fitted values.}
#' \item{residuals}{ residuals from the fitted model. That is 
#' observed minus fitted values.}
#' \item{x}{the original time series.}
#' 
#' @references
#' Haberman, S., & Renshaw, A. (2011). A comparative study of parametric 
#' mortality projection models. Insurance: Mathematics and Economics, 
#' 48(1), 35-55. 
#' 
#' @export
mrwd <- function(x) {  
  x <- as.matrix(x)
  if (ncol(x) == 1L) 
    x <- t(x)
  nYear <- ncol(x)
  N <- nrow(x)
  d <- t(colMeans(diff(t(x)), na.rm = TRUE))
  fits <- cbind(array(NA, c(N, 1)), x[, -nYear] + array(d, c(N, nYear - 1)))
  res <- x - fits
  dimnames(fits) <- dimnames(res) <- dimnames(x)
  sigma <- cov(t(res), use = "complete.obs")
  structure(list(drift = d, sigma = sigma, fitted = fits, residuals = res, 
                 x = x), class = "mrwd")  
}

#' Forecast a Multivariate Random Walk with Drift
#' 
#' Returns forecasts and other information for a Multivariate 
#' Random Walk with Drift model.
#' 
#' @param object an object of class \code{"mrwd"}.
#' @param h Number of periods for forecasting.
#' @param level confidence level for prediction intervals.
#' @param fan if \code{TRUE}, level is set to \code{seq(50, 99, by = 1)}. 
#' This is suitable for fan plots.
#' @param ... other arguments.
#' 
#' @return An object of class \code{"mrwdForecast"} with components:
#' \item{model}{a list containing information about the fitted model.}
#' \item{mean}{ array with the central forecast.}
#' \item{lower}{ three dimensional array with lower limits for prediction 
#'  intervals.}
#' \item{upper}{ three dimensional array with upper limits for prediction 
#'  intervals.}
#'  \item{level}{ the confidence values associated with the prediction 
#'  intervals.}
#'  @export
forecast.mrwd <- function(object, h = 10, level = c(80,95), fan = FALSE, ...) {
  
  x <- object$x
  nn <- 1:h
  nYear <- ncol(x)
  N <- nrow(x)  
  yearsFor <- (as.numeric(colnames(x)[nYear]) + 1):(as.numeric(colnames(x)[nYear]) + h)
  
  mean <- x[, nYear] + t(array(nn, c(h, N))) * array(object$drift, c(N,h))
  rownames(mean) <- rownames(x)
  colnames(mean) <- yearsFor
  
  if (fan) 
    level <- seq(51, 99, by = 3)
  else {
    if (min(level) > 0 & max(level) < 1) 
      level <- 100 * level
    else if (min(level) < 0 | max(level) > 99.99) 
      stop("Confidence limit out of range")
  }
  nn <- 1:h
  se <- sqrt(t(array(nn, c(h, N))) * array(diag(object$sigma), c(N, h)))
  nconf <- length(level)
  z <- qnorm(0.5 + level / 200)
  lower <- upper <- array(NA, c(N, h, nconf), 
                          dimnames = list(rownames(x), yearsFor, 
                                          paste(level, "%", sep = "")))
  for (i in 1:nconf) {
    lower[, , i] <- mean - z[i] * se
    upper[, , i] <- mean + z[i] * se
  }
      
  structure(list(model = object, level = level, mean = mean, lower = lower, 
                 upper = upper), class = "mrwdForecast")    
}

#' Simulate a Multivariate Random Walk with Drift
#' 
#' Returns one simulated path of the Multivariate 
#' Random Walk with Drift model in \code{object}.
#' 
#' @param object An object of class \code{"mrwd"}.
#' @param nsim number of periods for the simulated series.
#' @param seed either \code{NULL} or an integer that will be used in a 
#' call to \code{\link{set.seed}} before simulating the time series. 
#' The default, \code{NULL} will not change the random generator state.
#' @param ... other arguments.
#' 
#' @export
simulate.mrwd <- function(object, nsim = 10, seed = NULL, ...) {
  
  if (!exists(".Random.seed", envir = .GlobalEnv)) 
    runif(1)
  if (is.null(seed)) 
    RNGstate <- .Random.seed
  else {
    R.seed <- .Random.seed
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  
  #generate innovations
  x <- object$x
  nn <- 1:nsim
  nYear <- ncol(x)
  N <- nrow(x)  
  u <- MASS::mvrnorm(nsim, rep(0, N), Sigma = object$sigma)
  su <- t(apply(u, 2, cumsum))
  sim <- x[, nYear] + t(array(nn, c(nsim, N))) * array(object$drift, 
                                                       c(N, nsim)) + su
  
  yearsSim <- (as.numeric(colnames(x)[nYear]) + 1):(as.numeric(colnames(x)[nYear]) + nsim)
  rownames(sim) <- rownames(x)
  colnames(sim) <- yearsSim
  sim  
}
