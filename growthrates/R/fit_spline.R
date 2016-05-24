#' Fit Exponential Growth Model with Smoothing Spline
#'
#' Determine maximum growth rates from the first derivative of a smoothing spline.
#'
#'
#' @param time vector of independent variable.
#' @param y vector of dependent variable (concentration of organisms).
#' @param optgrid number of steps on the x-axis used for the optimum search .
#'  algorithm. The default should work in most cases, as long as the data are equally spaced.
#'  A smaller number may lead to non-detectable speed-up, but has the risk that
#'  the search gets trapped in a local minimum.
#' @param \dots other parameters passed to \code{\link{smooth.spline}}, see details.
#'
#' @return object with parameters of the fit
#'
#' @details The method was inspired by an algorithm of Kahm et al. (2010),
#'   with different settings and assumptions. In the moment, spline fitting
#'   is always done with log-transformed data, assuming exponential growth
#'   at the time point of the maximum of the first derivative of the spline fit.
#'
#'   All the hard work is done by function \code{\link{smooth.spline}} from package
#'   \pkg{stats}, that is highly user configurable. Normally, smoothness is
#'   automatically determined via cross-validation. This works well in many cases,
#'   whereas manual adjustment is required otherwise, e.g. by setting \code{spar}
#'   to a fixed value \eqn{[0, 1]} that also disables cross-validation.
#'
#' @references
#'
#' Kahm, M., Hasenbrink, G., Lichtenberg-Frate, H., Ludwig, J., Kschischo, M.
#' 2010. grofit: Fitting Biological Growth Curves with R.
#' Journal of Statistical Software, 33(7), 1-21. URL
#' \url{http://www.jstatsoft.org/v33/i07/}
#'
#' @family fitting functions
#'
#' @examples
#'
#' data(bactgrowth)
#' splitted.data <- multisplit(bactgrowth, c("strain", "conc", "replicate"))
#'
#' dat <- splitted.data[[2]]
#' time <- dat$time
#' y    <- dat$value
#'
#' ## automatic smoothing with cv
#' res <- fit_spline(time, y)
#'
#' plot(res, log="y")
#' plot(res)
#' coef(res)
#'
#' ## a more difficult data set
#' dat <- splitted.data[[56]]
#' time <- dat$time
#' y <- dat$value
#'
#' ## default parameters
#' res <- fit_spline(time, y)
#' plot(res, log="y")
#'
#' ## small optgrid, trapped in local minimum
#' res <- fit_spline(time, y, optgrid=5)
#' plot(res, log="y")
#'
#' ## manually selected smoothing parameter
#' res <- fit_spline(time, y, spar=.5)
#' plot(res, log="y")
#' plot(res, ylim=c(0.005, 0.03))
#'
#'
#' @rdname fit_spline
#' @export fit_spline
#'
fit_spline <- function(time, y, optgrid = length(time), ...) {

  #if (any(duplicated(time))) stop("x variable must not contain duplicated values")

  obs <- data.frame(time = time, y = y)
  ylog <- log(obs$y)

  ## fit smoothing spline with cross validation by default
  spl <- smooth.spline(time, ylog, ...)

  xnew <- seq(min(time), max(time), length.out = optgrid)
  delta.xnew <- diff(xnew)[1]

  ynew <- predict(spl, x=xnew)$y
  ## 1st derivative for n = optgrid candidate values
  dspl <- predict(spl, x=xnew, deriv=1)

  ## find approximate maximum of 1st derivative
  tmax <- dspl$x[which.max(dspl$y)]

  ## post-optimize maximum
  optmax <- function(x) predict(spl, x = x, deriv = 1)$y

  ## set search interval to tmax + 2 * delta.xnew
  ## but limited by [min, max] of time
  interval <- tmax + c(-2, 2) * delta.xnew
  interval[1] <- max(min(time), interval[1])
  interval[2] <- min(max(time), interval[2])

  tmax2 <- optimize(f = optmax, interval = interval,
                    maximum = TRUE)$maximum

  ## predict x, y and derivative at the maximum
  xy <- predict(spl, x = tmax2)
  dy <- predict(spl, x = tmax2, deriv = 1L)$y

  px <- xy$x; py <- xy$y

  ## Note: r2 applies to log transformed data
  r2 <- 1 - var(residuals(spl))/var(ylog)

  new("smooth.spline_fit", fit = spl,
      FUN = grow_exponential,
      par = c(y0 = exp(py)/exp((dy * px)), mumax = dy),
      xy = c(px, exp(py)), obs=obs, rsquared = r2)
}
