#' Fit Exponential Growth Model with a Heuristic Linear Method
#'
#' Determine maximum growth rates from the log-linear part of a growth curve using
#' a heuristic approach similar to the ``growth rates made easy''-method of
#' Hall et al. (2013).
#'
#' The algorithm works as follows:
#' \enumerate{
#'   \item Fit linear regressions to all subsets of \code{h} consecutive data
#'     points. If for example \eqn{h=5}, fit a linear regression to points
#'     1 \dots 5, 2 \dots 6, 3\dots 7 and so on. The method seeks the highest
#'     rate of exponential growth, so the dependent variable is of course
#'     log-transformed.
#'   \item Find the subset with the highest slope \eqn{$b_max$}{b_max} and
#'     include also the data points of adjacent subsets that have a slope of
#'     at least \eqn{$quota \cdot b_max$}{quota * b_max},
#'     e.g. all data sets that have at least 95\% of the maximum slope.
#'   \item Fit a new linear model to the extended data window identified in step 2.
#' }
#'
#' @param time vector of independent variable.
#' @param y vector of dependent variable (concentration of organisms).
#' @param h width of the window (number of data).
#' @param quota part of window fits considered for the overall linear fit
#'   (relative to max. growth rate)
#'
#' @return object with parameters of the fit.
#'
#' @references Hall, B. G., H. Acar and M. Barlow 2013. Growth Rates Made Easy.
#'   Mol. Biol. Evol. 31: 232-238 doi:10.1093/molbev/mst197
#'
#' @family fitting functions
#'
#' @examples
#' data(bactgrowth)
#'
#' splitted.data <- multisplit(bactgrowth, c("strain", "conc", "replicate"))
#' dat <- splitted.data[[1]]
#'
#' plot(value ~ time, data=dat)
#' fit <- fit_easylinear(dat$time, dat$value)
#'
#' plot(fit)
#' plot(fit, log="y")
#' plot(fit, which="diagnostics")
#'
#' fitx <- fit_easylinear(dat$time, dat$value, h=8, quota=0.95)
#'
#' plot(fit, log="y")
#' lines(fitx, pch="+", col="blue")
#'
#' plot(fit)
#' lines(fitx, pch="+", col="blue")
#'
#'
#' @rdname fit_easylinear
#' @export fit_easylinear
#'
fit_easylinear <- function(time, y, h=5, quota=0.95) {

  if (any(duplicated(time))) stop("time variable must not contain duplicated values")

  obs <- data.frame(time, y)
  obs$ylog <- log(obs$y)

  ## number of values
  N <- nrow(obs)

  ## repeat for all windows and save results in 'ret'
  ret <- matrix(0, nrow = N - h, ncol = 6)
  for(i in 1:(N - h)) {
    ret[i, ] <- c(i, with(obs, lm_parms(lm_window(time, ylog, i0 = i, h = h))))
  }

  ## indices of windows with high growth rate
  slope.quota <- quota * max(ret[ ,3])   # "3" is slope (b) because "i" is 1st
  candidates <- which(ret[ ,3] >= slope.quota)

  if(length(candidates) > 0) {
    tp <- seq(min(candidates), max(candidates) + h-1)
    m <- lm_window(obs$time, obs$ylog, min(tp), length(tp))
    p  <- c(lm_parms(m), n=length(tp))
  } else {
    p <- c(a=0, b=0, se=0, r2=0, cv=0, n=0)
  }
  # cat(strain, conc, replicate, p[2], "\n")
  #return(list(p=p, ndx=tp))
  obj <- new("easylinear_fit", FUN=grow_exponential, fit=m,
             par = c(y0 = unname(exp(coef(m)[1])), mumax = unname(coef(m)[2])), ndx = tp,
             obs = data.frame(time = obs$time, y = obs$y), rsquared = p["r2"])
  invisible(obj)
}
