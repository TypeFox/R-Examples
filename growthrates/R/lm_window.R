#' Fit Exponential Growth Model with a Heuristic Linear Method
#'
#' Helper functions for handling linear fits.
#'
#' The functions are used by a heuristic linear approach, similar to the
#' ``growth rates made easy''-method of Hall et al. (2013).
#'
#'
#' @param x vector of independent variable (e.g. time).
#' @param y vector of dependent variable (concentration of organisms).
#' @param i0 index of first value used for a window.
#' @param h with of the window (number of data).
#'
#' @return linear model object (\code{lm_window}
#'   resp. vector with parameters of the fit (lm_parms).
#'
#' @references Hall, B. G., H. Acar and M. Barlow 2013. Growth Rates Made Easy.
#'   Mol. Biol. Evol. 31: 232-238 doi:10.1093/molbev/mst197
#'
#' @keywords internal
#'
#' @export lm_window
#'
lm_window <- function(x, y, i0, h=5) {
  ## function that fits a linear model
  ## to a selected window of data from a series, e.g. 5 points
  #cat(h, "\n")
  x <- x[i0 - 1 + (1:h)]
  y <- y[i0 - 1 + (1:h)]
  m <- lm(y ~ x)
  return(m)
}

#' @param m linear model (\code{lm}) object
#'
#' @rdname lm_window
#' @export lm_parms
#'
#' @keywords internal
#'
lm_parms <- function(m) {
  sm   <- summary(m)
  a    <- sm$coefficients[1, 1] # intercept
  b    <- sm$coefficients[2, 1] # slope
  b.se <- sm$coefficients[2, 2] # SE of slope
  r2   <- sm$r.squared
  c(a=a, b=b, se=b.se, r2=r2, cv=b.se/b)
}

