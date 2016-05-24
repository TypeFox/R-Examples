#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# FILE:         fit_exponential.R
#
# AUTHOR:       Philippe Massicotte
#
# DESCRIPTION:  Fit an exponential curve to CDOM data.
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

#' Fit an exponential model to CDOM data.
#'
#' @details \deqn{y = a0 + e^{(-S(x - \lambda_0))} + K}
#'
#' @param wl The wavelength vector.
#' @param absorbance The absorbance vector.
#' @param wl0 The reference wavelength (ex.: 350).
#' @param startwl The starting wavelength (ex.: 240).
#' @param endwl The ending wavelength (ex.: 600).
#'
#' @return A list containing:
#' \describe{
#'   \item{params}{A data frame with values of fitted parameters.}
#'   \item{r2}{R2 of the nls model.}
#'   \item{data}{A data frame with fitted (predicted) values of the model.}
#' }
#'
#' The function will return \code{NULL} if the model did not converged.
#' @export
#' @import minpack.lm
#' @importFrom broom tidy augment
#' @importFrom stats splinefun predict var coef lm na.omit
#'
#' @examples
#' # Fit an exponential model using the reference wavelength 350 between 190 and 900 nm.
#'
#' data(spectra)
#'
#' fit <- cdom_fit_exponential(spectra$wavelength, spectra$spc1, 350, 190, 900)
#' str(fit)
#'
#' plot(spectra$wavelength, spectra$spc1)
#' lines(spectra$wavelength, fit$data$.fitted, col = "red")

cdom_fit_exponential <- function(wl, absorbance, wl0 = 350, startwl, endwl){

  stopifnot(length(wl) == length(absorbance),
            is.numeric(absorbance),
            is.numeric(wl),
            is.vector(wl),
            is.vector(absorbance),
            is.numeric(wl0),
            is.numeric(startwl),
            is.numeric(endwl))

  if(missing(startwl)){startwl = min(wl)}
  if(missing(endwl)){endwl = max(wl)}

  #--------------------------------------------
  # Get a0 value.
  #--------------------------------------------
  sf <- splinefun(wl, absorbance)
  a0 <- sf(wl0)

  #--------------------------------------------
  # Extract CDOM data based on user inputs.
  #--------------------------------------------
  x <- wl[which(wl >= startwl & wl <= endwl)]
  y <- absorbance[which(wl >= startwl & wl <= endwl)]

  #--------------------------------------------
  # Fit the data.
  #--------------------------------------------
  control <- list(minFactor = 1e-10,
                  warnOnly = FALSE,
                  maxiter = 1024,
                  maxfev = 600)

  tryCatch(
    {
      fit <- nlsLM(y ~ a0 * exp(-S * (x - wl0)) + K,
                   start = c(S = 0.02, K = 0.01, a0 = a0),
                   lower = c(S = 0, K = -Inf, a0 = 0),
                   upper = c(S = 1, K = Inf, a0 = max(y)),
                   control = control)

      r2 <- 1 - sum((y - predict(fit))^2) / (length(y) * var(y))

      return(list(params = tidy(fit), r2 = r2, data = augment(fit)))

    }, error = function(cond) {

      message("Error in fit_exponential() when trying to fit. Check your data.")
      message(cond)

      # Choose a return value in case of error
      return(NULL)

    },warning = function(cond) {

      message(cond)
      message("Warning in fit_exponential().")

      # Choose a return value in case of warning
      return(NULL)

    },finally = {

    }
  )
}
