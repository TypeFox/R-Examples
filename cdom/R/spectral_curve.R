#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# FILE:         spectral_curve.R
#
# AUTHOR:       Philippe Massicotte
#
# DESCRIPTION:  Calculate the spectral curve of CDOM spectra has proposed by
#               Loiselle et al. 2009
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

#'Calculate the spectral curve of CDOM spectra.
#'
#'@description Calculate the spectral curve of CDOM spectra has proposed by
#'  Loiselle et al. 2009.
#'
#'@inheritParams cdom_fit_exponential
#'
#'@param interval The interval used to claculate each slope (default = 21 nm).
#'
#'@param r2threshold The r2 threshold that determines if a slope is "valide".
#'  The default value is 0.8 meaning that the determination coefficient of the
#'  regression between log-transformed data and wavelength should be >= 0.8.
#'
#'@references \url{http://doi.wiley.com/10.4319/lo.2009.54.2.0590}
#'
#'@return A dataframe containing the centered wavelength, the calculated slope
#'  and the determination coefficient of the linear regression used to claculate
#'  the slope.
#'@export
#' @examples
#'data(spectra)
#'
#'res <- cdom_spectral_curve(spectra$wavelength, spectra$spc2)
#'plot(res$wl, res$s, type = "l")

cdom_spectral_curve <- function(wl, absorbance, interval = 21, r2threshold = 0.8) {

  stopifnot(length(wl) == length(absorbance),
            is.numeric(absorbance),
            is.numeric(wl),
            is.vector(wl),
            is.vector(absorbance),
            is.numeric(interval),
            is.numeric(r2threshold))

  #--------------------------------------------
  # Resample data by 1 nm increment.
  #--------------------------------------------
  sf <- splinefun(wl, absorbance)

  xx <- seq(from = min(wl), to = max(wl), by = 1)
  yy <- sf(xx)

  ## Adjust the offset of the spectral (we do not want negative values).
  if(min(yy) < 0){
    yy <- yy - min(yy)
  }

  #--------------------------------------------
  # Calculate the spectral slope.
  #--------------------------------------------
  n <- max(xx) - min(xx) - interval + 1

  res <- data.frame(wl = rep(NA, n), s = rep(NA, n), r2 = rep(NA, n))

  yy_log <- log(yy)
  yy_log[yy_log == -Inf] <- NA

  for(i in min(xx):(max(xx) - interval)){

    index <- xx >= i & xx <= (i + interval)

    fit <- lm(yy_log[index] ~ xx[index])

    mysummary <- summary(fit) ## Get the r2

    res$wl[i - min(xx) + 1] <- i + (interval / 2)
    res$s[i - min(xx) + 1] <- -coef(fit)[2]
    res$r2[i - min(xx) + 1] <- mysummary$r.squared
  }

  ## Filter data to keep only regression with r2 >= r2threshold
  res <- res[res$r2 >= r2threshold, ]
  res <- na.omit(res)

  return(res)

}
