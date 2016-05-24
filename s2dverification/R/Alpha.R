Alpha <- function(xdata, detrend = FALSE, filter = FALSE) {
  tmp <- xdata

  if (detrend == TRUE) {
    reg <- lm(xdata[is.na(xdata) == FALSE] ~ c(1:length(xdata[is.na(xdata) == FALSE])))
    fitted <- array(0, dim = length(xdata[is.na(xdata) == FALSE]))
    if (confint(reg)[2, 1] > 0 & confint(reg)[2, 2] > 0) {
      fitted <- c(1:length(xdata[is.na(xdata) == FALSE])) * min(confint(reg)[2, 2],
                  confint(reg)[2, 1])
    }
    if (confint(reg)[2, 1] < 0 & confint(reg)[2, 2] < 0) {
      fitted <- c(1:length(xdata[is.na(xdata) == FALSE])) * max(confint(reg)[2, 2],
                  confint(reg)[2, 1])
    }
    tmp[is.na(xdata) == FALSE] <- xdata[is.na(xdata) == FALSE] - fitted
  }

  if (filter == TRUE) {
    spec_z <- Spectrum(xdata)
    for (jlen in 1:dim(spec_z)[1]) {
      if (spec_z[jlen, 2] > spec_z[jlen, 4]) {
        tmp <- Filter(tmp, spec_z[jlen, 1])
      }
    }
  }

  estacf <- acf(tmp[is.na(xdata) == FALSE], plot = FALSE)$acf
  ##nacf <- length(estacf)
  ##estacf <- estacf[1:min(3, nacf)]
  ##alpha <- FitAutocor(estacf, c(0, 1), 0.001)
  alpha <- FitAcfCoef(max(estacf[2], 0), max(estacf[3], 0))

  alpha
}
