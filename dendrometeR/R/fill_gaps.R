#' Fill gaps in dendrometer series
#'
#' @description The function fills gaps in a \code{data.frame} with dendrometer series using an ARMA model (cf. Deslauriers et al. 2011), and is designed for single growing seasons. The function is able to fill gaps of short duration (i.e. several hours), but cannot sensibly handle long gaps.
#'
#' @usage fill_gaps(dm.data, Hz = 0.01, season = FALSE)
#'
#' @param dm.data a \code{data.frame} with a timestamp (\code{\%Y-\%m-\%d \%H:\%M:\%S} format) as row names, and dendrometer series in columns. Output as created using code from the \code{Import dendrometer data} vignette.
#' @param Hz a \code{numeric} specifying the parameter for smoothing with ARMA gap-filling. A higher value means rougher smoothing. Defaults to 0.01.
#' @param season a \code{logical} indicating whether \code{\link{auto.arima}} should check seasonal models; can be very slow. Defaults to FALSE, i.e. search restricted to non-seasonal models.
#'
#' @details The function uses \code{\link{auto.arima}} to fill missing records. The non-seasonal part of the model is specified by the three integer components: the AR order \emph{p}, the degree of differencing \emph{d}, and the MA order \emph{q}. For the seasonal part of the model, the period parameter is set equal to the number of daily measurements observed in the dendrometer data. The output of the ARMA model is smoothed using \code{\link{smooth.Pspline}}. The smoothing parameter Hz can be adjusted; defaults to 0.01.
#'
#' The function is designed for single growing seasons, amongst others because ARMA-based gap-filling routines will then perform best (i.e. ARMA parameters might be distinct for individual growing seasons). To allow the usage of \code{\link{fill_gaps}} for datasets from the Southern Hemisphere, the input data may contain two consecutive calendar years.
#'
#' @return The function returns a \code{data.frame} with gap-filled dendrometer series.
#'
#' @author Olivier Bouriaud, Ernst van der Maaten, Marieke van der Maaten-Theunissen and Marko Smiljanic.
#'
#' @references Deslauriers, A., Rossi, S., Turcotte, A., Morin, H. and Krause, C. (2011) A three-step procedure in SAS to analyze the time series from automatic dendrometers. \emph{Dendrochronologia} 29: 151-161.
#'
#' @examples
#' \dontrun{
#'
#' data(dmCD)
#' ## creating some artificial gaps (for demonstration purposes):
#' dmCD[c(873:877,985:990),1] <- NA
#' # slow, as also seasonal models are checked, but best possible gap-filling:
#' dm.gpf <- fill_gaps(dmCD, Hz = 0.01, season = TRUE)
#' }
#'
#' @import zoo
#' @import pspline
#' @import stats
#' @import forecast
#'
#' @export fill_gaps
#'
fill_gaps <- function(dm.data, Hz = 0.01, season = FALSE)
{
  nm <- deparse(substitute(dm.data))

  if(!is.dendro(dm.data)) {
    stop(paste("'", nm, "' is not in the required format", sep = ""))
  }
  if(!any(is.na(dm.data))) {
    stop(paste("'", nm, "' does not contain gaps, so no need for gap-filling", sep = ""))
  }

  gpf <- data.frame(date = rownames(dm.data))
  year <- as.POSIXlt(gpf$date, format = "%Y-%m-%d %H:%M:%S")$year
  if(length(unique(year)) > 2) {
    stop(paste("'fill.gaps' is designed to fill gaps for single growing seasons. '", nm,"' contains more than two years", sep = ""))
  }

  resolution <- dendro.resolution(dm.data, "hours")
  period.var <- 24 / resolution

  input <- dm.data

  for(i in 1:ncol(input)) {
    indices <- as.POSIXct(rownames(input), tz = "GMT")
    sensor <- zoo(input[,i], indices)
    smoothint <- na.spline(sensor, na.rm = FALSE)

    name.smoothint <- paste(names(input)[i], "smoothint", sep = ".")
    name.armafit <- paste(names(input)[i], "armafit", sep = ".")
    name.corrected <- paste(names(input)[i], "corrected", sep = ".")
    name.dendro2 <- paste(names(input)[i], "dendro2", sep = ".")

    input[name.smoothint] <- ts(smoothint, frequency = period.var)
    arma.model <- auto.arima(input[[name.smoothint]], stepwise = TRUE, seasonal = season)

    arma.var <- arma.model$residuals
    input[[name.armafit]] <- input[[name.smoothint]] + arma.var
    input[[name.corrected]] <- sensor
    input[[name.corrected]][is.na(sensor) == TRUE] <- input[[name.armafit]][is.na(sensor) == TRUE]
    input[[name.dendro2]] <- smooth.Pspline(seq(1:length(input[[name.corrected]])), input[[name.corrected]], spar = Hz)$ysmth
    gpf <- cbind(gpf, input[[name.dendro2]])
    names(gpf)[i+1] <- colnames(input[i])
  }

  rownames(gpf) <- gpf$date
  gpf$date <- NULL

  dmGpf <- dm.data
  for(i in 1:ncol(dm.data)) {
    tst.na <- which(is.na(dm.data[,i]) )
    dmGpf[tst.na,i] <- gpf[tst.na,i]
  }
  return(dmGpf)
}
