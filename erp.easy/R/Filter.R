#' Filter ERP data using a Butterworth filter from the "signal" package
#'
#' \code{erp.filter} Using the "signal" package, erp.easy data frames can be
#'    filtered with a variety of Butterworth filters.
#'
#' @param data A data frame in the format returned from \code{\link{load.data}}
#' @param freq A value between 0 and 1 specifying the cutoff frequency (the erp.easy
#'    adapted function supports digital filters only, where 1 is the Nyquist
#'    frequency). For example, for a sampling rate of 1000 Hz, a cutoff value
#'    of 12 Hz would be designated as 0.012 (12 / 1000).  A cutoff value of 8
#'    Hz would be 0.008 (8 / 1000).
#' @param ftype The type of filter to be used: values include: high-pass
#'    ("high"), low-pass ("low"), stop-band ("stop"), and pass-band ("pass").
#'    The default is set to "low."
#' @param order The complexity, or "order" of the Butterworth filter.  The default
#'    value is a 4th order filter (24 dB/octave).
#'
#' @details See \href{https://cran.r-project.org/package=signal}{signal}
#'    package reference manual for detailed instructions for using the \code{butter}
#'    function.
#'
#' @return A data frame identical to the original data frame, but with filtered amplitudes.
#'
#' @examples
#' # Low pass filter ERP data (sampled at 250 Hz) with a cutoff frequency of 8 Hz
#' filtered <- erp.filter(ERPdata, freq = 0.032, ftype = "low", order = 4)
#' grandaverage(filtered, "V78", window = c(1000, 1500))
#'
#' @author Travis Moore

erp.filter <- function(data, freq, ftype = "low", order = 4) {
  butterworth <- signal::butter(n = order, W = freq, type = ftype)
  filtered1 <- apply(data[ , 4:ncol(data)], 2, function(x) signal::filtfilt(butterworth, x))
  filtered2 <- data.frame(filtered1)
  filtered3 <- cbind.data.frame(data[ , 1:3], filtered2)
return(filtered3)
}
