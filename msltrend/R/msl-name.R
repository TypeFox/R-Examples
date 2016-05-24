#' msltrend: A package providing improved techniques to estimate trend, velocity
#' and acceleration from sea level records.
#'
#' The 'msltrend' package provides improved estimates of trend (mean sea
#' level) and associated real-time velocities and accelerations from long (minimum
#' 80 years), individual, annual average ocean water level data records. Improved
#' trend estimates are based on Singular Spectrum Analysis methods. Various gap-filling
#' options are included to accommodate incomplete time series records. The package
#' also contains a forecasting module to consider the implication of user defined
#' quantum of sea level rise between the end of the available historical record
#' and the year 2100. A wide range of screen and pdf plotting options are available
#' within the package.
#'
#' @section msltrend functions:
#' The \code{\link{msl.trend}} function is the key entry point to the package
#' deconstructing annual average time series data records into a trend and
#' associated real-time velocities and accelerations, filling necessary internal
#' structures which facilitate all functions in this package (Refer Watson 2016a,b
#' for more detail).
#'
#' The \code{\link{msl.forecast}} function enables a user defined quantum of sea
#' level rise to be added from the end of the deconstructed historical record to
#' the year 2100. Similalrly, this function estimates real-time velocities and
#' accelerations from the start of the available historical record to the year 2100.
#'
#' @references Watson, P.J., 2016a. Identifying the best performing time series
#' analytics for sea-level research. In: \emph{Time Series Analysis and
#' Forecasting, Contributions to Statistics}, ISBN 978-3-319-28725-6, Springer
#' International Publishing (in press).
#'
#' Watson, P.J., 2016b. How to improve estimates of real-time acceleration in
#' the mean sea level signal. In: Vila-Concejo, A., Bruce, E., Kennedy, D.M.,
#' and McCarroll, R.J. (eds.), Proceedings of the 14th International Coastal
#' Symposium (Sydney, Australia). \emph{Journal of Coastal Research},
#' Special Issue, No. 75. Coconut Creek (Florida), ISSN 0749-0208 (in press).
#'
#' @docType package
#'
#' @name msltrend
NULL
