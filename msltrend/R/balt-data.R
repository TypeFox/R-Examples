#' Ocean water level data for Baltimore, USA
#'
#' Annual average ocean water level data from Permanent Service for Mean Sea
#' Level (UK).
#'
#' @docType data
#'
#' @usage data(Balt)
#'
#' @format Time series data file with the first column the year and the second
#' column the corresponding annual average ocean water level (mm). File contains
#' 112 records spanning the period from 1904 to 2014 with a single missing value
#' in 1990.
#'
#' @details The raw (*.csv) form of this data set is used extensively in the
#' examples throughout this manual.
#'
#' @references Holgate, S.J., Matthews, A., Woodworth, P.L., Rickards, L.J.,
#' Tamisiea, M.E., Bradshaw, E., Foden, P.R., Gordon, K.M., Jevrejeva, S. and
#' Pugh, J., 2013. New data systems and products at the Permanent Service for
#' Mean Sea Level. \emph{Journal of Coastal Research}, 29(3), pp. 493-504.
#'
#' @source \href{http://www.psmsl.org/data/obtaining/map.html}{Permanent Service
#' for Mean Sea Level (2015)}
#'
#' @seealso \code{\link{msl.trend}}, \code{\link{msl.forecast}},
#' \code{\link{msl.plot}}, \code{\link{msl.pdf}}, \code{\link{summary}}.
#'
#' @examples
#' data(Balt)
#' plot(Balt, type = "l", xlab = "Year", ylab = "Annual Average Mean Sea Level (mm)",
#' main = 'BALTIMORE, USA')
#' str(Balt) # check structure of data file
"Balt"
