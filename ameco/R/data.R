#' European Commission Annual macro-economic database (AMECO)
#'
#' The dataset contains the tntire annual macrro-economic database provided by
#' the European Commission. Last update: 4 February 2015.
#'
#' @format A data frame with eight variables:
#' \describe{
#' \item{\code{code}}{Code that uniquely identifies a series.}
#' \item{\code{country}}{Country.}
#' \item{\code{sub.chapter}}{Groups series into larger categories to facilitate
#' finding a series of interest.}
#' \item{\code{title}}{Human-readable title of a series.}
#' \item{\code{unit}}{Unit of observation.}
#' \item{\code{cntry}}{Country ISO code.}
#' \item{\code{year}}{Observation year.}
#' \item{\code{value}}{Value of the given series.}
#' }
#'
#' For further details, see
#' \url{http://ec.europa.eu/economy_finance/db_indicators/ameco/index_en.htm}
#'
"ameco"