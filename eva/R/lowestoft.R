#' Top Ten Annual Sea Levels: Lowestoft, UK (1964 - 2014)
#'
#'
#' Top ten annual sea levels at the LoweStoft Station tide gauge from 1964 - 2014.
#' From 1964 - 1992, raw data is collected in hour intervals; from 1993 - present,
#' raw data is collected in fifteen minute intervals. Data is pre-processed here to
#' account for storm length - see reference for details.
#'
#'
#'
#' @docType data
#'
#' @usage data(lowestoft)
#' @name lowestoft
#'
#' @format A data matrix with 51 observations. Each year is considered an observation, with the top ten annual sea level events.
#'
#' @keywords datasets
#'
#' @references Bader B., Yan J., & Zhang X. (2015). Automated Selection of r for the r Largest Order Statistics Approach with Adjustment for Sequential Testing. Department of Statistics, University of Connecticut.
#'
#' @source UK Tide Gauge Network (Lowestoft Station): https://www.bodc.ac.uk/data/online_delivery/ntslf/processed/
#'
#' @examples
#' data(lowestoft)
#' gevrSeqTests(lowestoft, method = "ed")
#' ## Not run
#' ## Look at the difference in confidence intervals between r = 1 and r = 10
#' # z1 <- gevrFit(lowestoft[, 1])
#' # z2 <- gevrFit(lowestoft)
#' # gevrRl(z1, 50, method = "profile")
#' # gevrRl(z2, 50, method = "profile")
NULL
