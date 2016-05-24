#' Temperature cycling dataset from Wang & Dillon NCC 2014. doi:10.1038/nclimate2378
#'
#' \pkg{tempcyclesdata} is dataset containing metadata, linear, and cycling data,
#' to be used with the \pkg{tempcycles} package. Only data passing all checks in included.
#'
#' @format A data frame with 77181 rows and 27 variables:
#' \describe{
#'   \item{id}{station id, USAF-WBAN}
#'   \item{name}{station name}
#'   \item{lat}{latitude, negative values indicate South}
#'   \item{lon}{longitude, negative values indicate West}
#'   \item{el}{elevation, in meters}
#'   \item{period}{time period. "all": all data for the station,
#'                              "stdr": standard reference period,
#'                              or middle of five year window.}
#'   \item{region}{geographical zone}
#'   \item{shore_dist_km}{Distance to shoreline (GSHHG 2)}
#'   \item{start_date}{start of data window}
#'   \item{end_date}{end of data window}
#'   \item{num_samp}{number of observations}
#'   \item{Ta_mean}{mean temperature, C}
#'   \item{Ta_min}{minimum temperature, C}
#'   \item{Ta_max}{maximum temperature, C}
#'   \item{Ta_var}{temperature variance}
#'   \item{Ta_slope}{linear slope of record}
#'   \item{Ta_int}{intercept of linear model, for detrending}
#'   \item{DTC}{Daily temperature cycling range, (2 * amplitude), C}
#'   \item{ATC}{Annual temperature cycling range, (2 * amplitude), C}
#'   \item{DTC_red}{Redfit AR1 corrected DTC, C}
#'   \item{ATC_red}{Redfit AR1 corrected ATC, C}
#'   \item{day_tau}{tau lag for day signal}
#'   \item{year_tau}{tau lag for year signal}
#'   \item{day_phase}{DTC phase}
#'   \item{year_phase}{ATC phase}
#'   \item{lnDA}{\emph{ln}DTC / ATC}
#'   \item{lnDA_red}{\emph{ln} DTC_red / ATC_red}
#'   \item{mean_resid}{mean per-sample residual, C}
#'   \item{mean_resid_red}{mean per-sample residual using redfit corrected values, C}
#' }
#' @source Wang & Dillon NCC 2014. doi:10.1038/nclimate2378
#' @examples
#' summary(tempcyclesdata)
#' if (require("dplyr")) {
#' tempcyclesdata %>%
#'   filter(period == "stdr") %>%
#'   group_by(region) %>%
#'   summarise(mean_DTC = mean(DTC))
#' }
"tempcyclesdata"
