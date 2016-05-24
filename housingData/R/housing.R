#' Housing data
#'
#' Monthly median home listing, sale price per square foot, and number of units sold for 2984 counties in the contiguous United States From 2008 to January 2016, harvested from Quandl's Zillow Housing Data \url{https://www.quandl.com/collections/housing}.  Provided from Quandl with rights to distribute without restriction.
#'
#' List price data is the most complete.  Many counties are missing sale price and number of units sold.  Only counties with at least one non-NA listing price entry are included in the data.  Also, the data only represents counties in the contiguous United States (Alaska, Hawaii, Puerto Rico, and Virgin Islands are not included).  A few major outliers in sale prices in small counties have been removed.
#'
#' @section Variables:
#'
#' \itemize{
#'  \item \code{fips}: FIPS county code
#'  \item \code{county}: county name
#'  \item \code{state}: state abbreviation
#'  \item \code{time}: time of class "Date" - each month's values are placed at the first day of the month
#'  \item \code{nSold}: number of units sold
#'  \item \code{medListPriceSqft}: median list price per square foot
#'  \item \code{medSoldPriceSqft}: median sold price per square foot
#' }
#' @docType data
#' @name housing
#' @usage housing
#' @format A data frame with 247082 rows and 7 columns.
#' @examples
#' head(housing)
NULL

