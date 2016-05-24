#' Lifetime Duration Data for Customer Retention (Chapter 4)
#' @name customerRetentionLifetimeDuration
#' @docType data
#' @usage customerRetentionLifetimeDuration
#' @format Data frame with the following 8 variables
#' \describe{
#'   \item{\code{customer}}{customer number (from 1 to 500)}
#'   \item{\code{x}}{The number of transactions by a given customer over all time periods. Here we assume that it is the sum of the variable Purchase where customers at most made 1 purchase per quarter.}
#'   \item{\code{tx}}{time of the last transaction, i.e. the last quarter where purchase == 1}
#'   \item{\code{T}}{total time between the first purchase and the end of the observation window, i.e. 12 quarters for all customers}
#' }
#' @seealso customerRetentionTransactions
#' @examples data(customerRetentionLifetimeDuration)
#'   str(customerRetentionLifetimeDuration)
NULL
