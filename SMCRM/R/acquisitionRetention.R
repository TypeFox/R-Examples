#' Acquisition-Retention Data from Chapter 5
#' @name acquisitionRetention
#' @docType data
#' @usage acquisitionRetention
#' @format Data frame with the following 15 variables
#' \describe{
#'   \item{\code{customer}}{customer number (from 1 to 500)}
#'   \item{\code{acquisition}}{1 if the prospect was acquired, 0 otherwise}
#'   \item{\code{duration}}{number of days the customer was a customer of the firm, 0 if acquisition == 0}
#'   \item{\code{profit}}{customer lifetime value (CLV) of a given customer, -(Acq_Exp) if the customer is not acquired}
#'   \item{\code{acq_exp}}{total dollars spent on trying to acquire this prospect}
#'   \item{\code{ret_exp}}{total dollars spent on trying to retain this customer}
#'   \item{\code{acq_exp_sq}}{square of the total dollars spent on trying to acquire this prospect}
#'   \item{\code{ret_exp_sq}}{square of the total dollars spent on trying to retain this customer}
#'   \item{\code{freq}}{number of purchases the customer made during that customer's lifetime with the firm, 0 if acquisition == 0}
#'   \item{\code{freq_sq}}{square of the number of purchases the customer made during that customer's lifetime with the firm}
#'   \item{\code{crossbuy}}{number of product categories the customer purchased from during that customer's lifetime with the firm, 0 if acquisition = 0}
#'   \item{\code{sow}}{Share-of-Wallet; percentage of purchases the customer makes from the given firm given the total amount of purchases across all firms in that category}
#'   \item{\code{industry}}{1 if the customer is in the B2B industry, 0 otherwise}
#'   \item{\code{revenue}}{annual sales revenue of the prospect's firm (in millions of dollar)}
#'   \item{\code{employees}}{number of employees in the prospect's firm}
#' }
#' @examples data(acquisitionRetention)
#'   str(acquisitionRetention)
NULL
