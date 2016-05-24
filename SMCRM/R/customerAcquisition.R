#' Customer Acquisition Data from Chapter 3
#' @name customerAcquisition
#' @docType data
#' @usage customerAcquisition
#' @format Data frame with the following 17 variables
#' \describe{
#'   \item{\code{customer}}{customer number (from 1 to 500)}
#'   \item{\code{acquisition}}{1 if the prospect was acquired, 0 otherwise}
#'   \item{\code{first_purchase}}{dollar value of the first purchase (0 if the customer was not acquired)}
#'   \item{\code{clv}}{the predicted customer lifetime value score. It is
#'     0 if the prospect was not acquired or has already churned 
#'     from the firm.}
#'   \item{\code{duration}}{time in days that the acquired prospect has been 
#'     or was a customer, right-censored at 730 days}
#'   \item{\code{censor}}{1 if the customer was still a customer at the end
#'     of the observation window, 0 otherwise}
#'   \item{\code{acq_expense}}{dollars spent on marketing efforts to try and acquire that prospect}
#'   \item{\code{acq_expense_sq}}{square of dollars spent on marketing efforts to try and acquire that
#'     prospect}
#'   \item{\code{industry}}{1 if the customer is in the B2B industry, 0 otherwise}
#'   \item{\code{revenue}}{annual sales revenue of the prospect's firm (in millions of dollar)}
#'   \item{\code{employees}}{number of employees in the prospect's firm}
#'   \item{\code{ret_expense}}{dollars spent on marketing efforts to try and retain that customer}
#'   \item{\code{ret_expense_sq}}{square of dollars spent on marketing efforts to try and retain that customer}
#'   \item{\code{crossbuy}}{the number of categories the customer has purchased}
#'   \item{\code{frequency}}{the number of times the customer purchased during the observation window}
#'   \item{\code{frequency_sq}}{the square of the number of times the customer purchased during the observation window}
#' }
#' @examples data(customerAcquisition)
#'   str(customerAcquisition)
NULL
