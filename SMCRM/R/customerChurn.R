#' Customer Churn Data from Chapter 6
#' @name customerChurn
#' @docType data
#' @usage customerChurn
#' @format Data frame with the following 11 variables
#' \describe{
#'   \item{\code{customer}}{customer number (from 1 to 500)}
#'   \item{\code{duration}}{time in days that the acquired prospect has been or was a customer, right-censored at 730 days}
#'   \item{\code{censor}}{1 if the customer was still a customer at the end of the observation window, 0 otherwise}
#'   \item{\code{avg_ret_exp}}{average number of dollars spent on marketing efforts to try and retain that customer per month}
#'   \item{\code{avg_ret_exp_sq}}{square of the average number of dollars spent on marketing efforts to try and retain that customer per month}
#'   \item{\code{total_crossbuy}}{total number of categories the customer has purchased during the customer's lifetime}
#'   \item{\code{total_freq}}{total number of purchase occasions the customer had with the firm in the customer's lifetime}
#'   \item{\code{total_freq_sq}}{square of the total number of purchase occasions the customer had with the firm in the customer's lifetime}
#'   \item{\code{industry}}{1 if the customer is in the B2B industry, 0 otherwise}
#'   \item{\code{revenue}}{annual sales revenue of the prospect's firm (in millions of dollar)}
#'   \item{\code{employees}}{number of employees in the prospect's firm}
#' }
#' @examples data(customerChurn)
#'   str(customerChurn)
NULL
