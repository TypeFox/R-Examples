#' Customer Win-Back from Chapter 7
#' @name customerWinBack
#' @docType data
#' @usage customerWinBack
#' @format Data frame with the following 10 variables
#' \describe{
#'   \item{\code{customer}}{customer number (from 1 to 500)}
#'   \item{\code{reacquire}}{1 if the customer is reacquired, 0 if not}
#'   \item{\code{duration_2}}{time in days of the customer's second lifecycle with the company, 0 if not reacquired}
#'   \item{\code{slcv}}{CLV of the customer in the second lifecycle}
#'   \item{\code{duration_1}}{time in days of the customer's first lifecycle with the company}
#'   \item{\code{offer}}{value of the offer provided to the customer for reacquisition}
#'   \item{\code{duration_lapse}}{time in days since the customer was lost to when the offer to reacquire was given}
#'   \item{\code{price_change}}{increase (or decrease) in price of the subscription the customer received between the first lifecycle and the second lifecycle, 0 if not reacquired}
#'   \item{\code{gender}}{1 if male, 0 if female}
#'   \item{\code{age}}{age in years of the customer at the time of the attempt to reacquire}
#' }
#' @examples data(customerWinBack)
#' str(customerWinBack)
NULL
