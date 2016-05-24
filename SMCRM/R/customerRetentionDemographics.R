#' Demographics Data for Customer Retention (Chapter 4)
#' @name customerRetentionDemographics
#' @docType data
#' @usage customerRetentionDemographics
#' @format Data frame with the following 8 variables
#' \describe{
#'   \item{\code{customer}}{customer number (from 1 to 500)}
#'   \item{\code{gender}}{1 if the customer is male, 0 if the customer is female}
#'   \item{\code{married}}{1 if the customer is married, 0 if the customer is not married}
#'   \item{\code{income}}{1 if income < \$30,000
#'                 2 if \$30,001 < income < \$45,000
#'                 3 if \$45,001 < income < \$60,000
#'                 4 if \$60,001 < income < \$75,000
#'                 5 if \$75,001 < income < \$90,000
#'                 6 if income > \$90,001}
#'   \item{\code{first_purchase}}{value of the first purchase made by the customer in quarter 1}
#'   \item{\code{loyalty}}{1 if the customer is a member of the loyalty program, 0 if not}
#'   \item{\code{sow}}{share-of-wallet; the percentage of purchases the customer makes from the
#'              given firm given the total amount of purchases across all firms in that
#'              category}
#'   \item{\code{clv}}{discounted value of all expected future profits, or customer lifetime value}
#' }
#' @examples data(customerRetentionDemographics)
#'   str(customerRetentionDemographics)
NULL
