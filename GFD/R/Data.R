#' Startup Costs of five different companies
#'
#' A dataset containing the startup costs (in thousands of dollars) of five companies.
#'
#' @format A data frame with 60 rows and 2 variables:
#' \describe{
#'   \item{Costs}{price, in thousands of dollars}
#'   \item{company}{company, a factor with levels "pets", "pizza", "gifts", "shoes" and "bakery"}
#' }
#' 
#' @usage data(startup)
#' 
#' @source \url{http://college.cengage.com/mathematics/brase/understandable_statistics/7e/students/datasets/owan/frames/frame.html}
"startup"


#' Pizza delivery times
#'
#' A dataset containing the delivery times for pizza (in minutes) under different conditions.
#'
#' @format A data frame with 16 rows and 6 variables:
#' \describe{
#'   \item{Crust}{a factor with levels "thick" and "thin"}
#'   \item{Coke}{whether or not Coke was ordered with the pizza ("yes" or "no")}
#'   \item{Bread}{whether or not garlic bread was ordered with the pizza ("yes" or "no")}
#'   \item{Driver}{the sex of the driver, a factor with levels "M" and "F"}
#'   \item{Hour}{time of order in hours after midnight}
#'   \item{Delivery}{Delivery time in minutes}
#' }
#' 
#' @usage data(pizza)
#' 
#' @source \url{http://www.statsci.org/data/oz/pizza.html}
"pizza"


#' Curdies river data set
#'
#' A dataset containing the number of flatworms (dugesia) sampled in two seasons at different sites in the Curdies River in Western Victoria.
#'
#' @format A data frame with 36 rows and 3 variables:
#' \describe{
#'   \item{season}{a factor with levels "SUMMER" and "WINTER"}
#'   \item{site}{a factor with levels 1 to 6, nested within "season"}
#'   \item{dugesia}{number of flatworms counted on a particular stone (in no./dm^2)}
#' }
#' 
#' @usage data(curdies)
#' 
#' @source \url{http://users.monash.edu.au/~murray/AIMS-R-users/ws/ws7.html}
"curdies"