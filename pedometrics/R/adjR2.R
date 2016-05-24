#' Adjusted coefficient of determination
#' 
#' Calculates the adjusted coefficient of determination of a multiple linear 
#' regression model.
#' 
#' @param r2 Numeric vector with the coefficient of determination to be 
#' adjusted.
#' @param n Numeric vector providing the number of observations used to fit the
#' multiple linear regression model.
#' @param p Numeric vector providing the number of parameters included in the
#' multiple linear regression model.
#' 
#' @details Details will be added later.
#' 
#' @return A numeric vector with the adjusted coefficient of determination.
#' 
#' @author Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' 
#' @references Coefficient of determination. Wikipedia, The Free Encyclopedia.
#' Available at \url{http://en.wikipedia.org/wiki/Coefficient_of_determination}.
#' [Online; accessed 31-July-2014].
#' 
#' @keywords misc
#' @export
#' @examples
#' adjR2(r2 = 0.95, n = 100, p = 80)
# FUNCTION #####################################################################
adjR2 <- 
  function (r2, n, p) {
    r2 <- 1 - (1 - r2) * ((n - 1) / (n - p - 1))
    return(r2)
  }
