#' dydt.Linear returns the derivative of the linear function.
#' @title Return the derivative of the linear function
#' @author Marc Girondot
#' @return A list with the derivative
#' @param t The time in any unit
#' @param size The current size
#' @param parms A vector with alpha being c(alpha=x1, K=x2). Only alpha is used.
#' @description Return the derivative of the linear function\cr
#' dydt.linear(t, size, parms)
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' formated <- FormatNests(nest)
#' # The initial parameters value can be:
#' # "T12H", "DHA",  "DHH", "Rho25"
#' # Or
#' # "T12L", "DT", "DHA",  "DHH", "DHL", "Rho25"
#' x <- structure(c(306.174998729436, 333.708348843241, 299.856306141849,  
#' 	149.046870203155), .Names = c("DHA", "DHH", "T12H", "Rho25"))
#' # K or rK are not used for dydt.linear or dydt.exponential
#' resultNest_4p_linear <- searchR(parameters=x, fixed.parameters=NULL,  
#' 	temperatures=formated, derivate=dydt.linear, M0=1.7,  
#' 	test=c(Mean=39.33, SD=1.92))
#' }
#' @export

dydt.linear <- function(t, size, parms) {
		dy.dt <- parms["alpha"]
		list(dy.dt)
		}
