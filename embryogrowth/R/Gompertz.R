#' dydt.Gompertz returns the derivative of the Gompertz function.
#' @title Return the derivative of the Gompertz function
#' @author Marc Girondot
#' @return A list with the derivative
#' @param t The time in any unit
#' @param size The current size
#' @param parms A vector with alpha and K values being c(alpha=x1, K=x2)
#' @description Return the derivative of the Gompertz function\cr
#' dydt.Gompertz(t, size, parms)
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' formated <- FormatNests(nest)
#' # The initial parameters value can be:
#' # "T12H", "DHA",  "DHH", "Rho25"
#' # Or
#' # "T12L", "DT", "DHA",  "DHH", "DHL", "Rho25"
#' x <- structure(c(118.768297442004, 475.750095909406, 306.243694918151, 
#' 116.055824800264), .Names = c("DHA", "DHH", "T12H", "Rho25"))
#' # pfixed <- c(K=82.33) or rK=82.33/39.33
#' pfixed <- c(rK=2.093313)
#' # K or rK are not used for dydt.linear or dydt.exponential
#' resultNest_4p <- searchR(parameters=x, fixed.parameters=pfixed,  
#' 	temperatures=formated, derivate=dydt.Gompertz, M0=1.7,  
#' 	test=c(Mean=39.33, SD=1.92))
#' data(resultNest_4p)
#' }
#' @export

dydt.Gompertz <- function(t, size, parms) {
		dy.dt <- parms["alpha"]*log(parms["K"]/size)*size
		list(dy.dt)
		}
