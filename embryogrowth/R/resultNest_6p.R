#' Result of the fit using the nest database
#' @title Fit using the nest database
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @docType data
#' @name resultNest_6p
#' @description Fit using the nest database
#' @references Girondot, M., & Kaska, Y. (2014). A model to predict 
#'             the thermal reaction norm for the embryo growth rate 
#'             from field data. Journal of Thermal Biology, 45, 96-102. 
#'             doi: 10.1016/j.jtherbio.2014.08.005
#' @keywords datasets
#' @usage resultNest_6p
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' formated <- FormatNests(nest)
#' # The initial parameters value can be:
#' # "T12H", "DHA",  "DHH", "Rho25"
#' # Or
#' # "T12L", "DT", "DHA",  "DHH", "DHL", "Rho25"
#' x <- structure(c(115.770032186805, 428.649144244891, 503.716700005443, 
#' 12.2641939339981, 306.305498098871, 116.380086027401), 
#' .Names = c("DHA", "DHH", "DHL", "DT", "T12L", "Rho25"))
#' pfixed <- c(rK=2.093313)
#' resultNest_6p <- searchR(parameters=x, fixed.parameters=pfixed, 
#' 	temperatures=formated, derivate=dydt.Gompertz, M0=1.7, 
#' 	test=c(Mean=39.33, SD=1.92))
#' }
#' @format A list with fitted information about data(nest)
NULL
