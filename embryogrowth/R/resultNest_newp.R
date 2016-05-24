#' Result of the fit using the nest database with anchored parameters
#' @title Fit using the nest database with anchored parameters
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @docType data
#' @name resultNest_newp
#' @description Fit using the nest database with anchored parameters
#' @references Girondot, M., & Kaska, Y. (2014). A model to predict 
#'             the thermal reaction norm for the embryo growth rate 
#'             from field data. Journal of Thermal Biology, 45, 96-102. 
#'             doi: 10.1016/j.jtherbio.2014.08.005
#' @keywords datasets
#' @usage resultNest_newp
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' formated <- FormatNests(nest)
#' newp <- GenerateAnchor(nests=formated, number.anchors=7)
#' pfixed <- c(rK=2.093313)
#' resultNest_newp <- searchR(parameters=newp, fixed.parameters=pfixed, 
#'   temperatures=formated, derivate=dydt.Gompertz, M0=1.7, 
#' 	test=c(Mean=39.33, SD=1.92))
#' data(resultNest_newp)
#' plotR(resultNest_newp)
#' }
#' @format A list with fitted information from data(nest) with anchored parameters
NULL
