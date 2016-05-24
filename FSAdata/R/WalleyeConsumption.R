#' @title Consumption of prey by Walleye.
#' 
#' @description Consumption of prey by Walleye (\emph{Sander vitreus}) at different prey densities.
#' 
#' @name WalleyeConsumption
#' 
#' @docType data
#' 
#' @format A data frame of 23 observations on the following 2 variables:
#'  \describe{
#'    \item{PreyDensity}{Density of prey (mg per g per day).}
#'    \item{FoodConsump}{Food consumption by predator (mg per cubic meter)} 
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Nonlinear modeling 
#'  }
#'  
#' @source From Figure 3 in Madenjian, C.P., and S.R. Carpenter. 1991.  Individual-based model for growth of young-of-the-year walleye: A piece of the recruitment puzzle.  Ecological Applications.  1:268-279.  Data were originally from Swenson, W. A. 1977. Food consumptions of walleye (\emph{Stizostedion vitreum vitreum}) and sauger (\emph{S. canadense}) in relation to food availability and physical conditions in Lake of the Woods, Minnesota, Shagawa Lake, and western Lake Superior.  Journal of the Fisheries Research Board of Canada 34:1643-1654.  [Madenjian et al. (1991) was (is?) from http://www.esajournals.org/doi/abs/10.2307/1941756.  Swenson (1977) was (is?) from http://www.nrcresearchpress.com/doi/abs/10.1139/f77-229?journalCode=jfrbc.]
#' 
#' @keywords datasets
#' 
#' @examples
#' data(WalleyeConsumption)
#' str(WalleyeConsumption)
#' head(WalleyeConsumption)
#' plot(FoodConsump~PreyDensity,data=WalleyeConsumption,pch=16)
#' 
NULL