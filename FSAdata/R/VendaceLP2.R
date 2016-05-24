#' @title Stock and recruitment data for Vendace from Lake Pyhajarvi.
#' 
#' @description Vendace (\emph{Coregonus albula}) recruitment in Lake Pyhajarvi.
#' 
#' @note Original authors fit an exponential curve to the fecundity-recruits relationship.
#' 
#' @name VendaceLP2
#' 
#' @docType data
#' 
#' @format A data frame of 9 observations on the following 2 variables:
#'  \describe{
#'    \item{fecundity}{Total fecundity (10^9 eggs) of spawning stock}
#'    \item{recruits}{Number of recruits (10^6 fish) in Autumn after hatching}
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Stock-Recruit
#'    \item Recruitment
#'  }
#' 
#' @concept 'Stock-Recruit' Recruitment
#' 
#' @source From (approximately) Figure 6 in Helminen, H. and J. Sarvala. 1994. Population regulation of vendance (\emph{Coregonus albula}) in Lake Pyhajarvi, southwest Finland.  Journal of Fish Biology 45:387-400.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(VendaceLP2)
#' str(VendaceLP2)
#' head(VendaceLP2)
#' plot(recruits~fecundity,data=VendaceLP2)
#' 
NULL