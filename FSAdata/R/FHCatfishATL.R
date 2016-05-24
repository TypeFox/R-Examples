#' @title Catch-at-age of Flathead Catfish from three Atlantic rivers.
#' 
#' @description Catch-at-age of Flathead Catfish (\emph{Pylodictis olivaris}) from three populations of Atlantic rivers -- Lumber River, Northeast Cape Fear River (NCF), and Neuse River.
#' 
#' @name FHCatfishATL
#' 
#' @docType data
#' 
#' @format A data frame of 44 observations on the following 3 variables:
#'  \describe{
#'    \item{river}{Collection river (\code{Lumber}, \code{NCF}, and \code{Neuse}).} 
#'    \item{age}{Age (yrs) assessed by otolith.}
#'    \item{number}{Number of captured fish.} 
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Mortality 
#'    \item Catch curve
#'  }
#'  
#' @concept Mortality 'Catch Curve'
#' 
#' @source From (approximately) Figure 2 in Kwak, T.J., W.E. Pine III, and D.S. Waters.  2006.  Age, growth, and mortality of introduced flathead catfish in Atlantic rivers and a review of other populations.  North American Journal of Fisheries Management 26:73-87.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(FHCatfishATL)
#' str(FHCatfishATL)
#' head(FHCatfishATL)
#' op <- par(mfrow=c(2,2),pch=19)
#' plot(log(number)~age,data=FHCatfishATL,subset=river=="Lumber",main="Lumber")
#' plot(log(number)~age,data=FHCatfishATL,subset=river=="NCF",main="NCF")
#' plot(log(number)~age,data=FHCatfishATL,subset=river=="Neuse",main="Neuse")
#' par(op)
#' 
NULL