#' @title Catch-at-age of Flathead Catfish from three southeastern rivers.
#' 
#' @description Numbers of Flathead Catfish (\emph{Pylodictis olivaris}) captured by electrofishing in three rivers -- Coosa River, AL; Ocmulgee River, GA; and Satilla River, GA.
#' 
#' @name FHCatfish
#' 
#' @docType data
#' 
#' @format A data frame of 39 observations on the following 3 variables:
#'  \describe{
#'    \item{river}{Location of collection (\code{Coosa}, \code{Ocmulgee}, and \code{Satilla})}
#'    \item{age}{Age (years) assigned from otolith}
#'    \item{abundance}{Number of captured fish with boat electrofishing}
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
#' @source From (approximately) Figure 3 in Sakaris, P.C., E.R. Irwin, J.C. Jolley, and D. Harrison.  2006.  Comparison of Native and Introduced Flathead Catfish Populations in Alabama and Georgia: Growth, Mortality, and Management.  North American Journal of Fisheries Management 26:867-874.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(FHCatfish)
#' str(FHCatfish)
#' head(FHCatfish)
#' op <- par(mfrow=c(2,2),pch=19)
#' plot(log(abundance)~age,data=FHCatfish,subset=river=="Coosa",main="Coosa")
#' plot(log(abundance)~age,data=FHCatfish,subset=river=="Ocmulgee",main="Ocmulgee")
#' plot(log(abundance)~age,data=FHCatfish,subset=river=="Satilla",main="Satilla")
#' par(op)
#' 
NULL