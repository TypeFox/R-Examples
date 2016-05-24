#' @title Back-calculated lengths-at-age for Walleye from Lake Mille Lacs, 2000-2011.
#' 
#' @description Back-calculated lengths-at-age for Walleye (\emph{Sander vitreus})  from Lake Mille Lacs.  Walleye were captured by Minnesota Department of Natural Resources personnel in fishery-independent gillnets (five multifilament nylon mesh panels with each panel measuring 15.2 m wide and 1.8 m high; bar-measure mesh sizes of the panels were 19.1, 25.4, 31.7, 38.1, and 50.8 mm) set in the fall (mid September to early October) from 2000 to 2011.
#' 
#' @name WalleyeML
#' @docType data
#' @format A data frame of 14583 observations on the following 9 variables:
#' \describe{
#'  \item{ID}{A unique fish identification number.} 
#'  \item{Year}{Year of data.} 
#'  \item{Sex}{Sex (female, male).} 
#'  \item{Est.Age}{Estimated (from otoliths) age (yrs) at capture.}
#'  \item{TL}{Total length (mm).}
#'  \item{Scale.Rad}{Total scale radius (mm) at capture.} 
#'  \item{Dist.Ann}{Scale radisus (mm) to annulus given in \code{BC.Age}.} 
#'  \item{BC.Age}{Annulus or previous age.} 
#'  \item{BC.Len}{Back-calculated length at \code{BC.Age}.  Lengths were back-calculated using the \emph{Scale-Proportional Hypothesis} method.} 
#' }
#' @section Topic(s): \itemize{
#'  \item Growth
#'  \item von Bertalanffy
#'  \item Back-calculation
#' }
#' @concept 'Growth' 'von Bertalanffy' 'Back-calculation'
#' @source These unpublished data are from the Minnesota Department of Natural Resources, Section of Fisheries (via Melissa Treml).  \bold{Do not use for other than educational purposes without permission from the source.}
#' 
#' @keywords datasets
#' 
#' @examples
#' data(WalleyeML)
#' str(WalleyeML)
#' head(WalleyeML)
#' xtabs(~Year+Est.Age+Sex,data=WalleyeML)
#' 
NULL