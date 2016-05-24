#' @title Ages, lengths, and sexes of Atlantic Croaker by sex.
#' 
#' @description Assigned ages (by otoliths), total lengths, and sexes of Atlantic Croaker (\emph{Micropogonias undulatus}).
#' 
#' @name Croaker2
#' 
#' @docType data
#' 
#' @format A data frame of 318 observations on the following 3 variables:
#'  \describe{
#'    \item{age}{Otolith age-at-capture (years).}
#'    \item{tl}{Total length (nearest mm) at capture.} 
#'    \item{sex}{Sex of the fish (\code{M}=male and \code{F}=female).} 
#'  }
#' 
#' @section Topic(s):
#'  \itemize{
#'    \item Growth
#'    \item von Bertalanffy
#'  }
#'  
#' @concept Growth 'von Bertalanffy'
#' 
#' @source From Figure 4 in Chapter 1 (Atlantic Croaker) of the VMRC Final Report on Finfish Ageing, 1999 by the Center for Quantitative Fisheries Ecology at Old Dominion University.
#'
#' @keywords datasets
#' 
#' @examples
#' data(Croaker2)
#' str(Croaker2)
#' head(Croaker2)
#' op <- par(mfrow=c(1,2),pch=19)
#' plot(tl~age,data=Croaker2,subset=sex=="F",main="Female")
#' plot(tl~age,data=Croaker2,subset=sex=="M",main="Male")
#' par(op)
#' 
NULL