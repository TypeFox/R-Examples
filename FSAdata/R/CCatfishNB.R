#' @title Catch-at-age of Channel Catfish from two sections of the Platte River, NB.
#' 
#' @description Catch-at-age of Channel Catfish (\emph{Ictalurs puncatatus}) from two sections of the Platte River, NB, in 2007 and 2008.
#' 
#' @name CCatfishNB
#' 
#' @docType data
#' 
#' @format A data frame of 26 observations on the following 3 variables:
#'  \describe{
#'    \item{age}{Age (years) assigned from pectoral spines}
#'    \item{catch}{Number of captured fish with baited hoopnets and electrofishing}
#'    \item{loc}{Location of collection (\code{Central} and \code{Lower})}
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
#' @note Used in the \href{http://derekogle.com/IFAR/}{Introductory Fisheries Analyses with R} book.
#' 
#' @source From (approximately) Figure 3-14 in Barada, T.J. 2009. Catfish population dynamics in the Platte River, Nebraska.  Master's thesis, University of Nebraska, Lincoln, NE.  [Was (is?) from http://nlc1.nlc.state.ne.us/epubs/U1500/B013-2009.pdf.]  
#' 
#' @keywords datasets
#' 
#' @examples
#' data(CCatfishNB)
#' str(CCatfishNB)
#' head(CCatfishNB)
#' op <- par(mfrow=c(1,2),pch=19)
#' plot(log(catch)~age,data=CCatfishNB,subset=loc=="Central",main="Central")
#' plot(log(catch)~age,data=CCatfishNB,subset=loc=="Lower",main="Lower")
#' par(op)
#' 
NULL
