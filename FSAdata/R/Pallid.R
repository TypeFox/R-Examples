#' @title Lengths and weights for Pallid Sturgeon from four locations in the Missouri River.
#' 
#' @description Lengths (standard, fork, and total) and weights for Pallid Sturgion (\emph{Scaphirhynchus albus}) collected at four locations in the Missouri River drainage.
#' 
#' @name Pallid
#' 
#' @docType data
#' 
#' @format A data frame with 30 observations on the following 7 variables:
#'  \describe{
#'    \item{date}{Date of collection} 
#'    \item{sl}{Standard length (mm)}
#'    \item{fl}{Fork length (mm)} 
#'    \item{tl}{Total length (mm)} 
#'    \item{w}{Weight (g)}
#'    \item{status}{Living status of fish at time of collection (\code{Frozen}, \code{Live}, \code{Dead}).} 
#'    \item{loc}{Location of fish collection (\code{NB}=Nebraska, \code{SD}=South Dakota, \code{ND}=North Dakota, \code{MT}=Montana)} 
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Weight-Length 
#'    \item Length Conversion 
#'  }
#'  
#' @concept 'Weight-Length' 'Length Conversion'
#' 
#' @source From Table 1 in Keenlyne, K.D. and S.J. Maxwell.  1993.  Length conversions and length-weight relations for pallid sturgeon.  North American Journal of Fisheries Management. 13:395-397.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(Pallid)
#' str(Pallid)
#' head(Pallid)
#' op <- par(mfrow=c(2,2),pch=19)
#' ## four (of many possible) examples
#' plot(w~tl,data=Pallid,subset=loc=="NB",main="Nebraska")
#' plot(w~tl,data=Pallid,subset=loc=="SD",main="South Dakota")
#' plot(w~tl,data=Pallid,subset=loc=="ND",main="North Dakota")
#' plot(w~tl,data=Pallid,subset=loc=="MT",main="Montana")
#' par(op)
#' 
NULL