#' @title Lengths and weights for Yellow Perch from Trout Lake, WI.
#' 
#' @description Lengths and weights for Yellow Perch \emph{Perca flavescens} from Trout Lake, WI, 1981-2006.  Fish were collected with a variety of gears.
#' 
#' @name YPerchTL
#' 
#' @docType data
#' 
#' @format A data frame with 7238 observations on the following 7 variables:
#'  \describe{
#'    \item{lakeid}{Lake name (all \code{TR}=Trout Lake)}
#'    \item{year4}{Year of capture} 
#'    \item{sampledate}{Date of capture}
#'    \item{gearid}{Capture gear type}
#'    \item{spname}{Species name (all \code{YELLOWPERCH})}
#'    \item{length}{Total Length (nearest mm) at capture}
#'    \item{weight}{Weight (nearest 0.1 or 1 g) at capture} 
#'  }
#' 
#' @section Topic(s):
#'  \itemize{
#'    \item Weight-Length
#'    \item Length Frequency
#'    \item Size Structure
#'    \item PSD
#'  }
#'  
#' @concept 'Weight-Length' 'Length Frequency' 'Size Structure' PSD
#' 
#' @source From a query to the North Temperate Lakes Long Term Ecological Research, Fish Lengths and Weights Database which was (is?) at https://lter.limnology.wisc.edu/dataset/north-temperate-lakes-lter-fish-lengths-and-weights-1981-current.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(YPerchTL)
#' str(YPerchTL)
#' head(YPerchTL)
#' plot(weight~length,data=YPerchTL)
#' 
NULL