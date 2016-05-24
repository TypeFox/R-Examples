#' @title Lengths and weights for fish captured in Inch Lake
#' 
#' @description Total lengths and weights for a subsample of fish captured in Inch Lake, WI in May, 2007 and May, 2008.
#' 
#' @name InchLake2
#' 
#' @docType data
#' 
#' @format A data frame of 516 observations on the following 6 variables:
#'  \describe{
#'    \item{netID}{A unique identifier for the sampling event}
#'    \item{fishID}{A unique identifier for the individual fish}
#'    \item{species}{Species name}
#'    \item{length}{Total length (inches to nearest 0.1)}
#'    \item{weight}{Wet weight (grams to nearest 0.1)} 
#'    \item{year}{Year of capture}
#'   }
#'   
#' @section Topic(s):
#'  \itemize{
#'    \item Weight-Length
#'    \item Condition
#'    \item Length Frequency
#'  }
#'  
#' @concept 'Weight-Length' Condition 'Length Frequency'
#' 
#' @seealso See \code{\link{InchLake1}} for the entire sample, but without weights.
#' 
#' @source Derek H. Ogle, personal collection
#' 
#' @keywords datasets
#' 
#' @examples
#' data(InchLake2)
#' str(InchLake2)
#' head(InchLake2)
#' 
#' ## Isolate just Bluegills
#' bg.il <- subset(InchLake2,species=="Bluegill")
#' 
#' ## Isolate just largemouth bass from 2007
#' lmb7.il <- subset(InchLake2,species=="Largemouth Bass" & year==2007)
#' 
NULL