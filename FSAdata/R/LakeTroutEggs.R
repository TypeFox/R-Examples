#' @title Length and egg deposition of Lake Superior Lake Trout.
#' 
#' @description Length and egg deposition of Lake Superior Lake Trout (\emph{Salvelinus namaycush}).
#' 
#' @name LakeTroutEggs
#' 
#' @docType data
#' 
#' @format A data frame of 101 observations on the following 2 variables:
#'  \describe{
#'    \item{tl}{Total length (mm) of an individual fish.}
#'    \item{eggs}{Estimated number of eggs.} 
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Other
#'  }
#'
#' @concept Other
#'   
#' @source From (approximately) Figure 2 of Schram, S.T. 1993.  Fecundity and egg deposition of a wild Lake Superior Lake Trout stock.  Wisconsin Department of Natural Resources, Fisheries Management Report no. 149.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(LakeTroutEggs)
#' str(LakeTroutEggs)
#' head(LakeTroutEggs)
#' plot(eggs~tl,data=LakeTroutEggs)
#' 
NULL