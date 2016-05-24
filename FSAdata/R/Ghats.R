#' @title Species accumulation data for fish of the Western Ghats of India.
#' 
#' @description Species accumulation data for fish of the Western Ghats of India derived from nine random samples of publications.
#' 
#' @name Ghats
#' 
#' @docType data
#' 
#' @format A data frame with 350 observations on the following 2 variables.
#'  \describe{
#'    \item{unit}{a manuscript that was reviewed.}
#'    \item{cumspec}{cumulative number of species described in the reviewed manuscripts.} 
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Other
#'  }
#'
#' @concept Other
#'   
#' @source From (approximately) Figure 1 in Dahanukar, N., R. Raut,  and A. Bhat.  2004.  Distribution, endemism and threat status of freshwater fishes in the Western Ghats of India.  Journal of Biogeography 31:123-126.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(Ghats)
#' str(Ghats)
#' head(Ghats)
#' plot(cumspec~unit,data=Ghats)
#' 
NULL