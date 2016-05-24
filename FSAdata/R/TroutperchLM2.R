#' @title Lengths for Troutperch from Lake Michigan, 1977.
#' 
#' @description Troutperch (\emph{Percopsis omiscomaycus}) fork lengths from near Grand Haven, Lake Michigan, 1977.
#' 
#' @name TroutperchLM2
#' 
#' @docType data
#' 
#' @format A data frame of 3385 observations on the following 1 variable:
#'  \describe{
#'    \item{fl}{fork length (mm)} 
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Length Frequency 
#'    \item Size Structure
#'  }
#'  
#' @concept 'Length Frequency' 'Size Structure'
#' 
#' @source From Brandt, S.B., J.J. Magnuson, and L.B. Crowder.  1980.  Thermal habitat partitioning by fishes in Lake Michigan.  Canadian Journal of Fisheries and Aquatic Sciences.  37:1557-1564.  Data was simulated (uniform distribution of values within length bin) from summarized length frequencies in http://fishbase.org/.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(TroutperchLM2)
#' str(TroutperchLM2)
#' head(TroutperchLM2)
#' hist(TroutperchLM2$fl,main="")
#' 
NULL