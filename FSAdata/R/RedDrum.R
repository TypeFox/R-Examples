#' @title Ages and lengths for Red Drum from the Atlantic Coast.
#' 
#' @description Assigned ages (from otoliths) and fork lengths of Red Drum (\emph{Sciaenops ocellatus}) from various areas of the Atlantic Coast, 1981-1988.
#' 
#' @name RedDrum
#' 
#' @docType data
#' 
#' @format A data frame with 393 observations on the following 2 variables.
#'  \describe{
#'    \item{age}{Age (from otoliths to the nearest years but recorded at half-years)} 
#'    \item{fl}{Fork length (mm)} 
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
#' @source From (approximately) Figure 27 in Vaughan, D.S. and T.E. Helser.  1990.  Status of the red drum stock of the Atlantic Coast: Stock assessment report for 1989.  NOAA Technical Memorandum, NMFS-SEFC-263.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(RedDrum)
#' str(RedDrum)
#' head(RedDrum)
#' plot(fl~age,data=RedDrum)
#' 
NULL