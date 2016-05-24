#' @title Ages and lengths of Lake Ontario Rock Bass.
#' 
#' @description Assigned ages (from scales) and measured total lengths for each of 1288 Rock Bass (\emph{Ambloplites rupestris}) from Lake Ontario.
#' 
#' @name RockBassLO1
#' 
#' @docType data
#' 
#' @format A data frame with 1288 observations on the following 2 variables:
#'  \describe{
#'    \item{age}{Assigned ages (from scales).} 
#'    \item{tl}{Measured total lengths (mm).} 
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
#' @seealso \code{\link{RockBassLO2}}.
#' 
#' @source Simulated from Table 1 of Wolfert, D.R.  1980.  Age and growth of rock bass in Eastern Lake Ontario.  New York Fish and Game Journal, 27:88:90.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(RockBassLO1)
#' str(RockBassLO1)
#' head(RockBassLO1)
#' plot(tl~age,data=RockBassLO1)
#' 
NULL