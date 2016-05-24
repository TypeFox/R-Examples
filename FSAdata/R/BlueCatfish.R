#' @title Ages and lengths of Blue Catfish.
#' 
#' @description Ages and total lengths of Blue Catfish (\emph{Ictalurus furcatus}) collected form the Wilson Reservoir on the Tennessee River, AL.
#' 
#' @name BlueCatfish
#' 
#' @docType data
#' 
#' @format A data frame with 119 observations on the following 2 variables.
#'  \describe{
#'    \item{age}{Age (from otoliths)}
#'    \item{tl}{Total length (mm)}
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
#' @source From (approximately) Figure 2 of Maceina, M.J.  2007.  Use of piecewise nonlinear models to estimate variable size-related mortality rates.  North American Journal of Fisheries Management, 27:971-977.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(BlueCatfish)
#' str(BlueCatfish)
#' head(BlueCatfish)
#' plot(tl~age,data=BlueCatfish)
#' 
NULL
