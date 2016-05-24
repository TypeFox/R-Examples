#' @title Ages and lengths of Spotted Sucker.
#'
#' @description Ages and total lengths of Spotted Sucker (\emph{Minytrema melanops}) collected from the Apalachicola River, Florida.
#'
#' @name SpottedSucker1
#' 
#' @docType data
#' 
#' @format A data frame with 96 observations on the following 2 variables.
#'   \describe{
#'     \item{tl}{Total length (mm).}
#'     \item{age}{Age (from scales).} 
#'   }
#'   
#' @section Topic(s):
#'   \itemize{
#'     \item Growth
#'     \item von Bertalanffy
#'   }
#'   
#' @concept Growth 'von Bertalanffy'
#'
#' @source From Box 5.4 in Iseley, J.J. and T.B. Grabowski.  2007.  Age and Growth in Guy, C.S. and M.L. Brown, editors.  Analysis and Interpretation of Freshwater Fisheries Data.  American Fisheries Society.  Likely originally from Grabowski, T.B., S.P. Young, J.J. Isely, and P.C. Ely.  Age, growth, and reproductive biology of three catostomids from the Apalachicola River, Florida.  Journal of Fish and Wildlife Management 3:223-237.  [Was (is?) from http://www.fwspubs.org/doi/pdf/10.3996/012012-JFWM-008.]
#'
#' @keywords datasets
#'
#' @examples
#' data(SpottedSucker1)
#' str(SpottedSucker1)
#' head(SpottedSucker1)
#' plot(tl~age,data=SpottedSucker1)
#'
NULL