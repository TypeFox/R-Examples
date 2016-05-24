#' @title Ages of Lake Huron Alewife assigned from otoliths and scales.
#' 
#' @description Ages of Alewife (\emph{Alosa pseudoharengus}) assessed from otoliths and scales.
#' 
#' @name AlewifeLH
#' 
#' @docType data
#' 
#' @format A data frame of 104 observations on the following 2 variables:
#'  \describe{
#'    \item{otoliths}{Age assigned from examination of otoliths}
#'    \item{scales}{Age assigned from examination of scales}
#'  }
#' 
#' @section Topic(s):
#'  \itemize{
#'    \item Age Comparison
#'    \item Age Precision 
#'    \item Age Bias
#'    \item Ageing Error
#'  }
#' 
#' @concept Age Precision Bias 'Age Comparison'
#' 
#' @source From Table 1A of Hoenig, J.M., M.J. Morgan, and C.A. Brown. 1995.  Analysing differences between two age determination methods by tests of symmetry.  Canadian Journal of Fisheries And Aquatic Systems, 52:364-368.  Originally from O'Gorman, R., D.H. Barwick, and C.A. Bowen.  1987.  Discrepancies between ages determined from scales and otoliths for alewives from the Great Lakes.  pp. 203-210 In.  Summerfelt, R.C. and G.E. Hall (Eds.)  Age and Growth of Fish.  Iowa State University Press, Ames, IA.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(AlewifeLH)
#' str(AlewifeLH)
#' head(AlewifeLH)
#' plot(scales~otoliths,data=AlewifeLH)
#' xtabs(~otoliths+scales,data=AlewifeLH)
#' 
NULL
