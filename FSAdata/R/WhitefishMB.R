#' @title Ages of Lake Whitefish from four lakes assigned from scales and fin-rays.
#' 
#' @description Ages of Lake Whitefish (\emph{Coregonus clupeaformis}) from four lakes as determined by scales and fin-rays.
#' 
#' @name WhitefishMB
#' 
#' @docType data
#' 
#' @format A data frame with 859 observations on the following 3 variables:
#'  \describe{
#'    \item{fin}{Ages assigned from fin-ray sections} 
#'    \item{scale}{Ages assigned from scales} 
#'    \item{lake}{Lake from which the fish was captured (\code{L122}, \code{L226}, \code{Huron}, or \code{Dezadeash})} 
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
#' @source From (approximately) Figure 1 of Mills, K.H., and R.J. Beamish.  1980.  Comparison of fin-ray and scale age determinations for lake whitefish (\emph{Coregonus clupeaformis}) and their implications for estimates of growth and annual survival.  Canadian Journal of Fisheries and Aquatic Sciences, 37:534-544.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(WhitefishMB)
#' str(WhitefishMB)
#' head(WhitefishMB)
#' op <- par(mfrow=c(2,2),pch=19)
#' plot(scale~fin,data=WhitefishMB,subset=lake=="L122",main="Lake L122")
#' plot(scale~fin,data=WhitefishMB,subset=lake=="L226",main="Lake L226")
#' plot(scale~fin,data=WhitefishMB,subset=lake=="Huron",main="Huron")
#' plot(scale~fin,data=WhitefishMB,subset=lake=="Dezadeash",main="Dezadeash")
#' par(op)
#' 
NULL