#' @title Ages of Yellowtail Flounder assigned from scales and otoliths.
#' 
#' @description Ages of commercially caught Georges Bank Yellowtail Flounder (\emph{Limanda ferruginea}) as determined by scales, whole otoliths, or otolith cross-sections.
#' 
#' @name YTFlounder
#' 
#' @docType data
#' 
#' @format A data frame with 27 paired observations on the following 3 variables.
#'  \describe{
#'    \item{scale}{Ages assigned from scales}
#'    \item{whole}{Ages assigned from whole otoliths} 
#'    \item{cross}{Ages assigned from cross-sections of otoliths} 
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
#' @source From tables in Annex 3 of Walsh, S.J. and J. Burnett. 2002.  The Canada-United States yellowtail flounder age reading workshop: 28-30 November2000, St. John's, Newfoundland.  North Atlantic Fisheries Organization.  Scientific Council Studies 35:1-59.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(YTFlounder)
#' str(YTFlounder)
#' head(YTFlounder)
#' op <- par(mfrow=c(2,2),pch=19)
#' plot(scale~whole,data=YTFlounder)
#' plot(scale~cross,data=YTFlounder)
#' plot(whole~cross,data=YTFlounder)
#' par(op)
#' 
NULL