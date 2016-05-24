#' @title Ages of Walleye assigned from otoliths, scales, and spines.
#' 
#' @description Age of Pymatuning Sanctuary (PA) Walleye (\emph{Sander vitreus}) assessed from three calcified structures -- sectioned otoliths, sectioned dorsal spines, and scale impressions.
#' 
#' @note Relationships between otoliths and spines and otoliths and scales are exact according to Figure 2.  Relationship between spines and scales is approximate as Figure 2 did not show this exact relationship.
#' 
#' @name WalleyePS
#' 
#' @docType data
#' 
#' @format A data frame with 61 observations on the following 4 variables:
#'  \describe{
#'    \item{otolith}{Age (years) assigned from broken, ground, and polished otolith sections} 
#'    \item{spine}{Age (years) assigned from dorsal spine sections}
#'    \item{scale}{Age (years) assigned from scale impressions}
#'    \item{sex}{Sex of fish (\code{female} and \code{male})}
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
#' @source From Figure 2 in Kocovsky, P.M., and R.M. Carline.  2000.  A comparison of methods for estimating ages of unexploited walleyes.  North American Journal of Fisheries Management 20:1044-1048.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(WalleyePS)
#' str(WalleyePS)
#' head(WalleyePS)
#' op <- par(mfrow=c(3,2),pch=19)
#' plot(scale~otolith,data=WalleyePS,subset=sex=="female",main="Female")
#' plot(scale~otolith,data=WalleyePS,subset=sex=="male",main="Male")
#' plot(scale~spine,data=WalleyePS,subset=sex=="female",main="Female")
#' plot(scale~spine,data=WalleyePS,subset=sex=="male",main="Male")
#' plot(spine~otolith,data=WalleyePS,subset=sex=="female",main="Female")
#' plot(spine~otolith,data=WalleyePS,subset=sex=="male",main="Male")
#' par(op)
#' 
NULL