#' @title Ages, lengths, and sexes of Troutperch.
#' 
#' @description The assigned ages (by scales), total lengths (mm), and sexes of Troutperch (\emph{Percopsis omsicomaycus}) captured in southeastern Lake Michigan.
#' 
#' @name TroutperchLM1
#' 
#' @docType data
#' 
#' @format A data frame with 431 observations on the following 3 variables:
#'  \describe{
#'    \item{age}{Assigned ages (by scales).} 
#'    \item{tl}{Measured total length (mm).} 
#'    \item{sex}{Sex (\code{f}=female and \code{m}=male).} 
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
#' @source Simulated from the age-length data provided in Table 1 of House, R., and L. Wells.  1973.  Age, growth, spawning season, and fecundity of the trout-perch (\emph{Percopsis omsicomaycus}) in southeastern Lake Michigan.  Journal of the Fisheries Research Board of Canada.  30:1221-1225.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(TroutperchLM1)
#' str(TroutperchLM1)
#' head(TroutperchLM1)
#' op <- par(mfrow=c(1,2),pch=19)
#' plot(tl~age,data=TroutperchLM1,subset=sex=="f",main="female")
#' plot(tl~age,data=TroutperchLM1,subset=sex=="m",main="male")
#' par(op)
#' 
NULL