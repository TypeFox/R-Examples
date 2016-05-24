#' @title Lengths for Largemouth Bass from Boomer Lake, OK.
#' 
#' @description Length measurements of Largemouth Bass (\emph{Micropterus salmoides}) from Boomer Lake, Oklahoma, in 1974.
#' 
#' @name LMBassBL
#' 
#' @docType data
#' 
#' @format A data frame of 447 observations on the following 1 variable:
#'  \describe{
#'    \item{tl}{measured total length (cm)} 
#'  }
#' 
#' @section Topic(s):
#'  \itemize{
#'    \item Length Frequency
#'    \item Size Structure
#'    \item PSD
#'  }
#' 
#' @concept 'Length Frequency' 'Size Structure' PSD
#' 
#' @source From McNew, R.W. and R.C. Summerfelt.  1978.  Evaluation of a maximum-likelihood estimator for analysis of length-frequency distributions.  Transactions of the American Fisheries Society 107:730-736.  Data was simulated (uniform distribution of values within length bin) from summarized length frequencies in http://fishbase.org/.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(LMBassBL)
#' str(LMBassBL)
#' head(LMBassBL)
#' hist(LMBassBL$tl,main="")
#' 
NULL