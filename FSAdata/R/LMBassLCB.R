#' @title Lengths for Largemouth Bass from Lake Carl Blackwell, OK.
#' 
#' @description Lengthsfor Largemouth Bass (\emph{Micropterus salmoides}) from Lake Carl Blackwell, Oklahoma, in 1973.
#' 
#' @name LMBassLCB
#' 
#' @docType data
#' 
#' @format A data frame of 289 observations on the following variable:
#'  \describe{
#'    \item{tl}{Measured total length (cm).} 
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
#' @source From McNew, R.W. and R.C. Summerfelt.  1978.  Evaluation of a maximum-likelihood estimator for analysis of length-frequency distributions. Transactions of the American Fisheries Society 107:730-736.  Data was simulated (uniform distribution of values within length bin) from summarized length frequencies in http://fishbase.org/.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(LMBassLCB)
#' str(LMBassLCB)
#' head(LMBassLCB)
#' hist(LMBassLCB$tl,main="")
#' 
NULL