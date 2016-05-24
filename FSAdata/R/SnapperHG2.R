#' @title Ages (subsample) and lengths (all fish) for Snapper.
#' 
#' @description A large sample (approximately fixed sample size per length interval) of Snapper (\emph{Pagrus auratus}) were aged, with the remainder of the fish just measured for length.  Note that age-16 is actually age 16+ and length 60 is for 60-64 cm and 65 if for 65+ cm.
#' 
#' @name SnapperHG2
#' 
#' @docType data
#' 
#' @format A data frame of 6724 observations on the following 2 variables:
#'  \describe{
#'    \item{len}{Measured lengths (cm)} 
#'    \item{age}{Ages assigned}
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Age-Length Key
#'  }
#'  
#' @concept 'Age-Length Key'
#' 
#' @source Recreated from summarized results in Table 8.3 of Quinn, T. J. and R. B. Deriso. 1999. Quantitative Fish Dynamics. Oxford University Press, New York, NY. 542 p.
#' 
#' @seealso See the same data in summarized format as \code{alkdata} in \pkg{fishmethods}.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(SnapperHG2)
#' str(SnapperHG2)
#' head(SnapperHG2)
#' 
#' ## Extract the aged sample
#' sn2.aged <- subset(SnapperHG2,!is.na(age))
#' str(sn2.aged)
#' 
#' ## Extract the length sample
#' sn2.length <- subset(SnapperHG2,is.na(age))
#' str(sn2.length)
#' 
NULL