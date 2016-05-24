#' Ages (subsample) and lengths (all fish) for Striped Bass.
#' 
#' As many as 10 fish per 1-inch total length intervals from the \code{StripedBass2} data frame were obtained for age assignment.  The remaining fish in the file were only measured for length (i.e., the ages were deleted).  This data file can be used to demonstrate the use of age-length keys.
#' 
#' @name StripedBass3
#'
#' @docType data
#' 
#' @format A data frame of 1201 observations on the following 2 variables:
#'  \describe{
#'    \item{tl}{Measured total lengths (in inches).} 
#'    \item{age}{Ages assigned from examination of otoliths.} 
#'  }
#' 
#' @section Topic(s):
#'  \itemize{
#'    \item Age-Length Key 
#'  }
#' 
#' @concept 'Age-Length Key'
#'
#' @seealso \code{\link{StripedBass2}}.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(StripedBass3)
#' str(StripedBass3)
#' head(StripedBass3)
#' 
#' ## Extract the aged sample
#' sb.aged <- subset(StripedBass3,!is.na(age))
#' str(sb.aged)
#' 
#' ## Extract the length sample
#' sb.length <- subset(StripedBass3,is.na(age))
#' str(sb.length)
#' 
NULL