#' Sequential numbering of cases
#'
#' R implementation of the SPSS \code{$CASENUM} system variable. xpssCasenum counts the number of cases within a variable, or respectively the number of observations in the dataset.
#' 
#'  xpssCasenum fit well as ID-Variable.
#'
#' @usage xpssCasenum(x)
#' @param x a (non-empty) data.frame or input data of class \code{"xpssFrame"}. 
#' @return Returns a sequential atomic numeric or numeric vector.
#' @author Bastian Wiessner
#' @seealso Related Functions \code{\link{xpssCount}} , \code{\link{xpssDate}} , \code{\link{xpssDate11}}
#' @importFrom data.table is.data.table
#' @examples
#' data(fromXPSS)
#' 
#' fromXPSS$id <- xpssCasenum(fromXPSS)
#' @export

xpssCasenum <- function(x)
{
  stopifnot(is.data.frame(x) | is.data.table(x) | "xpssFrame" %in% class(x))
  
  casenum <- 1:length(x[[1]])
  return(casenum)
}
