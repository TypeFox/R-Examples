#' Temporary modification of the data
#'
#'  R implementation of the SPSS \code{TEMPORARY} Function. 
#' 
#' xpssTemporary signals the beginning of temporary transformation. Only the following data-management procedure takes affect on the data: , \code{xpssCount}, \code{xpssDoIf}, \code{xpssFilter}, \code{xpssMissingValues}, \code{xpssNumeric}, \code{xpssNofCases}, \code{xpssRecode}, \code{xpssSelectIf}, \code{xpssSample},   \code{xpssString}, \code{xpssVariableLabels}, \code{xpssValueLabels},.
#' \cr\cr All the changes that are made are temporary. After the next modification the data is restored. \cr
#' For example:  all created variables, e.g. numeric or string variables created while \code{xpssTemporary} is in effect are temporary variables! \cr
#' Any changes or modifications made to existing variables while the \code{xpssTemporary} command is in effect are also \strong{temporary}! \cr
#' Any variables which are created or modified after this procedure are again permanent. 
#' 
#'  \tabular{rlll}{
#' 
#' \tab \code{Function Status} \tab Effect on created variables \tab Effect on modified variables \cr 
#' \tab \code{Temporary ON} \tab Any created variable is temporary. \tab Any change on a variables is temporary. \cr 
#' \tab \code{Temporary OFF} \tab Any created variable is permanent. \tab Any change on a variables is permanent.} 
#' 
#' The xpssTemporary function allows analyses for subgroups without affecting the data and then repeat the analysis for the file as a whole. 
#'
#' @param x a (non-empty) data.frame or input data of class "xpssFrame". 
#' @author Andreas Wygrabek
#' @examples
#' data(fromXPSS)
#' obj <- xpssTemporary(fromXPSS)
#' @seealso Related Functions \code{\link{xpssDoIf}} \code{\link{xpssFilter}} \code{\link{xpssSample}} \code{\link{xpssSelectIf}}
#' @export
xpssTemporary <- function(x){
  
  functiontype <- "ME"
  x <- applyMetaCheck(x)
  
  attr(x, "ORIGIN") <- x

    attributes(x)$TEMPORARY <- TRUE
    
    return(x)
}
