#' Creates a DO IF - END IF subset
#'
#' R implementation of the SPSS \code{DO IF} argument. 
#'
#' xpssDoIf selects cases for analysis based on one or more logical conditions. The conditions to select cases are specified in a logical expression. These logical expressions can contain relational operators, logical operators and arithmetic operations. xpssDoIf creates a subset of the actual dataset, without deleting the excluded variables, respectively without deleting the excluded values.\cr \cr
#' 
#' If its needed to modify more than one subset, every following subset get selected by \code{\link{xpssElseIf}}. \code{xpssDoIf} is working similar to \code{xpssElseIf}, the only difference is that \code{xpssDoIf} has to initiate the subsetting.
#' 
#' The data is subsetted until \code{\link{xpssEndIf}} restores the excluded data. All changes made at the subsetted data will be taken over, the excluded data will remain untouched! As noted before, those cases are not actually deleted and will be available after \code{xpssEndIf} restores the excluded data. \cr \cr
#'
#'  In a different way to SPSS. Not only data management functions like \code{xpssRecode} can be used within \code{xpssDoIf}, it is possible to use statistical and descriptiv functions like \code{xpssFrequencies} too. \cr \cr
#'  \strong{NOTE:} For temporary case selection, specify \code{\link{xpssTemporary}} before \code{xpssDoIf}.

#' @usage xpssDoIf(x, cond = NULL)
#' @param x a (non-empty) data.frame or input data of class \code{"xpssFrame"}. 
#' @param cond logical expression indicating the condition for subsetting.
#' @return Output is a subset of the actual dataset under the condition of the logical expression.
#' @author Andreas Wygrabek
#' @seealso Related Functions \code{\link{xpssElseIf}}, \code{\link{xpssEndIf}}, \code{\link{xpssFilter}}, \code{\link{xpssSelectIf}}, \code{\link{xpssTemporary}}
#' @examples
#' data(fromXPSS)
#' 
#' temp <- xpssDoIf(x=fromXPSS, cond = "V3 == 1")
#' 
#' temp <- xpssRecode(x=temp,variables="V5",rec="lo:78 = 1; else = 2")
#' 
#' temp <- xpssEndIf(x=temp)
#' @export
xpssDoIf <- function(x, cond = NULL){
  
  
  
  #backup attributes
  attr_backup <- attributesBackup(x)
  
  functiontype <- "ME"
  x <- applyMetaCheck(x)
  
    attr(x, "DO_IF_INVERSE") <- x
    
    attributes(x)$DO_IF <- cond
    
    attBack <- attributesBackup(x)
    
    x <- subset(x, subset = eval(parse(text = cond)))
    
    
    x <- applyAttributes(x, attBack)
    
    logVec <-  !is.element(rownames(attributes(x)$DO_IF_INVERSE), rownames(x))
    attributes(x)$DO_IF_INVERSE <- attributes(x)$DO_IF_INVERSE[logVec,]
    
    attributes(x)$DO_IF_INVERSE <- applyAttributes(attributes(x)$DO_IF_INVERSE, attBack)
  
    
    return(x)
}


