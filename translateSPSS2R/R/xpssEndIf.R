#' Ends a DO IF - END IF subset
#'
#' R implementation of the SPSS \code{END IF} argument.
#' 
#'   \code{xpssEndIf} determines the end of the analysis based on logical conditions via \code{\link{xpssDoIf}}.  \code{xpssEndIf} merge the excluded data with the actual dataset, after \code{xpssDoIf} subsetted the data. All changes which were made until \code{xpssEndIf} will be taken over, the excluded data will remain untouched!\cr \cr
#' 
#'  \strong{NOTE:} For temporary case selection, specify \code{xpssTemporary} before \code{xpssDoIf}.
#'
#' @usage xpssEndIf(x)
#' @param x a (non-empty) data.frame or input data of class \code{"xpssFrame"}. 
#' @return Output is the original dataset.
#' @author Andreas Wygrabek
#' @examples
#' data(fromXPSS)
#' 
#' temp <- xpssDoIf(x=fromXPSS, cond = "V3 == 1")
#' 
#' temp <- xpssRecode(x=temp,variables="V5",rec="lo:78 = 1; else = 2")
#' 
#' temp <- xpssEndIf(x=temp)
#' 
#' @export
xpssEndIf <- function(x){
    
  functiontype <- "ME"
  x <- applyMetaCheck(x)
    
    attributes(x)$DO_IF <- FALSE
    
    attBack <- attributesBackup(x)
    
  if(length(x) != length(attributes(x)$DO_IF_INVERSE)){
    cols <- which(!(names(x) %in% names(attributes(x)$DO_IF_INVERSE)))
    attributes(x)$DO_IF_INVERSE[names(x[cols])] <- NA
  }
  
    x <- rbind(x,attributes(x)$DO_IF_INVERSE)
    
    x <- x[order(as.numeric(rownames(x))),] 

    x <- applyAttributes(x, attBack)
    
    attributes(x)$DO_IF_INVERSE <- NULL
    
    return(x)
}
