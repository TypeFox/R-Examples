#' Moves cases of variables
#'
#' computeLag shifts the dataset forward or backward by a given number of observations.
#' 
#' Creates shifted data dependent upon the direction in which the data got moved. \cr A positive indicator refers a shift to the right side, a negative indicator refers a shift of the data to the left side. Empty cases get filled with NA.
#'
#' @usage computeLag(x, move = 0, value = NA)
#' @param x a (non-empty) data.frame, data.table object or input data of class "xpssFrame". 
#' @param move atomic integer, either positive or negative that defines the cases to move. The algebraic sign indicates the direction. 
#' @param value atomic numeric or atomic character value that replaces the skipped cases.
#' @return Output is the shifted, respectively "lagged" input vector. \cr Length of the new lagged vector is identical with the length of the input vector.
#' @author Andreas Wygrabek
#' @examples 
#' 
#' data(fromXPSS)
#' computeLag(x=fromXPSS$V6, move = 2, value = NA)
#' @keywords internal
#' @export
computeLag <- function(x, move = 0, value = NA){ 
    
  #functiontype  <- "SB"
  
    dfcheck <- deparse(substitute(x))
    dfcheck <- strsplit(dfcheck, "\\$")[[1]][1]
  
  absMove <- abs(move)   
    
    if (move < 0 ){       
        OBJ <- c(tail(x,-absMove),rep(value,absMove))
    }else if (move > 0 ){        
        OBJ <- c(rep(value,absMove), head(x,-absMove))
    }else {       
        OBJ <- x
       }   
    return(OBJ)
}
