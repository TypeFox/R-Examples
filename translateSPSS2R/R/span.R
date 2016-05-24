#' Indicates range of varlist
#'
#' Creates a list of variables within a specific range.
#'
#' @usage span(x, from = NULL, to = NULL, addDF = FALSE)
#' @param x a (non-empty) data.frame or input data of class "xpssFrame". 
#' @param from  the variable that opens the span.
#' @param to the variable that closes the span.
#' @param addDF if the name of the input data should be used?
#' @return Returns a varlist with the name of the variables which are within the range of the span indicator. 
#' @author Andreas Wygrabek
#' @examples 
#' data(fromXPSS)
#' span(x=fromXPSS,from="V1",to="V5",addDF=FALSE)
#' span(x=fromXPSS,from="V3",to="V1",addDF=TRUE)
#' @export
span <- function(x = NULL, from = NULL, to = NULL, addDF = FALSE){

    if(!(is.character(from) & is.character(to))){
        stop("from/to has to be character")
    } else {
    
    
    startindex <- which(colnames(x) == from)
    endindex <- which(colnames(x) == to)
    
    if(addDF != TRUE && addDF !=  FALSE) {
      stop("only TRUE or FALSE are allowed as statement")
    }
    
    if(addDF == FALSE){
        varlist <- colnames(x)[ startindex : endindex ]
        return(varlist)
    } else {
    varlist <- eval(paste(substitute(x),"$", colnames(x)[ startindex : endindex ], sep=""))
    return(varlist)
        }
    }
}





