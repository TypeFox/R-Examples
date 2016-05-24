#' Print a parameters Object in the console
#' @name print.parameters
#' @aliases print.parameters
#' @title Print a parameters Object.
#' @method print parameters
#' @param x A parameters Object
#' @param \dots unused
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
#' @seealso \code{\link{parameterDef}}

print.parameters <-
function(x,...){
    y<-x[sapply(x,length)!=1]
    if (length(y)!=0){
        for (i in 1:length(y)){
            cat(names(y)[i]," : \n")
            print(y[[i]])
            cat("\n")
        }
    }
    
    y<-x[sapply(x,length)==1]
    cat(paste(names(y),y,sep='=',collapse=', '))
    cat("\n")
}
