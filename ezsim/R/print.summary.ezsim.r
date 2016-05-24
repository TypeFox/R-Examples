#' Print a summary.ezsim Object in the console. See \code{\link{summary.ezsim}} for details
#' @name print.summary.ezsim
#' @aliases print.summary.ezsim
#' @title Print a summary.ezsim Object.
#' @usage \method{print}{summary.ezsim}(x,digits=4,...)
#' @param x A summary.ezsim Object
#' @param digits Number of digits the data will be rounded to.
#' @param \dots unused
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export 
#' @seealso \code{\link{summary.ezsim}}

print.summary.ezsim <-
function(x,digits=4,...){

    for( i in 1:length(x) ) 
        if (is.numeric(x[[i]]))
            x[[i]]<-round(x[[i]],digits=4)
            

    print.data.frame(x)
}
