##' timeop function
##'
##' A function to time an operation in R
##'
##' @param expr an expression to evaluate
##' @return The time it took to evaluate the expression in seconds
##' @export

timeop <- function(expr){
    s<-Sys.time()
    eval(expr,envir=parent.frame())
    e<-Sys.time()
    return(difftime(e,s,units="secs"))
}
