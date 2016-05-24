"getInitial" <- function(object)
{
    initval <- object$"start"
    names(initval) <- object$"parNames"[[2]]
    
    initval
}