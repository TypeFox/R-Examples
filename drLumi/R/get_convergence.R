## Get convergence. The criteria from the model is not enough
get_convergence <- function(x) {
    convergence <- x$convInfo$isConv
    ans <- 2
    if(convergence==TRUE) ans <- 1
    return(ans)
}


