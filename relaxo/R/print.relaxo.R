`print.relaxo` <-
function(x, ...) {
    cat("Relaxed Lasso ('relaxo'): object\n")
    ## This is very cheap: 
    str(x, ...)
    invisible(x)
}

