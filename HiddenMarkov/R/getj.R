getj <- function(x, j){
    #   get the jth "row" from a list
    if (is.null(x)) return(NULL)
    n <- length(x)
    for (i in 1:n)
        x[[i]] <- x[[i]][j]
    return(x)
}

