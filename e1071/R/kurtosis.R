kurtosis <-
function(x, na.rm = FALSE, type = 3)
{
    if(any(ina <- is.na(x))) {
        if(na.rm) 
            x <- x[!ina]
        else
            return(NA)
    }
    
    if(!(type %in% (1 : 3)))
       stop("Invalid 'type' argument.")
    
    n <- length(x)
    x <- x - mean(x)
    r <- n * sum(x ^ 4) / (sum(x ^ 2) ^ 2)
    y <- if(type == 1)
        r - 3
    else if(type == 2) {
        if(n < 4)
            stop("Need at least 4 complete observations.")
        ((n + 1) * (r - 3) + 6) * (n - 1) / ((n - 2) * (n - 3))
    }
    else
        r * (1 - 1 / n) ^ 2 - 3

    y
}
