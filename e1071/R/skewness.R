skewness <-
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
    y <- sqrt(n) * sum(x ^ 3) / (sum(x ^ 2) ^ (3/2))
    if(type == 2) {
        if(n < 3)
            stop("Need at least 3 complete observations.")
        y <- y * sqrt(n * (n - 1)) / (n - 2)
    } else if(type == 3)
        y <- y * ((1 - 1 / n)) ^ (3/2)

    y
}

