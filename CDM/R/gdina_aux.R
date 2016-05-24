
#*********************************
# copied squeeze function from mice package
squeeze.cdm <- function (x, bounds = c(min(x[r]), max(x[r])), r = rep(TRUE, 
    length(x))) 
{
    if (length(r) != length(x)) 
        stop("Different length of vectors x and r")
    x[x < bounds[1]] <- bounds[1]
    x[x > bounds[2]] <- bounds[2]
    return(x)
}