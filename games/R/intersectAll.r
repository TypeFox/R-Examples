##
## INPUT:
## ...: vectors of any standard class
##
## RETURN:
## vector of elements contained in all vectors in '...'
##
intersectAll <- function(...)
{
    x <- list(...)
    ans <- x[[1]]
    for (i in 1:length(x)) ans <- intersect(ans, x[[i]])

    return(ans)
}
