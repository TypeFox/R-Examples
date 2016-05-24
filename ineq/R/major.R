major <- function(x,y)
{
    x <- sort(as.numeric(x))
    y <- sort(as.numeric(y))
    n <- length(x)
    if((length(y)==n)&(sum(x)==sum(y)))
        all((cumsum(x)-cumsum(y))<=0)
    else
        stop("incomparable arguments")
}
