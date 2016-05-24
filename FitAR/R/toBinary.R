toBinary <-
function(n, k=ceiling(logb(n+1, base=2)))
{
if (n < 0)
    stop("n must be non-negative integer")
if (!is.loaded("asBinary"))
    stop("asBinary not loaded")
m <- numeric(k)
z<-.C("asBinary", as.integer(n), as.integer(m), as.integer(k))[[2]]
z   
}

