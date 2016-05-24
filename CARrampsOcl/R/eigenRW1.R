eigenRW1 <- function(n)
{
# compute eigenvals and eigenvects of RW1

w <-( pi/n) * ((n-1):0)
evals <-  2 * (1-cos(w) )

evects <- sin( outer( 1:n, w ) ) - sin( outer( (0:(n-1)), w))
evects[,n] <- rep(1,n)

evects <- apply( evects, 2, function(v) v / sqrt( sum(v^2) ) )

list( values = evals, vectors = evects)

}

