subset.dist <- function (dist,subset) 
{
    if (!inherits(dist,'dist')) stop('You must pass an object of class dist')

    dist <- as.dist(as.matrix(dist)[subset,subset])
    dist
}
