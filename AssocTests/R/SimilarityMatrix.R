## eigenstrat Correlation Matrix
SimilarityMatrix <- function(x, k)
{
    x <- x/2
    if( k>1 ) 
    {
        z <- apply(x, 1, ModifyNormalization)
    }    
    else 
    {
        z <- ModifyNormalization(x)
        z <- matrix(z, ncol=1)
    }

    z %*% t(z)
}
