
cv.split <- function( y, fold )
{
    n <- length(y)
    group <- table(y)
    x <- c()
    for ( i in 1:length(group) )
    {
        x.group <- c(1:n)[ y==names(group)[i] ] 
        x <- c( x, sample(x.group) )
    }
    foldi <- split( x, rep(1:fold, length = n) )
    
    return( foldi )
}
