
# Univariate Soft Thresholding Estimator

"ust" <-
function( b, eta )
{
    b.ust <- matrix( 0, length(b), 1 )
    if ( eta < 1 )
    {
        valb <- abs(b) - eta * max( abs(b) )
        b.ust[ valb>=0 ] <- valb[ valb>=0 ] * (sign(b))[ valb>=0 ]
    }
    return(b.ust)
}
