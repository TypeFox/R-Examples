
# Correct the coefficients by setting zero
# if the confidence interval contains zero

"correct.spls" <-
function( object, plot.it=TRUE )
{
    betahat <- object$betahat
    lbmat <- object$lbmat
    ubmat <- object$ubmat

    # The plot of the original coefficient

    if ( plot.it==TRUE )
    {
        heatmap.spls( mat=betahat, coln=16,
        main='Original Coefficient Estimates',
        as='i', xlab='Predictors', ylab='Responses' )
    }

    # correct coefficients by setting b=0
    # if CI includes zero

    cimat <- lbmat*ubmat
    betahat2 <- betahat
    betahat2[ cimat<0 ] <- 0

    # The plot of the corrected coefficient

    if ( plot.it==TRUE )
    {
        dev.new()
        heatmap.spls( mat=betahat2, coln=16,
        main='Corrected Coefficient Estimates',
        as='i', xlab='Predictors', ylab='Responses' )
    }

    invisible(betahat2)
}
