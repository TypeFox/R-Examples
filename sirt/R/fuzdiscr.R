fuzdiscr <- function( X , theta0 = NULL , maxiter =200 , conv=.0001 ){
    if ( is.null(theta0) ){ theta0 <- rep( 1/ ncol(X) , ncol(X) ) }
    theta <- theta0
    iter <- 0
    change <- 1000
    while( ( iter < maxiter ) & ( change > conv ) ){
        # update xsi
        thetaM <- matrix( theta ,  nrow=nrow(X) , ncol=ncol(X) , byrow=TRUE)
        xsi <- thetaM * X  / rowSums( thetaM * X )
        # update theta
        theta <- colMeans( xsi)
        change <- max( abs( theta - theta0 ) )
        theta0 <- theta
        iter <- iter + 1
                }
    return(theta )
        }
