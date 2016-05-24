optimizelogpost <-
function( alpha,beta,D, y, By,  k, func )
{

# Note: this function needs to be made generic to accept different values
#   of F; right now it works only if F = 3.

F <- ncol( D )  + 1

# F = number of precision parms (meas err precision and for all Q's)

optimouts <- numeric()

# We need to run optimization from different starting values to increase
# chance of finding global max.  

if( F ==5 ) {
        optimout1 <- optim( c( 6/10, 1/10, 1/10, 1/10 ), func, alpha=alpha, beta = beta,
        D = D,   y = y,By = By,  k = k,              
        control = list(fnscale = -1) )$value

        optimout2 <- optim( c( 1/10, 1/10, 1/10, 1/10 ), func, alpha=alpha, beta = beta,
        D = D, y = y, By = By, k = k,
        control = list(fnscale = -1) )$value

        optimout3 <- optim( c( 1/10, 1/10, 1/10, 6/10 ), func, alpha=alpha, beta = beta,
        D = D,  y = y, By = By, k = k,              
        control = list(fnscale = -1) )$value
        optimout4 <- optim( c( 1/10, 1/10, 6/10, 1/10 ), func, alpha=alpha, beta = beta,
        D = D,  y = y, By = By, k = k,              
        control = list(fnscale = -1) )$value
        optimout5 <- optim( c( 1/10, 6/10, 1/10, 1/10 ), func, alpha=alpha, beta = beta,
        D = D,  y = y, By = By, k = k,              
        control = list(fnscale = -1) )$value

        gmax <- max( optimout1, optimout2, optimout3, optimout4, optimout5 )
}
else if( F ==4 ) {
        optimout1 <- optim( c( 1/16, 3/16, 5/16 ), func, alpha=alpha, beta = beta,
        D = D,   y = y,By = By,  k = k,              
        control = list(fnscale = -1) )$value

        optimout2 <- optim( c( 3/16, 5/16, 7/16 ), func, alpha=alpha, beta = beta,
        D = D, y = y, By = By, k = k,
        control = list(fnscale = -1) )$value

        optimout3 <- optim( c( 5/16, 7/16, 1/16 ), func, alpha=alpha, beta = beta,
        D = D,  y = y, By = By, k = k,              
        control = list(fnscale = -1) )$value
        optimout4 <- optim( c( 7/16, 1/16, 3/16 ), func, alpha=alpha, beta = beta,
        D = D,  y = y, By = By, k = k,              
        control = list(fnscale = -1) )$value

        gmax <- max( optimout1, optimout2, optimout3, optimout4 )
}
else if( F ==3 ) {
        #optimout1 <- optim( c( 1/9 , 3/9 ), func, alpha=alpha, beta = beta,
        optimout1 <- optim( c( 1/6 , 1/2 ), func, alpha=alpha, beta = beta,
        D = D,   y = y,By = By,  k = k,              
        control = list(fnscale = -1) )$value

        optimout2 <- optim( c( 1/2 , 1/3 ), func, alpha=alpha, beta = beta,
        D = D, y = y, By = By, k = k,
        control = list(fnscale = -1) )$value

        optimout3 <- optim( c( 1/3 , 1/6 ), func, alpha=alpha, beta = beta,
        D = D,  y = y, By = By, k = k,              
        control = list(fnscale = -1) )$value

        gmax <- max( optimout1, optimout2, optimout3 )
}

else
    if (F==2) {
        max1 <- optimize( func,c(0,0.5),alpha=alpha, beta = beta,
        D = D, y = y, By = By, k = k,       maximum=T)$objective
        max2 <- optimize( func,c(0.5,1.0), alpha=alpha, beta = beta,
        D = D, y = y, By = By, k = k,       maximum=T)$objective

        gmax<- max(max1,max2)

}

gmax


}

