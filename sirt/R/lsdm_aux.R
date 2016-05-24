

#...........................................................................................##
# Function for calculating logistic functions and probability quantiles,                    ##
# especially Item Response Curves                                                           ##
#
# probcurves <- data ; quantiles <- quant.list ; est.icc <- TRUE
#
##NS # export(est.logist.quant)
est.logist.quant <- function( probcurves , theta , quantiles , est.icc = TRUE ){
    # INPUT:
    # probcurves    ... ( I x N ) matrix of Item response curves
    # theta         ... vector of length with discrete points of latent trait
    # quant         ... quantiles of probability curves which has to be estimated
    # estimate parameters of attribute response probcurves
    I <- nrow(probcurves)
    if (est.icc ){
            pars.probcurves <- matrix( 0 , nrow= I, ncol= 5 )
            colnames(pars.probcurves) <- c("b.2PL" , "a.2PL" , "sigma.2PL" , "b.1PL" , "sigma.1PL")
            rownames(pars.probcurves) <- rownames(probcurves)
            for (kk in 1:I ){
                pars.probcurves[kk,1:3] <- .est.logist( y = probcurves[kk,] , theta = theta )
                pars.probcurves[kk,4:5] <- .est.logist.rasch( y = probcurves[kk,] , theta = theta )
                    }
            }
    # quantiles of Item Response Curves (Logistic Functions)
    probcurves.quant <- sapply( quantiles , FUN = function( ql ){ 
            sapply( 1:I , FUN = function(kk){ 
                    .extract.probquantile( vec = probcurves[kk,] , theta = theta  , quant = ql )
                            } )
            } )
    probcurves.quant <- as.data.frame( probcurves.quant)
    colnames(probcurves.quant) <- paste( "Q" , 100*quantiles , sep="")
    rownames(probcurves.quant) <- rownames(probcurves)
    if (est.icc == TRUE ){    pars.probcurves <- cbind( probcurves.quant , pars.probcurves ) } else 
                            { pars.probcurves <- probcurves.quant }
    for (vv in 1:( length(quantiles) )    ){   pars.probcurves[,vv] <- as.numeric( pars.probcurves[,vv] ) }
    return( pars.probcurves )
    }
#...........................................................................................##  
#...........................................................................................##
# auxiliary function for estimating a logistic item response curve with                     ##
# slope and intercept parameter                                                             ##
.est.logist <- function( y , theta ){
    # INPUT:
    # y ... vector of y values (probabilities)
    # theta ... vector of theta values
    y <- as.numeric(y)
    objfct <- function( x ){ 
        mean( ( y - (1 + exp(  - x[2] * ( theta - x[1] )) )^(-1) )^2 )
        }
    res <- stats::optim( c(0,1) , objfct )
    return( c(  res$par , sqrt( res$value ) ) )
            }
#...........................................................................................##
#...........................................................................................##
# auxiliary function for estimating a logistic item response curve with                     ##
# intercept parameter (Rasch model)                                                         ##
.est.logist.rasch <- function( y , theta ){
    # INPUT:
    # y ... vector of y values (probabilities)
    # theta ... vector of theta values
    y <- as.numeric(y)
    objfct <- function( x ){ 
        mean( ( y - (1 + exp(  - 1 * ( theta - x )) )^(-1) )^2 )
        }
    res <- stats::optimize( objfct , lower = -10 , upper = 10 )
    return( c(  res$minimum , sqrt( res$objective ) ) )
            }
#...........................................................................................##
#.........................................................................................
# auxiliary function for extracting quantiles of curves
.extract.probquantile <- function( vec , theta , quant ){
    # INPUT:
    # vec   ... vector (probability function)
    # theta ... Grid of theta values
    # quant ... Quantile which has to be calculated
    x2 <- theta[ vec >= quant ][1]
    x1 <- sort( theta[ vec < quant ] , decreasing=T)[1]
    value <- - Inf 
    if ( (1-is.na(x2)) * (1-is.na(x1)) ==1 ){
        y1 <- vec[ theta == x1]
        y2 <- vec[ theta == x2]
        value <- x1 + ( quant - y1 ) * ( x2 - x1 ) / ( y2 - y1 )
        }
    if ( is.na(x2) ){ value <- Inf }
    return( value )
    }
#.........................................................................................


