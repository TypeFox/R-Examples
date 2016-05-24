
#............................................................
# generalized logistic distribution
##NS export(pgenlogis)
pgenlogis <- function( x , alpha1 = 0 , alpha2 = 0 ){
    x <- as.vector(x)
    xp <- x[ x >= 0 ]
    xn <- x[ x < 0 ]
    indp <- which( x >= 0 )
    indn <- which( x < 0 )
    # positive part of x
    if( ( alpha1 > 0 ) ){ y <-  ( ( exp( alpha1 * xp ) - 1  ) / alpha1  ) }
    if( ( alpha1 == 0 ) ) { y <-  xp  }
    if( ( alpha1 < 0 ) ) { y <- ( - log( 1 - alpha1 * xp )/alpha1 ) }
    # negative part of x
    if( ( alpha2 > 0 ) ){ y1 <-  - ( ( exp( alpha2 * abs(xn) ) - 1  ) / alpha2  ) }
    if( ( alpha2 == 0 ) ) { y1 <-  xn  }
    if( ( alpha2 < 0 ) ) { y1 <- (  log( 1 - alpha2 * abs(xn) )/alpha2 ) }
    # rearrange entries
    w <- rep( 0 , length(x) )
    if ( length(indp) >  0 ){ w[ indp ] <- y }
    if ( length(indn) >  0 ){ w[ indn ] <- y1    }
    # transform to inverse logistic probability
    y <- stats::plogis(w)
    return(y)
        }
#............................................................
# moments of generalized logistic distribution
genlogis.moments <- function( alpha1 , alpha2){
    x0 <- seq(-30,30 , len= 30000 )
    y0 <- pgenlogis( x=x0 , alpha1 = alpha1 , alpha2 = alpha2 )
    wgt <- y0[-1] - y0[ - length(y0) ] 
    wgt <- wgt / sum(wgt)
    out <- ( x0[ -1 ] + x0[ - length(x0) ] ) / 2
    M <- sum( wgt * out )
    SD <- sqrt( sum( wgt*out^2 ) - M^2 )
    moments <- c(M,SD,SD^2 )
    names(moments) <- c("M" , "SD" , "Var" )
    return(moments)
        }
