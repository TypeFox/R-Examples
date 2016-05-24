marginal.truescore.reliability <-
function(  b , a=1+0*b , c =0*b , d=1+0*b ,
    mean.trait=0 , sd.trait=1 , theta.k = seq( -6 , 6 , len=200) ){

    TT <- length(theta.k)
    I <- length(b)
    phi.k <- stats::dnorm( theta.k , mean=mean.trait , sd= sd.trait )
    phi.k <- phi.k / sum( phi.k )
    phi.kM <- matrix( phi.k , nrow=TT , ncol=I )
    
    aM <- matrix( a , nrow=TT , ncol= I , byrow=T )
    bM <- matrix( b , nrow=TT , ncol= I , byrow=T )
    cM <- matrix( c , nrow=TT , ncol= I , byrow=T )
    dM <- matrix( d , nrow=TT , ncol= I , byrow=T )
    
    # item response functions
    theta.kM <- matrix( theta.k , nrow=TT , ncol=I )
    icc.theta <- cM + (dM-cM)* stats::plogis( aM*(theta.kM - bM ) )
    
    # calculate pi_i (predicted probabilities)
    pi.i <- colSums( icc.theta * phi.kM )
    
    #  expected number correct
    mu <- sum(pi.i)
    
    # compute error variance
    sig2.error <- colSums( icc.theta * ( 1 - icc.theta ) * phi.k )
    
    # compute total error variance
    sig2.total.error <- sum( sig2.error )
    
    # compute total true score variance
    iccmean <- rowMeans( icc.theta )
    sig2.total.tau <- sum( ( I*iccmean )^2 * phi.k ) - ( sum( I*iccmean*phi.k ) )^2 
    
    # item level true score variance
    sig2.tau <- pi.i * ( 1 - pi.i ) - sig2.error
    
    # collect all formulas for item
    item <- data.frame( "item" = 1:I , "pi" = pi.i ,
            "sig2.tau" = sig2.tau , "sig2.error" = sig2.error )
    item$rel.item <- item$sig2.tau /  ( item$sig2.tau + item$sig2.error )
    
    rel.test <- sig2.total.tau / ( sig2.total.tau + sig2.total.error )
#    rel.test
    cat("Reliability=" , round( rel.test ,3 ),"\n")
    # Formula (15)
    # sum( sqrt( outer( sig2.tau , sig2.tau ) ) )
    res <- list("rel.test" = rel.test , "item" = item , "pi"=mu/I , 
                "sig2.tau" = sig2.total.tau , "sig2.error" = sig2.total.error )
    invisible(res)
        }
