###############################################
# conversion of the partial credit parametrization 
pcm.conversion <- function( b ){
    maxK <- rowSums( 1 - is.na(b) )
    I <- nrow(b)
    # define delta 
    delta <- b[ cbind( 1:I , maxK ) ] / maxK
    # define tau
    tau <- matrix( NA , nrow=I , ncol= ncol(b) )
    for (kk in 1:I){
    # kk <- 1
    if ( maxK[kk] == 1 ){
    tau[kk,1] <- 0
            }
    if ( maxK[kk] > 1){   
    # kk <- 4     
        for (vv in 1:maxK[kk] ){        
        # vv <- 1
        t1 <- 0
        if (vv > 1){ t1 <- sum( tau[kk,seq(1,vv-1)] )  }
        tau[kk,vv] <- b[kk,vv] - vv*delta[kk] - t1
                }
            }
    }
    res <- list( "delta"=delta , "tau" = tau )
    return(res)
        }
