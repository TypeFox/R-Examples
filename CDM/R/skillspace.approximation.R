###########################################
# skill space approximation
skillspace.approximation <- function( L , K , nmax=5000 ){
    n <- nmax
	ndim <- K
    res <- sfsmisc::QUnif (n, p=ndim, leap = 409)
    res <- 1*(res>.5)
    res <- rbind( rep( 0,ndim) , rep(1,ndim) , res )                
    v1 <- paste0("P" , res[,1] )
    for (vv in 2:ndim){ v1 <- paste0( v1 , res[,vv] ) }
    rownames(res) <- v1
    res <- res[ ! duplicated(v1) , ]
    res <- res[ 1:L , ]
    res <- res[ order( rownames(res) ) , ]        
    return(res)
        }
###########################################