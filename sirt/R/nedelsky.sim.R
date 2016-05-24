

################################################################
# simulate data according to the Nedelsky model
nedelsky.sim <- function( theta , b , a=NULL , tau=NULL ){
    Theta <- matrix( theta , ncol=1 )
    N <- length(theta)
    # generate all combinations
    nodes <- c(0,1)
    ndim <- K <- ncol(b)
    I <- nrow(b)
    if ( is.null(a) ){
        a <- rep(1,I)
                }
    combis <- as.matrix( expand.grid( as.data.frame( matrix( rep(nodes, ndim) , 
            ncol = ndim ) ) ) )
    if ( is.null(tau) ){
        tau <- matrix( 1 , nrow=I , ncol=K+1 )
                        }
    # dataset
    dat <- matrix( NA , nrow=N , ncol=I)
    for (ii in 1:I){
        # ii <- 1
        b0 <- as.vector(b[ii,])
        a0 <- a[ii]
        thdim <- 1
        tau0 <- as.vector(tau[ii,])
        probs <- nedelsky.irf( Theta , K , b=b0 , a=a0 , tau=tau0 , combis  )$probs
        cprobs <- rowCumsums.sirt(matr=probs)
        rn1 <- stats::runif(N)
        dat[,ii] <- rowIntervalIndex.sirt(matr=cprobs,rn=rn1) - 1
                    }
	colnames(dat) <- paste0( "I" , 100+1:I )
    return(dat)
        }
################################################################

