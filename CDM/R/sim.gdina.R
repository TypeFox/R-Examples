


################################################################################
# Simulation of the GDINA model
sim.gdina <- function( n , q.matrix , delta , link = "identity" , 
				thresh.alpha=NULL , cov.alpha=NULL , alpha=NULL , 
                Mj , Aj , necc.attr ){
    I <- length(delta)
    # simulate alpha
	if ( is.null(alpha) ){
#		library(MASS)
		alpha <- 1* ( MASS::mvrnorm( n , mu = thresh.alpha , Sigma = cov.alpha ) > 0 )
				}
    dat <- matrix( NA , n , I )
    for (ii in 1:I){
        # ii <- 1
        na.ii <- necc.attr[[ii]]
        Aj.ii <- Aj[[ii]]
        Mj.ii <- Mj[[ii]][[1]]
        delta.ii <- delta[[ii]]
        alpha.ii <- alpha[ , na.ii , drop=FALSE ] 
        # calculate probability for every attribute pattern
        patt.prob <- rowSums( Mj.ii * outer( rep(1,nrow(Mj.ii) ) , delta.ii  ) )
        # create patterns for alpha.ii and Aj.ii
        l1.Aj <- l1.al <- "P"
        for (vv in 1:( ncol(alpha.ii) ) ){
            l1.Aj <- paste( l1.Aj , Aj.ii[,vv] , sep="") 
            l1.al <- paste( l1.al , alpha.ii[,vv] , sep="") 
                            }
        resp.ii <- patt.prob[ match( l1.al , l1.Aj ) ]
        if ( link=="logit"){ resp.ii <- stats::plogis( resp.ii ) }
        if ( link=="log"){ resp.ii <- exp( resp.ii ) }
        dat[ , ii ] <- 1 * ( stats::runif( n ) < resp.ii )
                }    
    res <- list( "data" = dat , "alpha" = alpha , "q.matrix" = q.matrix ,
                        "delta" = delta , "Aj" = Aj , "Mj" = Mj , "link" = link )
    return(res)                    
    }
######################################################################################
# Function for preparation of GDINA simulation
sim.gdina.prepare <- function( q.matrix ){ 
    I <- nrow(q.matrix)             # number of items
    rsqm <- rowSums(q.matrix)       # row sums in Q matrix
    necc.attr <- delta <- Mj <- Aj <- as.list( rep(1,I) )
    for (ii in 1:I){
        # ii <- 1
        necc.attr[[ii]] <- which( q.matrix[ii,] > 0 )
        Aj[[ii]] <- .create.Aj( nq = rsqm[ii] )
        Mj[[ii]] <- .create.Mj( Aj = Aj[[ii]] , rule = "GDINA" )
        delta[[ii]] <- rep( 0 , ncol( Mj[[ii]][[1]] ) )
                    }
    res <- list( delta=delta , necc.attr = necc.attr , Aj = Aj ,  Mj = Mj )
                                        }
######################################################################################
