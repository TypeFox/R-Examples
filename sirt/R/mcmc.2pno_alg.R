

######################
# draw latent responses Z
.draw.Z.2pl <- function( aM , bM, theta , N , I , threshlow , threshupp ){
    # calculate means
    mij <- aM * theta - bM
    # simulate uniform data
    rij <- matrix( stats::runif( N*I ) , nrow=N , ncol=I )
    # calculate corresponding value
    pl <- stats::pnorm( threshlow , mean=mij) 
    pu <- stats::pnorm( threshupp , mean=mij)
    pij <- pl + (pu-pl)*rij
    # simulate Z
    Zij <- stats::qnorm( pij , mean = mij )
    return(Zij)
        }
		
########################################
# draw theta
.draw.theta.2pl <- function( aM , bM , N , I , Z ){
    vtheta <- 1 / ( rowSums( aM^2 ) + 1  )
    mtheta <- rowSums( aM * ( Z + bM ) ) * vtheta
    theta <- stats::rnorm( N , mean=mtheta , sd = sqrt( vtheta ) )
#	v1 <- var(mtheta)
#	EAP.rel <- v1 / ( mean(vtheta) + v1 )	
#	res <- list("theta" = theta , "EAP.rel" = EAP.rel )
    res <- list("theta" = theta )
    return(res)
            }

##########################################
# draw item parameters
.draw.itempars.2pl <- function( theta , Z , I , N , weights){
	#--------------
	# sampling without weights
	Xast <- as.matrix( cbind( theta , -1 ) )
    if ( is.null(weights) ){
#		Sigma <- solve( t(Xast) %*% Xast )
		Sigma <- solve( crossprod(Xast) )
		# calculate mean
#		mj <- Sigma %*% t(Xast) %*% Z
		mj <- Sigma %*% crossprod(Xast , Z )
		mj <- as.matrix( t(mj))		    	
						}
	#--------------
    # sampling with weights						
	if ( ! is.null( weights ) ){
		# compute elements of Xast
		Xast11 <- sum( theta^2 * weights )
		Xast12 <- - sum( theta * weights )
		Xast22 <- sum( weights )
		# compute inverse of Xast
		Xastdet <- Xast11*Xast22 - Xast12^2 
		Xastinv11 <- Xast22 / Xastdet
		Xastinv22 <- Xast11 / Xastdet	
		Xastinv12 <- - Xast12 / Xastdet
		Sigma <- matrix( c(Xastinv11 , Xastinv12 , Xastinv12 , Xastinv22) , 2 ,2 )
		# compute t(Xast) %*% Z (weighted)
#		mj <- Sigma %*% t( Xast * weights ) %*% Z
		mj <- Sigma %*% crossprod( Xast * weights , Z )
	    mj <- as.matrix( t(mj))	
			}		
    #--------------							
    # draw item parameters
    ipars <- mvtnorm::rmvnorm( I , sigma=Sigma ) + mj
    a <- ipars[,1]
    b <- ipars[,2]
    return( list( "a"=a , "b"=b) )
            }

			
#################################################
# compute deviance
.mcmc.deviance.2pl <- function( aM , bM , theta , dat , dat.resp , 
				weights , eps ){
	pij <- stats::pnorm( aM * theta - bM )
	llij <- log( dat.resp * ( dat*pij + ( 1-dat )*(1-pij) ) + eps )
	if ( is.null( weights ) ){ deviance <- -2*sum( llij ) }
	if ( ! is.null( weights ) ){ 
		deviance <- -2*sum( rowSums(llij) * weights ) 
					}
	return(deviance)
				}
				
######################################################				
# subfunction for calculating the DIC
.mcmc.ic.2pl <- function( a.chain , b.chain , theta.chain , N , I , 
		dat , dat.resp , weights , eps ,
		deviance.chain ){
	#************
	aM <- matrix( colMeans( a.chain ) , nrow=N , ncol=I , byrow=TRUE )
	bM <- matrix( colMeans( b.chain ) , nrow=N , ncol=I , byrow=TRUE )	
	theta <- colMeans( theta.chain )
	Dhat <- .mcmc.deviance.2pl( aM , bM , theta , dat , dat.resp , 
				weights , eps )
	Dbar <- mean( deviance.chain )
	pD <- Dbar - Dhat
	ic <- list( "Dhat"=Dhat , "Dbar"=Dbar , "pD"=pD , "DIC" = Dhat + 2*pD )
	return(ic)
		}
#######################################################				
				
				
				
####################################################################
# calculate EAP reliability, person parameter EAPs
# and corresponding posterior SDs
.mcmc.person.2pno <- function( theta.chain , weights ){	
	###################
	# EAP reliability
		v1 <- stats::var( colMeans( theta.chain ) )
		if ( is.null(weights) ){ EAP.rel <- v1 / 1 }
		if ( ! is.null(weights) ){ 
				w1 <- weights / sum(weights )
				m1 <- colMeans( theta.chain )
				v1 <- sum( m1^2 * w1 ) - ( sum( m1*w1 ) )^2
				wM <- matrix( w1 , nrow=nrow(theta.chain) , ncol=ncol(theta.chain) , byrow=TRUE )			
				h1 <- rowSums( wM * theta.chain^2 ) - ( rowSums( wM * theta.chain ) )^2
				h1 <- mean(h1)
				EAP.rel <- v1 / h1			
						}	
	###################
	# person parameter estimates
	person <- data.frame( "EAP" = colMeans( theta.chain ) ,
						"SD" = colSds(theta.chain) )
	# output
	res <- list( "EAP.rel" = EAP.rel , "person" = person )
	return(res)
			}
###############################################################