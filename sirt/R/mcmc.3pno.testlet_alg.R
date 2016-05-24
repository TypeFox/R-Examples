
##############################################
# draw latent response data W
.draw.W.3pno.testlet <- function( aM , bM, theta , gamma.testlet ,
				N , I , threshlow , threshupp , testletgroups , param , 
				dat , guess , tau.ni ){	
	# P(Y=1) = c_i + ( 1 - c_i ) Phi ( tau )				
	#		 = Phi(tau) + c_i ( 1 - Phi(tau) )
    # calculate probabilities P(W=1|Y=1)
    p1 <- stats::pnorm( tau.ni )
	# calculate P(W=0|Y=1)
	p0 <- matrix( guess , N , I , byrow=TRUE ) * ( 1 - p1 )
	p1 <- p1 / ( p1 + p0 )
	# simulate random draw
    rij <- matrix( stats::runif( N*I ) , nrow=N , ncol=I )
    # calculate corresponding value
	W <- 1 * ( p1 > rij  ) * dat
    return(W)
        }

###############################################
# draw latent responses Z
.draw.Z.3pno.testlet <- function( aM , bM, theta , gamma.testlet ,
				N , I , threshlow , threshupp , testletgroups , 
				param , W , dat.resp , tau.ni){				
	# calculate indicator testletgroups
	gamma.testletM <- gamma.testlet[ , testletgroups ]	
    # calculate means
#    mij <- aM * theta + gamma.testletM + bM
	mij <- tau.ni
    # simulate uniform data
    rij <- matrix( stats::runif( N*I ) , nrow=N , ncol=I )
	# calculate thresholds if they should calculated
		# define lower and upper thresholds 
		ZZ <- 1000
		threshlow <- -ZZ + ZZ*W
		threshlow <- threshlow - (1-dat.resp)*ZZ
		threshupp <- ZZ*W
		threshupp <- threshupp + (1-dat.resp)*ZZ
    # calculate corresponding value
    pl <- stats::pnorm( threshlow , mean=mij) 
    pu <- stats::pnorm( threshupp , mean=mij)
    pij <- pl + (pu-pl)*rij
    # simulate Z
    Zij <- stats::qnorm( pij , mean = mij )
    return(Zij)
        }
##################################################		
		
#############################################
# draw theta from 3PNO testlet model
.draw.theta.3pno.testlet <- function( aM , bM , N , I , Z , param ,
		gamma.testletM , sigma , a.testletM ){
    if (param==1){ V <- Z - gamma.testletM }
    if (param==2){ V <- Z - aM*gamma.testletM }	
    if (param==3){ V <- Z - a.testletM*gamma.testletM }		
    vtheta <- 1 / ( rowSums( aM^2 ) + 1 / sigma^2  )
    mtheta <- rowSums( aM * V ) * vtheta
    theta <- stats::rnorm( N , mean=mtheta , sd = sqrt( vtheta ) )
	theta <- theta - mean( theta )
	if (param==3){ theta <- theta / stats::sd( theta ) }
    return(theta)
            }
#############################################
	
#############################################
# draw gamma from 3PNO testlet model
.draw.gamma.3pno.testlet <- function( aM , bM , N , I , Z , param ,
		theta , testletgroups , TT , N.items , gamma.testlet , 
		sigma.testlet , a.testletM){
	Zast <- Z - aM * theta - bM	
	for (tt in 1:TT){		
#		tt <- 1
		ind.tt <- which( testletgroups==tt)
		#*** param=1
		if (param==1){
			Zast.tt <- Zast[ , ind.tt ]
			nu.tt <- 1 / N.items[tt]
			vtheta <- 1 / ( 1 / nu.tt + 1 / sigma.testlet[tt]^2 ) 
			mtheta <- rowMeans( Zast.tt ) * 1 / nu.tt * vtheta
					}
		#*** param=2
		if (param==2){
			Zast.tt <- Zast[ , ind.tt ]
			nu.tt <- rowSums( aM[,ind.tt]^2 )
			vtheta <- 1 / ( 1 / nu.tt + 1 / sigma.testlet[tt]^2 ) 
			mtheta <- rowSums( aM[,ind.tt] * Zast.tt ) * 1 / nu.tt * vtheta
					}					
		#*** param=3
		if (param==3){
			Zast.tt <- Zast[ , ind.tt ]
			nu.tt <- rowSums( a.testletM[,ind.tt]^2 )
			vtheta <- 1 / ( 1 / nu.tt + 1 / sigma.testlet[tt]^2 ) 
			mtheta <- rowSums( a.testletM[,ind.tt] * Zast.tt ) * 1 / nu.tt * vtheta
					}
		# sample gamma parameter					
		gamma.testlet[,tt] <- stats::rnorm(N , mean=mtheta , sd = sqrt(vtheta) )
		gamma.testlet[,tt] <- gamma.testlet[,tt] - mean( gamma.testlet[,tt] )
		
#		gamma.testlet[ , 1:TT] <- gamma.testlet[ , 1:TT ] - rowMeans( gamma.testlet[,1:TT] )		
		if (param==3){ 
     		gamma.testlet[,tt] <- gamma.testlet[,tt] / stats::sd( gamma.testlet[,tt] )		
					}			
			}	
    return(gamma.testlet)
            }
#############################################	

			
		
##########################################
# draw item parameters a and b 
.draw.itempars.3pno.testlet <- function( theta , Z , I , N , weights ,
	gamma.testlet , testletgroups , param , TT , a.testletM){
	# define adjusted Z values
	gamma.testletM <- gamma.testlet[ , testletgroups ]
# print( round( a.testletM[1,],2))	
	if (param==1){ Z <- Z - gamma.testletM }
	if (param==3){ Z <- Z - a.testletM*gamma.testletM }	
	if (param==2){ # Z <- Z 
	#	gamma.testletM <- gamma.testlet[ , testletgroups ]
#		theta <- theta + gamma.testletM		
		theta0 <- theta
		Z0 <- Z
				}
	# for parametrization 2, this function must be rewritten
	# because "the theta" is now item specific
	# loop over testlets tt=1,...,TT
	# maybe for TT+1 some adjustment has to be done
	#''''''''''''''''''''''''''''''''''''''''
	# parametrization param=1
	if ( (param==1) | (param==3) ){			
		#--------------
		# sampling without weights
		Xast <- as.matrix( cbind( theta , 1 ) )
		if ( is.null(weights) ){
#			Sigma <- solve( t(Xast) %*% Xast )
			Sigma <- solve( crossprod(Xast) )
			# calculate mean
#			mj <- Sigma %*% t(Xast) %*% Z
			mj <- Sigma %*% crossprod( Xast , Z )
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
#			mj <- Sigma %*% t( Xast * weights ) %*% Z
			mj <- Sigma %*% crossprod( Xast * weights , Z )
			mj <- as.matrix( t(mj))	
				}		
		#--------------							
		# draw item parameters
		ipars <- mvtnorm::rmvnorm( I , sigma=Sigma ) + mj
		a <- ipars[,1]
		b <- ipars[,2]
				}
	#''''''''''''''''''''''''''''''''''''''''
	# parametrization param=2
	if (param==2){
		a <- rep(NA,I)
		b <- rep(NA,I)
		TTT <- TT
		if ( sum( testletgroups== TT+1 ) > 0 ){
				TTT <- TT + 1 }
		for (tt in 1:TTT){
			#tt <- 1
			theta <- theta0
			Z <- Z0
			ind.tt <- which( testletgroups== tt)
			Itt <- length(ind.tt)		
			theta <- theta0 + gamma.testlet[ , tt]
			Z <- Z[ , ind.tt , drop=FALSE]
			#--------------
			# sampling without weights
			Xast <- as.matrix( cbind( theta , 1 ) )
			if ( is.null(weights) ){
#				Sigma <- solve( t(Xast) %*% Xast )
				Sigma <- solve( crossprod(Xast) )
				# calculate mean
#				mj <- Sigma %*% t(Xast) %*% Z
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
				mj <- Sigma %*% t( Xast * weights ) %*% Z
				mj <- as.matrix( t(mj))	
					}		
			#--------------							
			# draw item parameters
			ipars <- mvtnorm::rmvnorm( Itt , sigma=Sigma ) + mj
			a[ind.tt] <- ipars[,1]
			b[ind.tt] <- ipars[,2]
					}			# end testlet tt
				} # end param=2
	#******************************
    return( list( "a"=a , "b"=b) )
            }
############################################################

#############################################
# draw b from 3PNO testlet model if slopes are fixed at one
.draw.est.b.3pno.testlet <- function( aM , bM , N , I , Z , param ,
		gamma.testletM , sigma , weights , theta ){		
	#************************
	V <- Z - aM * theta - gamma.testletM			
	if ( is.null(weights) ){
		mb <- colMeans( V )
		vb <- 1 / N
			}
	if ( ! is.null(weights) ){
		mb <- colSums( V * weights ) / N
		vb <- sum( weights^2 ) / N^2
			}		
    b <- stats::rnorm( I , mean=mb , sd = sqrt( vb ) )
    return(b)
            }
#############################################


#############################################
# draw a.testlet from 3PNO testlet model
.draw.est.a.testlet.3pno.testlet <- function( aM , bM , N , I , Z , param ,
		gamma.testletM , sigma , weights , theta , testletgroups ,
		gamma.testlet , TT	){		
	#************************
	# Z = a * theta + a.testlet *gamma + bM + eps
	# Z - a*theta - bM = a.testlet * gamma.testletM + esp
	
	V <- Z - aM * theta - bM
	a.testlet <- rep(0,I)
	for (tt in 1:TT){
	# tt <- 1	
		ind.tt <- which( testletgroups == tt )
		Itt <- length(ind.tt)
		Vtt <- V[,ind.tt]			
		v.ik <- 1 / ( colSums( Vtt^2 ) )		
		if ( is.null(weights ) ){
			m.ik <- colSums(  gamma.testlet[,tt] * Vtt  )
						}
		if ( ! is.null(weights ) ){
			m.ik <- colSums( weights*gamma.testlet[,tt] * Vtt  )
						}						
		m.ik <- m.ik * v.ik	
		a.testlet[ind.tt] <- stats::rnorm( Itt , mean = m.ik , sd = sqrt(v.ik ) )
				}
# cat(".....\n")					
    return(a.testlet)
            }
#############################################


##########################################
# draw guessing parameter
.draw.guess.3pno.testlet <- function( guess.prior , W , 
		dat , dat.resp , weights , I ){
	if ( is.null(weights) ){
		ri <- colSums( (W==0) )
		si <- colSums( (W==0) * (dat==1) * dat.resp)
				}
	if ( ! is.null(weights) ){
		ri <- colSums( (W==0) * weights )
		si <- colSums( (W==0) * (dat==1) * dat.resp * weights)
				}	
	guess <- stats::rbeta(I , si + guess.prior[,1] , ri-si+guess.prior[,2] )		
    return( guess )
            }
#################################################			
# draw theta variance
.draw.theta.variance.3pno.testlet <- function( theta , weights , N ){
		if ( is.null(weights) ){  sig2 <- stats::var(theta) }
		if ( ! is.null(weights) ){
			sig2 <- ( sum( theta^2 * weights ) - ( sum( theta*weights )/N )^2 )/ N
						}				
		sigma <- sqrt( .mcmc.draw.variance( 1 , w0=1 , sig02=1 , n=N , sig2=sig2 ) )						
		return(sigma)
			}			
############################################################
# draw testlet variance
.draw.testlet.variance.3pno.testlet <- function( gamma.testlet , N ,
	sigma.testlet , testlet.variance.prior , weights , TT){
	for (tt in 1:TT){
		if ( is.null(weights)){ sig2 <- sum( gamma.testlet[,tt]^2 ) / N }
		if ( ! is.null(weights)){ sig2 <- sum( weights*gamma.testlet[,tt]^2 ) / N }		
		sigma.testlet[tt] <- sqrt( .mcmc.draw.variance( 1 , w0=testlet.variance.prior[1] , 
				sig02=testlet.variance.prior[2] , n=N , sig2=sig2 ) )
					}
	return(sigma.testlet)
			}
######################################################
	