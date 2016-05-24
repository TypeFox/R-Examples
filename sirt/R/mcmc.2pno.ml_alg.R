
#***********************************************
# parts of algorithm:
# * draw Z
# * draw theta2
# * draw theta1 (centering!)
# * draw sigma2
# * draw sigma1
# * draw b_class
# * draw b
# * draw variance(b_class) 
# * draw a_class
# * draw a
# * draw variance(a_class)
#********************************************

#############################################
# compute group mean by the rowsum function
.mcmc.groupmean <- function( matr , group , groupsize=NULL ){
    r1 <- rowsum( matr , group )
    if ( is.null(groupsize) ){ 
			groupsize <- rowsum( 1+0*matr[,1] , group )[,1] 
					}
    r1 / groupsize
        }
############################################

####################################################
# draw latent responses Z
.draw.Z.2pno.ml <- function( aM , bM, theta , N , I , threshlow , threshupp ){
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
#############################################

########################################
# draw theta = theta.L2 + theta.L1
# here, the "total" theta is sampled!
.draw.theta.2pno.ml <- function( aM , b , bM , N , I , Z ,
	sigma1 , sigma2 , sigma.res , link , theta2 , idgroup ){
	#**************************************
	# link = logit
	# Z = aM * theta - bM + eps
	# Z + bM = aM * theta + eps
#	bM <- matrix(b,N,I,byrow=TRUE )
	if (link=="logit"){
#		vtheta <- 1 / ( rowSums( aM^2 ) + ( sigma1^2+sigma2^2 )  )
#		vtheta <- 1 / ( rowSums( aM^2 ) + 1/( sigma1^2+sigma2^2 )  )
#		mtheta <- rowSums( aM * ( Z + bM ) ) * vtheta
		Zres <- Z + bM
		m1ij <- rowSums( aM * Zres )
		prec1 <- rowSums(aM^2 )
		m1ij <- m1ij / prec1
		m2ij <- theta2[ idgroup ]
#		prec2 <- 1 / ( sigma1^2 + sigma2^2 )	
		prec2 <- 1 / ( sigma1^2)	
					} # end logit
	#***************************************
	# link = normal
	if (link == "normal"){
# print(sigma.res)	
		sigma.resM <- matrix( sigma.res , N , I , byrow=TRUE)
		Zres <- Z + bM
		# correct precision and estimator
		# something seems to go wrong here
#		m1ij <- rowSums( aM * Zres * sigma.resM)
		m1ij <- rowSums( aM * Zres / sigma.resM^2 )
		m2ij <- theta2[ idgroup ]
#		prec1 <- rowSums( aM^2 * sigma.resM^2 )
		prec1 <- rowSums( aM^2 / sigma.resM^2 )
		m1ij <- m1ij / prec1
		prec2 <- 1 / ( sigma1^2 )
				} # end normal
	prectotal <- prec1 + prec2	
	mtheta <- ( prec1*m1ij + prec2 * m2ij ) / prectotal
	vtheta <- 1 / prectotal				
    theta <- stats::rnorm( N , mean=mtheta , sd = sqrt( vtheta ) )
	theta <- theta - mean(theta)
    return(theta)
            }
#########################################
# draw level 2 latent class mean
.draw.theta2.2pno.ml <- function( theta , idgroup , groupsize , 
		sigma1 , sigma2 , G ){
		# compute latent mean
		mij1 <- .mcmc.groupmean( theta , idgroup , groupsize)[,1]
		prec.ij1 <- 1 / ( sigma1^2 / groupsize )
		prec.ij2 <- 1 / sigma2^2
		prec.tot <- prec.ij1 + prec.ij2
		vtheta <- 1 / prec.tot
		mtheta <- mij1 * prec.ij1 / prec.tot
		theta2 <- stats::rnorm( G , mean=mtheta , sd = sqrt( vtheta ))
		theta2 <- theta2 - mean(theta2)
		return(theta2)
			}
######################################

############################################################
# draw level 1 and level 2 variances
.draw.sigma12.2pno.ml <- function( theta , theta2 , idgroup , N ,
		G , prior.sigma2 ){
		# level 1 variance
		sig2 <- sum( (theta - theta2[ idgroup ])^2 ) / N
		sigma1 <- sqrt( .mcmc.draw.variance( 1 , w0= 1 , sig02= .7 , n=N , sig2=sig2 ) )		
		# level 2 variance
		sig2 <- sum( theta2^2 ) / G
		sigma2 <- sqrt( .mcmc.draw.variance( 1 , 
					w0= prior.sigma2[1] , sig02= prior.sigma2[2]^2 , n=G , sig2=sig2 ) )				
		res <- list( "sigma1"=sigma1 , "sigma2" = sigma2 )
		return(res)
			}
######################################################
# b item parameter single level case
.draw.est.b.sl <- function(	Z , aM , theta , N , I , omega.b , mu.b ,
		sigma.res , link ){			
		#     Z = aM * theta - bM + eps
		# =>  bM = Z - aM*theta - eps
	Zres <- Z - aM * theta
	mij1 <- -colMeans(Zres)				
	#****
	# link = logit
	if (link=="logit"){
		prec1ij <- N
		prec2ij <- 1 / omega.b^2
		prectotal <- prec1ij + prec2ij	
		mij <- ( mij1 * prec1ij + mu.b * prec2ij ) / prectotal
		vij <- 1 / prectotal
				}
	#***
	# link = normal
	if (link=="normal"){
		prec1ij <- N / sigma.res^2
		prec2ij <- 1 / omega.b^2
		prectotal <- prec1ij + prec2ij	
		mij <- ( mij1 * prec1ij + mu.b * prec2ij ) / prectotal
		vij <- 1 / prectotal	
				} # end link=normal
	#*** sample b parameter
	b <- stats::rnorm( I , mean=mij , sd = sqrt( vij ) )
	return(b)
			}
####################################################
# draw hyperparameters of b
.draw.est.b.hyperpars <- function( b, mu.b , omega.b , I , 
	prior.omega.b , est.b.M  ){		
		if ( est.b.M=="h"){
			# sample mu.b
			mij <- mean(b)
			vij <- omega.b^2 / I
			mu.b <- stats::rnorm(1 , mean=mij , sd = sqrt(vij) )
			# sample omega.b
			sig2 <- sum( ( b - mu.b )^2 ) / I
			omega.b <- sqrt( .mcmc.draw.variance( 1 , 
						w0= prior.omega.b[1] , sig02= prior.omega.b[2]^2 , 
						n=I , sig2=sig2 ) )				
						}
		res <- list( "mu.b"=mu.b , "omega.b" = omega.b )
		return(res)
			}

#####################################################
# sampling of b in case of multilevel DIF
.mcmc.est.b.2pno.ml.v2 <- function( N , Z , aM , theta , idgroup , groupsize , 
							b , bG , G , I , sigma.b , omega.b , mu.b , sigma.res , link ){			
			# sampling of b (fixed item difficulties)
			# Z = aM *  theta - b - bG + eps
			# Z - aM * theta + bG = - b + eps			
			Zres <- Z - aM * theta + bG[ idgroup , ]
#			Zres <- Z - aM * theta			
			mij1 <- -colMeans(Zres)				
			#****
			# link = logit
			if (link=="logit"){
				prec1ij <- N
				prec2ij <- 1 / omega.b^2
				prectotal <- prec1ij + prec2ij	
				mij <- ( mij1 * prec1ij + mu.b * prec2ij ) / prectotal
				vij <- 1 / prectotal
							}
			#***
			# link = normal
			if (link=="normal"){
				prec1ij <- N / sigma.res^2
				prec2ij <- 1 / omega.b^2
				prectotal <- prec1ij + prec2ij	
				mij <- ( mij1 * prec1ij + mu.b * prec2ij ) / prectotal
				vij <- 1 / prectotal	
						} # end link=normal
			#*** sample b parameter
			b <- stats::rnorm( I , mean=mij , sd = sqrt( vij ) )			
			return(b)
				}
####################################################################
			
			
#########################################################
# estimation of b for clusters
.mcmc.est.b.group.2pno.ml <- function( Z , aM , theta , idgroup , groupsize , 
		b , G , I , sigma.b , sigma.res , link  , N  ){
		#*****		
		# Z = a * theta - b - bG + eps
		# => - b =  Z - a * theta + bG + eps
		bM1 <- matrix( b , nrow=N , ncol=I , byrow=TRUE )
		Zres <- Z - aM * theta + bM1
		# compute means
		mij1 <- - .mcmc.groupmean( matr=Zres , group=idgroup , groupsize )		
		mij2 <- matrix( 0 , G , I , byrow=TRUE )
		# compute precisions
		prec1 <- matrix( groupsize  , G , I , byrow=FALSE) / 1	
		if (link == "normal"){						
			prec1 <- prec1 / matrix( sigma.res^2 , G , I , byrow=TRUE )							
							}
		prec2 <- 1 / matrix( sigma.b^2 , G , I , byrow=TRUE )						
		prectot <- prec1 + prec2
		# compute total means
		mtot <- ( mij1*prec1 + mij2*prec2 ) / prectot	
		vtot <- 1 / prectot
		# sampling of bG
		bG <- matrix( stats::rnorm( G*I , mean=mtot  , sd = sqrt(vtot) ) , G , I )
		# adjustment
#		bG1 <- rowMeans( bG )
#		bG1 <- bG1 - mean(bG1)
#		bG <- bG - bG1
#		bG1 <- colMeans( bG )
		bG <- as.matrix( base::scale( bG , scale=FALSE ) )
		return(bG)
			}
##################################################################
# estimation of hierarchical distribution
.mcmc.sigma.b.2pno.ml <- function( bG , mu.b , omega.b , G , I , 
			est.b.Var , prior.sigma.b , sigma.b ){
		#*****
		# draw item group standard deviations
		bresG <- bG
		if ( est.b.Var == "i"){ 
			sig2b <- colSums( bresG^2 ) / G	
			sigma.b <- sqrt( .mcmc.draw.variance( I , 
							w0= prior.sigma.b[1] , sig02= prior.sigma.b[2]^2 , 
							n=G , sig2=sig2b ) )
								}
		if ( est.b.Var == "j"){ 
			sig2b <- sum( bresG^2 ) / (G*I)
			sigma.b <- sqrt( .mcmc.draw.variance( 1 , 
							w0= prior.sigma.b[1] , sig02= prior.sigma.b[2]^2 , 
							n=G*I , sig2=sig2b ) )	
			sigma.b <- rep( sigma.b , I )							
								}
		return(sigma.b)
			}
##################################################################
# sampling of a parameters
.draw.est.a.sl <- function( Z , bM , theta , mu.a , omega.a , I ,
		sigma.res , link ){		
		# Z = a * theta - b + eps
		# => Z + b = a*theta + eps
		# a is obtained as a regression estimate
		eps <- .01
		Zres <- Z + bM		
		h1 <- sum( theta^2 )
		h2 <- colSums( Zres * theta )
		# calculate means
		m1ij <- h2 / h1
		m2ij <- mu.a
		# calculate precisions
		prec1 <- h1 * 1  # (X'X)^(-1) * sigma^2_{res}
		if ( link == "normal"){
			prec1 <- prec1 / sigma.res^2		
						}
		prec2 <- 1 / omega.a^2
		prectotal <- prec1 + prec2		
		# define mean and variance of posterior
		m1 <- ( m1ij * prec1 + m2ij * prec2 ) / prectotal
		# sampling of a
		a <- stats::rnorm( I , mean=m1 , sd = 1/sqrt(prectotal) )
#		a <- a - ( mean(a) - 1 )	
		a[ a < 0 ] <- eps
		a <- exp( log(a) - mean( log( a ) ) )
#		a <- a / prod(a)
		return(a)
			}
######################################################################			
# draw hyperparameters of a
.draw.est.a.hyperpars <- function( a, mu.a , omega.a , I , 
	prior.omega.a , est.a.M  ){		
		if ( est.a.M=="h"){
			# set mu.a to one
			mu.a <- 1
			# sample omega.a
			sig2 <- sum( ( a - mu.a )^2 ) / I
			omega.a <- sqrt( .mcmc.draw.variance( 1 , 
						w0= prior.omega.a[1] , sig02= prior.omega.a[2]^2 , 
						n=I , sig2=sig2 ) )				
						}
		res <- list( "mu.a"=mu.a , "omega.a" = omega.a )
		return(res)
			}
###################################################################		
# sampling of a parameters
.mcmc.a.est.a.2pno.ml <- function( Z , bM , aG , idgroup , theta , 
				mu.a , omega.a , I , link , sigma.res ){
		# Z = a * theta + aG*theta - bM + eps
		# => Z + bM - aG * theta = a*theta + eps
		# a is obtained as a regression estimate
		Zres <- Z + bM - aG[ idgroup , ] * theta		
		h1 <- sum( theta^2 )
		h2 <- colSums( Zres * theta )
		# calculate means
		m1ij <- h2 / h1
		m2ij <- mu.a
		# calculate precisions
		prec1 <- h1 * 1  # (X'X)^(-1) * sigma^2_{res}
		if (link=="normal"){
			prec1 <- prec1 / sigma.res^2
							}
		prec2 <- 1 / omega.a^2
		prectotal <- prec1 + prec2
		# define mean and variance of posterior
		m1 <- ( m1ij * prec1 + m2ij * prec2 ) / prectotal
		# sampling of a
		a <- stats::rnorm( I , mean=m1 , sd = 1/sqrt(prectotal) )
		a <- a - ( mean(a) - 1 )				
		return(a)
			}
###################################################################				
			
			
			
########################################################################
# draw a parameters group wise
.mcmc.est.aG.2pno.ml.v2 <- function( Z , bM , theta , idgroup , G , I ,
			a , sigma.a , N , link , sigma.res  ){
		#****
		# Z = a * theta + aG*theta - bM + eps
		# Z + bM - a*theta = aG * theta + eps		
		aM1 <- matrix( a , N , I , byrow=TRUE )
		Zres <- Z + bM	- aM1*theta
		# calculate means of a parameters
		theta2l <- rowsum( theta^2 , idgroup )[,1]
		Zrestheta <- rowsum( Zres*theta , idgroup )	
		m1ij <- Zrestheta / theta2l
		m2ij <- matrix( 0 , G , I , byrow=TRUE )		
		# calculate precisions
		prec1 <- matrix( theta2l , G , I )
		if (link=="normal"){
			prec1 <- prec1 / sigma.res^2 
							}		
		# take sigma.res into account!!
		prec2 <- matrix( 1 / sigma.a^2 , G , I , byrow=TRUE )
		prectotal <- prec1 + prec2 
		m1 <- ( m1ij*prec1 + m2ij * prec2 ) / prectotal 		
		aG <- matrix( stats::rnorm(G*I , mean = m1 , sd = sqrt( 1 / prectotal )) , G , I )
		# center aG parameters within each group
#		aG <- aG - ( rowMeans( aG ) - 0 )
		aG <- scale( aG , scale=FALSE)
		return(aG)
			}
#########################################################					
# sampling from hierarchical a distribution
.mcmc.a.grouphier.2pno.ml <- function( aG , mu.a , G , omega.a , I ,
		prior.sigma.a , est.a.Var , sigma.a ){
		#***
		# draw item group standard deviations
		aresG <- aG
		if ( est.a.Var == "i"){ 
			sig2b <- colSums( aresG^2 ) / G	
			sigma.a <- sqrt( .mcmc.draw.variance( I , 
							w0= prior.sigma.a[1] , sig02= prior.sigma.a[2]^2 , 
							n=G , sig2=sig2b ) )
								}								
		if ( est.a.Var == "j"){ 
			sig2b <- sum( aresG^2 ) / (G*I)
			sigma.a <- sqrt( .mcmc.draw.variance( 1 , 
							w0= prior.sigma.a[1] , sig02= prior.sigma.a[2]^2 , 
							n=G*I , sig2=sig2b ) )	
			sigma.a <- rep( sigma.a , I )							
								}
		return(sigma.a)
			}
####################################################################			
# draw residual standard deviations
.draw.sigma.res.2pno.ml <- function( Z , aM , bM , theta , N , I ){		
		# Z = a * theta - b + eps
		Zres <- Z - aM * theta + bM
		sig2 <- colSums( Zres^2 ) / N
		sigma.res <- sqrt( .mcmc.draw.variance( I , 
							w0= .001 , sig02= 1, 
							n=N , sig2=sig2 ) )
		return(sigma.res)
			}	
