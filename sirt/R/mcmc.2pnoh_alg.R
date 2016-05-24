

######################
# draw latent responses Z
.draw.Z.2pnoh <- function( aM , bM, theta , N , I , threshlow , threshupp ){
    # calculate means
    mij <- aM * ( theta - bM )
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
.draw.theta.2pnoh <- function( aM , bM , N , I , Z ){
    vtheta <- 1 / ( rowSums( aM^2 ) + 1  )
    mtheta <- rowSums( aM * ( Z + aM*bM ) ) * vtheta
    theta <- stats::rnorm( N , mean=mtheta , sd = sqrt( vtheta ) )
	v1 <- stats::var(mtheta)
    res <- list("theta" = theta )
    return(res)
            }

##########################################
# draw item parameters in the 2PNOH model
.draw.itempars.2pnoh <- function( theta , Z , I , N , a , b , 
		xi , omega , sig , nu , itemgroup , K , Ik ,  weights){

	aM <- matrix( a , nrow=N , ncol=I , byrow=TRUE )	

	#************************
	# sampling beta parameters
	v.ik <- 1 / ( N * a^2 + 1 / sig^2 )
	if ( is.null(weights ) ){
		m.ik <- a * colSums( aM * theta - Z  ) + xi[ itemgroup ] / sig^2
					}
	if ( ! is.null(weights ) ){
		m.ik <- a * colSums( weights*(aM * theta - Z)  ) + xi[ itemgroup ] / sig^2
					}					
    m.ik <- m.ik * v.ik
	b <- stats::rnorm( I , mean = m.ik , sd = sqrt(v.ik ) )
	bM <- matrix( b , nrow=N , ncol=I , byrow=TRUE )		

	#************************
	# sampling alpha parameters
	v.ik <- 1 / ( colSums( ( theta - bM  )^2 ) + 1 /nu^2 )		
	if ( is.null(weights ) ){
		m.ik <- colSums( ( theta - bM) *Z  ) + omega[ itemgroup ] / nu^2
					}
	if ( ! is.null(weights ) ){
		m.ik <- colSums( weights*( theta - bM) *Z  ) + omega[ itemgroup ] / nu^2
					}
    m.ik <- m.ik * v.ik
	a <- stats::rnorm( I , mean = m.ik , sd = sqrt(v.ik ) )
	
	#*********************************
	# sampling of xi parameter	
	m.xi <- stats::aggregate( b , list(itemgroup) , mean )[,2]
	xi <- stats::rnorm(K , mean=m.xi , sd = sqrt( sig^2 / Ik) )

	#*********************************
	# sampling of omega parameter	
	m.xi <- stats::aggregate( a , list(itemgroup) , mean )[,2]
	omega <- stats::rnorm(K , mean=m.xi , sd = sqrt( nu^2 / Ik) )
	
	# output
    return( list( "a"=a , "b"=b , "xi"=xi , "omega" = omega ) )
            }
#################################################
# compute deviance
.mcmc.deviance.2pnoh <- function( aM , bM , theta , dat , dat.resp , 
				weights , eps ){
	pij <- stats::pnorm( aM * theta - bM )
	llij <- log( dat.resp * ( dat*pij + ( 1-dat )*(1-pij) ) + eps )
	if ( is.null( weights ) ){ deviance <- -2*sum( llij ) }
	if ( ! is.null( weights ) ){ 
		deviance <- -2*sum( rowSums(llij) * weights ) 
					}
	return(deviance)
				}
####################################################
# draw item variances from inverse chi square distributions
.draw.itemvariances.2pnoh <- function(a,b,I , itemgroup , xi , omega , prior.variance ){		
		# inverse gamma distribution
		# n <- 20
		# sig2 <- 1
		# 1 / stats::rgamma( 1 , n/2 , n*sig2/2 )
		a0 <- prior.variance[1]
		b0 <- prior.variance[2]
		#****
		# sample sig2
#		tau2 <- ( sum( ( b - xi[ itemgroup ] )^2 ) + 1 )
#		p1 <- 1 / rchisq( 1 , df= I + 1 )
#		sig <- sqrt( tau2 * p1 )	
		sig <- sqrt( 1 / stats::rgamma( 1 , (I+a0)/2 , ( sum( ( b - xi[ itemgroup ] )^2 ) + b0 ) / 2 )	)	
		
		#****
		# sample sig2
#		tau2 <- ( sum( ( a - omega[ itemgroup ] )^2 ) + 1 )
#		p1 <- 1 / rchisq( 1 , df= I + 1 )
#		nu <- sqrt( tau2 * p1 )		
		nu <- sqrt( 1 / stats::rgamma( 1 , (I+a0)/2 , ( sum( ( a - omega[ itemgroup ] )^2 ) + b0 ) / 2 )	)
		res <- list("sig" = sig , "nu" = nu )
		return(res)
		}
#################################################
# calculation of mastery probabilities		
.mcmc.mastery.2pnoh <- function( xi , omega , N , K , weights , theta , prob.mastery ){
	h1 <- stats::pnorm( matrix( omega , nrow=N , ncol=K , byrow=TRUE) * 
			( theta - matrix( xi , nrow=N , ncol=K , byrow=TRUE) ) )
	if ( is.null(weights)){
		nonmastery <- colMeans( h1 < prob.mastery[1] )
		transition <- colMeans( ( h1 < prob.mastery[2] ) * (h1> prob.mastery[1] ) )
		mastery <- colMeans( h1 > prob.mastery[2] )
				}
	if ( ! is.null(weights)){
		nonmastery <- colMeans( weights*(h1 < prob.mastery[1] ) )
		transition <- colMeans( weights*( h1 < prob.mastery[2] ) * (h1> prob.mastery[1] ) )
		mastery <- colMeans( weights*(h1 > prob.mastery[2] ) )
				}			
	res <- list( "nonmastery" = as.vector(nonmastery) , "transition" = as.vector(transition) , 
			"mastery" = as.vector(mastery) )
	return(res)
		}		