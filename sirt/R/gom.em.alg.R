


###########################################################
# calculate theta grid
.gom.calc.theta <- function( K , problevels ){
	m1 <- problevels
	if ( ! is.matrix(problevels) ){
        PL <- length(problevels)	
		m1 <- matrix( problevels , PL , 1 )
		for (kk in 2:K){
			# kk <- 2
			NM <- nrow(m1)
			m1 <- cbind( m1[ rep( 1:NM , PL)  , ] , rep( problevels, each=NM)  )
			m1 <- m1[ rowSums(m1) <= 1 , ]
				}
			}
	m1 <- m1[ rowSums(m1) == 1 , ]
    return(m1)
        }
################################################################
# gom calcprobs
.gom.calcprobs <- function( lambda , theta.k , b=NULL , theta0.k=NULL ){
	if ( ! is.null( b ) ){
#		lambda <- t( plogis( outer( theta0.k , b , "-" )	) )
		lambda <- stats::plogis( - b + matrix( theta0.k , length(b) , length(theta0.k) ,
					byrow=TRUE )	)
					}
#    probs <- lambda %*% t(theta.k)
    probs <- base::tcrossprod( lambda ,theta.k)
	probsL <- array( 0 , dim=c( nrow(lambda) , 2 , nrow(theta.k) ) )
	probsL[,2,] <- probs
	probsL[,1,] <- 1-probs	
	res <- list("probs"=probs , "probsL"=probsL)
	return(res)
            }
#*************************************************************************************
# E Step GOM
.e.step.gom.em <- function( dat2 , dat2.resp , pi.k , probs , I , 
							theta.k  ){
    #...................................                    
		#***
		# array notation of probabilities
		pjk <- probs
		pjkL <- array( NA , dim=c(2 , nrow(pjk) , ncol(pjk) ) )
		pjkL[1,,] <- 1 - pjk
		pjkL[2,,] <- pjk	
		f.yi.qk <- matrix( 1 , nrow(dat2) , nrow(theta.k) )	
		for (ii in 1:ncol(dat2)){
		#	ii <- 1
			ind.ii <- which( dat2.resp[,ii] == 1 )
			f.yi.qk[ind.ii,] <- f.yi.qk[ind.ii,] * pjkL[ dat2[ind.ii,ii]+1 , ,ii]
						}
		#******
    f.qk.yi <- 0 * f.yi.qk
    pi.k <- matrix( pi.k , ncol=1 )
    f.qk.yi <- f.yi.qk * matrix( pi.k , nrow(f.yi.qk) , ncol(f.yi.qk) ,
				byrow=TRUE )
    f.qk.yi <- f.qk.yi / rowSums( f.qk.yi )		
	# expected counts at theta.k
    n.k <- matrix( 0 , length(theta.k) , 1 )
    r.jk <- n.jk <- array( 0 , dim=c( ncol(dat2) , nrow(theta.k) , 1) )
    ll <- rep(0,1)
    n.k <- colSums( f.qk.yi )
    # expected counts at theta.k and item j
#    n.jk[,,1] <-  t(dat2.resp)  %*% f.qk.yi
    n.jk[,,1] <-  base::crossprod(dat2.resp , f.qk.yi )
    # compute r.jk (expected counts for correct item responses at theta.k for item j
#    r.jk[,,1] <- t( dat2 * dat2.resp )  %*% f.qk.yi
    r.jk[,,1] <- base::crossprod( dat2 * dat2.resp , f.qk.yi )
	K <- 1
	# expected counts
	TP <- nrow(theta.k)
	n.ik <- array( 0 , dim =c(TP , I , K+1 ) )
	N.ik <- array( 0 , dim=c(TP , I ) )
	n.ik[,,1] <- t( n.jk[,,1] - r.jk[,,1] )
	n.ik[,,2] <- t( r.jk[,,1] )	
	for (kk in 1:(K+1) ){
#		# kk <- 1
#		n.ik[,,kk] <- t( t(dat2.resp*(dat2==(kk-1)) ) %*% f.qk.yi )
		N.ik <- N.ik + n.ik[,,kk]
						}	
    # compute log-Likelihood
    ll <- sum(  log( rowSums( f.yi.qk * matrix( pi.k , nrow(f.yi.qk) , ncol(f.yi.qk) ,
				byrow=TRUE ) ) ) )
    res <- list( "n.k" = n.k , "n.ik"=n.ik , "pi.k"= n.k / sum( n.k) , 
			"n.jk" = n.jk , "r.jk" = r.jk , "f.qk.yi" = f.qk.yi , "pjk" = pjk  ,
            "f.yi.qk" = f.yi.qk , "ll" = sum(ll) )
    return(res)
    }
#*************************************************************************************
# estimation of lambda parameters
.gom.est.lambda <- function( lambda , I , K , n.ik , 
		numdiff.parm=.001 , max.increment=1,theta.k , msteps ,
		mstepconv , eps = .001 , progress=TRUE){
    h <- numdiff.parm
	diffindex <- 1:I
	Q0 <- 0*lambda
	se.lambda <- Q0
	an.ik <- aperm( n.ik , c(2,3,1) )
	if (progress){ 
		cat("  M steps lambda parameter |")
					}
	it <- 0 ;	conv1 <- 1000
	while( ( it < msteps ) & ( conv1 > mstepconv ) ){	
		lambda0 <- lambda
		for (kk in 1:K){
	#		kk <- 1
			Q1 <- Q0
			Q1[,kk] <- 1			
			pjk <- .gom.calcprobs( lambda , theta.k )$probsL
			pjk1 <- .gom.calcprobs( lambda+h*Q1 , theta.k )$probsL
			pjk2 <- .gom.calcprobs( lambda-h*Q1 , theta.k )$probsL
			# numerical differentiation			
			res <- .gom.numdiff.index( pjk , pjk1 , pjk2 , an.ik , diffindex , 
                max.increment=max.increment , numdiff.parm )				
			increment <- Q1*matrix( res$increment , nrow=I , ncol=K)	
			lambda <- lambda + increment
			lambda[ lambda[,kk] < eps , kk ] <- eps
			lambda[ lambda[,kk] > 1-eps , kk ] <- 1 - eps			
			se.lambda[,kk] <- sqrt(-1/res$d2)		
					}
		conv1 <- max( abs( lambda - lambda0 ) )
		it <- it+1
		if (progress){ 
				cat("-") 
						}
			}
	if (progress){		
		cat(" " , it , "Step(s) \n")	
				}
	res <- list("lambda" = lambda , "se.lambda" = se.lambda , 
			"ll" = sum(res$ll0) )
    return(res)
			}
######################################################################
# estimation of b parameters
.gom.est.b <- function( lambda , I , K , n.ik , 
		b , theta0.k , numdiff.parm=.001 , max.increment,theta.k , msteps ,
		mstepconv , eps = .001 , progress=progress ){
    h <- numdiff.parm
	diffindex <- 1:I
#	Q0 <- 0*lambda
#	se.b <- Q0
    if (progress){
		cat("  M steps b parameter |")
				  }
	an.ik <- base::aperm( n.ik , c(2,3,1) )	
	it <- 0 ;	conv1 <- 1000
	while( ( it < msteps ) & ( conv1 > mstepconv ) ){	
		b0 <- b
		pjk <- .gom.calcprobs( lambda , theta.k , b , theta0.k )$probsL
		pjk1 <- .gom.calcprobs( lambda , theta.k , b+h , theta0.k )$probsL
		pjk2 <- .gom.calcprobs( lambda , theta.k , b-h , theta0.k)$probsL
		# numerical differentiation			
		res <- .gom.numdiff.index( pjk , pjk1 , pjk2 , an.ik , diffindex , 
               max.increment=max.increment , numdiff.parm )				
#		increment <- Q1*matrix( res$increment , nrow=I , ncol=K)	
		b <- b + res$increment
		#b <- b - mean(b)
		conv1 <- max( abs( b - b0 ) )
		it <- it+1
		if (progress){
			cat("-") 
				}
			}
	if (progress){ 
		cat(" " , it , "Step(s) \n")	
					}
	res <- list("b" = b , "se.b" = sqrt( - 1/res$d2 ) , 
			"ll" = sum(res$ll0) )
    return(res)
			}
######################################################################



####################################################################
# general function for numerical differentiation
# diffindex aggregates across super items
.gom.numdiff.index <- function( pjk , pjk1 , pjk2 , an.ik , diffindex , 
		max.increment , numdiff.parm , eps=10^(-20) ){					
	h <- numdiff.parm
    ll0 <- rowSums( an.ik * log(pjk+eps) )
    ll1 <- rowSums( an.ik * log(pjk1+eps) )
    ll2 <- rowSums( an.ik * log(pjk2+eps) )	
    d1 <- ( ll1 - ll2  ) / ( 2 * h )    # negative sign?
    # second order derivative
    # f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
    d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2
    # change in item difficulty
    d2[ abs(d2) < 10^(-10) ] <- 10^(-10)
    increment <- - d1 / d2
	ci <- ceiling( abs(increment) / ( abs( max.increment) + 10^(-10) ) )
    increment <- ifelse( abs( increment) > abs(max.increment)  , 
                                 increment/(2*ci) , increment )	
#    increment <- ifelse( abs( increment) > max.increment , 
#            max.increment*sign(increment) , increment ) 
	res <- list("increment"=increment , "d2"=d2 , "ll0"=ll0)
	return(res)
		}
#################################################################
# estimation of covariance
.gom.est.covariance <- function( f.qk.yi , Sigma , theta.kM , N  ){
		Sigma.cov <- Sigma
		delta.theta <- 1
		hwt <- f.qk.yi	
		theta.k <- theta.kM
		thetabar <- hwt%*%theta.k
		# calculation of mu
		mu <- colSums( thetabar ) / N
		mu[1] <- 0
		# calculation of the covariance matrix
		theta.k.adj <- theta.k - matrix( mu , nrow=nrow(theta.k) , 
									ncol=ncol(theta.k) , byrow=TRUE)
		D <- 2
		for (dd1 in 1:D){
			for (dd2 in dd1:D){
				tk <- theta.k.adj[,dd1]*theta.k.adj[,dd2]
				h1 <- ( hwt %*% tk ) * delta.theta
				Sigma.cov[dd1,dd2] <- sum( h1 ) / N
				if (dd1 < dd2 ){ Sigma.cov[dd2,dd1] <- Sigma.cov[dd1,dd2] }
									}
								}
		diag(Sigma.cov) <- diag(Sigma.cov) + 10^(-10)
		pi.k <- matrix( mvtnorm::dmvnorm( theta.k , mean = mu , sigma = Sigma.cov )	, 
						ncol=1 )		
		pi.k <- pi.k / sum( pi.k )			
		res <- list("mu"=mu , "Sigma"=Sigma.cov , "pi.k"= pi.k )				
		return(res)
					}
###########################################################