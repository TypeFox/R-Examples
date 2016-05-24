
############################################
# probability in noncompensatory model
## extern "C" {
## SEXP SMIRT_CALCPROB_NONCOMP( SEXP a, SEXP b, SEXP Q, SEXP thetak, SEXP cc, SEXP dd) ;
calcprob.noncomp <- function (a,b,Q,thetak,cc,dd){ 
	.Call("SMIRT_CALCPROB_NONCOMP", a,b,Q,thetak,cc,dd, PACKAGE = "sirt")
					}
#######################################
# calculation of posterior distribution
## extern "C" {
## SEXP SMIRT_CALCPOST( SEXP dat2, SEXP dat2resp, SEXP probs, SEXP dat2ind, SEXP pik, SEXP K) ;
calcpost <- function (dat2 , dat2resp , probs, dat2ind , pik , K){ 
	.Call("SMIRT_CALCPOST", dat2 , dat2resp , probs, dat2ind , pik , 
				K, PACKAGE = "sirt")
					}
######################################
#calcpost2 <- function (dat2 , dat2resp , probs){ 
#.Call("SMIRT_CALCLIKE", dat2 , dat2resp , probs , PACKAGE = "sirt")
# 					}

					
######################################
# estimation of covariance
.smirt.est.covariance <- function( f.qk.yi , Sigma , theta.k , N ,
		mu.fixed , variance.fixed , D , est.corr , irtmodel 	){
		Sigma.cov <- Sigma
		delta.theta <- 1
		hwt <- f.qk.yi	
#		if (qmcnodes){ 
#			hwt <- hwt / nrow(theta.k) 
#			hwt <- hwt / rowSums( hwt )
#				}	
		thetabar <- hwt%*%theta.k	
		# calculation of mu
		mu <- colSums( thetabar ) / N
		if ( ! is.null(mu.fixed ) ){
			  if (is.matrix(mu.fixed) ){	 
			    mu0 <- mu
				mu[ mu.fixed[,1] ] <- mu.fixed[,2]
									}
							}
		# calculation of the covariance matrix
		theta.k.adj <- theta.k - matrix( mu , nrow=nrow(theta.k) , 
									ncol=ncol(theta.k) , byrow=TRUE)							
		for (dd1 in 1:D){
			for (dd2 in dd1:D){
				tk <- theta.k.adj[,dd1]*theta.k.adj[,dd2]
				h1 <- ( hwt %*% tk ) * delta.theta
				Sigma.cov[dd1,dd2] <- sum( h1 ) / N
				if (dd1 < dd2 ){ Sigma.cov[dd2,dd1] <- Sigma.cov[dd1,dd2] }
									}
								}								
		if ( est.corr ){ Sigma.cov <- stats::cov2cor(Sigma.cov ) }					
		if ( ! is.null(variance.fixed ) ){
				Sigma.cov[ variance.fixed[,1:2,drop=FALSE] ] <- variance.fixed[,3]
				Sigma.cov[ variance.fixed[,c(2,1),drop=FALSE] ] <- variance.fixed[,3]		
#										
									}
		diag(Sigma.cov) <- diag(Sigma.cov) + 10^(-10)
		pi.k <- matrix( mvtnorm::dmvnorm( theta.k , mean = mu , sigma = Sigma.cov )	, ncol=1 )		
		pi.k <- pi.k / sum( pi.k )			
		res <- list("mu"=mu , "Sigma"=Sigma.cov , "pi.k"= pi.k )				
		return(res)
					}
#################################################################
# convert matrix with probabilities
problong2probarray <- function( probres , I , TP ){
		probs <- array(0 , dim=c(I,2,TP))
		probs[,1,] <- 1 - probres
		probs[,2,] <- probres
		return(probs)
		}

###########################################
# estimation of b
.smirt.est.b.noncomp <- function(   b , a , c , d , Qmatrix , est.b , theta.k , 
        n.ik , I , K , TP , D , numdiff.parm=.001 , max.increment=1,
		msteps ,  mstepconv  , increment.factor){
    h <- numdiff.parm
	diffindex <- est.b		# zeros are allowed!
	cat("  M steps b parameter   |")
	it <- 0 ;	conv1 <- 1000	
	Q2 <- Q1 <- 0*Qmatrix
	Q <- Qmatrix	
	se.b <- b
	b00 <- b
	while( ( it < msteps ) & ( conv1 > mstepconv ) ){	
		b0 <- b
		for (dd in 1:D){
#				dd <- 2
			Q2 <- Q1
			Q2[,dd] <- 1 * ( Qmatrix[,dd] != 0 )
			
			probres <- calcprob.noncomp( a , b, Q , thetak=theta.k , c , d )
			pjk <- problong2probarray( probres , I , TP )		

			probres <- calcprob.noncomp( a , b + h*Q2, Q , thetak=theta.k , c , d )
			pjk1 <- problong2probarray( probres , I , TP )		

			probres <- calcprob.noncomp( a , b - h*Q2, Q , thetak=theta.k , c , d )
			pjk2 <- problong2probarray( probres , I , TP )		
			
			# numerical differentiation			
			res <- .rm.numdiff.index( pjk , pjk1 , pjk2 , n.ik , diffindex[,dd] , 
					max.increment=max.increment , numdiff.parm )
			ind <- match( diffindex[,dd] , sort(unique( diffindex[,dd] )) )
			b[,dd] <- b[,dd] + (res$increment)[ind]
			se.b[,dd] <- (sqrt(  1 / abs(res$d2) ))[ind]	
						}   # end dd
		conv1 <- max( abs( b - b0 ) )
		it <- it+1
		cat("-") # ; flush.console()
			}
	cat(" " , it , "Step(s) \n")	#; flush.console()	
	if ( increment.factor > 1){
		b <- .adj.maxincrement.parameter( oldparm=b00 , newparm=b , 
					max.increment=max.increment )		
						}
    res <- list("b" = b , "se.b" = se.b , 
			"ll" = sum(res$ll0) )
    return(res)
			}			
		
###########################################
# estimation of a
.smirt.est.a.noncomp <- function(   b , a , c , d , Qmatrix , est.a , theta.k , 
        n.ik , I , K , TP , D , numdiff.parm=.001 , max.a.increment ,
		msteps ,  mstepconv , increment.factor){
    h <- numdiff.parm
	diffindex <- est.a		# zeros are allowed!
	cat("  M steps a parameter   |")
	it <- 0 ;	conv1 <- 1000	
	Q2 <- Q1 <- 0*Qmatrix
	Q <- Qmatrix	
	se.a <- a
	a00 <- a
	while( ( it < msteps ) & ( conv1 > mstepconv ) ){	
		a0 <- a		
		for (dd in 1:D){
#				dd <- 2
			Q2 <- Q1
			Q2[,dd] <- 1 * ( Qmatrix[,dd] != 0 )
			
			probres <- calcprob.noncomp( a , b, Q , thetak=theta.k , c , d )
			pjk <- problong2probarray( probres , I , TP )		

			probres <- calcprob.noncomp( a + h*Q2 , b , Q , thetak=theta.k , c , d )
			pjk1 <- problong2probarray( probres , I , TP )		

			probres <- calcprob.noncomp( a - h*Q2 , b , Q , thetak=theta.k , c , d )
			pjk2 <- problong2probarray( probres , I , TP )		
			
			# numerical differentiation			
			res <- .rm.numdiff.index( pjk , pjk1 , pjk2 , n.ik , diffindex[,dd] , 
					max.increment=max.a.increment[,dd] , numdiff.parm )
			ind <- match( diffindex[,dd] , sort(unique( diffindex[,dd] )) )
			a[,dd] <- a[,dd] + (res$increment)[ind]
#			max.a.increment[,dd] <- abs( (res$increment)[ind] )
			se.a[,dd] <- (sqrt(  1 / abs(res$d2) ))[ind]			
						}   # end dd
		conv1 <- max( abs( a - a0 ) )
		it <- it+1
		cat("-") # ; flush.console()
			}
	cat(" " , it , "Step(s) \n")	#; flush.console()	
	#****
	# post-processing of a parameters	
	if ( increment.factor > 1){
		a <- .adj.maxincrement.parameter( oldparm=a00 , newparm=a , 
					max.increment=max.a.increment )		
						}	
    res <- list("a" = a , "se.a" = se.a , "ll" = sum(res$ll0) )
    return(res)
			}				
		
			
			
###########################################
# estimation of c
.smirt.est.c.noncomp <- function(   b , a , c , d , Qmatrix , est.c , theta.k , 
        n.ik , I , K , TP , D , numdiff.parm=.001 , max.increment=max.increment,
		msteps ,  mstepconv  , increment.factor){
    h <- numdiff.parm
	diffindex <- est.c		# zeros are allowed!
	Q1 <- rep(1,I)
	Q1[ est.c == 0 ] <- 0
	Q <- Qmatrix	
	c00 <- c
	cat("  M steps c parameter   |")
	it <- 0 ;	conv1 <- 1000	
	while( ( it < msteps ) & ( conv1 > mstepconv ) ){	
		c0 <- c
			probres <- calcprob.noncomp( a , b, Q , thetak=theta.k , c , d )
			pjk <- problong2probarray( probres , I , TP )		

			probres <- calcprob.noncomp( a  , b , Q , thetak=theta.k , c+h*Q1 , d )
			pjk1 <- problong2probarray( probres , I , TP )		

			probres <- calcprob.noncomp( a  , b , Q , thetak=theta.k , c-h*Q1 , d )
			pjk2 <- problong2probarray( probres , I , TP )		
			
			# numerical differentiation			
			res <- .rm.numdiff.index( pjk , pjk1 , pjk2 , n.ik , diffindex , 
					max.increment=max.increment , numdiff.parm )
			ind <- match( diffindex , sort(unique( diffindex ))	)
			incr <- res$increment
			incr[ is.na(incr ) ] <- 0
			c <- c + incr[ind]
			c[ c< 0 ] <- .001
#			max.a.increment[,dd] <- abs( (res$increment)[ind] )
			se.c <- (sqrt(  1 / abs(res$d2) ))[ind]	
		conv1 <- max( abs( c - c0 ) )
		it <- it+1
		cat("-") # ; flush.console()
			}
	cat(" " , it , "Step(s) \n")	#; flush.console()	
	if ( increment.factor > 1){
		c <- .adj.maxincrement.parameter( oldparm=c00 , newparm=c , 
					max.increment=max.increment )		
						}	
    res <- list("c" = c , "se.c" = se.c , 
			"ll" = sum(res$ll0) )
    return(res)
			}				
			

###########################################
# estimation of c
.smirt.est.d.noncomp <- function(   b , a , c , d , Qmatrix , est.d , theta.k , 
        n.ik , I , K , TP , D , numdiff.parm=.001 , max.increment=max.increment,
		msteps ,  mstepconv , increment.factor){
    h <- numdiff.parm
	diffindex <- est.d		# zeros are allowed!
	Q1 <- rep(1,I)
	Q1[ est.d == 0 ] <- 0
	Q <- Qmatrix	
	d00 <- d
	cat("  M steps d parameter   |")
	it <- 0 ;	conv1 <- 1000	
	while( ( it < msteps ) & ( conv1 > mstepconv ) ){	
		d0 <- d
			probres <- calcprob.noncomp( a , b, Q , thetak=theta.k , c , d )
			pjk <- problong2probarray( probres , I , TP )		

			probres <- calcprob.noncomp( a  , b , Q , thetak=theta.k , c , d +h*Q1)
			pjk1 <- problong2probarray( probres , I , TP )		

			probres <- calcprob.noncomp( a  , b , Q , thetak=theta.k , c , d-h*Q1 )
			pjk2 <- problong2probarray( probres , I , TP )		
			
			# numerical differentiation			
			res <- .rm.numdiff.index( pjk , pjk1 , pjk2 , n.ik , diffindex , 
					max.increment=max.increment , numdiff.parm )
			ind <- match( diffindex , sort(unique( diffindex ))		)
			incr <- res$increment
			incr[ is.na(incr ) ] <- 0			
			d <- d + incr[ind]
			d[ d>1 ] <- .999
#			max.a.increment[,dd] <- abs( (res$increment)[ind] )
			se.d <- (sqrt(  1 / abs(res$d2) ))[ind]	
		conv1 <- max( abs( d - d0 ) )
		it <- it+1
		cat("-") # ; flush.console()
			}
	cat(" " , it , "Step(s) \n")	#; flush.console()	
	if ( increment.factor > 1){
		d <- .adj.maxincrement.parameter( oldparm=d00 , newparm=d , 
					max.increment=max.increment )		
						}		
    res <- list("d" = d , "se.d" = se.d , "ll" = sum(res$ll0) )
    return(res)
			}				
			
			