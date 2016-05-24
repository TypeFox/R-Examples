
#######################################
# calculate probability in the GDM		
.gdm.calc.prob <- function( a,b,thetaDes,Qmatrix,I,K,TP,TD){
		probs <- array( 0 , dim=c(I,K+1,TP) )	# categories 0 , ... , K
		for (kk in 1:K){
			l0 <- matrix( b[,kk] , nrow=I,ncol=TP)
			for (td in 1:TD){ 	# kk <- 1	# category 1
				# td <- 1
				l0 <- l0 + a[ , td , kk ] * Qmatrix[ , td, kk] * matrix( thetaDes[ , td ] , nrow=I,ncol=TP , byrow=T)
							}
				probs[,kk+1,] <- l0
						}
		probs <- exp( probs )
		probs1 <- probs[,1,]
		for (kk in 2:(K+1)){ probs1 <- probs1 + probs[,kk,] }
		for (kk in 1:(K+1)){ 
			probs[,kk,] <- probs[,kk,] / probs1 
					}
		return(probs)
			}
###############################################################
# calculation of posterior probabilities
.gdm.calc.post <- function(pi.k,group,p.xi.aj,weights,G,ind.group ,
		use.freqpatt ){
		# posterior probabilities  P( \alpha_l | X_i ) 		
	sel <- 1
	if ( use.freqpatt & (G>1) ){ sel <- 2 }
	#*****************
	if ( sel == 1 ){
		prior <- ( t( pi.k ) )[ group , ]
		p.aj.xi <- prior * p.xi.aj 
		p.aj.xi <- p.aj.xi / rowSums( p.aj.xi )
		# calculate pi.k
		for (gg in 1:G){ # gg <- 1
#			ind.gg <- which( group == gg )
			ind.gg <- ind.group[[gg]]
			wgg <- weights[ind.gg]
			if (G==1){ 	pi.k[,gg] <- colSums( p.aj.xi * wgg ) / sum( wgg ) }
			if (G>1){ 	pi.k[,gg] <- colSums( p.aj.xi[ ind.gg , ] * wgg ) / sum( wgg ) }			
					}
				}
	#***********************
	if ( sel == 2 ){  # if use.freqpatt == TRUE for multiple groups
		# calculate pi.k
		p.aj.xi <- list(1:G)
		for (gg in 1:G){ # gg <- 1
			wgg <- weights[,gg]
			ind.gg <- which( wgg > 0 )
			NP <- length(ind.gg)
			wgg <- wgg[ind.gg]
			prior <- ( t( pi.k[,gg] ) )[ rep(1,NP) , ]
			p.aj.xi.gg <- prior * p.xi.aj[ind.gg,]
			p.aj.xi.gg <- p.aj.xi.gg / rowSums( p.aj.xi.gg )		
			p.aj.xi[[gg]] <- p.aj.xi.gg			
			pi.k[,gg] <- colSums( p.aj.xi.gg * wgg ) / sum( wgg )
					}
				}
	#**********************	
		res <- list("pi.k"=pi.k , "p.aj.xi"=p.aj.xi )
		return(res)		
			}
			
################################################	
# calculation of expected counts
.gdm.calc.counts <- function(G, weights, dat.ind, dat, dat.resp,
			p.aj.xi, K, n.ik, TP,I,group , dat.ind2 , ind.group ,
			use.freqpatt ){
	# n.ik [ 1:TP , 1:I , 1:(K+1) , 1:G ]
	# N.ik [ 1:TP , 1:I ,  1:G ]
# z0 <- Sys.time()	
	N.ik <- array( 0 , dim=c(TP,I,G) )
    if (G==1){
	gg <- 1
		for (kk in 1:(K+1)){   #		kk <- 1	# category 0 ( -> 1 )
#			dkk <- (dat.ind[[kk]])
#			dkk2 <- dkk * dat.resp * weights
            dkk2 <- dat.ind2[[kk]][[gg]]
#			n.ik[,,kk,gg] <- t( p.aj.xi ) %*% dkk2
			n.ik[,,kk,gg] <- crossprod( p.aj.xi , dkk2 )
			N.ik[,,gg] <- N.ik[,,gg] + n.ik[,,kk,gg]
						}	
				}
	if (G>1){
		for (gg in 1:G){	# gg <- 1

		if ( ! use.freqpatt ){	
			ind.gg <- ind.group[[gg]]
			t.p.aj.xi.gg <- t( p.aj.xi[ind.gg,] )
#			p.aj.xi.gg <- p.aj.xi[ind.gg,] 
					} 
		if (  use.freqpatt ){	
			t.p.aj.xi.gg <- t( p.aj.xi[[gg]] )
#			p.aj.xi.gg <- p.aj.xi[[gg]]
					} 
					
			for (kk in 1:(K+1)){   #		kk <- 1	# category 0 ( -> 1 )
#				dkk <- (dat.ind[[kk]])[ ind.gg , ]
#				dkk2 <- dkk * dat.resp[ind.gg,] * weights[ind.gg] 
				dkk2 <- dat.ind2[[kk]][[gg]]
				if ( use.freqpatt ){
				  if (G>1){ dkk2 <- dkk2[ which(weights[,gg] > 0) , ] }
#				  if (G==1){ dkk2 <- dkk2[ which(weights[,gg] > 0) , ] }
									}
				n.ik[,,kk,gg] <- t.p.aj.xi.gg %*% dkk2
#				n.ik[,,kk,gg] <- crossprod( p.aj.xi.gg , dkk2 )
				N.ik[,,gg] <- N.ik[,,gg] + n.ik[,,kk,gg]
						}						
					}
				}
	res <- list("n.ik" = n.ik , "N.ik" = N.ik )					
	return( res)
	}
	
	
###########################################################################
# estimation of b parameters
.gdm.est.b <- function(probs, n.ik, N.ik, I, K, G,b,b.constraint,
	max.increment,a,thetaDes,Qmatrix,TP,TD,msteps,convM ,
	centerintercepts ){		
 	max.increment <- 1
	iter <- 1
	parchange <- 1
	b00 <- b
	
	while( ( iter <= msteps ) & ( parchange > convM)  ){
		b0 <- b
		probs <- .gdm.calc.prob( a,b,thetaDes,Qmatrix,I,K,TP,TD)								
		d2.b <- d1.b <- matrix( 0 , nrow=I,ncol=K)				
		for (kk in 2:(K+1)){
			for (gg in 1:G){
				d1.b[,kk-1] <- d1.b[,kk-1] - rowSums( t(n.ik[,,kk,gg]) - t(N.ik[,,gg]) * probs[,kk,] )
				d2.b[,kk-1] <- d2.b[,kk-1]  + rowSums( t(N.ik[,,gg]) * ( 1 - probs[,kk,] ) * probs[,kk,] )
							}
						}		
#		increment <- - d1.b / d2.b 
		increment <-  - d1.b / ( abs( d2.b + 10^(-10) ) )
		increment[ is.na(increment) ] <- 0		
		increment <- ifelse(abs(increment)> max.increment, 
					sign(increment)*max.increment , increment )						
# print(increment)					
#        increment.temp <- diff.temp*abs(1/( deriv.temp + 10^(-20) ) )  
#        ci <- ceiling( abs(increment) / ( abs( max.increment ) + 10^(-10) ) )
 #       increment <- ifelse( abs(increment)> max.increment  , 
 #                                   increment /(2*ci) , increment)										
		max.increment <- max(abs(increment)) / .98
		b <- b + increment
		se.b <- sqrt( 1 / abs( d2.b+10^(-10)) )
		if ( ! is.null( b.constraint) ){
			b[ b.constraint[,1:2,drop=FALSE] ] <- b.constraint[,3,drop=FALSE]
			se.b[ b.constraint[,1:2,drop=FALSE] ] <- 0		
				}
		# centerintercepts
		if ( centerintercepts) {
		   if (TD==1){
				b <- b - mean(b)		
						}
			if (TD > 1){		
				for (dd in 1:TD){
					ind.dd <- which( Qmatrix[,dd,1] > 0 )
					m1 <- sum( b[ind.dd,] ) / ( ncol(b) * length(ind.dd) )	
					b[ind.dd,] <- b[ind.dd,] - 	m1
							}
						  }
						}				
		iter <- iter + 1
		parchange <- max( abs(b0-b))
# cat(iter,parchange , "\n" )
			}
		max.increment <- max( abs( b - b00 ))
		res <- list("b" = b , "se.b" = se.b , "max.increment.b"=max.increment)
		}

		
###########################################
# estimation of a
# Q matrix [1:I , 1:TD , 1:K]
# thetaDes [TP,TD]
# n.ik [ TP , I , K+1 , G ]
# N.ik [ TP , I , G ]
# probs [I , K+1 , TP ]
.gdm.est.a <- function(probs, n.ik, N.ik, I, K, G,a,a.constraint,TD,
				Qmatrix,thetaDes,TP, max.increment ,
				b , msteps , convM , centerslopes ){
	iter <- 1
	parchange <- 1
	a00 <- a

	maxa <- max.increment + 0 * a
# Reval("print(maxa)")	
	while( ( iter <= msteps ) & ( parchange > convM )  ){
		a0 <- a
		probs <- .gdm.calc.prob( a,b,thetaDes,Qmatrix,I,K,TP,TD)
		# 1st derivative	
		d2.b <- d1.b <- array( 0 , dim=c(I , TD ) )
		for (td in 1:TD){
			for (gg in 1:G){
				for (kk in 2:(K+1)){		
					QM <- matrix( Qmatrix[  , td , kk-1 ]  , nrow=TP , ncol=I , byrow=TRUE )
					v1 <- colSums( n.ik[,,kk,gg] * QM * thetaDes[ , td ] )
					v2 <- N.ik[,,gg] * QM * thetaDes[,td] *  t( probs[,kk,] )
					v2 <- colSums(v2)

					d1.b[  , td] <- d1.b[,td] + v1 - v2
				#*****************					
				# old gdm version
				#					v1 <- colSums( n.ik[,,kk,gg] * Qmatrix[  , td , kk-1 ] * thetaDes[ , td ] )
				#					v2 <- N.ik[,,gg] * matrix( Qmatrix[,td,kk-1] , 
				#                     nrow=TP , ncol=I) * thetaDes[,td] * 
				#							t( probs[,kk,] )
				#					v2 <- colSums(v2)
				# cat( "v1=",v1 , "v2=" ,v2 , "\n")			
				#					d1.b[  , td] <- d1.b[,td] + v1 - v2					
				#******************					
							}
						}
						}
		# 2nd derivative
		for (td in 1:TD){
			for (ii in 1:I){
				v1 <- l0 <- 0
				for (gg in 1:G){
				for (kk in 2:(K+1) ){		# kk <- 2
					v1 <- v1 + N.ik[,ii,gg] * as.vector( ( Qmatrix[ii,td,kk-1] * 
							thetaDes[ , td ] )^2 * t( probs[ii,kk,] ) )
					l0 <- l0 + as.vector ( Qmatrix[ii,td,kk-1] * thetaDes[ , td ]  * t( probs[ii,kk,] ) )
			
								}
								
								
							}				
				d2.b[ii,td] <- sum(v1) - sum( l0^2 * N.ik[,ii,gg] )
					}
					}
			increment <-  d1.b / abs( d2.b + 10^(-10) )
			increment[ is.na(increment) ] <- 0		
			increment <- ifelse(abs(increment)> max.increment, 
						   sign(increment)*max.increment , increment )							   					  
#			increment <- ifelse(abs(increment)> max.increment / sqrt( iter)  , 
#						   sign(increment)*max.increment / sqrt(iter ) , increment )							   					  
#			ci <- ceiling( abs(increment) / ( abs( max.increment ) + 10^(-10) ) )
#	        increment <- ifelse( abs(increment)> max.increment  , 
#	                                  increment /(2*ci) , increment)					
			a[,,1] <- a[,,1] + increment	
			se.a <- sqrt( 1 / abs( d2.b + 10^(-10) ) )
			if (K>1){ for (kk in 2:K){ a[,,kk] <- a[,,1] }	 }	
			if ( ! is.null( a.constraint) ){
				a[ a.constraint[,1:3,drop=FALSE] ] <- a.constraint[,4,drop=FALSE]
				se.a[ a.constraint[,1:3,drop=FALSE] ] <- 0			
				increment[ a.constraint[,1:2,drop=FALSE] ] <- 0			
					}		
			if (centerslopes){
			  if (TD>1){
				m11 <- t( colSums( a[,,1] ) / colSums( Qmatrix ) )	
				a[,,1] <- a[,,1] / m11[ rep(1,I) , ]
						}
			  if (TD==1){
				m11 <- t( colSums( a ) / colSums( Qmatrix ) )	
				a <- a / m11[ rep(1,I) , ]
						}						
						}
			parchange <- max( abs(a-a0))
			iter <- iter + 1
			}	# end iter
#		max.increment <- max(abs(a-a00))
		max.increment <- max(abs(a-a00)) / 1.005
		res <- list( "a" = a , "se.a" = se.a , "max.increment.a" = max.increment)
		return(res)
		}		

###########################################
# estimation of a
# Q matrix [1:I , 1:TD , 1:K]
# thetaDes [TP,TD]
# n.ik [ TP , I , K+1 , G ]
# N.ik [ TP , I , G ]
# probs [I , K+1 , TP ]
.gdm.est.a.cat <- function(probs, n.ik, N.ik, I, K, G,a,a.constraint,TD,
				Qmatrix,thetaDes,TP, max.increment ,
				b , msteps , convM  ){
	iter <- 1
	parchange <- 1
	a00 <- a
	while( ( iter <= msteps ) & ( parchange > convM )  ){
		a0 <- a
		probs <- .gdm.calc.prob( a,b,thetaDes,Qmatrix,I,K,TP,TD)

		# 1st derivative
		d2.b <- d1.b <- array( 0 , dim=c(I , TD , K ) )
		for (td in 1:TD){
		for (kk in 2:(K+1)){	
			for (gg in 1:G){
					QM <- matrix( Qmatrix[,td,kk-1] , nrow=TP , ncol=I, byrow=TRUE)
					v1 <- colSums( n.ik[,,kk,gg] * QM * thetaDes[ , td ] )
					v2 <- N.ik[,,gg] * QM * thetaDes[,td] *  t( probs[,kk,] )
					v2 <- colSums(v2)
					d1.b[  , td , kk-1] <- d1.b[  , td , kk-1] + v1 - v2
							}
						}	
						}	
		# 2nd derivative
		for (td in 1:TD){
			for (ii in 1:I){
				v1 <- l0 <- 0
				for (kk in 2:(K+1) ){		# kk <- 2
				v1 <- l0 <- 0
				  for (gg in 1:G){			
					v1 <- N.ik[,ii,gg] * as.vector( ( Qmatrix[ii,td,kk-1] * 
							thetaDes[ , td ] )^2 * t( probs[ii,kk,] ) )
					l0 <- as.vector ( Qmatrix[ii,td,kk-1] * thetaDes[ , td ]  * t( probs[ii,kk,] ) )
					d2.b[ii,td,kk-1] <- d2.b[ii,td,kk-1] + sum(v1) - sum( l0^2 * N.ik[,ii,gg] )				
								}
							}				
					}
					}				
	#		increment <-  d1.b / d2.b 
			increment <-  d1.b / abs( d2.b + 10^(-10) )
			increment[ is.na(increment) ] <- 0		
			increment <- ifelse(abs(increment)> max.increment, 
						    sign(increment)*max.increment , increment )	
	#        ci <- ceiling( abs(increment) / ( abs( max.increment ) + 10^(-10) ) )
	#        increment <- ifelse( abs(increment)> max.increment  , 
	#                                  increment /(2*ci) , increment)					
			a <- a + increment
			se.a <- sqrt( 1 / abs( d2.b + 10^(-10) ) )
	#		if (K>1){ for (kk in 2:K){ a[,,kk] <- a[,,1] }	 }
			if ( ! is.null( a.constraint) ){
				a[ a.constraint[,1:3,drop=FALSE] ] <- a.constraint[,4,drop=FALSE]
				se.a[ a.constraint[,1:3,drop=FALSE] ] <- 0		
				increment[ a.constraint[,1:3,drop=FALSE] ] <- 0						
					}
			iter <- iter + 1
			parchange <- max( abs( a - a0 ))
				} # iter
		max.increment <- max(abs(a-a00)) / .95				
		res <- list( "a" = a , "se.a" = se.a , "max.increment.a" = max.increment)
		return(res)
		}	

		
###########################################################################
# reduced skillspace estimation
.gdm.est.skillspace <- function(Ngroup, pi.k , Z, G , delta , eps=10^(-10) ){		
		# gg <- 1
	covdelta <- as.list(1:G)
	for (gg in 1:G){
		ntheta <- Ngroup[gg] * pi.k[,gg]
		#*****
		# ARb 2014-01-14 inclusion
		ntheta <- ntheta / sum(ntheta )		
		lntheta <- log(ntheta+eps)
		mod <- stats::lm( lntheta ~ 0 + Z , weights = ntheta )
		covbeta <- vcov(mod)		
		beta <- coef(mod)		
		pi.k[,gg] <- exp( Z %*% beta ) / Ngroup[gg]
		pi.k[,gg] <- pi.k[,gg] / sum( pi.k[,gg] )
		delta[,gg] <- beta
		covdelta[[gg]] <- covbeta
			}
	res <- list( "pi.k"=pi.k , "delta"=delta , 
			"covdelta" = covdelta )			
			}
			
##############################################################
# estimation of skill distribution under normality
.gdm.est.normalskills <- function( pi.k , theta.k , irtmodel,G , D ,
	mean.constraint , Sigma.constraint , standardized.latent ,
	p.aj.xi , group , ind.group , weights , b , a ){
	# mean.constraint [ dimension , group , value ]
	# Sigma.constraint [ dimension1 , dimension2 , group , value ]	
   ####################################
   # unidimensional model
   if (D==1){
	for (gg in 1:G){
		# gg <- 1
#		weights.gg <- weights[ ind.group[[gg]] ]
#		weights.gg <- weights.gg / sum(weights.gg )
		mg <- sum( theta.k[,1] * pi.k[,gg] )
		sdg <- sqrt( sum( theta.k[,1]^2 * pi.k[,gg] ) - mg^2 )
# cat( mg , sdg , "\n" )			
	if ( (! is.null ( mean.constraint ))  ){
		i1 <- mean.constraint[ mean.constraint[,2] == gg , , drop=FALSE]	
#		  if ( ( nrow(i1) == 1 ) & (G>=1) ){ 	
		  if ( ( nrow(i1) == 1 ) & (G>1) ){ 	
				if ( ( gg==1 ) & (i1[,1]==1) & (i1[,2]==1) ){ 
					b <- b + ( mg - i1[3] ) 
						}
				mg <- i1[3] 
						}
		  if ( nrow(i1) > 0 ){				
				mg <- i1[,3]				
							}						
					}				
	if ( ( ! is.null ( Sigma.constraint ) )  ){
		i1 <- Sigma.constraint[ Sigma.constraint[,3] == gg , , drop=FALSE]
#		  if ( ( nrow(i1) == 1 ) & (G>=1) ){ 	
		  if ( ( nrow(i1) == 1 ) & (G>1) ){ 		  
			if ( ( gg==1 ) & (i1[,1]==1) & (i1[,2]==1) ){ 
	                a <- a * sdg / sqrt(i1[4])  
							}		
				sdg <- sqrt(i1[4]) 
						}
					}
# cat( mg , sdg , "\n" )						
#	if (standardized.latent){ mg <- 0 ; sdg <- 0 }					
#		if (gg==1){ mg <- 0 }			
		pi.k[,gg] <- stats::dnorm( theta.k[,1] ,mean=mg , sd=sdg)		
		pi.k[,gg] <- pi.k[,gg] / sum( pi.k[,gg] )
		
			}			
		}
	#####################################
    # multidimensional model	
	if (D>1){
	  for (gg in 1:G){
		# gg <- 1
		mean.gg <- rep(0,D)
		Sigma.gg <- diag(0,D)
		
		for (dd in 1:D){
			# dd <- 1
			mean.gg[dd] <- sum( pi.k[,gg] * theta.k[,dd] )
				}
		for (dd1 in 1:D){
			for (dd2 in dd1:D){
#		dd1 <- 1 ; 	dd2 <- 1
		Sigma.gg[dd1,dd2] <- sum( pi.k[,gg] * (theta.k[,dd1] - mean.gg[dd1] )*(theta.k[,dd2] - mean.gg[dd2] ) ) 
#		Sigma.gg[dd1,dd2] <- Sigma.gg[dd1,dd2] - mean.gg[dd1] * mean.gg[dd2]
		Sigma.gg[dd2,dd1] <- Sigma.gg[dd1,dd2]
							}
						}
		Sigma.gg <- Sigma.gg + diag(10^(-10) , D )
		m.gg <- mean.constraint[ mean.constraint[,2] == 1 , ]
		if ( ! is.null(mean.constraint)){
		if( dim(m.gg)[1] > 0 ){
			mean.gg[ m.gg[,1] ] <- m.gg[,3]
								}}
		s.gg <- Sigma.constraint[ Sigma.constraint[,3] == 1 , ]	
		
#		if ( standardized.latent & ( gg == 1 )){
#			Sigma.gg <- cov2cor( Sigma.gg )
#				}
		
		
		if ( ! is.null(Sigma.constraint)){
		if( dim(s.gg)[1] > 0 ){		
			c1 <- stats::cov2cor( Sigma.gg )
				d1 <- diag(Sigma.gg)
				s.gg1 <- s.gg[ s.gg[,1] == s.gg[,2] , ]
				if ( nrow(s.gg1) > 0 ){
					d1[ s.gg1[,1:2] ] <- s.gg[,4]
						}
				d1 <- outer( sqrt(d1) , sqrt(d1) )*c1			
				s.gg2 <- s.gg[ s.gg[,1] != s.gg[,2] , ]
				if ( nrow(s.gg1) > 0 ){
					d1[ s.gg1[,1:2] ] <- s.gg[,4]
					d1[ s.gg1[,c(2,1)] ] <- s.gg[,4]					
						}
			Sigma.gg <- d1
						}
								}
		pi.k[,gg] <- mvtnorm::dmvnorm( theta.k , mean=mean.gg , sigma = Sigma.gg )			
		pi.k[,gg] <- pi.k[,gg] / sum( pi.k[,gg] )					
					}
				}
	res <- list("pi.k"=pi.k , "b" = b , "a"=a )
	return(res)
				}
#*************************************************************

#####################################################
# estimation of skill space
.gdm.est.skillspace.traits <- function( n.ik , a , b , theta.k , Qmatrix , I , K , TP,
		TD , numdiff.parm , max.increment , msteps , convM ){
	# sum over groups
	n.ik0 <- apply( n.ik , c(1,2,3) , sum )
	h <- numdiff.parm
	parchange <- 1000
	iter <- 1
	se.theta.k <- 0 * theta.k
	Q1 <- matrix( 0 , nrow=TP , ncol=TD)
	while( ( iter <= msteps ) & (parchange > convM ) ){
		theta.k0 <- theta.k
		for ( dd in 1:TD){
#			dd <- 1		
			Q0 <- Q1
			Q0[,dd] <- 1
			# calculate log-likelihood
			pjk <- .gdm.calc.prob( a,b,thetaDes=theta.k,Qmatrix,I,K,TP,TD)
			theta.k1 <- theta.k + h*Q0
			pjk1 <- .gdm.calc.prob( a,b,thetaDes=theta.k1,Qmatrix,I,K,TP,TD)		
			theta.k2 <- theta.k - h*Q0
			pjk2 <- .gdm.calc.prob( a,b,thetaDes=theta.k2,Qmatrix,I,K,TP,TD)		
			res <- .gdm.numdiff.index( pjk , pjk1 , pjk2 , n.ik=n.ik0 , 
					max.increment , numdiff.parm , eps=10^(-80) )
			theta.k[,dd] <- theta.k[,dd] + res$increment		
			se.theta.k[,dd] <- 1 / sqrt( abs(res$d2) )
					}
		iter <- iter + 1 
		parchange <- max( abs( theta.k - theta.k0 ))
			}	
#	thetaDes <- theta.k
	res <- list( "theta.k" = theta.k , "se.theta.k" = se.theta.k )
	return(res)
	}
##########################################################

####################################################################
# general function for numerical differentiation
# diffindex aggregates across super items
.gdm.numdiff.index <- function( pjk , pjk1 , pjk2 , n.ik , 
		max.increment , numdiff.parm , eps=10^(-80) ){					
	h <- numdiff.parm
    an.ik <- aperm( n.ik , c(2,3,1) )
    ll0 <- colSums( colSums( an.ik * log(pjk+eps) ))
    ll1 <- colSums( colSums(an.ik * log(pjk1+eps) ) )
    ll2 <- colSums( colSums( an.ik * log(pjk2+eps) ) )
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
	res <- list("increment"=increment , "d2"=d2 , "ll0"=ll0)
	return(res)
		}