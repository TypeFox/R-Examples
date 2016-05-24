# probs <- .gdm.calc.prob( a,b,thetaDes,Qmatrix,I,K,TP,TD)
#############################################################
.slca.calc.prob0 <- function( Xdes , Xlambda , I,K,TP){
	# Xdes [ 1:I , 1:(K+1) , 1:TP , 1:Nlam ] 
	p1 <- probs <- array( 0 , dim=c(I,K+1,TP) )
	for (tt in 1:TP){
		tmp0 <- 0
		for (hh in 1:(K+1) ){
			tmp1 <- exp( Xdes[ , hh , tt , ] %*% Xlambda )
			tmp0 <- tmp0 + tmp1
			p1[ , hh , tt  ] <- tmp1 
								}
		for (hh in 1:(K+1) ){
			probs[,hh,tt] <- p1[ , hh , tt  ] / tmp0
								}
				}
	return(probs)
			}
#############################################################
# Rcpp function for calculating probabilities
.slca.calc.prob <- function( XdesM , dimXdes , Xlambda ){
#   res <- calc_slca_probs( XdesM , dimXdes , Xlambda )  # .Call
   res <- .Call( "calc_slca_probs" , XdesM , dimXdes , Xlambda , packlage="CDM")  # .Call
   I <- dimXdes[1]
   maxK <- dimXdes[2]
   TP <- dimXdes[3]
   probs <- array( res , dim=c( I , maxK , TP ))   
   return(probs)
		}
# .slca.calc.prob( XdesM , dimXdes , Xlambda )
## [ii,kk,jj] =Cpp= [ii+kk*I+jj*I*K]  (ii,kk,jj=0,...)

###########################################################################
# estimation of Xlambda parameters
.slca.est.Xlambda <- function(Xlambda , Xdes , probs, n.ik1, N.ik1, I, K, G,
	max.increment,TP,msteps,convM , Xlambda.fixed , XdesM , dimXdes , oldfac){		
 	max.increment <- 1
	iter <- 1
	eps <- 1e-8
	parchange <- 1
	Xlambda00 <- Xlambda
	Nlam <- length(Xlambda)
	n.ik <- aperm( n.ik1 , c(2,3,1) )
	N.ik <- aperm( N.ik1 , c(2,1) )
	maxK <- K+1
	while( ( iter <= msteps ) & ( parchange > convM)  ){
		Xlambda0 <- Xlambda
		probs <- .slca.calc.prob( XdesM , dimXdes , Xlambda )
		d2.b <- d1.b <- rep(eps,Nlam)	
		# probs  num [1:I, 1:maxK, 1:TP]
		# n.ik  num [1:I, 1:maxK, 1:TP]
		# N.ik  num [1:I,1:TP]
		# Xdes  num [1:I, 1:maxK, 1:TP, 1:Nlam] 	    
				
#		res <- calc_slca_deriv( XdesM , dimXdes , Xlambda , as.vector(probs) ,
#			as.vector(n.ik) , as.vector(N.ik) )   # .Call
		res <- .Call( "calc_slca_deriv" , XdesM , dimXdes , Xlambda , as.vector(probs) ,
			as.vector(n.ik) , as.vector(N.ik) , PACKAGE="CDM" )   
		res$d1b -> d1.b
		res$d2b -> d2.b							

#		increment <-  - d1.b / ( abs( d2.b + 10^(-10) ) )
		increment <-   d1.b / ( abs( d2.b + 10^(-10) ) )
		increment[ is.na(increment) ] <- 0		
		increment <- ifelse(abs(increment)> max.increment, 
					sign(increment)*max.increment , increment )						
				
		max.increment <- max(abs(increment)) / .98
		Xlambda <- Xlambda + increment
		se.Xlambda <- sqrt( 1 / abs( d2.b+10^(-10)) )
		if ( ! is.null( Xlambda.fixed) ){
			Xlambda[ Xlambda.fixed[,1] ] <- Xlambda.fixed[,2]
			se.Xlambda[ Xlambda.fixed[,1] ] <- 0		
				}
				
		iter <- iter + 1
		parchange <- max( abs(Xlambda0-Xlambda))
# cat(iter,parchange , "\n" )
			}
		if (oldfac > 0 ){
			Xlambda <- oldfac*Xlambda00 + ( 1 - oldfac ) *Xlambda
						}
		max.increment <- max( abs( Xlambda - Xlambda00 ))
		res <- list("Xlambda" = Xlambda , "se.Xlambda" = se.Xlambda , "max.increment.b"=max.increment)
		}

		
		
		
################################################	
# calculation of expected counts
.slca.calc.counts <- function(G, weights, dat.ind, dat, dat.resp,
			p.aj.xi, K, n.ik, TP,I,group , dat.ind2 , ind.group ,
			use.freqpatt ){
	# n.ik [ 1:TP , 1:I , 1:(K+1) , 1:G ]
	# N.ik [ 1:TP , 1:I ,  1:G ]
# z0 <- Sys.time()	
	N.ik <- array( 0 , dim=c(TP,I,G) )
	N.ik1 <- array( 0 , dim=c(TP,I) )
	n.ik1 <- array( 0 , dim=c(TP,I,K+1 ) ) 

    if (G==1){
	gg <- 1
		for (kk in 1:(K+1)){   #		kk <- 1	# category 0 ( -> 1 )
#			dkk <- (dat.ind[[kk]])
#			dkk2 <- dkk * dat.resp * weights
            dkk2 <- dat.ind2[[kk]][[gg]]
#			n.ik[,,kk,gg] <- t( p.aj.xi ) %*% dkk2
			n.ik[,,kk,gg] <- crossprod( p.aj.xi , dkk2 )
			n.ik1[,,kk] <- n.ik[,,kk,gg]
			N.ik[,,gg] <- N.ik[,,gg] + n.ik[,,kk,gg]
						}
            N.ik1 <- N.ik1 + N.ik[,,gg]						
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
				n.ik1[,,kk] <- n.ik1[,,kk] + n.ik[,,kk,gg]
#				n.ik[,,kk,gg] <- crossprod( p.aj.xi.gg , dkk2 )
				N.ik[,,gg] <- N.ik[,,gg] + n.ik[,,kk,gg]
						}						
				N.ik1 <- N.ik1 + N.ik[,,gg]
					}
				}
	res <- list("n.ik" = n.ik , "N.ik" = N.ik , "n.ik1" = n.ik1,
				"N.ik1" = N.ik1)					
	return( res)
	}

	

###########################################################################
# reduced skillspace estimation
.slca.est.skillspace <- function(Ngroup, pi.k , 
			delta.designmatrix , G , delta , delta.fixed ,
			eps=10^(-10) , oldfac , delta.linkfct){		
		# gg <- 1
	covdelta <- as.list(1:G)
	Z <- delta.designmatrix	
	delta0 <- delta
	ND <- length(delta)
	for (gg in 1:G){
	   if ( delta.linkfct == "log"){
		ntheta <- Ngroup[gg] * pi.k[,gg]
		#*****
		# ARb 2014-01-14 inclusion
		ntheta <- ntheta / sum(ntheta )		
		lntheta <- log(ntheta+eps)
		mod <- stats::lm( lntheta ~ 0 + Z , weights = ntheta )
		covbeta <- vcov(mod)		
		beta <- coef(mod)
							}
	   if ( delta.linkfct == "logit"){
			nj <- Ngroup[gg] * pi.k[,gg]
			pij <- stats::qlogis( pi.k[,gg] + eps )
			wj <- stats::plogis( delta.designmatrix %*% delta[,gg,drop=FALSE] )
			wj <- wj / sum(wj )
			wj <- wj[,1]
			n <- Ngroup[gg]	
			mod1 <- stats::lm( pij ~ 0 + delta.designmatrix ) 
			beta <- coef(mod1)		
			covbeta <- vcov(mod1)
				}						
		if ( ! is.null( delta.fixed ) ){
		# delta.fixed: 1st column: parameter index
		#              2nd column: group index
		#              3rd column: parameter value 
		    ind.gg <- which( delta.fixed[ ,2] == gg )
			if ( length(ind.gg) > 0 ){
				beta[ delta.fixed[ind.gg,1] ] <- delta.fixed[ind.gg,3]
									}
							}
		if ( delta.linkfct == "log"){
			pi.k[,gg] <- exp( Z %*% beta ) / Ngroup[gg]
									}
		if ( delta.linkfct == "logit"){
			pi.k[,gg] <- exp( delta.designmatrix %*% beta ) 
							}
		pi.k[,gg] <- pi.k[,gg] / sum( pi.k[,gg] )
		if ( oldfac > 0 ){
			beta <- oldfac*delta0[,gg] + ( 1 - oldfac)*beta
						}
		delta[,gg] <- beta
		covdelta[[gg]] <- covbeta
			}
	res <- list( "pi.k"=pi.k , "delta"=delta , 
			"covdelta" = covdelta )			
			}
	