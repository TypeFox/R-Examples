
####################################################################
# calculation of probabilities in the missing data IRT model
.calcprob.missing1 <- function( theta.k , b , beta , delta.miss , pjk ){	
	M1 <- t( stats::plogis( outer( theta.k[,1] , b , "-" ) ))
	# probability of a response for incorrect item responses
	M2a <- t( stats::plogis( outer( theta.k[,2] , beta , "-" ) ))
	# probability of a response for correct item responses
	M2b <- t( stats::plogis( outer( theta.k[,2] , beta + delta.miss , "-" ) ))
	pjk[,1,] <- 1 - M1
	pjk[,2,] <- M1
	# compute joint probability (R,Y)
	pjk[,1,] <- M2a*pjk[,1,]
	pjk[,2,] <- M2b*pjk[,2,]
	# P(R=0) = P(R=0|Y=0)P(Y=0)+P(R=0|Y=1)P(Y=1)
	pjk[,3,] <- (1-M2a)*(1-M1) + (1-M2b)*M1
	return(pjk)
	}

#######################################################################	
# E step missing data IRT model			
.e.step.missing1 <- function( dat2 , dat2.resp , theta.k , b , beta , delta.miss , I , CC ,
		TP , group_ , pi.k , pjk  , weights ){
# zz0 <- Sys.time()		
	#********************
	# calculate probability
	# probability correct response
	pjk <- .calcprob.missing1( theta.k , b , beta , delta.miss , pjk )	
#  cat("-- calcprob missing1") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1	
	
	#******
	# calculate likelihood
	# SEXP probs_pcm_groups_C( SEXP dat_, SEXP dat_resp_, SEXP group_, SEXP probs_, 
	#	 SEXP CC_, SEXP TP_ ){
	probs_ <- as.matrix( array( pjk , dim=c(I ,CC*TP) ) )
	#     	// probs ( ii , cc , tt , gg ) =   
	#     	// probs_C(ii ,  cc + tt*CC + gg * CC*TP ) 
	m1 <- .Call( "probs_pcm_groups_C" , dat2 , dat2.resp , group_ , probs_ , CC , TP , 
				package="sirt" )
	f.yi.qk <- m1$fyiqk	
# cat("-- likelihood") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1		
# likelihood purely written in R is much slower
	#****
	# calculate expected counts
	# SEXP calccounts_pcm_groups_C( SEXP dat_, SEXP dat_resp_, SEXP group_, SEXP fyiqk_, 
	# 	SEXP pik_, SEXP CC_, SEXP weights_ ){
	e1 <- .Call( "calccounts_pcm_groups_C" , 
					dat2 , dat2.resp , group_, f.yi.qk , pi.k , CC , weights ,
					package="sirt" )
# cat("-- posterior and counts") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1							
	e1$f.yi.qk <- f.yi.qk
	v1 <- array( e1$nik , dim=c(I,CC,TP) )
	e1$pjk <- pjk
	e1$n.k <- e1$count_pik
	e1$r.jk <- e1$n.jk <- NULL
	e1$n.ik <- v1
	e1$f.qk.yi <- e1$fqkyi
	return(e1)
		}
#########################################################################



####################################################################
# general function for numerical differentiation
# diffindex aggregates across super items
.mml2.numdiff.index <- function( pjk , pjk1 , pjk2 , n.ik , diffindex , 
		max.increment , numdiff.parm , eps=10^(-80) ){					
	h <- numdiff.parm
    an.ik <- n.ik
    ll0 <- rowSums( an.ik * log(pjk+eps) )
    ll1 <- rowSums( an.ik * log(pjk1+eps) )
    ll2 <- rowSums( an.ik * log(pjk2+eps) )
    ll0 <- stats::aggregate( ll0 , list(diffindex) , sum )[,2]
    ll1 <- stats::aggregate( ll1 , list(diffindex) , sum )[,2]
    ll2 <- stats::aggregate( ll2 , list(diffindex) , sum )[,2]
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
########################################################
########################################################
# M-step for missing data model			
.mstep.mml.missing1 <- function( theta.k , n.ik , mitermax , conv1 , 
		b , beta , delta.miss , pjk , numdiff.parm ,
		constraints , est.delta , min.beta , est_delta ){
    h <- numdiff.parm
	se.delta <- 0
	miter <- 0 
	diffindex <- seq( 1 , length(b) )
	
	
	diffindex1 <- est.delta
	
	# diffindex1 <- rep( 1 , length(b) )
	# diffindex1 <- seq( 1 , length(b) )
	
	# est_delta <- sum( 1 - is.na( est.delta) ) > 0
	
#	Q0 <- 0 * c.rater
	max_incr_b <- 1
	max_incr_beta <- 1
	max_incr_delta <- 1
	miterchange <- 1000
	while( ( miterchange > conv1 ) & ( miter < mitermax ) ){
	    #--- update b
		b0 <- b
		pjk <- .calcprob.missing1( theta.k , b , beta , delta.miss , pjk )				
		pjk1 <- .calcprob.missing1( theta.k , b+h , beta, delta.miss , pjk )				
		pjk2 <- .calcprob.missing1( theta.k , b-h , beta, delta.miss , pjk )				
		# numerical differentiation			
		res <- .mml2.numdiff.index( pjk , pjk1 , pjk2 , n.ik , diffindex , 
				max.increment= max_incr_b , numdiff.parm )			
		b <- b0 + res$increment
		max_incr_b <- max( abs( res$increment ) )		
		se.b <- sqrt( abs( 1 / res$d2 ) )
        if ( ! is.null( constraints) ){   
			b[ constraints[,1] ] <- constraints[,2] 
					}		

	    #--- update beta
		beta0 <- beta
		pjk <- .calcprob.missing1( theta.k , b , beta , delta.miss , pjk )				
		pjk1 <- .calcprob.missing1( theta.k , b , beta+h, delta.miss , pjk )				
		pjk2 <- .calcprob.missing1( theta.k , b , beta-h, delta.miss , pjk )				
		# numerical differentiation			
		res <- .mml2.numdiff.index( pjk , pjk1 , pjk2 , n.ik , diffindex , 
				       max.increment= max_incr_beta , numdiff.parm )									
		beta <- beta0 + res$increment		
		max_incr_beta <- max( abs( res$increment ) )
		se.beta <- sqrt( abs( 1 / res$d2 ) )

	    #--- update delta
		if (est_delta ){
			delta0 <- delta.miss
			pjk <- .calcprob.missing1( theta.k , b , beta , delta.miss , pjk )				
			pjk1 <- .calcprob.missing1( theta.k , b , beta, delta.miss+h , pjk )				
			pjk2 <- .calcprob.missing1( theta.k , b , beta, delta.miss-h , pjk )				
			# numerical differentiation			
			res <- .mml2.numdiff.index( pjk , pjk1 , pjk2 , n.ik , diffindex1 , 
					        max.increment= max_incr_delta , numdiff.parm )								
			incr <- res$increment
			d2 <- res$d2
			d2 <- d2[ est.delta ]
			d2[ is.na(incr) ] <- 1E10	
			incr <- incr[ est.delta ]
			incr[ is.na(incr) ] <- 0
					
			delta.miss <- delta0 + incr
			
			
			max_incr_delta <- max( abs( incr ) )
			se.delta <- sqrt( abs( 1 / d2 ) )		
				}
		
		miter <- miter + 1
				}
      res <- list("b"=b , "se.b"=b , "beta" = beta , "se.beta"=se.beta ,
			"delta.miss" = delta.miss , "se.delta" = se.delta )
	  return(res)
		}