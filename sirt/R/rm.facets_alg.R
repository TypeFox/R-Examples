
#######################################################
# parameters expanded dataset
rm.facets.itempar.expanded <- function( b.item , b.rater , Qmatrix , tau.item ,
        VV , K , I , TP , a.item , a.rater , item.index , rater.index ,
		theta.k , RR ){
	b <- tau.item[ item.index , ]
	b0 <- ( matrix( b.rater , nrow= RR , ncol=K) )[ rater.index , ] * 	Qmatrix[ item.index ,]	 		
	b <- b + b0
	# a parameter
	a <- a.item[ item.index ] * a.rater[ rater.index ]
	res <- list("a" = a , "b" = b )
 	return(res)
    }
#########################################################

################################################
# calculation of the likelihood
rm_calclike <- function (dat2,dat2resp,probs,K){ 
	.Call("RM_CALCPOST", dat2,dat2resp,probs,K, PACKAGE = "sirt")
					}

########################################################
# calculation of probabilities in the partial credit model
.rm.pcm.calcprobs <- function( a , b , Qmatrix , theta.k , I , K , TP ){
    probs <- array( 0 , dim=c(I,K+1,TP) )   # categories 0 , ... , K	
    for (kk in 1:K){
        l0 <- matrix( - b[,kk] , nrow=I,ncol=TP)
        l0 <- l0 + outer( a * Qmatrix[ , kk] , theta.k )
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
######################
sumtau <- function(tau.item){
	K <- ncol(tau.item)
	matr <- tau.item
	for (kk in 2:K){
		matr[,kk] <- rowSums( tau.item[,1:kk] )
				}
	return(matr)
			}
#############################################################################
# calculation of probabilities in the facet model
.rm.facets.calcprobs <- function( b.item , b.rater , Qmatrix , tau.item ,
        VV , K , I , TP , a.item , a.rater , item.index , rater.index ,
		theta.k , RR ){
    # b parameter
#	pcm.param <- FALSE
#	if (pcm.param){
#		b.item0 <- b.item
#		b <- matrix( b.item0 , nrow= VV , ncol=K ) * Qmatrix + sumtau(tau.item)
#		b <- b[ item.index , ]
#		b0 <- ( matrix( b.rater , nrow= RR , ncol=K) )[ rater.index , ] * 
#				Qmatrix[ item.index ,]
#	
#		b <- b + b0				
#				} else { 

		b <- tau.item[ item.index , ]
		b0 <- ( matrix( b.rater , nrow= RR , ncol=K) )[ rater.index , ] * 	Qmatrix[ item.index ,]	 		
		b <- b + b0
#			}
	# a parameter
	a <- a.item[ item.index ] * a.rater[ rater.index ]
    res <- .rm.pcm.calcprobs( a , b , Qmatrix=Qmatrix[item.index,] , theta.k , I , K , TP )
	return(res)
    }
#############################################################################	
# cpp implementation of calculation of facets probabilities
.rm.facets.calcprobs2 <- function( b.item , b.rater , Qmatrix , tau.item ,
        VV , K , I , TP , a.item , a.rater , item.index , rater.index ,
		theta.k , RR ){						
	## SEXP rm_facets_calcprobs_cpp( SEXP b_item_, SEXP b_rater_, SEXP Qmatrix_, 
	##	SEXP tau_item_, SEXP K_, SEXP I_, SEXP TP_, SEXP a_item_, SEXP a_rater_, 
	##	SEXP item_index_, SEXP rater_index_, SEXP theta_k_ ) ;
#	probs <- rm_facets_calcprobs_cpp( b.item , b.rater ,Qmatrix , tau.item ,
#         K , I , TP , a.item , a.rater , item.index-1 , rater.index -1 ,
#         theta.k )
	probs <- .Call("rm_facets_calcprobs_cpp", 
				b.item , b.rater ,Qmatrix , tau.item ,
				K , I , TP , a.item , a.rater , item.index-1 , rater.index -1 ,
				theta.k  , 
				PACKAGE = "sirt")		 
	probs <- array( probs , dim=c(I , K+1 , TP ) )
	return(probs)
	  }
	
#######################################################
# calculate posterior and counts
.rm.posterior <- function( dat2 , dat2.resp , TP , pi.k ,
	K, I , probs ){
	# calculate likelihood
#	f.yi.qk <- matrix( 1 , nrow(dat2) , TP )
#	for (ii in 1:ncol(dat2)){
		 #   ii <- 1
#		ind.ii <- which( dat2.resp[,ii] == 1 )
#		f.yi.qk[ind.ii,] <- f.yi.qk[ind.ii,] * probs[ ii , dat2[ind.ii,ii]+1 , ]
#							}
	probsM <- matrix( aperm( probs , c(2,1,3) ) , nrow=I*(K+1) , ncol=TP )
	f.yi.qk <- rm_calclike( dat2=dat2 , dat2resp = dat2.resp , 
					probs = probsM , K=K)$fyiqk
							
							
	####################							
	# calculate posterior and expected counts
	prior <- matrix( pi.k , nrow=nrow(dat2) , ncol=TP , byrow =TRUE )
	f.qk.yi <- f.yi.qk * prior
	f.qk.yi <- f.qk.yi / rowSums( f.qk.yi )
	# expected counts
	n.ik <- array( 0 , dim =c(TP , I , K+1 ) )
	N.ik <- array( 0 , dim=c(TP , I ) )
	for (kk in 1:(K+1) ){
		# kk <- 1
#		n.ik[,,kk] <- t( t(dat2.resp*(dat2==(kk-1)) ) %*% f.qk.yi )
		n.ik[,,kk] <- crossprod( f.qk.yi , dat2.resp*(dat2==(kk-1)) )
		N.ik <- N.ik + n.ik[,,kk]
						}
				
    ll <- sum( log( rowSums( f.yi.qk * outer( rep(1,nrow(f.yi.qk)) , pi.k) ) ) )
					
	# compute pi.k
	pi.k <- colMeans( f.qk.yi )
	res <- list( "f.yi.qk" = f.yi.qk , "f.qk.yi" = f.qk.yi , 
			"n.ik" = n.ik , "N.ik" = N.ik , "pi.k" = pi.k , "ll"=ll)
	return(res)
	}
#####################################################
# estimation of slope parameter for items
.rm.facets.est.a.item <- function( b.item , b.rater , Qmatrix , tau.item ,
        VV , K , I , TP , a.item , a.rater , item.index , rater.index ,
        n.ik , numdiff.parm=.001 , max.b.increment=1,theta.k , msteps ,
		mstepconv , a.item.center , a.item.fixed  ){
    h <- numdiff.parm
	diffindex <- item.index
	RR <- length(b.rater)
	cat("  M steps a.item parameter   |")
	it <- 0 ;	conv1 <- 1000	
	while( ( it < msteps ) & ( conv1 > mstepconv ) ){	
		a.item0 <- a.item	
		pjk <- .rm.facets.calcprobs2( b.item , b.rater , Qmatrix , tau.item ,
				VV , K , I , TP , a.item , a.rater , diffindex , rater.index , theta.k,RR)
		pjk1 <- .rm.facets.calcprobs2( b.item , b.rater , Qmatrix , tau.item ,
				VV , K , I , TP , a.item+h , a.rater , diffindex , rater.index , theta.k,RR)
		pjk2 <- .rm.facets.calcprobs2( b.item , b.rater , Qmatrix , tau.item ,
				VV , K , I , TP , a.item-h , a.rater , diffindex , rater.index , theta.k,RR)
		# numerical differentiation			
		res <- .rm.numdiff.index( pjk , pjk1 , pjk2 , n.ik , diffindex , 
				       max.increment=max.b.increment , numdiff.parm )					
		a.item <- a.item + res$increment
		a.item[ a.item < .05 ] <- .05
#		a.item <- a.item - mean(a.item ) + 1
		if ( ! is.null( a.item.fixed) ){		
			ind <- which( ! is.na( a.item.fixed  ) )
		    a.item[ind] <- a.item.fixed[ind]
            res$d2[ ind ] <- -1E10			
								}

		if ( a.item.center ){
			b1 <- mean( log( a.item ) )
			a.item <- a.item / exp( b1 )
						}
		conv1 <- max( abs( a.item - a.item0 ) )
		it <- it+1
		cat("-") 
			}
	cat(" " , it , "Step(s) \n")	
    res <- list("a.item" = a.item , "se.a.item" = sqrt( abs(-1/res$d2 )) , 
			"ll" = sum(res$ll0) )
    return(res)
			}			
#####################################################
# estimation of tau.item parameters
.rm.facets.est.tau.item <- function( b.item , b.rater , Qmatrix , tau.item ,
        VV , K , I , TP , a.item , a.rater , item.index , rater.index ,
        n.ik , numdiff.parm=.001 , max.b.increment=1,theta.k , msteps ,
		mstepconv , tau.item.fixed , tau.item.fixed_val ){
zz0 <- Sys.time()		
active <- TRUE
active <- FALSE
		
    h <- numdiff.parm
	diffindex <- item.index
	RR <- length(b.rater)	
	Q0 <- matrix(0,nrow=VV, ncol=K)
	se.tau.item <- Q0
	cat("  M steps tau.item parameter |")
	it <- 0 ;	conv1 <- 1000
	while( ( it < msteps ) & ( conv1 > mstepconv ) ){	
		tau.item0 <- tau.item
		for (kk in 1:K){
	#		kk <- 1
			Q1 <- Q0
			Q1[,kk] <- 1
			pjk <- .rm.facets.calcprobs2( b.item , b.rater , Qmatrix , tau.item ,
					VV , K , I , TP , a.item , a.rater , item.index , rater.index , theta.k,RR)
			pjk1 <- .rm.facets.calcprobs2( b.item , b.rater , Qmatrix , tau.item+h*Q1 ,
					VV , K , I , TP , a.item , a.rater , item.index , rater.index , theta.k,RR)
			pjk2 <- .rm.facets.calcprobs2( b.item , b.rater , Qmatrix , tau.item-h*Q1 ,
					VV , K , I , TP , a.item , a.rater , item.index , rater.index , theta.k,RR)
			# numerical differentiation			
			res <- .rm.numdiff.index( pjk , pjk1 , pjk2 , n.ik , diffindex , 
					max.increment=max.b.increment , numdiff.parm )					
			increment <- Q1*matrix( res$increment , nrow=VV , ncol=K)	
			tau.item <- tau.item + increment
			se.tau.item[,kk] <- sqrt(abs(-1/res$d2)	)
zz0 <- sirtcat( "  *** tau.item for one k and one iteration   " , zz0 , active )			
					}
					
		if ( ! is.null( tau.item.fixed_val ) ){
              MK <- ncol( tau.item.fixed_val )
			  for ( kk in 1:MK){
			  ind <- which( ! is.na( tau.item.fixed_val[,kk]) )
			  if ( length(ind) > 0 ){
				tau.item[ ind , kk] <- tau.item.fixed_val[ ind , kk]
								}
								}
							}

					
		conv1 <- max( abs( tau.item - tau.item0 ) )
		it <- it+1
		cat("-") 
		if (!is.null(tau.item.fixed)){
			tau.item[ tau.item.fixed[,1:2,drop=FALSE] ] <- tau.item.fixed[,3]
								}
			}
	cat(" " , it , "Step(s) \n")	
	res <- list("tau.item" = tau.item , "se.tau.item" = se.tau.item , 
			"ll" = sum(res$ll0) )
    return(res)
			}
			
			
#########################################
# estimation of rater severity
.rm.facets.est.b.rater <- function( b.item , b.rater , Qmatrix , tau.item ,
        VV , K , I , TP , a.item , a.rater , item.index , rater.index ,
        n.ik , numdiff.parm=.001 , max.b.increment=1 , theta.k , msteps ,
		mstepconv , b.rater.center , b.rater.fixed  ){
    h <- numdiff.parm
	diffindex <- rater.index
	RR <- length(b.rater)	
	cat("  M steps b.rater parameter  |")
	it <- 0 ;	conv1 <- 1000
	while( ( it < msteps ) & ( conv1 > mstepconv ) ){
		b0 <- b.rater
		pjk <- .rm.facets.calcprobs2( b.item , b.rater , Qmatrix , tau.item ,
				VV , K , I , TP , a.item , a.rater , item.index , rater.index , theta.k,RR)
		pjk1 <- .rm.facets.calcprobs2( b.item , b.rater+h , Qmatrix , tau.item ,
				VV , K , I , TP , a.item , a.rater , item.index , rater.index , theta.k,RR)
		pjk2 <- .rm.facets.calcprobs2( b.item , b.rater-h , Qmatrix , tau.item ,
				VV , K , I , TP , a.item , a.rater , item.index , rater.index , theta.k,RR)
		# numerical differentiation			
		res <- .rm.numdiff.index( pjk , pjk1 , pjk2 , n.ik , diffindex , 
				max.increment=max.b.increment , numdiff.parm )					
		increment <- res$increment
#		increment <- increment - mean(increment )
		b.rater <- b.rater + increment
		if ( ! is.null( b.rater.fixed) ){		
			ind <- which( ! is.na( b.rater.fixed  ) )
		    b.rater[ind] <- b.rater.fixed[ind]
            res$d2[ ind ] <- -1E10			
								}
		brc <- mean( b.rater )								
		# centering
		if ( b.rater.center){
			b.rater[RR] <- - sum( b.rater[-RR] )
						}
#		max.b.increment <- abs( b.rater - b0 )
		conv1 <- max( abs( b.rater - b0 ) )
		it <- it+1
		cat("-")  
			}
	cat(" " , it , "Step(s) \n")	
    res <- list("b.rater" = b.rater , "se.b.rater" = sqrt( abs(-1/res$d2 ) ) , 
			"ll" = sum(res$ll0)  , "brc" = brc 
				)
    return(res)
			}	

#####################################################
# estimation of slope parameter for rater
.rm.facets.est.a.rater <- function( b.item , b.rater , Qmatrix , tau.item ,
        VV , K , I , TP , a.item , a.rater , item.index , rater.index ,
        n.ik , numdiff.parm=.001 , max.b.increment=1,theta.k , msteps ,
		mstepconv , a.rater.center , a.rater.fixed ){
    h <- numdiff.parm
	diffindex <- rater.index
	RR <- length(b.rater)
	cat("  M steps a.rater parameter  |")
	it <- 0 ;	conv1 <- 1000	
	while( ( it < msteps ) & ( conv1 > mstepconv ) ){	
		a.rater0 <- a.rater	
		pjk <- .rm.facets.calcprobs2( b.item , b.rater , Qmatrix , tau.item ,
				VV , K , I , TP , a.item , a.rater , item.index , rater.index , theta.k,RR)
		pjk1 <- .rm.facets.calcprobs2( b.item , b.rater , Qmatrix , tau.item ,
				VV , K , I , TP , a.item , a.rater+h , item.index , rater.index , theta.k,RR)
		pjk2 <- .rm.facets.calcprobs2( b.item , b.rater , Qmatrix , tau.item ,
				VV , K , I , TP , a.item , a.rater-h , item.index , rater.index , theta.k,RR)
		# numerical differentiation			
		res <- .rm.numdiff.index( pjk , pjk1 , pjk2 , n.ik , diffindex , 
				max.increment=max.b.increment , numdiff.parm )					
		a.rater <- a.rater + res$increment
		a.rater[ a.rater < .05 ] <- .05
#		a.rater <- a.rater - mean(a.rater ) + 1
		if ( ! is.null( a.rater.fixed) ){		
			ind <- which( ! is.na( a.rater.fixed  ) )
		    a.rater[ind] <- a.rater.fixed[ind]
            res$d2[ ind ] <- -1E10			
								}

		if ( a.rater.center){
			b1 <- mean( log( a.rater ) )
			a.rater <- a.rater / exp( b1 )
							}
		conv1 <- max( abs( a.rater - a.rater0 ) )
		it <- it+1
		cat("-") 
			}
	cat(" " , it , "Step(s) \n")	
    res <- list("a.rater" = a.rater , "se.a.rater" = sqrt( abs(-1/res$d2 )) , 
			"ll" = sum(res$ll0) )
    return(res)
			}				
			
####################################################################
# general function for numerical differentiation
# diffindex aggregates across super items
.rm.numdiff.index <- function( pjk , pjk1 , pjk2 , n.ik , diffindex , 
		max.increment , numdiff.parm , eps=10^(-80) ){					
	h <- numdiff.parm
    an.ik <- aperm( n.ik , c(2,3,1) )
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

		
		
		
		