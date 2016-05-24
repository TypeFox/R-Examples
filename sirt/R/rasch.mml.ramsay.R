

######################################################
# probability ramsay qm
.prob.ramsay <- function( theta , b ,  fixed.K = 3 + 0*b , pow = 1){
	XXT <- matrix( theta , nrow=length(theta) , ncol=length(b))
	XXb <- matrix( b , nrow=length(theta) , ncol=length(b) , byrow=T )
	XX2 <- exp( (exp( XXT ))^pow / exp(XXb) )
	XX2[ XX2 > 10^15 ] <- 10^15
	KM <- matrix( fixed.K , nrow=length(theta) , ncol=length(b) , byrow=T )	
	pm <- XX2 / ( KM + XX2 )
	pm[ pm == Inf ] <- 1
    return(pm)
    }
#-----------------------------------------------------------------	


##################################################################
# simulate ramsay quotient model
##NS export(sim.qm.ramsay)
sim.qm.ramsay <- function( theta , b , K  ){ 
    N <- length(theta)
    I <- length(b)
    p1 <- exp( outer( theta , 1/b ) )
    p1 <- p1 / ( outer( rep(1,length(theta)) , K ) + p1 )
    dat <- 1* ( p1 > matrix( stats::runif( N*I ) , ncol=I) )
	colnames(dat) <- paste( "I" , substring(100+1:I,2) , sep="")	
    return(dat )
    }
##################################################################	



#*************************************************************************************
# E Step Ramsay QM 
.e.step.ramsay <- function( dat1 , dat2 , dat2.resp , theta.k , pi.k , I , n , b , fixed.K ,
			group , pow.qm , ind.ii.list , pjk = NULL){
    #...................................                    
    # arrange groups
    if ( is.null(group) ){ group <- rep( 1 , nrow(dat1)) }
    G <- length( unique( group) )    
    # probabilities of correct item at theta_k
	if ( is.null( pjk )){
		pjk <- .prob.ramsay( theta.k , b ,  fixed.K , pow = pow.qm )
					}
	TP <- dim(pjk)[1]
					
	#***
	# array notation of probabilities
#	pjkL <- array( NA , dim=c(2 , nrow(pjk) , ncol(pjk) ) )
#	pjkL[1,,] <- 1 - pjk
#	pjkL[2,,] <- pjk
#	f.yi.qk <- matrix( 1 , nrow(dat2) , length(theta.k) )
#	for (ii in 1:ncol(dat2)){
#	#	ii <- 1
#		ind.ii <- ind.ii.list[[ii]]
#		if ( length(ind.ii) == nrow(dat2) ){ 
#			f.yi.qk <- f.yi.qk * pjkL[ dat2[,ii]+1 , ,ii]
#					} else { 
#			f.yi.qk[ind.ii,] <- f.yi.qk[ind.ii,] * pjkL[ dat2[ind.ii,ii]+1 , ,ii]
#							}
#					}
        pjkt <- t(pjk)
		pjkL <- array( NA , dim=c( I , 2 , TP  ) )
		pjkL[,1,] <- 1 - pjkt
		pjkL[,2,] <- pjkt	
		probsM <- matrix( aperm( pjkL , c(2,1,3) ) , nrow=I*2 , ncol=TP )
		f.yi.qk <- mml_calc_like( dat2=dat2 , dat2resp = dat2.resp , 
					probs = probsM )$fyiqk
#cat("   calc likelihood") ; vv1 <- Sys.time(); print(vv1-vv0) ; vv0 <- vv1				

	#******	
    f.qk.yi <- 0 * f.yi.qk
    if ( G==1 ){ pi.k <- matrix( pi.k , ncol=1 ) }
    for ( gg in 1:G){ 
        f.qk.yi[ group == gg , ] <- f.yi.qk[ group == gg , ] * outer( rep( 1 , nrow(dat2[ group==gg,]) ) , pi.k[,gg] )
                    }
    f.qk.yi <- f.qk.yi / rowSums( f.qk.yi )
    # expected counts at theta.k
    n.k <- matrix( 0 , TP , G )
    r.jk <- n.jk <- array( 0 , dim=c( ncol(dat2) , TP , G) )
    ll <- rep(0,G)
    for (gg in 1:G){ 
		ind.gg <- which( group == gg )
	    res <- mml_raschtype_counts( dat2=dat2[ind.gg,] , dat2resp=dat2.resp[ind.gg,] , 
					dat1=dat1[ind.gg,2] , fqkyi=f.qk.yi[ind.gg,] ,
					pik=pi.k[,gg] , fyiqk=f.yi.qk[ind.gg,]  )
		n.k[,gg] <- res$nk
		n.jk[,,gg] <- res$njk
        r.jk[,,gg] <- res$rjk
		ll[gg] <- res$ll

        }
    res <- list( "n.k" = n.k , "n.jk" = n.jk , "r.jk" = r.jk , "f.qk.yi" = f.qk.yi , "pjk" = pjk  ,
            "f.yi.qk" = f.yi.qk , "ll" = sum(ll) )
    return(res)
    }
#*************************************************************************************



#*************************************************************************************
# M Step Ramsay's quotient model                                                               
.m.step.ramsay <- function( theta.k , b , n.k , n , n.jk , r.jk , I , conv1 , constraints ,
			mitermax , pure.rasch ,  trait.weights , fixed.K , 
			designmatrix = designmatrix , group = group ,
            numdiff.parm=numdiff.parm , pow.qm = 1){
    abs.change <- 1
    miter <- 0
    # group estimation
    G <- ncol(n.k)
    # number of subjects within groups
    n <- colSums(n.k)
    # update pi.k 
    pi.k <- matrix( 0 , nrow=length(theta.k) , ncol=G)
    for (gg in 1:G){ 
        if (is.null( trait.weights) | G > 1 ){ 
                pi.k[,gg] <- n.k[,gg] / n[gg]  } else  {
                 pi.k <- stats::dnorm( theta.k)
                 pi.k <- pi.k/sum(pi.k) 
                 pi.k <- matrix( pi.k , ncol=1 )
                        }
               }
    #*****
    # begin loop
    eps <- numdiff.parm 
    # eps <- .00005
    h <- h1 <- eps
    h2 <- 1 + 2*eps
    while( abs.change > conv1 & miter < mitermax ){
        b0 <- b     
            pjk <- .prob.ramsay( theta.k , b , fixed.K , pow=pow.qm )
            pjk <- ( pjk + h1 ) / h2
            pjk.M <- t(pjk)
            qjk.M <- 1 - pjk.M
            pjk1 <- .prob.ramsay( theta.k , b + h , fixed.K , pow=pow.qm )
            pjk1 <- ( pjk1 + h1 ) / h2
            pjk1.M <- t(pjk1)
            qjk1.M <- 1 - pjk1.M
            pjk2 <- .prob.ramsay( theta.k , b - h , fixed.K, pow=pow.qm )
            pjk2 <- ( pjk2 + h1 ) / h2
            pjk2.M <- t(pjk2)
            qjk2.M <- 1 - pjk2.M
            # update item parameter
            # first order derivative
            # f(x+h) - f(x-h) = 2* f'(x) * h
            ll0 <- ll1 <- ll2 <- matrix( 0 , nrow=nrow(n.jk) , ncol=G)
            for ( gg in 1:G){
                ll0[,gg] <- rowSums( r.jk[,,gg] * log( pjk.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk.M ) )
                ll1[,gg] <- rowSums( r.jk[,,gg] * log( pjk1.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk1.M ) )
                ll2[,gg] <- rowSums( r.jk[,,gg] * log( pjk2.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk2.M ) )
                            }
            # sum across all groups
            ll0 <- rowSums(ll0)
            ll1 <- rowSums(ll1)            
            ll2 <- rowSums(ll2)            
            d1 <- ( ll1 - ll2  ) / ( 2 * h )    # negative sign?
            # second order derivative
            # f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
            d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2
        # change in item difficulty
		eps <- 10^(-15)
        b.change <- - d1 / d2
		b.change[ abs(d2) < eps ] <- 0
		bstep <- .5
		b.change <- ifelse( abs( b.change ) > bstep , sign(b.change)*bstep , b.change )
        b <- b + b.change
        # linear parameter constraints
        if ( ! is.null( designmatrix )){ 
            mod <- stats::lm( b ~ 0 + designmatrix  )
            b <- stats::fitted(mod)
                }
        # last item is the sum of all other item difficulties
        center <- is.null(constraints) 
        if ( !is.null( constraints) ){   b[ constraints[,1] ] <- constraints[,2] }

        abs.change <- max( abs( b0 - b ) )
        miter <- miter+1
        }   #*** end loop
       pjk <- .prob.ramsay( theta.k , b , fixed.K )
        # calculate log likelihood
       ll <- sapply( 1:G , FUN = function(gg){ 
                    sum( rowSums( r.jk[,,gg]  * t( log( pjk ) ) + ( n.jk[,,gg] - r.jk[,,gg] ) * ( t( log(1 - pjk)) ) ) )  
                            } )
       ll <- sum(ll)
       list( "b" = b , "pi.k" = pi.k , "ll" = ll , "miter" = miter , "center" = center , "G" = G ,
			"se.b" = sqrt(1/abs(d2)) )
    }
#*************************************************************************************
# estimation of K parameter in Ramsay's QM
.mml.ramsay.est.K <- function( theta.k , b , fixed.a , fixed.c , fixed.d ,
					fixed.K , pjk , alpha1 , alpha2 , h , G , I , r.jk , n.jk , est.K ,
					min.K , max.K , iter , pow.qm ){
					
					#.prob.ramsay <- function( theta , b ,  fixed.K = 3 + 0*b)			
					#****
					# kk
					pjk <- .prob.ramsay( theta.k , b , fixed.K , pow = pow.qm )
					pjk <- ( pjk + .000000005 ) / 1.00000001 
					pjk.M <- t(pjk) ;  qjk.M <- 1 - pjk.M
					pjk1 <- .prob.ramsay( theta.k , b , fixed.K+h , pow = pow.qm  )
					pjk1 <- ( pjk1 + .000000005 ) / 1.00000001 
					pjk1.M <- t(pjk1) ; qjk1.M <- 1 - pjk1.M
					pjk2 <- .prob.ramsay( theta.k , b , fixed.K-h , pow = pow.qm  )
					pjk2 <- ( pjk2 + .000000005 ) / 1.00000001 
					pjk2.M <- t(pjk2) ; qjk2.M <- 1 - pjk2.M

					ll0 <- ll1 <- ll2 <- matrix(0, I , G)
				# was ist hier das G? => G ist hier Anzahl der Gruppen
					for (gg in 1:G){ 
						ll0[,gg] <- rowSums( r.jk[,,gg] * log( pjk.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk.M ) )
						ll1[,gg] <- rowSums( r.jk[,,gg] * log( pjk1.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk1.M ) )
						ll2[,gg] <- rowSums( r.jk[,,gg] * log( pjk2.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk2.M ) 	)		
								   }  
					ll0 <- rowSums(ll0)
					ll1 <- rowSums(ll1)
					ll2 <- rowSums(ll2)
					# aggregate with respect to estimation of a
					a1 <- stats::aggregate( cbind( ll0 , ll1 , ll2 ) , list(est.K) , sum , na.rm=T)		
					a1 <- a1[ a1[,1] > 0 , ]										
					ll0 <- a1[,2]
					ll1 <- a1[,3]
					ll2 <- a1[,4]			
					d1 <- ( ll1 - ll2  ) / ( 2 * h )    # negative sign?
					# second order derivative
					# f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
					d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2
					# change in item difficulty
					K.change <- - d1 / d2
#					a.change <- ifelse( abs( a.change ) > .2 , .2*sign(a.change) , a.change )              
					# dampening parameter as in tam
#					ci <- ceiling( abs(K.change) / ( abs( old_increment.d) + 10^(-10) ) )
					K.change <- ifelse( abs( K.change) > 1 , 1 , K.change )
					K.change <- K.change[ match( est.K , a1[,1] ) ]						
					if ( any( est.K == 0 ) ){
						K.change[ est.K == 0 ] <- 0
									}											
					fixed.K <- fixed.K + K.change
					fixed.K[ fixed.K < min.K ] <- min.K
					fixed.K[ fixed.K > max.K ] <- max.K
				    res <- list("fixed.K" = fixed.K , "se.K" = sqrt( 1 /abs(d2) ) )
					return( res )						
											}
#####################################################################




