
########################################
# tam.mml.3pl estimate item slopes
.mml.3pl.est.slopes <- function( max.increment , np , 
			Msteps , nitems , A , AXsi , B , xsi , guess , theta , nnodes , maxK ,
			progress ,ItemScore , fac.oldxsi , rprobs , xsi.fixed , convM , rprobs0 ,
			n.ik , N.ik , gammaslope , E , FdesM , dimFdes ,
			gammaslope.fixed , gammaslope.prior , maxgamma = 9.99 , Edes ){
	  	  
	  if (progress){ cat("\nM Step Slopes       |"); utils::flush.console() }
		eps <- 1e-10
		Nlam <- length(gammaslope)
		Xlambda00 <- Xlambda0 <- Xlambda <- gammaslope
        Xlambda.fixed <- gammaslope.fixed		
      msteps <- Msteps
	  parchange <- 1
	  oldfac <- fac.oldxsi
	  iter <- 1
	  	  

		  
	while( ( iter <= msteps ) & ( parchange > convM)  ){		
	    Xlambda0 <- gammaslope <- Xlambda		
# a0 <- Sys.time()		
#		B <- .mml.3pl.computeB( E , gammaslope )
		B <- .mml.3pl.computeB.v2( Edes , gammaslope , E )
# cat(" +++ compute B") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1		

 
		# calculate probabilities
		res <- .mml.3pl.calc_prob.v5(iIndex=1:nitems, A, AXsi, B, xsi, theta, 
					nnodes, maxK, recalc=TRUE , guess=guess )
		rprobs <- res$rprobs
		rprobs0 <- res$rprobs0	
# cat(" +++ calc prob") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1				
		# init derivatives		
		d2.b <- d1.b <- rep(eps,Nlam)			
		res <- .Call("mml3_slca_deriv" ,FdesM , dimFdes , gammaslope , as.vector(rprobs) ,
			as.vector(n.ik) , as.vector(N.ik) , guess , as.vector(rprobs0)  ,
			PACKAGE="TAM")   
		d1.b <- res$d1b
		d2.b <- res$d2b	
# cat(" +++ calc slca deriv") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1					
		# prior distribution for gammaslope
		  if ( ! is.null(gammaslope.prior) ){
			  h <- .0001	
			  if ( ncol(gammaslope.prior) == 2 ){
				  d0  <- log( stats::dnorm( Xlambda , mean=gammaslope.prior[,1] , 
								   sd=gammaslope.prior[,2] ) + eps)
				  d0p  <- log( stats::dnorm( Xlambda + h , mean=gammaslope.prior[,1] , 
								   sd=gammaslope.prior[,2] ) + eps)
				  d0m  <- log( stats::dnorm( Xlambda - h , mean=gammaslope.prior[,1] , 
								   sd=gammaslope.prior[,2] ) + eps)
											}
			  if ( ncol(gammaslope.prior) == 4 ){
				  d0  <- log( msm::dtnorm( Xlambda , mean=gammaslope.prior[,1] , 
								   sd=gammaslope.prior[,2] , lower=gammaslope.prior[,3] ,
								   upper=gammaslope.prior[,4] ) + eps)
				  d0p  <- log( msm::dtnorm( Xlambda + h , mean=gammaslope.prior[,1] , 
								   sd=gammaslope.prior[,2] , lower=gammaslope.prior[,3] ,
								   upper=gammaslope.prior[,4] ) + eps)
				  d0m  <- log( msm::dtnorm( Xlambda - h , mean=gammaslope.prior[,1] , 
								   sd=gammaslope.prior[,2] , lower=gammaslope.prior[,3] ,
								   upper=gammaslope.prior[,4] ) + eps)
											}			  
			  
			  d1 <- ( d0p - d0 ) / h
			  d2 <- ( ( d0p - d0 ) - ( d0 - d0m ) ) / h^2		
              d1.b <- d1.b + d1
              d2.b <- d2.b + d2			  
								}						
		increment <-   d1.b / ( abs( d2.b + 10^(-10) ) )
		increment[ is.na(increment) ] <- 0		
		increment <- ifelse(abs(increment)> max.increment, 
					sign(increment)*max.increment , increment )											
		max.increment <- max(abs(increment)) / .98
		Xlambda <- Xlambda + increment
		se.Xlambda <- sqrt( 1 / abs( d2.b+10^(-10)) )

		Xlambda <- ifelse( Xlambda > maxgamma , maxgamma , Xlambda )
		Xlambda <- ifelse( Xlambda < - maxgamma , - maxgamma , Xlambda )
		if ( ! is.null(gammaslope.prior) ){
			if ( ncol(gammaslope.prior) == 4 ){
					  Xlambda <- ifelse( Xlambda < gammaslope.prior[,3] ,
										gammaslope.prior[,3] + 1.3* h , Xlambda )
					  Xlambda <- ifelse( Xlambda > gammaslope.prior[,4] ,
										gammaslope.prior[,4] - 1.3* h , Xlambda )
								}
					}
				
		if ( ! is.null( Xlambda.fixed) ){
			Xlambda[ Xlambda.fixed[,1] ] <- Xlambda.fixed[,2]
			se.Xlambda[ Xlambda.fixed[,1] ] <- 0		
				}
				
		if (progress){ cat("-") ; utils::flush.console() }		
		iter <- iter + 1
		parchange <- max( abs(Xlambda0-Xlambda))
# cat(" +++ calc rest") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1							
			}
		#********* end algorithm
		if (oldfac > 0 ){
			Xlambda <- oldfac*Xlambda00 + ( 1 - oldfac ) *Xlambda
						}				
		max.increment <- max( abs( Xlambda - Xlambda00 ))
		res <- list("gammaslope" = Xlambda , "se.gammaslope" = se.Xlambda , 
		         "max.increment.b"=max.increment , 
				 "gammachange"= max( abs( Xlambda00 - Xlambda) ) )
			}
################################################################################			
# centering gammaslope vector
.mml.3pl.gammaslope.center <- function( gammaslope , gammaslope.center.index  ,
				gammaslope.center.value  ){														
		if ( ! is.null( gammaslope.center.index ) ){
			M <- max( gammaslope.center.index )
		for (mm in 1:M){
			ind.mm <- which( gammaslope.center.index == mm )
			IM <- length(ind.mm)
			rmm <- gammaslope[ ind.mm ] - mean( gammaslope[ind.mm] )
			gammaslope[ ind.mm ] <- gammaslope.center.value[mm] / IM + rmm
						}
						}
		return(gammaslope)
				}
###########################################################