

##########################################################################
# estimation of item intercepts
.mml.3pl.est.intercepts <- function( max.increment , np , est.xsi.index0 , 
		Msteps , nitems , A , AXsi , B , xsi , guess , theta , nnodes , maxK ,
		progress , itemwt , indexIP.no , indexIP.list2 ,
		ItemScore , fac.oldxsi , rprobs , xsi.fixed , convM , rprobs0 ,
		n.ik , N.ik  , xsi.prior ){
      converge <- FALSE
      Miter <- 1
	  eps <- 1e-10
	  oldxsi <- xsi
      old_increment <- rep( max.increment , np )
      est.xsi.index <- est.xsi.index0
      while ( !converge & ( Miter <= Msteps ) ) {	
        #      xbar2 <- xxf <- xbar <- rep(0,np)
        # Only compute probabilities for items contributing to param p
        if (Miter > 1){ 
          res.p <- .mml.3pl.calc_prob.v5( iIndex=1:nitems , A=A , AXsi=AXsi , B=B , 
                                 xsi=xsi , theta=theta , nnodes=nnodes, maxK=maxK , 
								 guess=guess )					
          rprobs <- res.p[["rprobs"]]
          rprobs0 <- res.p$rprobs0		  
        }

        res <- .mml.3pl.calc_exp_TK3( rprobs , A , np , est.xsi.index , itemwt ,
                             indexIP.no , indexIP.list2 , rprobs0 , guess ,
							 n.ik , N.ik )
        xbar <- res$xbar
        xbar2 <- res$xbar2
        xxf <- res$xxf
        ItemScore <- res$iscore		
        
        # Compute the difference between sufficient statistic and expectation
        diff <- as.vector(ItemScore) - xbar
        #Compute the Newton-Raphson derivative for the equation to be solved
        deriv <- xbar2 - xxf 	

		#***********************
		# xsi prior
		  if ( ! is.null(xsi.prior) ){
			  h <- .0001	
			  d0  <- log( stats::dnorm( xsi , mean=xsi.prior[,1] , 
						       sd=xsi.prior[,2] ) + eps)
			  d0p  <- log( stats::dnorm( xsi + h , mean=xsi.prior[,1] , 
						       sd=xsi.prior[,2] ) + eps)
			  d0m  <- log( stats::dnorm( xsi - h , mean=xsi.prior[,1] , 
						       sd=xsi.prior[,2] ) + eps)
			  d1 <- ( d0p - d0 ) / h
			  d2 <- ( ( d0p - d0 ) - ( d0 - d0m ) ) / h^2		
              diff <- diff + d1
              deriv <- deriv + d2			  
								}			
		#************************
		
        increment <- diff*abs(1/( deriv + 10^(-20) ) )
        if ( !is.null( xsi.fixed) ){ increment[ xsi.fixed[,1] ] <- 0 } 
        #!!!	  necessary to include statement to control increment?
        ci <- ceiling( abs(increment) / ( abs( old_increment) + 10^(-10) ) )
        increment <- ifelse( abs( increment) > abs(old_increment)  , 
                             increment/(2*ci) , 
                             increment )
        
        old_increment <- increment
        
        
        ##**SE
        se.xsi <- sqrt( 1 / abs(deriv) )
        if ( ! is.null( xsi.fixed) ){ se.xsi[ xsi.fixed[,1] ] <- 0 } 
        ##**

        xsi <- xsi+increment   # update parameter p
        #	  est.xsi.index <- which( abs(increment) > convM )
        if ( max(abs(increment)) < convM ) { converge <- TRUE }
        Miter <- Miter + 1						
        
        # stabilizing the algorithm | ARb 2013-09-10
        if (fac.oldxsi > 0 ){
          xsi <-  (1-fac.oldxsi) * xsi + fac.oldxsi *oldxsi
        }	  	  
        
        # progress bar
        if (progress){ 
          #        cat( paste( rep("-" , sum( mpr == p ) ) , collapse="" ) )
          cat("-") ; utils::flush.console()
        }
      } # end of all parameters loop
	  
	  res <- list( "xsi" = xsi , "se.xsi" = se.xsi )
	  return(res)
		}
#######################################################################


###########################################################
# faster version of calc_exp_TP3
.mml.3pl.calc_exp_TK3 <- function( rprobs , A , np , est.xsi.index , itemwt ,
	indexIP.no , indexIP.list2 , rprobs0 , guess , n.ik , N.ik){
	CC <- dim(rprobs)[2]
	TP <- dim(rprobs)[3]
	NXSI <- dim(A)[3]
	NI <- dim(A)[1]
	# restructure rprobs and AL	
	AL <- matrix( A , nrow=NI*CC , ncol=NXSI )	
	AL[ is.na(AL) ] <- 0
	rprobsL <- matrix( rprobs , nrow=NI*CC , ncol=TP )
	rprobsL[ is.na(rprobsL) ] <- 0
	rprobsL0 <- matrix( rprobs0 , nrow=NI*CC , ncol=TP )
	rprobsL0[ is.na(rprobsL0) ] <- 0
	# use nik
	nik <- as.vector(n.ik)
	ni <- as.vector(N.ik)

	#******
	# Call .mml.3pl function  (define CALL)	
	res <- .Call("mml3pl_tam_calcexp" ,  np , rprobsL , AL ,	indexIP.no , 
			    indexIP.list2 , est.xsi.index , CC , itemwt , rprobsL0 , 
			    guess , nik , ni , PACKAGE="TAM")
	return(res)
		}
		