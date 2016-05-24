
  
###########################################################################  
tam.jml.xsi <-
  function ( resp , resp.ind, A, B, nstud, nitems, maxK, convM, 
             ItemScore, theta, xsi, Msteps, pweightsM,
             est.xsi.index
  ){
    
    #Update item parameters
    
    AXsi <- matrix(0, nrow=nitems, ncol=maxK) 
    r <- matrix(0,nrow=nitems,ncol=maxK) 
    rr <- array(0,dim=c(nitems,maxK,maxK))
    AA <- array (0, dim=c(nitems,maxK,maxK))
    
    A.0 <- A
    A.0[ is.na(A.0) ] <- 0  	
	
    maxChangeP <- 0
    errorP <- rep(0, max(est.xsi.index))
    convergeAllP <- FALSE
    p_loop <- est.xsi.index
	PP1 <- dim(A)[3]
	nonest.xsi.index <- setdiff( seq(1,PP1) , est.xsi.index )
    convergeP <- rep(FALSE,max(est.xsi.index))
    
    # begin loop 
    iterP <- 1
	old_xsi <- xsi
    old_increment <- rep(5,max(p_loop))
    cat(" Item parameter estimation |")
    while (!convergeAllP & ( iterP <= Msteps ) ) {  
      res.p <- calc_prob.v5( iIndex = 1:nitems , A , AXsi , 
                             B , xsi , theta , nstud, maxK , TRUE )        	
      rprobs <- res.p[["rprobs"]]               
      
      #compute probability weights, summed over students, so that there is no cycling through students for parameter estimation (p loop)
      for (k1 in 1:maxK) {
        r[,k1] <- colSums(t(rprobs[,k1,]) * resp.ind * pweightsM, na.rm=TRUE)
        for (k2 in 1:maxK) {
          rr[,k1,k2] <- colSums(t(rprobs[,k1,]*rprobs[,k2,]) * resp.ind * pweightsM, na.rm=TRUE)
        }
      }


		  A_Sq <- AA_bari <- A_bari <- matrix( 0 , PP1 , nitems )
		  
		  for (kk in 1:maxK){ 
			A_bari <- A_bari + t( A.0[ , kk , ] * r[ , kk ] )
			AA_bari <- AA_bari + t( A.0[ , kk , ]^2 * r[ , kk ] )		
		  }
		  for (kk1 in 1:maxK){ 
			for (kk2 in 1:maxK){ 
			  A_Sq <- A_Sq + t( A.0[,kk1,] * A.0[,kk2,] * rr[ , kk1 , kk2 ] )	
			}
		  }		


        
#        expected <- sum (A_bari, na.rm=TRUE) # sum over items 
#        err <- sum (AA_bari - A_Sq, na.rm=TRUE)   #sum over the items  
		  expected <- rowSums (A_bari, na.rm=TRUE) # sum over items
		  err <- rowSums(AA_bari - A_Sq, na.rm=TRUE)   #sum over the items		
		
		
		  err_inv <- abs (1/( abs(err) + 10^(-10) ))		  
		  scores <- ItemScore * ( ! convergeP ) - expected
		  increment <-  err_inv*scores
		  ci <- ceiling( abs(increment) / ( abs( old_increment) + 10^(-10) ) )
		  increment <- ifelse( abs( increment) > abs(old_increment)  , 
							   increment/(2*ci) , 
							   increment )
		  increment[ nonest.xsi.index ] <- 0
		  xsi <- xsi + increment
		  old_increment <- increment    
		  
		  # xsi[ nonest.xsi.index ] <- old_xsi[ nonest.xsi.index ]
		  
		  errorP <- sqrt(err_inv)
		  convergeP[ abs(increment) < convM ] <- TRUE
		  flush.console()
		  iterP <- iterP + 1 
		  p_loop <- est.xsi.index[convergeP[est.xsi.index]==FALSE]
		  convergeAllP <- (sum(convergeP[est.xsi.index]) == length(est.xsi.index))  	  
		  cat("-")  
		  
    } # end of all parameters convergence
    
    res <- list( "xsi" = xsi , "errorP" = errorP, "maxChangeP" = max(abs( xsi - old_xsi ) ) )
    return (res)  
  }

######################################################################  
