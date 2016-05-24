

################################################################
################################################################
################################################################
tam.jml.WLE <-
  function ( tamobj, resp , resp.ind, A, B, nstud, nitems, maxK, convM, 
             PersonScores, theta, xsi, Msteps, WLE=FALSE ,
			 theta.fixed=NULL){
    
    AXsi <- matrix(0, nrow=nitems, ncol=maxK) 
    B1 <- B[,,1]
    BB <- array (0, dim=c(nitems,maxK))
    BBB <- array (0, dim=c(nitems,maxK)) 
    B_bari <- array(0,dim=c(nstud, nitems))
    BB_bari <- array(0, dim=c(nstud, nitems))
    BBB_bari <- array(0,dim=c(nstud, nitems))
    B_Sq <- array(0,dim=c(nstud, nitems))
    B2_B <- array(0,dim=c(nstud, nitems))
    B_Cube <- array(0,dim=nstud)
    
    #Calculate Axsi. Only need to do this once for ability estimates.
    for (i in 1:nitems) {
      for (k in 1:maxK){
        AXsi[i,k] <- ( A[i,k,] %*% xsi )
      }
    }
    cat("\n MLE/WLE estimation        |")
    
    #Compute WLE
    #similar to the M step in the tam function, but each student's theta is now one node.
    convergeWLE <- FALSE
    iterWLE <- 0
    BB <- B1^2
    BBB <- BB * B1  
    # BB	[ nitems , maxK ]
    # BBB	[ nitems , maxK ]
    maxChangeWLE <- 0
    thetaOld <- theta
    while (!convergeWLE & ( iterWLE <= Msteps ) ) {  
      resWLE <- calc_prob.v5(iIndex = 1:nitems , A , AXsi , 
                             B , xsi , theta , nstud, maxK , recalc=FALSE )      	
      rprobsWLE <- resWLE[["rprobs"]] 
      rprobsWLE[ is.na( rprobsWLE ) ] <- 0
      
      
      B_bari <- B1[,1] * rprobsWLE[ , 1 , ]
      BB_bari <- BB[,1] * rprobsWLE[ , 1 , ]
      for (kk in 2:maxK){ 
        B_bari <- B_bari + B1[,kk]*rprobsWLE[,kk,] 
        BB_bari <- BB_bari + BB[,kk] * rprobsWLE[ , kk , ]			
      }
      B_bari <- t(B_bari) * resp.ind 
      
      # B_bari.OLD <- sapply(1:nitems, function(i) colSums(B1[i,] * rprobsWLE[i,,] , na.rm = TRUE)) * resp.ind
      # B1		[ nitems , maxK ]
      # rprobsWLE [ nitems , maxK , nstud ]
      # resp.ind  [ nstud , nitems ]	
      # colSums(B1[i,] * rprobsWLE[i,,] , na.rm = TRUE))
      #	-> colSums( [ maxK , nstud ] ) = [nstud]
      # B_bari	[ nstud , nitems ]
      
      BB_bari <- t(BB_bari ) * resp.ind
      B_Sq <- B_bari^2
      expected <- rowSums(B_bari, na.rm=TRUE)
      err <- rowSums(BB_bari, na.rm=TRUE) - rowSums(B_Sq, na.rm=TRUE)  #sum over the items
      err_inv <- abs(1/err)
      scores <- PersonScores - expected    
      
      if (WLE) {
        BBB_bari <- BBB[,1] * rprobsWLE[ , 1 , ]
        for (kk in 2:maxK){ 
          BBB_bari <- BBB_bari + BBB[,kk] * rprobsWLE[ , kk , ]			
        }	  
        BBB_bari <- t(BBB_bari ) * resp.ind				
        B2_B <- BB_bari*B_bari   
        B_Cube <- B_Sq*B_bari
        warm <- -3*B2_B + 2*B_Cube + BBB_bari
        warmadd <- rowSums(warm, na.rm=TRUE)                 #sum over the items
        warmaddon <- err_inv*warmadd
        scores <- scores + warmaddon/2.0      
      }     
      increment <-  err_inv*scores
      
      if (maxChangeWLE < max(abs(increment))) {
        maxChangeWLE <- max(abs(increment))
      }
      
      increment[abs(increment)>3] <- sign(increment[abs(increment)>3])*3
      
	  if ( ! is.null(theta.fixed ) ){
	      increment[ theta.fixed[,1] ] <- 0
						}
	  
      theta <- theta + increment
      if ( max(abs(increment)) < convM ) {
        convergeWLE <- TRUE
      }
      iterWLE <- iterWLE + 1 
      cat( "-"  ) 
      utils::flush.console()
    }  # end of Newton-Raphson   
    cat("\n")
    meanChangeWLE <- mean(theta - thetaOld)
    #standard errors of theta estimates
    errorWLE <- sqrt(err_inv)
	  if ( ! is.null(theta.fixed ) ){
	      errorWLE[ theta.fixed[,1] ] <- 0
						}
    
    res <- list( "theta" = theta , "errorWLE" = errorWLE, "meanChangeWLE" = meanChangeWLE)
    return (res)
  }
