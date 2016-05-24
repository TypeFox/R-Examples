tam.wle <- function( tamobj, ... ){
  CALL <- match.call()
  if(class(tamobj) == "tam.mml"){
    res <- tam.mml.wle2( tamobj, ...)
  }
  if(class(tamobj) == "tamaan"){
    res <- tam.mml.wle2( tamobj, ...)
  }
  if(class(tamobj) == "tam.jml"){
    res <- tam.jml.WLE( tamobj, ...)
  }
  attr(res,"call") <- CALL
  return( res )
}

################################################################
################################################################
################################################################
tam.mml.wle <-
  function( tamobj, score.resp=NULL , WLE=TRUE , adj=.3 , Msteps=20 , 
            convM = .0001 , progress=TRUE ,
			output.prob=FALSE ){
    #########################################################
    # INPUT:
    # tamobj ... result from tam analysis
    # (WLE = TRUE) will produce WLE. Otherwise it will be MLE
    # 
    #########################################################
    #  adj <- 0.3
    #  Msteps <- 20
    #  convM <- .0001
	CALL <- match.call()
    B <- tamobj$B
    A <- tamobj$A
    nitems <- tamobj$nitems
    xsi <- ( tamobj$xsi )[,1]
    nstud <- tamobj$nstud
    AXsi <- tamobj$AXsi
    ndim <- tamobj$ndim
    maxK <- tamobj$maxK
    resp <- tamobj$resp
	pweights <- tamobj$pweights
    if ( ! is.null( score.resp) ){
      resp <- score.resp
      nstud <- nrow(resp)
      tamobj$resp.ind <- 1 - is.na(resp)
      tamobj$pid <- 1:nstud	
	  pweights <- 1 + 0 * tamobj$pid	  
    }  
    resp[ is.na(resp) ] <- 0  
    resp.ind <- tamobj$resp.ind  
    col.index <- rep( 1:nitems , each = maxK )
    cResp <- resp[ , col.index  ]*resp.ind[ , col.index ]
    cResp <- 1 * t( t(cResp) == rep(0:(maxK-1), nitems) )
    cB <- t( matrix( aperm( B , c(2,1,3) ) , nrow = dim(B)[3] , byrow = TRUE ) )
    cB[is.na(cB)] <- 0
    
    #Compute person sufficient statistics (total score on each dimension)
    PersonScores <- cResp %*% cB
    
    #Compute possible maximum score for each item on each dimension
    maxBi <- apply(B , 3 , rowMaxs , na.rm = TRUE)
    
    #Compute possible maximum score for each person on each dimension
    PersonMax <- resp.ind %*% maxBi
    PersonMax[ PersonMax == 0 ] <- 2 * adj
    
    #Adjust perfect scores for each person on each dimension
    PersonScores[PersonScores==PersonMax] <- PersonScores[PersonScores==PersonMax] - adj
    
    #Adjust zero scores for each person on each dimension
    PersonScores[PersonScores==0] <- PersonScores[PersonScores==0] + adj
    
    #Calculate Axsi. Only need to do this once.
    for (i in 1:nitems) {
      for (k in 1:maxK){
        AXsi[i,k] <- ( A[i,k,] %*% xsi )
      }
    }
    
    #Initialise theta (WLE) values for all students
    theta <- log(PersonScores/(PersonMax-PersonScores)) #log of odds ratio of raw score
    
    ######################################
    #Compute WLE
    #similar to the M step in the tam function, but each student's theta vector is now one node.
    converge <- FALSE
    Miter <- 0
    BB <- array (0, dim=c(nitems,maxK,ndim,ndim))
    BBB <- array (0, dim=c(nitems,maxK,ndim)) 
    for (i in 1:nitems) {
      for (k in 1:maxK) {
        BB[i,k,,] <- B[i,k,] %*% t(B[i,k,])
        BBB[i,k,] <- BB[i,k,,] %*% B[i,k,]
      }
    }
    increment <- array(0, dim=c(nstud,ndim))
    old_increment <- 3 + increment
    
    
    # Begin iterations
    while (!converge & ( Miter <= Msteps ) ) {  
      resWLE <- calc_prob.v5(iIndex = 1:nitems , A , AXsi , 
                             B , xsi , theta , nstud, maxK , recalc=FALSE )      	
      rprobsWLE <- resWLE[["rprobs"]] 
      B_bari <- array(0,dim=c(nstud, nitems, ndim))
      BB_bari <- array(0, dim=c(nstud, nitems, ndim, ndim))
      BBB_bari <- array(0,dim=c(nstud, nitems, ndim))
      for (d1 in 1:ndim) {
        B_bari[,,d1] <- sapply(1:nitems, function(i) colSums(B[i,,d1] * rprobsWLE[i,,] , na.rm = TRUE)) * resp.ind
        for (d2 in 1:ndim) {
          BB_bari[,,d1,d2] <- sapply(1:nitems, function(i) colSums(BB[i,,d1,d2] * rprobsWLE[i,,] , na.rm = TRUE)) *resp.ind
        }
        BBB_bari[,,d1] <- sapply(1:nitems, function(i) colSums(BBB[i,,d1] * rprobsWLE[i,,] , na.rm = TRUE)) *resp.ind  
      }
      
      B_Sq <- array(0,dim=c(nstud, nitems, ndim, ndim))
      B2_B <- array(0,dim=c(nstud, nitems, ndim))
      B_Cube <- array(0,dim=c(nstud, nitems, ndim))
      for (d1 in 1:ndim) {      
        B2_B[,,d1] <- 0
        B_Cube[,,d1] <- 0
        for (d2 in 1:ndim) {
          B_Sq[,,d1,d2] <- B_bari[,,d1]*B_bari[,,d2]
          B2_B[,,d1] <- B2_B[,,d1] + BB_bari[,,d1,d2]*B_bari[,,d2]   
          B_Cube[,,d1] <- B_Cube[,,d1] + B_Sq[,,d1,d2]*B_bari[,,d2]
        }
      }
      expected <- colSums(aperm(B_bari,c(2,1,3)))
      err <- colSums(aperm(BB_bari,c(2,1,3,4))) - colSums(aperm(B_Sq, c(2,1,3,4)))  #sum over the items
      if (ndim == 1) {
        #      err_inv <- apply(err,1,function(x) 1/x )
        err_inv <- 1 / err
      } else {
        #err_inv <- aperm(apply(err,1,solve),c(2,1))   
        err_inv <- aperm(apply(err,1,function(ee){
          ee1 <- ee		
          diag(ee1) <- diag(ee1) + 10^(-15)
          solve(ee1)	
        }
        ),c(2,1))            
      }
      err_inv <- array(abs(err_inv),dim=c(nstud,ndim,ndim))
      warm <- -3*B2_B + 2*B_Cube + BBB_bari
      warmadd <- colSums(aperm(warm,c(2,1,3)))  #sum over the items
      scores <- PersonScores - expected
      if (WLE) {
        warmaddon <- array(0,dim=c(nstud,ndim))
        for (d1 in 1:ndim) {
          warmaddon[,d1] <- 0
          for (d2 in 1:ndim) {
            warmaddon[,d1] <- warmaddon[,d1] + err_inv[,d1,d2]*warmadd[,d2]
          }
        }
        scores <- scores + warmaddon/2.0      
      }
      increment <- array(0, dim=c(nstud,ndim))
      for (d1 in 1:ndim) {
        increment[,d1] <- 0
        for (d2 in 1:ndim) {
          increment[,d1] <- increment[,d1] + err_inv[,d1,d2]*scores[,d2]
        }
      }
      # dampening the increment
      for ( d1 in 1:ndim){ 
        #	   increment[,d1] <- ifelse( abs(increment[,d1]) > 3 , sign( increment[,d1] )*3 , increment[,d1] )
        ci <- ceiling( abs(increment[,d1]) / ( abs( old_increment[,d1]) + 10^(-10) ) )	   
        increment[,d1] <- ifelse( abs( increment[,d1]) > abs(old_increment[,d1])  , 
                                  increment[,d1]/(2*ci) , 
                                  increment[,d1] )	   
        old_increment[,d1] <- increment[,d1] 
        #***
        # avoid NaNs in increment
        increment[,d1] <- ifelse( is.na(increment[,d1] ) , 0 , increment[,d1] )
        # increment[abs(increment)>3] <- sign(increment[abs(increment)>3])*3	
      }
      theta <- theta + increment
      if ( max(abs(increment)) < convM ) {
        converge <- TRUE
      }
      Miter <- Miter + 1 
	  if (progress){
      cat( paste( "Iteration in WLE/MLE estimation ", Miter, 
                  "  | Maximal change " , round( max(abs(increment)) , 4) , "\n" )  ) 
      utils::flush.console()
			}
    }  # end of Newton-Raphson
    
    #standard errors of theta estimates
    if (ndim == 1) {
      error <- apply(err_inv,1,sqrt) 
    } else {    
      error <- aperm(apply(sqrt(err_inv),1,diag), c(2,1))
    }
    
    # The output contains 
    #   Person Scores on the test, by dimension
    #   Person possible maximum score, by dimension (Each person could take 
    #    different items, so possible maximum could vary)
    #   WLE or MLE estimate, by dimension
    #   Standard errors of WLE/MLE estimates, by dimension
    
    if ( ndim> 1){
      colnames(error) <- paste0("error.Dim" , substring( 100+1:ndim , 2) )
    }
    res <- data.frame( "pid" = tamobj$pid , 
                       "N.items" = rowSums(resp.ind) , 
                       "PersonScores" = PersonScores, 
                       "PersonMax" = PersonMax, "theta" = theta , error )
    
    if (ndim==1){ colnames(res)[4:5] <- c("PersonMax" , "theta") }
    if (ndim>1){  
      colnames(res)[ 1:ndim + 2] <- paste0("PersonScores.Dim" , substring( 100+1:ndim , 2) )	
      ind <- grep( "theta" , colnames(res) )	
      colnames(res)[ind] <- 	paste0("theta.Dim" , substring( 100+1:ndim , 2) )	
    }
    ####################
    # correct personMax set theta and standard error to missing		
    # if there are no observations on one dimension
    ind1 <- grep("PersonMax" , colnames(res))
    check1 <- ( res[ , ind1 , drop=FALSE] == 2*adj )
    ind2 <- grep("theta" , colnames(res))
    D <- length(ind1)
    for (ii in 1:D){
      res[ check1[,ii] , ind2[ii] ] <- NA
    }
    ind2 <- grep("error" , colnames(res))
    for (ii in 1:D){
      res[ check1[,ii] , ind2[ii] ] <- NA
    }
    #***
    # WLE reliability
    if ( ndim==1 ){
      ind <- which( res$N.items > 0 )	  
	  WLE.rel <- WLErel(theta=theta[ind] , error = error[ind] , w = pweights[ind] )	  
	  if (WLE){ w1 <- "WLE" } else { w1 <- "MLE" }
	  if (progress){
		cat("----\n" , w1 ,"Reliability =" , round(WLE.rel,3) ,"\n" )
					}
      res$WLE.rel <- rep( WLE.rel , nrow(res) )
    }
    if ( ndim>1 ){
      cat("\n-------\n")
      for (dd in 1:ndim){
        #	dd <- 1
		ind1 <- paste0("theta.Dim" , substring( 100+1:ndim , 2))[dd]
        # v1 <- stats::var( res[, ind1 ] , na.rm=TRUE)
		ind2 <- paste0("error.Dim" , substring( 100+1:ndim , 2))[dd] 
        # v2 <- mean( res[, ind2]^2 , na.rm=TRUE)				
        #		v2 <- mean( error^2 )
        # res[ ,paste0("WLE.rel.Dim" , substring( 100+ dd , 2)) ]	<- h1 <- 1 - v2 / v1
		h1 <- WLErel( theta=res[,ind1] , error=res[,ind2] , w = pweights )
		res[ ,paste0("WLE.rel.Dim" , substring( 100+ dd , 2)) ]	<- h1
		if (WLE){ w1 <- "WLE" } else { w1 <- "MLE" }
        cat(paste0(w1 , " Reliability (Dimension" , dd , ") = " , round(h1,3) ) , "\n" )
        #	  res$WLE.rel <- rep( WLE.rel , nrow(res) )
      }
    }				
    #  res <- list( "PersonScores" = PersonScores, "PersonMax" = PersonMax, "theta" = theta , "error" =  error )
					
	
	attr(res,"ndim") <- ndim
	attr(res,"nobs") <- nrow(res)
	#*** collect reliabilities
	i1 <- grep( "WLE.rel" , colnames(res), fixed = TRUE )
    if (ndim==1){
		attr(res,"WLE.rel") <- res[[i1]][1]
				} else {
		v1 <- as.numeric(res[1,i1])
		names(v1) <- colnames(res)[i1]
		attr(res,"WLE.rel") <- v1	
				}	
	#	attr(res,"WLE.rel") <- res[1,i1]
	attr(res,"call") <- CALL
	class(res) <- c("tam.wle","data.frame")		
	if (output.prob){
		 res <- as.list(res)
	     res$probs <- rprobsWLE
					}	
	
	
    return(res)
  }

################################################################
################################################################
################################################################
tam.mml.wle2 <-
  function( tamobj, score.resp=NULL , WLE=TRUE , adj=.3 , Msteps=20 , 
            convM = .0001 , progress=TRUE , output.prob=FALSE ){
    #########################################################
    # INPUT:
    # tamobj ... result from tam analysis
    # (WLE = TRUE) will produce WLE. Otherwise it will be MLE
    # 
    #########################################################
    #  adj <- 0.3
    #  Msteps <- 20
    #  convM <- .0001
	CALL <- match.call()
	
    B <- tamobj$B
    A <- tamobj$A
    nitems <- tamobj$nitems
    xsi <- ( tamobj$xsi )[,1]
    nstud <- tamobj$nstud
    AXsi <- tamobj$AXsi
    ndim <- tamobj$ndim
    maxK <- tamobj$maxK
    resp <- tamobj$resp
    if ( ! is.null( score.resp) ){
      resp <- score.resp
      nstud <- nrow(resp)
      tamobj$resp.ind <- 1 - is.na(resp)
      tamobj$pid <- 1:nstud
	  tamobj$pweights <- 1+0*tamobj$pid
    }  
    resp[ is.na(resp) ] <- 0  
    resp.ind <- tamobj$resp.ind  
    col.index <- rep( 1:nitems , each = maxK )
    cResp <- resp[ , col.index  ]*resp.ind[ , col.index ]
    cResp <- 1 * t( t(cResp) == rep(0:(maxK-1), nitems) )
    cB <- t( matrix( aperm( B , c(2,1,3) ) , nrow = dim(B)[3] , byrow = TRUE ) )
    cB[is.na(cB)] <- 0
    #Compute person sufficient statistics (total score on each dimension)
    PersonScores <- cResp %*% cB
    
    #Compute possible maximum score for each item on each dimension
    maxBi <- apply(B , 3 , rowMaxs , na.rm = TRUE)
    
    #Compute possible maximum score for each person on each dimension
    PersonMax <- resp.ind %*% maxBi
    PersonMax[ PersonMax == 0 ] <- 2 * adj
    
    #Adjust perfect scores for each person on each dimension
    PersonScores[PersonScores==PersonMax] <- PersonScores[PersonScores==PersonMax] - adj
    
    #Adjust zero scores for each person on each dimension
    PersonScores[PersonScores==0] <- PersonScores[PersonScores==0] + adj
    
    #Calculate Axsi. Only need to do this once.
    for (i in 1:nitems) {
      for (k in 1:maxK){
        AXsi[i,k] <- ( A[i,k,] %*% xsi )
      }
    }
    
    #Initialise theta (WLE) values for all students
    theta <- log(PersonScores/(PersonMax-PersonScores)) #log of odds ratio of raw score
    
    ######################################
    #Compute WLE
    #similar to the M step in the tam function, but each student's theta vector is now one node.
    converge <- FALSE
    Miter <- 0
    BB <- array (0, dim=c(nitems,maxK,ndim,ndim))
    BBB <- array (0, dim=c(nitems,maxK,ndim)) 
    for (i in 1:nitems) {
      for (k in 1:maxK) {
        BB[i,k,,] <- B[i,k,] %*% t(B[i,k,])
        BBB[i,k,] <- BB[i,k,,] %*% B[i,k,]
      }
    }
    BL <- matrix(B, nitems*maxK, ndim)
    BBL <- matrix(BB, nitems*maxK, ndim*ndim)
    BBBL <- matrix(BBB, nitems*maxK, ndim)
    
    increment <- array(0, dim=c(nstud,ndim))
    old_increment <- 3 + increment
    
    
    # Begin iterations
    while (!converge & ( Miter <= Msteps ) ) {
      
      resWLE <- calc_prob.v5(iIndex = 1:nitems , A , AXsi , 
                             B , xsi , theta , nstud, maxK , recalc=FALSE )      	
      rprobsWLE <- resWLE[["rprobs"]] 
      rprobsWLEL <- matrix(rprobsWLE, nitems*maxK, nstud )
      
      rprobsWLEL[is.na(rprobsWLEL)] <- 0
      resB <- .Call("tam_wle_Bs", rprobsWLEL, resp.ind, BL, BBL, BBBL, 
                        ndim, nitems, maxK, nstud, PACKAGE="TAM")
      B_bari <- array(resB$B_bari, dim=c(nstud, nitems,ndim))
      BB_bari <- array(resB$BB_bari, dim=c(nstud, nitems, ndim, ndim))
      BBB_bari <- array(resB$BBB_bari, dim=c(nstud, nitems, ndim))
      
      B_Sq <- array(resB$B_Sq,dim=c(nstud, nitems, ndim, ndim))
      B2_B <- array(resB$B2_B,dim=c(nstud, nitems, ndim))
      B_Cube <- array(resB$B_Cube,dim=c(nstud, nitems, ndim))
       expected <- colSums(aperm(B_bari,c(2,1,3)))
      err <- colSums(aperm(BB_bari,c(2,1,3,4))) - colSums(aperm(B_Sq, c(2,1,3,4)))  #sum over the items
   
      if (ndim == 1) {
        # err_inv <- apply(err,1,function(x) 1/x )
        err_inv <- 1 / err
      } else {
        ## diag err_i= forall i
        diag_ind <- cbind(rep(1:nstud, each=ndim), 1:ndim, 1:ndim)
        err[diag_ind] <- err[diag_ind]+10^-15    
        errl <- matrix(err, nstud, ndim*ndim)
      
        err_inv <- .Call( "tam_wle_errinv" , errl, ndim, nstud , PACKAGE="TAM" )
      }
        
      err_inv <- array(abs(err_inv),dim=c(nstud,ndim,ndim))
      warm <- -3*B2_B + 2*B_Cube + BBB_bari
      warmadd <- colSums(aperm(warm,c(2,1,3)))  #sum over the items
      scores <- PersonScores - expected
      if (WLE) {
        warmaddon <- array(0,dim=c(nstud,ndim))
        for (d1 in 1:ndim) {
          warmaddon[,d1] <- 0
          for (d2 in 1:ndim) {
            warmaddon[,d1] <- warmaddon[,d1] + err_inv[,d1,d2]*warmadd[,d2]
          }
        }
        scores <- scores + warmaddon/2.0      
      }
      
      increment <- array(0, dim=c(nstud,ndim))
      for (d1 in 1:ndim) {
        increment[,d1] <- 0
        for (d2 in 1:ndim) {
          increment[,d1] <- increment[,d1] + err_inv[,d1,d2]*scores[,d2]
        }
      }
      
      # dampening the increment
      for ( d1 in 1:ndim){ 
        #	   increment[,d1] <- ifelse( abs(increment[,d1]) > 3 , sign( increment[,d1] )*3 , increment[,d1] )
        ci <- ceiling( abs(increment[,d1]) / ( abs( old_increment[,d1]) + 10^(-10) ) )	   
        increment[,d1] <- ifelse( abs( increment[,d1]) > abs(old_increment[,d1])  , 
                                  increment[,d1]/(2*ci) , 
                                  increment[,d1] )	   
#        old_increment[,d1] <- increment[,d1] 
		old_increment[,d1] <- .95 * old_increment[,d1]
        #***
        # avoid NaNs in increment
        increment[,d1] <- ifelse( is.na(increment[,d1] ) , 0 , increment[,d1] )
        # increment[abs(increment)>3] <- sign(increment[abs(increment)>3])*3	
      }
      
      theta <- theta + increment
      if ( max(abs(increment)) < convM ) {
        converge <- TRUE
      }
      
      Miter <- Miter + 1 

      if (progress){
		  cat( paste( "Iteration in WLE/MLE estimation ", Miter, 
					  "  | Maximal change " , round( max(abs(increment)) , 4) , "\n" )  ) 
		  utils::flush.console()
					}
    }  # end of Newton-Raphson

    #standard errors of theta estimates
    if (ndim == 1) {
      error <- apply(err_inv,1,sqrt) 
    } else {    
      error <- aperm(apply(sqrt(err_inv),1,diag), c(2,1))
    }
    
    # The output contains 
    #   Person Scores on the test, by dimension
    #   Person possible maximum score, by dimension (Each person could take 
    #    different items, so possible maximum could vary)
    #   WLE or MLE estimate, by dimension
    #   Standard errors of WLE/MLE estimates, by dimension

    
    if ( ndim> 1){
      colnames(error) <- paste0("error.Dim" , substring( 100+1:ndim , 2) )
    }
    res <- data.frame( "pid" = tamobj$pid , 
                       "N.items" = rowSums(resp.ind) , 
                       "PersonScores" = PersonScores, 
                       "PersonMax" = PersonMax, "theta" = theta , error )
    
    if (ndim==1){ colnames(res)[4:5] <- c("PersonMax" , "theta") }
    if (ndim>1){  
      colnames(res)[ 1:ndim + 2] <- paste0("PersonScores.Dim" , substring( 100+1:ndim , 2) )	
      ind <- grep( "theta" , colnames(res) )	
      colnames(res)[ind] <- 	paste0("theta.Dim" , substring( 100+1:ndim , 2) )	
    }
    ####################
    # correct personMax set theta and standard error to missing		
    # if there are no observations on one dimension
    ind1 <- grep("PersonMax" , colnames(res))
    check1 <- ( res[ , ind1 , drop=FALSE] == 2*adj )
    ind2 <- grep("theta" , colnames(res))
    D <- length(ind1)
    for (ii in 1:D){
      res[ check1[,ii] , ind2[ii] ] <- NA
    }
    ind2 <- grep("error" , colnames(res))
    for (ii in 1:D){
      res[ check1[,ii] , ind2[ii] ] <- NA
    }
    #***
    # WLE reliability
	pweights <- tamobj$pweights
    if ( ndim==1 ){
      ind <- which( res$N.items > 0 )  
	  WLE.rel <- WLErel(theta=theta[ind] , error = error[ind] , w = pweights[ind] )	  	  	  	 
	  if (WLE){ w1 <- "WLE" } else { w1 <- "MLE" }
	  if (progress){
		cat("----\n" , w1 ,"Reliability =" , round(WLE.rel,3) ,"\n" )	  
					}
#      cat("----\nWLE Reliability =" , round(WLE.rel,3) ,"\n" )
      res$WLE.rel <- rep( WLE.rel , nrow(res) )
    }
    if ( ndim>1 ){
      cat("\n-------\n")
      for (dd in 1:ndim){
        #	dd <- 1
		ind1 <- paste0("theta.Dim" , substring( 100+1:ndim , 2))[dd]
		ind2 <- paste0("error.Dim" , substring( 100+1:ndim , 2))[dd] 		
        h1 <- WLErel( theta=res[,ind1] , error=res[,ind2] , w = pweights )
        res[ ,paste0("WLE.rel.Dim" , substring( 100+ dd , 2)) ]	<- h1  #  <- 1 - v2 / v1
		if (WLE){ w1 <- "WLE" } else { w1 <- "MLE" }
        cat(paste0(w1 , " Reliability (Dimension" , dd , ") = " , round(h1,3) ) , "\n" )
        #	  res$WLE.rel <- rep( WLE.rel , nrow(res) )
      }
    }	

    
	#*******************************
	# check identifiability
	dimB <- dim(B)
    D <- dimB[3]
	I <- dimB[1]
	K <- dimB[2]
    identM <- matrix( 0 , nrow=I , ncol=D )
	for ( kk in 1:K){
		identM <- identM + 1 * ( B[,kk,] != 0 )
					}	
	betweenload <- 0*identM
	sumloads <- rowSums( identM > 0 )
	ind <- which( sumloads == 1 )
	if (length(ind) > 0 ){
		betweenload[ ind , ] <- identM[ ind , ]
						}
	loaddim <- colSums(betweenload)
	if ( min(loaddim) == 0 ){
		cat("\n * Not all dimensions do have items with simple loading structure.\n")
		cat(" * Maybe the WLE is not identified (i.e. estimable).\n")
		cat(" * Please proceed with caution.\n")		
						}
					
	res <- as.data.frame(res)				
	attr(res,"ndim") <- ndim
	attr(res,"nobs") <- nrow(res)	
	#*** collect reliabilities
	i1 <- grep( "WLE.rel" , colnames(res), fixed = TRUE )
    if (ndim==1){
		attr(res,"WLE.rel") <- res[[i1]][1]
				} else {
		v1 <- as.numeric(res[1,i1])
		names(v1) <- colnames(res)[i1]
		attr(res,"WLE.rel") <- v1	
				}
	attr(res,"call") <- CALL
	class(res) <- c("tam.wle","data.frame")
	if (output.prob){
		 res <- as.list(res)
	     res$probs <- rprobsWLE
					}	
    #  res <- list( "PersonScores" = PersonScores, "PersonMax" = PersonMax, "theta" = theta , "error" =  error )
    return(res)
  }

