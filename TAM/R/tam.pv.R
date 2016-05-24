tam.pv <- function( tamobj , nplausible = 10 , 
			ntheta = 2000 , 
			normal.approx = FALSE , 
            samp.regr = FALSE , theta.model = FALSE , np.adj = 8 ){
    #####################################################
    # INPUT:
    # tamobj ... result from tam analysis
    # nplausible ... number of plausible values
    # ntheta ... number of simulated theta values
    # samp.regr ... sample regression coefficients?
	#        when sampling regression coefficients,
	#        plausible values are used for recalculating
	#        regression coefficients
	#        (sampling of regression coefficients only
	#          works in the unidimensional case)
	# normal.approx ... use normal distribution as an 
	#					approximation of the posterior
    ####################################################

a0 <- Sys.time()	
    type <- "nonparm"		# there is no type='normal' up to now implemented
	latreg <- FALSE
	if ( class(tamobj) == "tam.latreg" ){
		theta.model <- TRUE
		latreg <- TRUE
		like <- tamobj$like
				}	
	if ( ! latreg ){
		if (class(tamobj)!= "tam.mml.3pl"){
			guess <- rep( 0 , dim(tamobj$B)[1] )		
			} else { guess <- tamobj$guess }				
		B <- tamobj$B
		A <- tamobj$A
		AXsi <- tamobj$AXsi		
		xsi <- ( tamobj$xsi )[,1]		
		maxK <- tamobj$maxK		
				}
    Y <- tamobj$Y
	YSD <- tamobj$YSD
    nitems <- tamobj$nitems
	snodes <- tamobj$control$snodes 
	
    beta <- tamobj$beta
    variance <- tamobj$variance
    nstud <- tamobj$nstud

		
	if ( theta.model ){
        ntheta <- nrow(tamobj$theta)			
			}	
	
	nthetal <- rep( 1 , ntheta )
	nnodes <- ntheta    
	ndim <- tamobj$ndim
	pweights <- tamobj$pweights

	#***************************
    # define theta grid
	#--- dim = 1
	if ( ndim == 1 ){
		MEAP <- mean( tamobj$person$EAP )
		SDEAP <- sqrt( stats::var( tamobj$person$EAP ) + mean( tamobj$person$SD.EAP^2 ) )
				}
	#--- dim > 1
	if ( ndim > 1 ){
		tp1 <- tamobj$person
		ind <- grep("EAP\\.Dim" , colnames(tp1) )
		ind <- ind[ seq( 1 , length(ind) , 2 ) ]
		dat1 <- tp1[ ,  ind ]
		mu1 <- as.vector( colMeans( dat1 ) )
		var1 <- apply( dat1 , 2 , stats::var ) / tamobj$EAP.rel
		Sigma1 <- stats::cov2cor(variance)
		Sigma1 <- np.adj * diag( sqrt( var1) ) %*% Sigma1 %*% diag( sqrt( var1 ))
					}
										
    # create pv matrix (uni- and multidimensional case)
    pv <- matrix( 0 , nrow=nstud , ncol= nplausible*ndim)     
    NPV <- nplausible
    pp <- 1
	cat("|")
	cat( paste( rep("*" , nplausible ) , collapse="") )
	cat("|\n|") ; utils::flush.console()
	

	###################################################
	# routine for drawing plausible values
	while ( pp <= NPV ){
	
	#*****************
	#***** sampling of theta values. These values can also be left fixed.
	#*****************
	if ( ! theta.model ){
	   #***************************
       # 1-dimensional PV imputation	   
	   if (ndim == 1){
			# unidimensional theta simulation
			if ( ! normal.approx){			
				theta <- matrix( stats::rnorm( ntheta , mean = MEAP , sd = np.adj*SDEAP )  , ncol= 1)	
				theta <- theta[ order( theta[,1] ) , , drop=FALSE]
								} else {
				theta <- matrix( SDEAP * seq( - 5 , 5 , len=ntheta ) + MEAP , ncol=1 )		 
							}
						}
	   #*****************************
       # multidimensional PV imputation									
	   if( ndim > 1 ){
	      theta <- MASS::mvrnorm( ntheta , mu = mu1 , Sigma = Sigma1 )
					}
				}
	 if ( theta.model ){
	      theta <- tamobj$theta 
						}

				
# cat("start prob") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1	
 
      if ( ! latreg ){				
 	    res <- .mml.3pl.calc_prob.v5( iIndex=1:nitems , A=A , AXsi=AXsi , B=B , xsi=xsi , theta=theta , 
 	                         nnodes=nnodes, maxK=maxK , recalc=TRUE , guess=guess)
		rprobs <- res[["rprobs"]]
		AXsi <- res[["AXsi"]]
					}
# cat("calc prob") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1			

		# calculate student's prior distribution    	
		gwt <- stud_prior.v2( theta=theta , Y=Y , beta=beta , variance=variance , nstud=nstud , 
                          nnodes=nnodes , ndim=ndim , YSD=YSD , unidim_simplify=FALSE,
						  snodes = snodes )
# cat("stud prior") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1							  
   
   
        #**** 
		# posterior distribution
		if ( ! latreg ){		
			hwt <- calc_posterior.v2( rprobs=rprobs , gwt=gwt , resp=tamobj$resp , nitems=nitems , 
		                          resp.ind.list=tamobj$resp.ind.list , normalization=TRUE , 
		                          thetasamp.density=NULL , snodes=0 )$hwt
						}
		if (latreg){
   		    hwt <- like * gwt
			hwt <- hwt / rowSums(hwt)	 		
					}
					
   	   hwt1 <- hwt			
# cat("posterior") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1							  		
		hwt1 <- rowCumsums.TAM(hwt1) # include this function in later versions!!
# cat("rowcumsums TAM") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1							  		

		#---------------------------------------------------------
		# sampling of regression coefficients
		if ( samp.regr ){
				
			#*****
			# no normal approximation
			if ( ! normal.approx){	
				rn1 <- stats::runif( nstud )
				ind <- rowSums( hwt1 < outer( rn1 , nthetal ) ) +1				
				if (ndim==1){
					pv[,pp] <- theta1 <- theta[ind]	
							} else {
					theta1 <- pv[ , (pp-1)*(ndim) + 1:ndim ] <- theta[ ind , ]
									}
								}
			#*****
			# normal approximation in unidimensional case
			if ( normal.approx & ( ndim == 1 ) ){
				thetaM <- matrix( theta[,1] , nstud , ntheta , byrow=TRUE )
				EAP <- rowSums( thetaM * hwt )
				SDPost <- sqrt( rowSums( thetaM^2 * hwt ) - EAP^2 )
				pv[,pp] <- theta1 <- stats::rnorm( nstud , mean = EAP , sd = SDPost )		
					           }

			 #------
			 # normal approximation (ndim > 1)
			 if (  normal.approx  & ( ndim > 1 ) ){		
			    N <- nrow(hwt)			
				MEAP <- matrix( 0 , nrow=N , ncol=ndim)	
				SDEAP <- matrix( 0 , nrow=N , ncol=ndim)
				nstudl <- rep(1,N)	
				for ( dd in 1:ndim ){
					MEAP[,dd] <- rowSums( hwt * outer( nstudl , theta[,dd] ) )
					SDEAP[,dd] <- sqrt(rowSums( hwt * outer( nstudl , theta[,dd]^2 ) ) -
											MEAP[,dd]^2)					
								}												 				
				thetaPV <- matrix( stats::rnorm( N * ndim ) , nrow=N , ncol=ndim )
				thetaPV <- SDEAP * thetaPV + MEAP
				theta1 <- pv[ , (pp-1)*(ndim) + 1:ndim ] <- thetaPV			 
									}
									
			pp <- pp + 1
						
			# bootstrap sample of persons to get sampled beta coefficients
			if ( ndim > 1 ){
				N <- nrow(theta1)			
				ind <- base::sample( 1:N , N , replace=TRUE)
				theta1 <- theta1[ ind , ]				
				Y1 <- Y[ ind , , drop=FALSE ]
						} else {
				N <- length(theta1)
				ind <- base::sample( 1:N , N , replace=TRUE)
				theta1 <- theta1[ ind ]				
				Y1 <- Y[ ind , , drop=FALSE ]				
								}
			
			modlm <- stats::lm( theta1  ~  -1 + as.matrix(Y1) , weights = pweights)
			beta <- modlm$coef			# sampled regression coefficients
			if ( ndim == 1 ){
				beta <- matrix( beta , ncol=1 )
							}
								}


		#---------------------------------------------------------
		# no sampling of regression cofficients
   		if ( ! samp.regr ){
			for ( pp in 1:nplausible ){
			 #------
			 # no normal approximation
			 if (  ! normal.approx  ){
				rn1 <- stats::runif( nstud )				
				ind <- interval_index( hwt1 , rn1 )
				ind <- ind - 1				
                pv[ , (pp-1)*(ndim) + 1:ndim ] <- theta[ind , ]			
									}

			 #------
			 # normal approximation (ndim > 1)
			 if (  normal.approx  & ( ndim > 1 ) ){		
			    N <- nrow(hwt)
				MEAP <- matrix( 0 , nrow=N , ncol=ndim)	
				SDEAP <- matrix( 0 , nrow=N , ncol=ndim)
				nstudl <- rep(1,N)	
				for ( dd in 1:ndim ){
					MEAP[,dd] <- rowSums( hwt * outer( nstudl , theta[,dd] ) )
					SDEAP[,dd] <- sqrt(rowSums( hwt * outer( nstudl , theta[,dd]^2 ) ) - MEAP[,dd]^2)					
								}												 				
				thetaPV <- matrix( stats::rnorm( N * ndim ) , nrow=N , ncol=ndim )
				thetaPV <- SDEAP * thetaPV + MEAP
				pv[ , (pp-1)*(ndim) + 1:ndim ] <- thetaPV			 
									}
									
									
			#-------						 
			# dim =1 and normal approximation
			if ( normal.approx & ( ndim == 1) ){
				thetaM <- matrix( theta[,1] , nstud , ntheta , byrow=TRUE )
				EAP <- rowSums( thetaM * hwt )
				SDPost <- sqrt( rowSums( thetaM^2 * hwt ) - EAP^2 )
				pv[,pp] <- stats::rnorm( nstud , mean = EAP , sd = SDPost )					  
									}
																		
            if (pp != nplausible){ 
					cat("-") ; utils::flush.console() 
									}
								}
			NPV <- nplausible / 2
						}   # end no plausible
		#--------------------------		
		cat("-" ) ; utils::flush.console()
					}  # end while
	##################################################	
			cat("|\n")
# cat("rest") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1							  			
	# label the pv matrix	
	colnames(pv) <- paste("PV" , rep(1:nplausible,each=ndim) , 
					".Dim" , rep(1:ndim,nplausible) , sep="")   
    pv <- data.frame( "pid" =tamobj$pid , pv )					
    res <- list( "pv" = pv , "hwt" = hwt , "hwt1" = hwt1 ,
                "theta" = theta , "ndim" = ndim , "nplausible" = nplausible ,
				"pid" = tamobj$pid , "pweights" = tamobj$pweights )
	class(res) <- "tam.pv"
    return(res)
    }
##################################################################
##################################################################