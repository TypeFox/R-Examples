
###################################################################
# latent regression
tam.latreg <- function( like , theta=NULL , Y=NULL , group=NULL , 
				formulaY = NULL , dataY = NULL , 
				beta.fixed = NULL , beta.inits = NULL , 
				variance.fixed = NULL , variance.inits = NULL , 
				est.variance = TRUE , pweights = NULL , pid=NULL , 
				userfct.variance = NULL , variance.Npars = NULL , 
				control = list() 
  ){
       
    s1 <- Sys.time()
	CALL <- match.call()
    # display
    disp <- "....................................................\n"    
    increment.factor <- progress <- nodes <- snodes <- ridge <- xsi.start0 <- QMC <- NULL
    maxiter <- conv <- convD <- min.variance <- max.increment <- Msteps <- convM <- NULL 
    pweightsM <- R <- NULL
    
    # attach control elements
    e1 <- environment()
    con <- list( convD = .001 ,conv = .0001 , snodes=0 , convM = .0001 , Msteps = 4 ,            
                 maxiter = 1000 ,min.variance = .001 , progress = TRUE , ridge=0 ,
                 seed = NULL )  	
    #a0 <- Sys.time()			   
    con[ names(control) ] <- control  
    Lcon <- length(con)
    con1a <- con1 <- con ; 
    names(con1) <- NULL
    for (cc in 1:Lcon ){
      assign( names(con)[cc] , con1[[cc]] , envir = e1 ) 
    }
    
	if ( is.null(theta) ){
	   theta <- attr( like , "theta" )
						}

						
    nodes <- theta 
	ndim <- ncol(theta)		 

	
    if (progress){ 
      cat(disp)	
      cat("Processing Data     ", paste(Sys.time()) , "\n") ; 
	  utils::flush.console()
    }  
    
    if ( ! is.null(group) ){ 
      con1a$QMC <- QMC <- FALSE
      con1a$snodes <- snodes <- 0
    }

    
    
    if ( !is.null(con$seed)){ set.seed( con$seed )	 }
    nullY <- is.null(Y)
    
    nstud <- nrow(like)        # number of students
    if ( is.null( pweights) ){
      pweights <- rep(1,nstud) # weights of response pattern
    }
    if (progress){ 
      cat("    * Response Data:" , nstud , "Persons \n" )  ;
      utils::flush.console()	  
    }  	  
    
    #!! check dim of person ID pid
    if ( is.null(pid) ){ 
			pid <- seq(1,nstud) 
				} else { 
			pid <- unname(c(unlist(pid))) 
						}
       
    # normalize person weights to sum up to nstud
    pweights <- nstud * pweights / sum(pweights)
        
    betaConv <- FALSE         #flag of regression coefficient convergence
    varConv <- FALSE          #flag of variance convergence
    # nnodes <- length(nodes)^ndim
	nnodes <- nrow(nodes)

	
	if ( snodes > 0 ){ nnodes <- snodes }
    
    #****
    # display number of nodes
    if (progress ){   
      l1 <- paste0( "    * ")
      if (snodes==0){ l1 <- paste0(l1 , "Numerical integration with ")}
      else{ 
        if (QMC){ 
          l1 <- paste0(l1 , "Quasi Monte Carlo integration with ")
        } else {
          l1 <- paste0(l1 , "Monte Carlo integration with ")					
        }
      }
      cat( paste0( l1 , nnodes , " nodes\n") )
      if (nnodes > 8000){
        cat("      @ Are you sure that you want so many nodes?\n")
        cat("      @ Maybe you want to use Quasi Monte Carlo integration with fewer nodes.\n")		
      }
    }
	
    #*********
    # variance inits  
    # initialise conditional variance 
    if ( !is.null( variance.inits ) ){
      variance <- variance.inits
    } else variance <- diag( ndim ) 
    if ( !is.null(variance.fixed) ){
      variance[ variance.fixed[,1:2 ,drop=FALSE] ] <- variance.fixed[,3]
      variance[ variance.fixed[,c(2,1) ,drop=FALSE] ] <- variance.fixed[,3]	
    }
    # group indicators for variance matrix
    if ( ! is.null(group) ){ 
      groups <- sort(unique(group))
      G <- length(groups)
      group <- match( group , groups )
      var.indices <- rep(1,G)
      for (gg in 1:G){
        var.indices[gg] <- which( group == gg )[1]				
      }
    } else { 
      G <- 1 
      groups <- NULL
    }  
    # beta inits
    # (Y'Y)
    if ( ! is.null( formulaY ) ){
      formulaY <- stats::as.formula( formulaY )
      Y <- stats::model.matrix( formulaY , dataY )[,-1]   # remove intercept
      nullY <- FALSE
    }
    
    
    #  if ( ! is.null(Y) ){ 
    if (! nullY){
      Y <- as.matrix(Y)
      nreg <- ncol(Y)
      if ( is.null( colnames(Y) ) ){
        colnames(Y) <- paste("Y" , 1:nreg , sep="")
      }
      if ( ! nullY ){ 		
        Y <- cbind(1,Y)          #add a "1" column for the Intercept
        colnames(Y)[1] <- "Intercept"
      }
    } else 
    {
      Y <- matrix( 1 , nrow=nstud , ncol=1 ) 
      nreg <- 0
    }
    if ( G > 1 & nullY ){	
      Y <- matrix( 0 , nstud , G )
      #		colnames(Y) <- paste("group" , 1:G , sep="")
      colnames(Y) <- paste("group" , groups , sep="")
      for (gg in 1:G){ Y[,gg] <- 1*(group==gg) }
      nreg <- G - 1
    }
    
    # W <- t(Y * pweights) %*% Y
	W <- crossprod(Y * pweights ,  Y )
    if (ridge > 0){ diag(W) <- diag(W) + ridge }
    YYinv <- solve( W )
    
    #initialise regressors
#    if ( is.null(beta.fixed)  ){
#      beta.fixed <- matrix( c(1,1,0) , nrow= 1) 
#      if (  ndim > 1){ 
#        for ( dd in 2:ndim){
#          beta.fixed <- rbind( beta.fixed , c( 1 , dd , 0 ) )
#        }}}
    
    #****
    if( ! is.matrix(beta.fixed) ){
      if ( ! is.null(beta.fixed) ){
        if ( ! beta.fixed   ){ beta.fixed <- NULL }
      }
    }
    #****
    
    beta <- matrix(0, nrow = nreg+1 , ncol = ndim)  
    if ( ! is.null( beta.inits ) ){ 
      beta[ beta.inits[,1:2] ] <- beta.inits[,3]
    }	
	   
    beta.min.deviance <- beta
    variance.min.deviance <- variance
    
    # cat("b200"); a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1  									  
    
    # nodes
    if ( snodes == 0 ){ 
#      theta <- as.matrix( expand.grid( as.data.frame( matrix( rep(nodes, ndim) , ncol = ndim ) ) ) )
      #we need this to compute sumsig2 for the variance
      #    theta2 <- matrix(theta.sq(theta), nrow=nrow(theta),ncol=ncol(theta)^2)            
      theta2 <- matrix(theta.sq2(theta), nrow=nrow(theta),ncol=ncol(theta)^2)            
      # grid width for calculating the deviance
      thetawidth <- diff(theta[,1] )
      thetawidth <- ( ( thetawidth[ thetawidth > 0 ])[1] )^ndim 
      thetasamp.density <- NULL
    } else {
      # sampled theta values
      if (QMC){						
        r1 <- sfsmisc::QUnif(n=snodes, min = 0, max = 1, n.min = 1, p=ndim, leap = 409)						
        theta0.samp <- stats::qnorm( r1 )
      } else {
        theta0.samp <- matrix( MASS::mvrnorm( snodes , mu = rep(0,ndim) , 
                                        Sigma = diag(1,ndim ) )	,
                               nrow= snodes , ncol=ndim )			
      }
      thetawidth <- NULL
    }
    
    
    deviance <- 0  
    deviance.history <- matrix( 0 , nrow=maxiter , ncol = 2)
    colnames(deviance.history) <- c("iter" , "deviance")
    deviance.history[,1] <- 1:maxiter
    
    iter <- 0 
    a02 <- a1 <- 999	# item parameter change
    a4 <- 0
    
    YSD <- max( apply( Y , 2 , stats::sd ) )
    if (YSD > 10^(-15) ){ YSD <- TRUE } else { YSD <- FALSE }
    
    # define progress bar for M step
#    mpr <- round( seq( 1 , np , len = 10 ) )
    
    hwt.min <- 0
    deviance.min <- 1E100
    itemwt.min <- 0

	nomiss <- TRUE
	Variance.fixed <- variance.fixed
	res.hwt <- list()
    
    ##############################################################   
    #Start EM loop here
    while ( ( (!betaConv | !varConv)  | ((a1 > conv) | (a4 > conv) | (a02 > convD)) )  & (iter < maxiter) ) { 
      
      iter <- iter + 1
      if (progress){ 
        cat(disp)	
        cat("Iteration" , iter , "   " , paste( Sys.time() ) )
        cat("\nE Step\n") ; 
		utils::flush.console()
      }
      # calculate nodes for Monte Carlo integration	
      if ( snodes > 0){
        #      theta <- beta[ rep(1,snodes) , ] +  t ( t(chol(variance)) %*% t(theta0.samp) )
        theta <- beta[ rep(1,snodes) , ] + theta0.samp %*% chol(variance) 
        # calculate density for all nodes
        thetasamp.density <- mvtnorm::dmvnorm( theta , mean = as.vector(beta[1,]) , sigma = variance )
        # recalculate theta^2
        #      theta2 <- matrix( theta.sq(theta) , nrow=nrow(theta) , ncol=ncol(theta)^2 )   
        theta2 <- matrix( theta.sq2(theta) , nrow=nrow(theta) , ncol=ncol(theta)^2 )   
      }			
      olddeviance <- deviance
# a0 <- Sys.time()	
      
      #***	
      # print(AXsi)	
      # AXsi[ is.na(AXsi) ] <- 0
      # print(AXsi)	
      
      # cat("calc_prob") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1
      
      # calculate student's prior distribution
      gwt <- stud_prior.v2(theta=theta , Y=Y , beta=beta , variance=variance , nstud=nstud , 
                           nnodes=nnodes , ndim=ndim,YSD=YSD, unidim_simplify=FALSE)
	  # compute posterior	  
	  hwt <- like * gwt
	  res.hwt$rfx <- rowSums(hwt)
	  hwt <- hwt / rowSums(hwt)	 
      # cat("calc_posterior") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1						 
      
     
      if (progress){ 
	       cat("M Step Intercepts   |")
		   utils::flush.console() 
			   }
      oldbeta <- beta
      oldvariance <- variance 
  
#	  cat("before mstep regression") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1						           
      # M step: estimation of beta and variance
	  resr <- latreg.mstep.regression( hwt , 
			pweights , pweightsM , Y , theta , theta2 , YYinv , ndim , 
			nstud , beta.fixed , variance , Variance.fixed , group , G , 
			snodes = snodes , thetasamp.density=thetasamp.density , nomiss=FALSE)
  
      # cat("mstep regression") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1						       
      beta <- resr$beta
      
      variance <- resr$variance	
      if( ndim == 1 ){  # prevent negative variance
        variance[ variance < min.variance ] <- min.variance 
      }
      itemwt <- resr$itemwt
      
      # constraint cases (the design matrix A has no constraint on items)
      if ( max(abs(beta-oldbeta)) < conv){    
        betaConv <- TRUE       # not include the constant as it is constrained
      }

      
      if (G == 1){
        diag(variance) <- diag(variance) + 10^(-10)
			}
      
	    
	  # function for reducing the variance	  
	  if ( ! is.null( userfct.variance ) ){  
			variance <- do.call( userfct.variance , list(variance ) )			
				}	        
      if (max(abs(variance-oldvariance)) < conv) varConv <- TRUE
      

      
      # calculate deviance
      if ( snodes == 0 ){ 
        deviance <- - 2 * sum( pweights * log( res.hwt$rfx * thetawidth ) )
#        deviance <- - 2 * sum( pweights * log( res.hwt$like * thetawidth ) )
      } else {
        #       deviance <- - 2 * sum( pweights * log( res.hwt$rfx ) )
        deviance <- - 2 * sum( pweights * log( rowMeans( res.hwt$swt ) ) )	
      }
      deviance.history[iter,2] <- deviance
      a01 <- abs( ( deviance - olddeviance ) / deviance  )
      a02 <- abs( ( deviance - olddeviance )  )	
      
      if( deviance > deviance.min ){ 	 
        beta.min.deviance <- beta.min.deviance
        variance.min.deviance <- variance.min.deviance
        hwt.min <- hwt.min
        deviance.min <- deviance.min
      }   else { 
        beta.min.deviance <- beta
        variance.min.deviance <- variance	
        hwt.min <- hwt	
        deviance.min <- deviance
      }
      
      a1 <- 0
      a2 <- max( abs( beta - oldbeta ))	
      a3 <- max( abs( variance - oldvariance ))
      if (progress){ 
        cat( paste( "\n  Deviance =" , round( deviance , 4 ) ))
        devch <- -( deviance - olddeviance )
        cat( " | Deviance change:", round( devch  , 4 ) )
        if ( devch < 0 & iter > 1 ){ 
          cat ("\n!!! Deviance increases!                                        !!!!") 
          cat ("\n!!! Choose maybe fac.oldxsi > 0 and/or increment.factor > 1    !!!!") 			
        }
        
        
        cat( "\n  Maximum regression parameter change:" , round( a2 , 6 ) )  
        if ( G == 1 ){ 
          cat( "\n  Variance: " , round( variance[ ! lower.tri(variance)] , 4 ) , " | Maximum change:" , round( a3 , 6 ) )  
        } else {
          cat( "\n  Variance: " , round( variance[var.indices] , 4 ) ,
               " | Maximum change:" , round( a3 , 6 ) )  		
        }					
        cat( "\n  beta ",round(beta,4)  )
        cat( "\n" )
        utils::flush.console()
      }
      
      # cat("rest") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1	
      
      
    } # end of EM loop
    #******************************************************

    beta.min.deviance -> beta
    variance.min.deviance -> variance	
    hwt.min -> hwt	
    deviance.min -> deviance

    ##*** Information criteria
	ic <- latreg_TAM.ic( nstud , deviance , 
				beta , beta.fixed , ndim , variance.fixed , G , 
				est.variance , variance.Npars=NULL , group )

			   
    #####################################################
    # post ... posterior distribution	
    # create a data frame person	
    person <- data.frame( "pid"=pid , "case" = 1:nstud , "pweight" = pweights )
    # person$score <- rowSums( resp * resp.ind )		
    # use maxKi here; from "design object"
    nstudl <- rep(1,nstud)
#    person$max <- rowSums( outer( nstudl , apply( resp ,2 , max , na.rm=TRUE) ) * resp.ind )
    # calculate EAP
    # EAPs are only computed in the unidimensional case for now,
    # but can be easily adapted to the multidimensional case
    if ( snodes == 0 ){ 
      hwtE <- hwt 
    } else { 	
      # hwtE <- hwt / snodes 
      hwtE <- hwt
    }
    if ( ndim == 1 ){
      person$EAP <- rowSums( hwtE * outer( nstudl , theta[,1] ) )
      person$SD.EAP <- sqrt( rowSums( hwtE * outer( nstudl , theta[,1]^2 ) ) - person$EAP^2)
      #***
      # calculate EAP reliability
      # EAP variance
      EAP.variance <- stats::weighted.mean( person$EAP^2 , pweights ) - ( stats::weighted.mean( person$EAP , pweights ) )^2
      EAP.error <- stats::weighted.mean( person$SD.EAP^2 , pweights )
      EAP.rel <- EAP.variance / ( EAP.variance + EAP.error )	
    } else { 
      EAP.rel <- rep(0,ndim)
      names(EAP.rel) <- paste("Dim",1:ndim , sep="")
      for ( dd in 1:ndim ){
        #	dd <- 1  # dimension
        person$EAP <- rowSums( hwtE * outer( nstudl , theta[,dd] ) )
        person$SD.EAP <- sqrt(rowSums( hwtE * outer( nstudl , theta[,dd]^2 ) ) - person$EAP^2)	
        #***
        # calculate EAP reliability
        # EAP variance
        EAP.variance <- stats::weighted.mean( person$EAP^2 , pweights ) - ( stats::weighted.mean( person$EAP , pweights ) )^2
        EAP.error <- stats::weighted.mean( person$SD.EAP^2 , pweights )
        EAP.rel[dd] <- EAP.variance / ( EAP.variance + EAP.error )	
        colnames(person)[ which( colnames(person) == "EAP" ) ] <- paste("EAP.Dim" , dd , sep="")
        colnames(person)[ which( colnames(person) == "SD.EAP" ) ] <- paste("SD.EAP.Dim" , dd , sep="")				
      }
#      person <- data.frame( "pid" = pid , person )
    }
    #cat("person parameters") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1				  
    ############################################################
    s2 <- Sys.time()
	
	

    if (progress){
      cat(disp)
      cat("Regression Coefficients\n")
      print( beta , 4  )
      cat("\nVariance:\n" ) # , round( varianceM , 4 ))
      if (G==1 ){ 
        varianceM <- matrix( variance , nrow=ndim , ncol=ndim ) 
        print( varianceM , 4 )	
      } else { 
        print( variance[ var.indices] , 4 )	}
      if ( ndim > 1){ 
        cat("\nCorrelation Matrix:\n" ) # , round( varianceM , 4 ))	
        print( cov2cor(varianceM) , 4 )	
      }
      cat("\n\nEAP Reliability:\n")
      print( round (EAP.rel,3) )
      cat("\n-----------------------------")
      devmin <- which.min( deviance.history[,2] )
      if ( devmin < iter ){
        cat(paste("\n\nMinimal deviance at iteration " , devmin , 
                  " with deviance " , round(deviance.history[ devmin , 2 ],3) , sep="") , "\n")
        cat("The corresponding estimates are\n")
        cat("  xsi.min.deviance\n  beta.min.deviance \n  variance.min.deviance\n\n")
      }
      cat( "\nStart: " , paste(s1))
      cat( "\nEnd: " , paste(s2),"\n")
      print(s2-s1)
      cat( "\n" )
    }
    
    # Output list
    deviance.history <- deviance.history[ 1:iter , ]
    res <- list( "beta" = beta , "variance" = variance ,
                 "person" = person , pid = pid , "EAP.rel" = EAP.rel , 
                 "post" = hwt , "theta" = theta , 
                 "Y" = Y ,  "group" = group , 
                 "G" = if ( is.null(group)){1} else { length(unique( group ) )} , 
                 "groups" = if ( is.null(group)){1} else { groups } , 			   
                 "formulaY" = formulaY , "dataY" = dataY , 
                 "pweights" = pweights , 
                 "time" = c(s1,s2,s2-s1) , 
                 "nstud" = nstud , 
				 "hwt" = hwt ,  "like" = like , 
				 "ndim" = ndim ,
                 "beta.fixed" = beta.fixed , 
                 "variance.fixed" = variance.fixed ,
                 "nnodes" = nnodes , "deviance" = deviance ,
                 "ic" = ic , 
                 "deviance.history" = deviance.history ,
                 "control" = con1a ,    "iter" = iter ,
                 "YSD"=YSD , CALL = CALL 

    )
    class(res) <- "tam.latreg"
    return(res)
  }


# tam.mml.output <- function(){
# 	}
