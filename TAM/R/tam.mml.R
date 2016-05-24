tam.mml <-
  function( resp , Y=NULL , group = NULL ,  irtmodel ="1PL" ,
            formulaY = NULL , dataY = NULL , 
            ndim = 1 , pid = NULL ,
            xsi.fixed=NULL ,  xsi.inits = NULL , 
            beta.fixed = NULL , beta.inits = NULL , 
            variance.fixed = NULL , variance.inits = NULL , 
            est.variance = FALSE , constraint="cases" , 
            A=NULL , B=NULL , B.fixed = NULL , 
            Q=NULL , est.slopegroups=NULL , E = NULL , 
            pweights = NULL , 
			userfct.variance = NULL , variance.Npars = NULL , 
			item.elim = TRUE , 
			control = list() 
            # control can be specified by the user 
  ){
    
    #------------------------------------
    # INPUT:
    # Y ... matrix of regression predictors
    # group ... group indicators
    # formulaY ... a formula object for creating regression indicators
    # dataY ... corresponding data frame for formula object of Y
    # pid ... person ID
    # xsi.fixed  matrix L * 2 where L<=np
    #		1st column: xsi parameter label
    #		2nd column: fixed item parameter value
    # beta.fixed ... L*3 matrix
    # 	1st column: beta parameter integer predictor in Y
    #   2nd column: dimension
    #   3rd column: which parameter should be fixed?
    # B.fixed ... fixed B entries
    #			first three columns are B entries, fourth entry is the value to be fixed
    # est.variance ... should variance be estimated? Relevant for 2PL model
    #  2PL ... default FALSE
    # control:
    #		control = list( nodes = seq(-6,6,len=15) , 
    #		          			convD = .001 ,conv = .0001 , convM = .0001 , Msteps = 30 ,            
    #                   maxiter = 1000 , progress = TRUE) 
    # progress ... if TRUE, then display progress
    #-------------------------------------
    
    s1 <- Sys.time()
	CALL <- match.call()
    # display
    disp <- "....................................................\n"
    
    increment.factor <- progress <- nodes <- snodes <- ridge <- xsi.start0 <- QMC <- NULL
    maxiter <- conv <- convD <- min.variance <- max.increment <- Msteps <- convM <- NULL 
    R <- NULL
    
    if ( is.null(A)){ printxsi <- FALSE  } else { printxsi <- TRUE }
    # attach control elements
    e1 <- environment()
    con <- list( nodes = seq(-6,6,len=21) , snodes = 0 , QMC=TRUE , 
                 convD = .001 ,conv = .0001 , convM = .0001 , Msteps = 4 ,            
                 maxiter = 1000 , max.increment = 1 , 
                 min.variance = .001 , progress = TRUE , ridge=0 ,
                 seed = NULL , xsi.start0=FALSE , increment.factor=1 , fac.oldxsi=0 ,
				 acceleration="none" , dev_crit = "absolute" )  
	#@@@@AAAA@@@@@				 
    #a0 <- Sys.time()			   
    con[ names(control) ] <- control  
    Lcon <- length(con)
    con1a <- con1 <- con ; 
    names(con1) <- NULL
    for (cc in 1:Lcon ){
      assign( names(con)[cc] , con1[[cc]] , envir = e1 ) 
    }
    fac.oldxsi <- max( 0 , min( c( fac.oldxsi , .95 ) ) )  
    acceleration <- con$acceleration
	resp <- as.matrix(resp)
	resp0 <- resp <- add.colnames.resp(resp)
	
	
    if (progress){ 
      cat(disp)	
      cat("Processing Data     ", paste(Sys.time()) , "\n") ; flush.console()
    }  
    
    if ( ! is.null(group) ){ 
      con1a$QMC <- QMC <- FALSE
      con1a$snodes <- snodes <- 0
    }
    # define design matrix in case of PCM2
#    if (( irtmodel=="PCM2" ) & (is.null(Q)) & ( is.null(A)) ){ 
#      A <- .A.PCM2( resp ) 
#    }  
	#**** constraints
	if ( constraint == "items" ){
		irtmodel <- "PCM2"
						}
						
						
    if (( irtmodel=="PCM2" ) & ( is.null(A)) ){ 
      A <- .A.PCM2( resp , constraint=constraint , Q=Q  ) 	    
    }  
    
    
    if ( !is.null(con$seed)){ set.seed( con$seed )	 }
    
    nullY <- is.null(Y)
    
    nitems <- ncol(resp)       # number of items
    if (is.null(colnames(resp))){
      colnames(resp) <- paste0( "I" , 100+1:nitems )
    }
    nstud <- nrow(resp)        # number of students
	#*****
	nstud1 <- sum(1*( rowSums( 1 - is.na(resp) ) > 0 ))
	
    if ( is.null( pweights) ){
      pweights <- rep(1,nstud) # weights of response pattern
    }
    if (progress){ 
      cat("    * Response Data:" , nstud , "Persons and " , 
          nitems , "Items \n" )  ;
      flush.console()	  
    }  	  
    
    #!! check dim of person ID pid
    if ( is.null(pid) ){ pid <- seq(1,nstud) }else{ pid <- unname(c(unlist(pid))) }
    
    # print( colSums( is.na(resp)) )
    
    # normalize person weights to sum up to nstud
	pweights0 <- pweights	
	pweights <- nstud * pweights / sum(pweights)		
	#@@--

				
	# a matrix version of person weights
    pweightsM <- outer( pweights , rep(1,nitems) )
    
    # calculate ndim if only B or Q are supplied
    if ( ! is.null(B) ){ ndim <- dim(B)[3] } 
    if ( ! is.null(Q) ){ ndim <- dim(Q)[2] }
    
    betaConv <- FALSE         #flag of regression coefficient convergence
    varConv <- FALSE          #flag of variance convergence
    nnodes <- length(nodes)^ndim
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
    
    # maximum no. of categories per item. Assuming dichotomous
    maxK <- max( resp , na.rm=TRUE ) + 1 
    
    # create design matrices
    modeltype <- "PCM"
    if (irtmodel=="RSM"){  modeltype <- "RSM" }
	#****
	# ARb 2015-12-08
	maxKi <- NULL
	if ( ! (item.elim ) ){
		maxKi <- rep( maxK - 1 , ncol(resp) )		
				}
    #*** 
    design <- designMatrices( modeltype= modeltype , maxKi=maxKi , resp=resp , 
                              A=A , B=B , Q=Q , R=R, ndim=ndim  , constraint=constraint )
    A <- design$A
    B <- design$B
    cA <- design$flatA
    cA[is.na(cA)] <- 0
    if (progress){ 
      cat("    * Created Design Matrices   (", 
          paste(Sys.time()) , ")\n") ; flush.console()	  
    }    
    
    design <- NULL				
   
    #   #---2PL---
    #   B_orig <- B  #keep a record of generated B before estimating it in 2PL model 
    #   #---end 2PL---
    
    ################################
    # number of parameters
    np <- dim(A)[[3]]
    # xsi inits
    if ( ! is.null(xsi.inits) ){
      #    xsi <- xsi.inits 
      xsi <- rep(0,np)
      xsi[ xsi.inits[,1] ] <- xsi.inits[,2]	
    } else { xsi <- rep(0,np)   } 
    
    
    if ( ! is.null( xsi.fixed ) ){
      xsi[ xsi.fixed[,1] ] <- xsi.fixed[,2]
      est.xsi.index <- setdiff( 1:np , xsi.fixed[,1] )
				} else { 
		est.xsi.index <- 1:np 
						}
    est.xsi.index -> est.xsi.index0
    
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
      # user must label groups from 1, ... , G
      #    if ( length( setdiff( 1:G , groups)  ) > 0 ){
      #      stop("Label groups from 1, ...,G\n")
      #				}							
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
    
    W <- t(Y * pweights) %*% Y
    if (ridge > 0){ diag(W) <- diag(W) + ridge }
    YYinv <- solve( W )
    
    #initialise regressors
    if ( is.null(beta.fixed) & (  is.null(xsi.fixed) ) & ( constraint=="cases") ){
      beta.fixed <- matrix( c(1,1,0) , nrow= 1) 
      if (  ndim > 1){ 
        for ( dd in 2:ndim){
          beta.fixed <- rbind( beta.fixed , c( 1 , dd , 0 ) )
        }}}
    
    #****
    # ARb 2013-08-20: Handling of no beta constraints	
    # ARb 2013-08-24: correction	
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
    #	if ( ! is.null( beta.inits ) ){ 
    #		beta <- beta.inits 
    #  } else {
    #		beta <- matrix(0, nrow = nreg+1 , ncol = ndim)        
    #			}
    
    
    # define response indicator matrix for missings
    resp.ind <- 1 - is.na(resp)
    nomiss <- sum( is.na(resp) ) == 0
    #*** included nomiss in M step regression
    resp.ind.list <- list( 1:nitems )
    for (i in 1:nitems){ resp.ind.list[[i]] <- which( resp.ind[,i] == 1)  }
    resp[ is.na(resp) ] <- 0 	# set all missings to zero
    #stop("here")  
    #@@ ARb:Include Starting values for xsi??
    #       xsi <- - qnorm( colMeans( resp ) )
    AXsi <- matrix(0,nrow=nitems,ncol=maxK )  #A times xsi
    
    # Create an index linking items and parameters
    indexIP <- colSums(aperm(A, c(2,1,3)) != 0, na.rm = TRUE)
    # define list of elements for item parameters
    indexIP.list <- list( 1:np )
    for ( kk in 1:np ){
      indexIP.list[[kk]] <- which( indexIP[,kk] > 0 )
    }
    lipl <- cumsum( sapply( indexIP.list , FUN = function(ll){ length(ll) } ) )
    indexIP.list2 <- unlist(indexIP.list)
    indexIP.no <- as.matrix( cbind( c(1 , lipl[-length(lipl)]+1 ) , lipl ) )
    
    
    # These sufficient statistics must be changed
    # to make it more general
    # First extension:  pweights and dependent on A; needs to be further extended (e.g., different number of categories)
    # Second extension: multiple category option       -> resp \in 0:maxKi (see method definition calc_posterior_TK)
    #                                                  -> length(ItemScore) = np (see diff computation in M Step)
    #                   multiple category option Bugfix
    #                                                  -> dim(cResp) = (nstud, nitems*maxK)
    #                                                  -> adapt dim(A) to dim(cResp) for 
    #							sufficient statistic (cf. print.designMatrices)
    #**************************************************	
    # ARb 2013-04-29
    # more efficient calculation of sufficient statistics
    col.index <- rep( 1:nitems , each = maxK )
    #  cResp <- (resp[ , col.index  ]+1) *resp.ind[ , col.index ]
    cResp <- (resp +1) *resp.ind
    
    cResp <- cResp[ , col.index  ]
    #   cResp <- 1 * t( t(cResp) == rep(1:(maxK), nitems) )
    cResp <- 1 * ( cResp == matrix( rep(1:(maxK), nitems) , nrow(cResp) , 
                                    ncol(cResp) , byrow=TRUE ) )
    if ( sd(pweights) > 0 ){ 
      ItemScore <- as.vector( t( colSums( cResp * pweights ) ) %*% cA )
    } else { 
      ItemScore <- as.vector( t( colSums( cResp) ) %*% cA )			
    }
    if (progress){ 
      cat("    * Calculated Sufficient Statistics   (", 
          paste(Sys.time()) , ")\n") ; flush.console()	  
    }    				   
    
    #**************************************************
    # starting values for xsi
    maxAi <-  - (apply(-(A) , 3 , rowMaxs , na.rm=TRUE))  
    personMaxA <- resp.ind %*% maxAi
    #ItemMax <- personMaxA %t*% pweights  
	ItemMax <- crossprod( personMaxA , pweights ) 
    xsi[est.xsi.index] <- - log(abs(( ItemScore[est.xsi.index]+.5)/
                                      (ItemMax[est.xsi.index]-ItemScore[est.xsi.index]+.5) ) )
    # starting values of zero
    if( xsi.start0 ){ xsi <- 0*xsi }
    
    #log of odds ratio of raw scores  
    if ( ! is.null(xsi.inits) ){   
      #		xsi <- xsi.inits  
      xsi[ xsi.inits[,1] ] <- xsi.inits[,2]
    }
    if ( ! is.null( xsi.fixed ) ){   xsi[ xsi.fixed[,1] ] <- xsi.fixed[,2] }
    
    xsi.min.deviance <- xsi
    beta.min.deviance <- beta
    variance.min.deviance <- variance
    
# cat("b200"); a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1  									  
    
    # nodes
    if ( snodes == 0 ){ 
      theta <- as.matrix( expand.grid( as.data.frame( matrix( rep(nodes, ndim) , ncol = ndim ) ) ) )
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
        fac <- 1
		# fac <- 2
        r1 <- sfsmisc::QUnif(n=snodes, min = 0, max = 1, 
						n.min = 1, p=ndim, leap = 409)
        theta0.samp <- fac * stats::qnorm( r1 )
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
    #---2PL---
    #   if (irtmodel %in% c("2PL","GPCM" , "GPCM.design", "2PL.groups") ){
    # 		a4 <- 999   
    # 		basispar <- NULL
    # 		} else{  a4 <- 0 }
    
    #   if (irtmodel == "GPCM.design" ){
    # 		ES <- ncol(E)
    # 		basispar <- matrix( 1 , ES , ndim )	
    # 		basispar1 <- solve( t(E) %*% E	, t(E) %*% rep(1,nitems))
    # 		for ( dd in 1:ndim){
    # 			basispar[,dd] <- basispar1
    # 						}
    # 			   }
    # 	
    #    if ( irtmodel %in% c("2PL","GPCM" , "GPCM.design","2PL.groups") ){
    # 	if ( ! is.null(B.fixed) ){
    # 			B[ B.fixed[,1:3] ] <- B.fixed[,4]	
    # 			B_orig[ B.fixed[,1:3] ] <- 0
    # 						}
    # 						}
    #   #---end 2PL---
    
    ##**SE
    se.xsi <- 0*xsi
    se.B <- 0*B
    
    YSD <- max( apply( Y , 2 , sd ) )
    if (YSD > 10^(-15) ){ YSD <- TRUE } else { YSD <- FALSE }
    
	#*****
	#@@@@ speed gains, further auxiliary objects, 2015-06-26
	Avector <- as.vector(A)
	Avector[ is.na(Avector) ] <- 0
	unidim_simplify <- TRUE
	if (G > 1){ unidim_simplify <- FALSE }
	if ( YSD){ unidim_simplify <- FALSE }	
	if (  is.null(beta.fixed) ){ unidim_simplify <- FALSE }
	
	# unidim_simplify <- FALSE
	#@@@@
	
    # define progress bar for M step
    mpr <- round( seq( 1 , np , len = 10 ) )
    
    hwt.min <- 0
    rprobs.min <- 0
    AXsi.min <- 0
    B.min <- 0
    deviance.min <- 1E100
    itemwt.min <- 0
    se.xsi.min <- se.xsi
    se.B.min <- se.B
    
	#@@@@AAAA@@@@@
	# acceleration
	xsi_acceleration <- list( "acceleration" = acceleration , "w" = .35 ,
							"w_max" = .95 , 
							parm_history = cbind( xsi, xsi , xsi ) ,
							"beta_new" = 0 ,
							"beta_old" = 0 
						)
	v1 <- as.vector(variance)					
	variance_acceleration <- list( "acceleration" = acceleration , "w" = .35 ,
							"w_max" = .95 , 
							parm_history = cbind( v1, v1 , v1) ,
							"beta_new" = 0 ,
							"beta_old" = 0 
						)
	if (G>1){  variance_acceleration$acceleration <- "none" }						
	#@@@@AAAA@@@@@			

	
	#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    ##############################################################   
    #Start EM loop here
    while ( ( (!betaConv | !varConv)  | ((a1 > conv) | (a4 > conv) | (a02 > convD)) )  & (iter < maxiter) ) { 
      	  	  
      iter <- iter + 1
      if (progress){ 
        cat(disp)	
        cat("Iteration" , iter , "   " , paste( Sys.time() ) )
        cat("\nE Step\n") ; flush.console()
      }
      # calculate nodes for Monte Carlo integration	
      if ( snodes > 0){
        #      theta <- beta[ rep(1,snodes) , ] +  t ( t(chol(variance)) %*% t(theta0.samp) )
		
		# theta <- beta[ rep(1,snodes) , ] + fac0* (theta0.samp %*% chol(variance) )
		theta <- beta[ rep(1,snodes) , ] + theta0.samp %*% chol(variance) 
		
        # calculate density for all nodes
        thetasamp.density <- mvtnorm::dmvnorm( theta , mean = as.vector(beta[1,]) , sigma = variance )
        # recalculate theta^2
        #      theta2 <- matrix( theta.sq(theta) , nrow=nrow(theta) , ncol=ncol(theta)^2 )   
        theta2 <- matrix( theta.sq2(theta) , nrow=nrow(theta) , ncol=ncol(theta)^2 )   
      }			
      olddeviance <- deviance
# a0 <- Sys.time()	
      # calculation of probabilities
      res <- calc_prob.v5(iIndex=1:nitems , A=A , AXsi=AXsi , B=B , xsi=xsi , theta=theta , 
                          nnodes=nnodes , maxK=maxK , recalc=TRUE )	
      rprobs <- res[["rprobs"]]
      AXsi <- res[["AXsi"]]
	  
      #***	      
      # AXsi[ is.na(AXsi) ] <- 0
      
#cat("calc_prob") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1
      
      # calculate student's prior distribution
      gwt <- stud_prior.v2(theta=theta , Y=Y , beta=beta , variance=variance , nstud=nstud , 
                           nnodes=nnodes , ndim=ndim,YSD=YSD , unidim_simplify=unidim_simplify ,
						   snodes = snodes )
						   
# cat("stud_prior") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1						 

      # calculate student's likelihood
      res.hwt <- calc_posterior.v2(rprobs=rprobs , gwt=gwt , resp=resp , nitems=nitems , 
                                   resp.ind.list=resp.ind.list , normalization=TRUE , 
                                   thetasamp.density=thetasamp.density , snodes=snodes ,
                                   resp.ind=resp.ind  ,  avoid.zerosum=TRUE	)	
      hwt <- res.hwt[["hwt"]]   
# cat("calc_posterior") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1						 
      
      
      if (progress){ cat("M Step Intercepts   |"); flush.console() }
      # collect old values for convergence indication
      oldxsi <- xsi
      oldbeta <- beta
      oldvariance <- variance 
# cat("before mstep regression") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1						           

      # M step: estimation of beta and variance
      resr <- mstep.regression( resp=resp , hwt=hwt , resp.ind=resp.ind , pweights=pweights , 
                                pweightsM=pweightsM , Y=Y , theta=theta , theta2=theta2 , YYinv=YYinv , 
                                ndim=ndim , nstud=nstud , beta.fixed=beta.fixed , variance=variance , 
                                Variance.fixed=variance.fixed , group=group ,  G=G , snodes = snodes ,
                                nomiss=nomiss ,  thetasamp.density= thetasamp.density )
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

        #@@@@AAAA@@@@@
		# variance acceleration
		if ( variance_acceleration$acceleration != "none" ){		
			variance_acceleration <- accelerate_parameters( xsi_acceleration=variance_acceleration , 
							xsi=as.vector(variance) , iter=iter , itermin=3)
			variance <- matrix( variance_acceleration$parm , nrow= nrow(variance) , ncol=ncol(variance) )
								}
	    #@@@@AAAA@@@@@

				
      
      if (max(abs(variance-oldvariance)) < conv) varConv <- TRUE
      
      ######################################
      #M-step
      converge <- FALSE
      Miter <- 1
      old_increment <- rep( max.increment , np )
      est.xsi.index <- est.xsi.index0
      while (!converge & ( Miter <= Msteps ) ) {	
        # xbar2 <- xxf <- xbar <- rep(0,np)
        # Only compute probabilities for items contributing to param p
# z0 <- Sys.time()

        
        if (Miter > 1){ 
          res.p <- calc_prob.v5( iIndex=1:nitems , A=A , AXsi=AXsi , B=B , 
                                 xsi=xsi , theta=theta , nnodes=nnodes, maxK=maxK)					
          rprobs <- res.p[["rprobs"]]            
        }
# cat("calc prob") ; z1 <- Sys.time(); print(z1-z0) ; z0 <- z1		

	
        res <- calc_exp_TK3( rprobs , A , np , est.xsi.index , itemwt ,
                             indexIP.no , indexIP.list2 , Avector )
# cat("calc_exp") ; z1 <- Sys.time(); print(z1-z0) ; z0 <- z1									 
        xbar <- res$xbar
        xbar2 <- res$xbar2
        xxf <- res$xxf	
        
        # Compute the difference between sufficient statistic and expectation
        diff <- as.vector(ItemScore) - xbar
        #Compute the Newton-Raphson derivative for the equation to be solved
        deriv <- xbar2 - xxf 			
        increment <- diff*abs(1/( deriv + 10^(-20) ) )
        if ( !is.null( xsi.fixed) ){ increment[ xsi.fixed[,1] ] <- 0 } 
        #!!!	  necesessary to include statement to control increment?
        ci <- ceiling( abs(increment) / ( abs( old_increment) + 10^(-10) ) )
        increment <- ifelse( abs( increment) > abs(old_increment)  , 
                             increment/(2*ci) , 
                             increment )						
#        increment <- ifelse( abs( increment) > abs(old_increment)  , 
#                             sign(increment) * max.increment , increment )	

							 
        old_increment <- increment
        
        ##**SE
        se.xsi <- sqrt( 1 / abs(deriv) )
        if ( ! is.null( xsi.fixed) ){ se.xsi[ xsi.fixed[,1] ] <- 0 } 
        ##**	
        
        xsi <- xsi+increment   # update parameter p
        #	est.xsi.index <- which( abs(increment) > convM )
		
# ask for est.xsi.index
		
        if ( max(abs(increment)) < convM ) { converge <- TRUE }
        Miter <- Miter + 1						
        
        # stabilizing the algorithm | ARb 2013-09-10
        if (fac.oldxsi > 0 ){
          xsi <-  (1-fac.oldxsi) * xsi + fac.oldxsi *oldxsi
        }	
		

        
        # progress bar
        if (progress){ 
          #        cat( paste( rep("-" , sum( mpr == p ) ) , collapse="" ) )
          cat("-") ;    flush.console()
        }
# cat("rest") ; z1 <- Sys.time(); print(z1-z0) ; z0 <- z1	
		
      } # end of all parameters loop

        #@@@@AAAA@@@@@
		# acceleration
		if ( xsi_acceleration$acceleration != "none" ){		
			xsi_acceleration <- accelerate_parameters( xsi_acceleration=xsi_acceleration , 
							xsi=xsi , iter=iter , itermin=3)
			xsi <- xsi_acceleration$parm
								}
	    #@@@@AAAA@@@@@
	  
# cat("mstep item parameters") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1						 	

     
      #***
      # decrease increments in every iteration
      if( increment.factor > 1){max.increment <-  1 / increment.factor^iter }    
      
      # calculate deviance
      if ( snodes == 0 ){ 
        deviance <- - 2 * sum( pweights * log( res.hwt$rfx * thetawidth ) )
      } else {
        #       deviance <- - 2 * sum( pweights * log( res.hwt$rfx ) )

        # deviance <- - 2 * sum( pweights * log( rowMeans( res.hwt$swt ) ) )
		#deviance <- - 2 * sum( pweights * log( rowSums( res.hwt$swt ) ) )
#Revalpr("deviance")		
		# deviance <- - 2 * sum( pweights * log( res.hwt$rfx  ) )
		deviance <- - 2 * sum( pweights * log( res.hwt$rfx   ) )
		      }
      deviance.history[iter,2] <- deviance
      a01 <- abs( ( deviance - olddeviance ) / deviance  )
      a02 <- abs( ( deviance - olddeviance )  )	
	  	  
	  if (con$dev_crit == "relative" ){ a02 <- a01 }
      
      if( deviance > deviance.min ){ 	
#      if( ( deviance - olddeviance < 0 ) | ( iter == 1)  ){ 
        xsi.min.deviance <- xsi.min.deviance 
        beta.min.deviance <- beta.min.deviance
        variance.min.deviance <- variance.min.deviance
        hwt.min <- hwt.min
        rprobs.min <- rprobs.min
        AXsi.min <- AXsi.min
        B.min <- B.min
        deviance.min <- deviance.min
        itemwt.min <- itemwt.min
        se.xsi.min <- se.xsi.min
        se.B.min <- se.B.min
      }   else { 
        xsi.min.deviance <- xsi 
        beta.min.deviance <- beta
        variance.min.deviance <- variance	
        hwt.min <- hwt	
        AXsi.min <- AXsi	
        B.min <- B
        deviance.min <- deviance
        itemwt.min <- itemwt
        se.xsi.min <- se.xsi
        se.B.min <- se.B
      }
      
      
      
      a1 <- max( abs( xsi - oldxsi ))	
      a2 <- max( abs( beta - oldbeta ))	
      a3 <- max( abs( variance - oldvariance ))
      if (progress){ 
        cat( paste( "\n  Deviance =" , round( deviance , 4 ) ))
        devch <- -( deviance - olddeviance )
        cat( " | Deviance change:", round( devch  , 4 ) )
		cat( " | Relative deviance change:", round( a01  , 8 ) )
        if ( devch < 0 & iter > 1 ){ 
          cat ("\n!!! Deviance increases!                                        !!!!") 
          cat ("\n!!! Choose maybe fac.oldxsi > 0 and/or increment.factor > 1    !!!!") 			
        }
        
        
        cat( "\n  Maximum intercept parameter change:" , round( a1 , 6 ) )
        # 	  if (irtmodel %in% c("GPCM","2PL","2PL.group") ){
        # 		cat( "\n  Maximum slope parameter change:" , round( a4 , 6 ) )
        # 							}
        cat( "\n  Maximum regression parameter change:" , round( a2 , 6 ) )  
        if ( G == 1 ){ 
          cat( "\n  Variance: " , round( variance[ ! lower.tri(variance)] , 4 ) , " | Maximum change:" , round( a3 , 6 ) )  
        } else {
          cat( "\n  Variance: " , round( variance[var.indices] , 4 ) ,
               " | Maximum change:" , round( a3 , 6 ) )  		
        }					
        cat( "\n  beta ",round(beta,4)  )
        cat( "\n" )
        flush.console()
      }
      
# cat("rest") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1	

# Revalpr("variance[1,1]")
      
      
    } # end of EM loop
    #******************************************************
	
	
# stop()	
# a0 <- Sys.time()	
    xsi.min.deviance -> xsi 
    beta.min.deviance -> beta
    variance.min.deviance -> variance	
    hwt.min -> hwt	
    AXsi.min -> AXsi	
    B.min -> B
    deviance.min -> deviance
    itemwt.min -> itemwt
    se.xsi.min -> se.xsi	
    se.B.min -> se.B		
    #******
    # a0 <- Sys.time()  
    # cat("restructure output") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1	
    
    ##**SE  
    # standard errors of AXsi parameters
    # check for missing entries in A
    se.AXsi <- 0*AXsi
    A1 <- A
    A1[ is.na(A) ] <- 0
    if ( length( se.xsi) > 1){
      se.xsiD <- diag( se.xsi^2 )
    } else {
      se.xsiD <- matrix( se.xsi^2,1,1)
    }
		
	#******
	# generate input for fixed parameters
	xsi.fixed.estimated <- generate.xsi.fixed.estimated( xsi , A )
	B.fixed.estimated <- generate.B.fixed.estimated(B)
	
    #*******
    # A1 [ items ,categs , params ]  -> xsi design matrix
    # se.xsiD  [ params , params ]
    #@@ ARb 2013-08-27: Correction in case of one item
    for (kk in 1:maxK){  # kk <- 1
      #	se.AXsi[,kk] <- sqrt( diag( A1[,kk,] %*% se.xsiD %*% t( A1[,kk,]) ) )
      #**** bugfix
	  dim_A1 <- dim(A1)
      A1_kk <- A1[,kk,]
      if ( is.vector(A1_kk) ){
        A1_kk <- matrix( A1_kk , nrow=dim_A1[1] , ncol=dim_A1[3] )
      }
      se.AXsi[,kk] <- sqrt( diag( A1_kk %*% se.xsiD %*% t( A1_kk ) ) )	
      #****	
    }
#    cat("se.axsi") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1	
    
    ##*** Information criteria
    ic <- .TAM.ic( nstud=nstud1 , deviance , xsi , xsi.fixed ,
                   beta , beta.fixed , ndim , variance.fixed , G ,
                   irtmodel , B_orig=NULL , B.fixed , E , est.variance =TRUE ,
                   resp ,  est.slopegroups=NULL , 
				   variance.Npars=variance.Npars , group )
#    cat("TAM.ic") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1			

    #***
    # calculate counts
    res <- .tam.calc.counts( resp, theta , resp.ind , 
                             group , maxK , pweights , hwt )
    n.ik <- res$n.ik
    pi.k <- res$pi.k 
#    cat("calculate counts") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1				
    #****
    # collect item parameters
    item1 <- .TAM.itempartable( resp , maxK , AXsi , B , ndim ,
                                resp.ind , rprobs,n.ik,pi.k)											
#     cat("tam itempartable") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1								 
    #####################################################
    # post ... posterior distribution	
    # create a data frame person	
    person <- data.frame( "pid"=pid , "case" = 1:nstud , "pweight" = pweights )
    person$score <- rowSums( resp * resp.ind )
		
    # use maxKi here; from "design object"
    nstudl <- rep(1,nstud)
    person$max <- rowSums( outer( nstudl , apply( resp ,2 , max , na.rm=T) ) * resp.ind )
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
      EAP.variance <- weighted.mean( person$EAP^2 , pweights ) - ( weighted.mean( person$EAP , pweights ) )^2
      EAP.error <- weighted.mean( person$SD.EAP^2 , pweights )
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
        EAP.variance <- weighted.mean( person$EAP^2 , pweights ) - ( weighted.mean( person$EAP , pweights ) )^2
        EAP.error <- weighted.mean( person$SD.EAP^2 , pweights )
        EAP.rel[dd] <- EAP.variance / ( EAP.variance + EAP.error )	
        colnames(person)[ which( colnames(person) == "EAP" ) ] <- paste("EAP.Dim" , dd , sep="")
        colnames(person)[ which( colnames(person) == "SD.EAP" ) ] <- paste("SD.EAP.Dim" , dd , sep="")				
      }
#      person <- data.frame( "pid" = pid , person )
    }
#    cat("person parameters") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1				  
    ############################################################
    s2 <- Sys.time()
    if ( is.null( dimnames(A)[[3]] ) ){  
      dimnames(A)[[3]] <- paste0("Xsi" , 1:dim(A)[3] )
    }
    item <- data.frame( "xsi.index" = 1:np , 
                        "xsi.label" = dimnames(A)[[3]] , 
                        "est" = xsi )
    if (progress){
      cat(disp)
      cat("Item Parameters\n")
      item2 <- item
      item2[,"est"] <- round( item2[,"est"] , 4 )
      print(item2)
      cat("...................................\n")
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
    
    # collect xsi parameters
    obji <- data.frame( "xsi"=xsi , "se.xsi"=se.xsi ) 
    rownames(obji) <- dimnames(A)[[3]]	
    xsi <- obji
    
	#**** calculate individual likelihood
      res.hwt <- calc_posterior.v2(rprobs=rprobs , gwt=1+0*gwt , resp=resp , nitems=nitems , 
                                   resp.ind.list=resp.ind.list , normalization=FALSE , 
                                   thetasamp.density=thetasamp.density , snodes=snodes ,
                                   resp.ind=resp.ind	)	
      res.like <- res.hwt[["hwt"]] 	
	
	
    # Output list
    deviance.history <- deviance.history[ 1:iter , ]
    res <- list( "xsi" = xsi ,
                 "beta" = beta , "variance" = variance ,
                 "item" = item1 , 
                 "person" = person , pid = pid , "EAP.rel" = EAP.rel , 
                 "post" = hwt ,  "rprobs" = rprobs , "itemweight" = itemwt ,
                 "theta" = theta , 
                 "n.ik" = n.ik , "pi.k" = pi.k ,
                 "Y" = Y , "resp" = resp0 , 
                 "resp.ind" = resp.ind , "group" = group , 
                 "G" = if ( is.null(group)){1} else { length(unique( group ) )} , 
                 "groups" = if ( is.null(group)){1} else { groups } , 			   
                 "formulaY" = formulaY , "dataY" = dataY , 
                 "pweights" = pweights0 , 
                 "time" = c(s1,s2,s2-s1) , "A" = A , "B" = B  ,
                 "se.B" = se.B , 
                 "nitems" = nitems , "maxK" = maxK , "AXsi" = AXsi ,
                 "AXsi_" = - AXsi ,
                 "se.AXsi" = se.AXsi , 
                 "nstud" = nstud , "resp.ind.list" = resp.ind.list ,
                 "hwt" = hwt ,  "like" = res.like , 
				 "ndim" = ndim ,
                 "xsi.fixed" = xsi.fixed , 
				 "xsi.fixed.estimated" = xsi.fixed.estimated , 
				 "B.fixed.estimated" = B.fixed.estimated ,
				 "beta.fixed" = beta.fixed , "Q" = Q  ,
                 "variance.fixed" = variance.fixed ,
                 "nnodes" = nnodes , "deviance" = deviance ,
                 "ic" = ic , 
                 "deviance.history" = deviance.history ,
                 "control" = con1a , "irtmodel" = irtmodel ,
                 "iter" = iter ,
                 "printxsi" = printxsi ,
                 "YSD"=YSD , CALL =CALL
                 #			   "design"=design
                 #			   "xsi.min.deviance" = xsi.min.deviance ,
                 #			   "beta.min.deviance" = beta.min.deviance , 
                 # "variance.min.deviance" = variance.min.deviance 
    )
    class(res) <- "tam.mml"
    return(res)
  }


# tam.mml.output <- function(){
# 	}
