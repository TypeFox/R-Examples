tam.mml.3pl <-
  function( resp , Y=NULL , group = NULL ,  
            formulaY = NULL , dataY = NULL , 
            ndim = 1 , pid = NULL ,
            xsi.fixed=NULL ,  xsi.inits = NULL , xsi.prior = NULL , 
            beta.fixed = NULL , beta.inits = NULL , 
            variance.fixed = NULL , variance.inits = NULL , 
            est.variance = TRUE , 
            A=NULL , notA = FALSE , Q=NULL , 
			Q.fixed = NULL , 
			E = NULL , gammaslope.des = "2PL" , 
			gammaslope=NULL , gammaslope.fixed=NULL ,
            est.some.slopes=TRUE ,    			
			gammaslope.constr.V=NULL , gammaslope.constr.c = NULL , 
			gammaslope.center.index=NULL ,  gammaslope.center.value=NULL , 
			gammaslope.prior=NULL ,  userfct.gammaslope = NULL ,
			gammaslope.constr.Npars = 0 ,
			est.guess = NULL ,  guess = rep(0,ncol(resp)) , 
			guess.prior=NULL ,
            skillspace = "normal" , theta.k = NULL , 
            delta.designmatrix=NULL , delta.fixed=NULL , 
			delta.inits = NULL ,  pweights = NULL , 
			item.elim = TRUE , 
			control = list() ,	Edes = NULL  
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
    maxgamma <- maxiter <- conv <- convD <- min.variance <- max.increment <- Msteps <- convM <- NULL 
    delta <- R <- NULL
    B <- NULL ; B.fixed <- NULL ; theta <- NULL	
	irtmodel <- "2PL" 
	est.slopegroups <- NULL	
	init.gammaslope <- ( ! is.null( gammaslope ) )
	
	
	resp <- as.matrix(resp)
	resp0 <- resp <- add.colnames.resp(resp)
	
	#********************
	# create E design matrix from different input matrices
	res0 <- .mml.3pl.create.E( resp , E , Q , gammaslope.des , Q.fixed )		
	E <- res0$E
	if ( is.null(gammaslope.fixed ) ){
		gammaslope.fixed <- res0$gammaslope.fixed
								}
	
	# calculation of not A if requested
	if ( notA) {
		res <- .mml.3pl.create.notA( E , notA )	
		A <- res$A
		xsi.fixed <- res$xsi.fixed
				}
						
	#********************			
	#********************
	# compute B from E or an input statement
	 # starting values gammaslope
	  if ( is.null( gammaslope ) ){
		   Ngam <- dim(E)[4]
		   if ( est.some.slopes){
			gammaslope <- stats::runif( Ngam , .9 , 1.1 )
					} else { 
			gammaslope <- rep(1,Ngam ) 
					}
				}

	if ( is.null(Edes) ){			
		Edes <- .Call("mml_3pl_nonzero_entries", as.vector(E) , dim(E) ,
				   PACKAGE="TAM")$E_design
						}
											
	 # B <-.mml.3pl.computeB( E , gammaslope )
	 B <- .mml.3pl.computeB.v2( Edes , gammaslope , E )
	 
	 #*********************** 	
    if ( is.null(A)){ printxsi <- FALSE  } else { printxsi <- TRUE }  
    # attach control elements
    e1 <- environment()
    con <- list( nodes = seq(-6,6,len=21) , snodes = 0 ,QMC=TRUE,
                 convD = .001 ,conv = .0001 , convM = .0001 , Msteps = 10 ,            
                 maxiter = 1000 , max.increment = 1 , 
                 min.variance = .001 , progress = TRUE , ridge=0,seed=NULL,
                 xsi.start0=FALSE , increment.factor=1 , fac.oldxsi=0 ,
				 maxgamma = 9.99 , acceleration="none" , dev_crit = "absolute"  )  	
    con[ names(control) ] <- control  
    Lcon <- length(con)
    con1a <- con1 <- con ; 
    names(con1) <- NULL
    for (cc in 1:Lcon ){
      assign( names(con)[cc] , con1[[cc]] , envir = e1 ) 
    }
    if ( !is.null(con$seed)){ set.seed( con$seed )	 }
    fac.oldxsi <- max( 0 , min( c( fac.oldxsi , .95 ) ) )
    acceleration <- con$acceleration
	
    
    if (progress){ 
      cat(disp)	
      cat("Processing Data     ", paste(Sys.time()) , "\n") ; flush.console()
    }     
    if ( ! is.null(group) ){ 
      con1a$QMC <- QMC <- FALSE
      con1a$snodes <- snodes <- 0
    }   
    if ( is.null(group) ){ 
	  group1 <- rep(1,nrow(resp) )
    }   	
    
    nullY <- is.null(Y)
    
    # define design matrix in case of PCM2
    if (( irtmodel=="PCM2" ) & (is.null(Q)) & ( is.null(A)) ){ 
      A <- .A.PCM2( resp ) 
    }  

  if ( ! is.null( variance.fixed ) ){
			est.variance <- TRUE			
				}    
	
	
	# manage guessing parameters
	if ( is.null(guess) ){
	   guess <- rep( 0 , ncol(resp) )
						}
	
	
	if ( ! is.null(est.guess) ){
#	    guess <- rep(0,ncol(resp))
		h1 <- setdiff( unique(est.guess) , 0 )
		est.guess <- match( est.guess ,h1 )
		est.guess[ is.na( est.guess ) ] <- 0
					}
    est.some.guess <- sum( est.guess > 0 )
    
    nitems <- ncol(resp)       # number of items
    nstud <- nrow(resp)        # number of students
    
	#*****
	nstud100 <- sum(1*( rowSums( 1 - is.na(resp) ) > 0 ))
  
  
	if ( is.null( pweights) ){
      pweights <- rep(1,nstud) # weights of response pattern
    }
    
    if (progress){ 
      cat("    * Response Data:" , nstud , "Persons and " , 
          nitems , "Items \n" )  ;
      flush.console()	  
    }  	  
    
    #!! check dim of person ID pid
    if ( is.null(pid) ){ pid <- seq(1,nstud) }
        
    # normalize person weights to sum up to nstud
	pweights0 <- pweights
    pweights <- nstud * pweights / sum(pweights)
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
	if ( skillspace=="normal"){
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
	}
    #*********  
    
    # maximum no. of categories per item. Assuming dichotomous
    maxK <- max( resp , na.rm=TRUE ) + 1 
	
  
	#****
	# ARb 2015-12-15
	maxKi <- NULL
	if ( ! (item.elim ) ){
		maxKi <- rep( maxK - 1 , ncol(resp) )		
				}
    #***     	
	
    # create design matrices
    design <- designMatrices( modeltype="PCM" , maxKi=NULL , resp=resp , 
                              A=A , B=B , Q=Q , R=R, ndim=ndim )
    A <- design$A	
    B <- design$B
    cA <- design$flatA
    cA[is.na(cA)] <- 0
    if (progress){ 
      cat("    * Created Design Matrices   (", 
          paste(Sys.time()) , ")\n") ; flush.console()	  
    }    
    design <- NULL
    #---2PL---
    B_orig <- B  #keep a record of generated B before estimating it in 2PL model 
    #---end 2PL---
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
    } else { est.xsi.index <- 1:np }
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
	  group1 <- group
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
      colnames(Y) <- paste("group" , groups , sep="")
      for (gg in 1:G){ Y[,gg] <- 1*(group==gg) }
      nreg <- G - 1
    }
    
    W <- t(Y * pweights) %*% Y
    if (ridge > 0){ diag(W) <- diag(W) + ridge }
    YYinv <- solve( W )
    
    #initialise regressors
    if ( is.null(beta.fixed) & (  is.null(xsi.fixed) ) ){
      beta.fixed <- matrix( c(1,1,0) , nrow= 1) 
      if ( ndim > 1){ 
        for ( dd in 2:ndim){
          beta.fixed <- rbind( beta.fixed , c( 1 , dd , 0 ) )
        }}}  
    
    #****	
    if( ! is.matrix(beta.fixed) ){
      if ( ! is.null(beta.fixed) ){
        if ( ! beta.fixed   ){ beta.fixed <- NULL }
      }
    }
    #*****
    beta <- matrix(0, nrow = nreg+1 , ncol = ndim)  
    if ( ! is.null( beta.inits ) ){ 
      beta[ beta.inits[,1:2] ] <- beta.inits[,3]
    }
    # define response indicator matrix for missings
    resp.ind <- 1 - is.na(resp)
    #  nomiss <- sum( is.na(resp) == 0 )  	#*** included nomiss in M step regression  
    nomiss <- sum( is.na(resp) ) == 0
    resp.ind.list <- list( 1:nitems )
    for (i in 1:nitems){ resp.ind.list[[i]] <- which( resp.ind[,i] == 1)  }
    resp[ is.na(resp) ] <- 0 	# set all missings to zero
    
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
    
    #**************************************************	
    col.index <- rep( 1:nitems , each = maxK )
    cResp <- (resp +1) *resp.ind
    cResp <- cResp[ , col.index  ]
    cResp <- 1 * ( cResp == matrix( rep(1:(maxK), nitems) , nrow(cResp) , 
                                    ncol(cResp) , byrow=TRUE ) )  
    if ( stats::sd(pweights) > 0 ){ 
      ItemScore <- as.vector( t( colSums( cResp * pweights ) ) %*% cA )
    } else { 
      ItemScore <- as.vector( t( colSums( cResp) ) %*% cA )			
    }
    #**************************************************
    
    if (progress){ 
      cat("    * Calculated Sufficient Statistics   (", 
          paste(Sys.time()) , ")\n") ; flush.console()	  
    }  


    # starting values for xsi
    maxAi <-  - (apply(-(A) , 3 , rowMaxs , na.rm=TRUE))  
    personMaxA <- resp.ind %*% maxAi
    # ItemMax <- personMaxA %t*% pweights  
	ItemMax <- crossprod( personMaxA , pweights )
	
	
    # maximum score in resp, equal categories?  
    maxscore.resp <- apply( resp , 2 , max , na.rm=TRUE)
    if ( ncol(resp)>1){ 
      sd.maxscore.resp <- stats::sd(maxscore.resp)
    } else { sd.maxscore.resp <- 0 }

    
    equal.categ <- if( sd.maxscore.resp > .00001 ){ FALSE } else { TRUE  }
    #  xsi[est.xsi.index] <- - log(abs(ItemScore[est.xsi.index]/(ItemMax[est.xsi.index]-
    #      ItemScore[est.xsi.index])))  #log of odds ratio of raw scores  
    xsi[est.xsi.index] <- - log(abs(( ItemScore[est.xsi.index]+.5)/
                                      (ItemMax[est.xsi.index]-ItemScore[est.xsi.index]+.5) ) )
    # starting values of zero
    if( xsi.start0 ){ xsi <- 0*xsi }				
    if ( ! is.null(xsi.inits) ){   
      #		xsi <- xsi.inits  
      xsi[ xsi.inits[,1] ] <- xsi.inits[,2]
    }
    if ( ! is.null( xsi.fixed ) ){   xsi[ xsi.fixed[,1] ] <- xsi.fixed[,2] }

 
    xsi.min.deviance <- xsi
    beta.min.deviance <- beta
    variance.min.deviance <- variance
    # nodes
    if ( snodes == 0 ){
     if ( skillspace=="normal"){	
      theta <- as.matrix( expand.grid( as.data.frame( matrix( rep(nodes, ndim) , ncol = ndim ) ) ) )
								}
	  if ( ( skillspace != "normal") & ( ! is.null(theta.k) ) ){	  
				theta <- as.matrix( theta.k )
				nnodes <- nrow(theta)
								}		
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
	  theta <- theta0.samp
    }
    deviance <- 0  
    deviance.history <- matrix( 0 , nrow=maxiter , ncol = 2)
    colnames(deviance.history) <- c("iter" , "deviance")
    deviance.history[,1] <- 1:maxiter
    iter <- 0 
    a02 <- a1 <- 999	# item parameter change
    #---2PL---
    # a4 <- 999   
    basispar <- NULL
    #		} else{  a4 <- 0 }
    
    
	#*************************
	# skill space
	if ( ! is.null(theta ) ){ 
			ntheta <- nrow(theta)
					} 															
	fulldesign <- FALSE
	if ( skillspace != "normal" ){	    
		gwt <- hwt <- matrix( 1/ntheta , nrow=nstud , ncol=ntheta)				 						
		covdelta <- group1.list <- list(1:G)
		Ngroup <- rep(0,G)
		for (gg in 1:G){ 
		    ind.gg <- which( group1 == gg )
			Ngroup[gg] <- sum( pweights[ind.gg] )
			group1.list[[gg]] <- ind.gg
						}
		pi.k <- matrix( 1/ntheta , nrow=ntheta , ncol=G)
			
		if ( ! is.null( delta.designmatrix) ){
			if ( ncol(delta.designmatrix) == ntheta ){
				fulldesign <- TRUE 
							}		
						}
		
		# design matrix		
		if ( is.null( delta.designmatrix) ){
		    delta.designmatrix <- diag( ntheta )
			fulldesign <- TRUE
										}
        delta <- matrix( 0 , nrow= ncol(delta.designmatrix) , ncol=G)	
				}
	gwt1 <- matrix( 1 , nrow=nstud , ncol=ntheta )	

	#***** inits for delta
	if ( ! is.null(delta.inits) ){
		delta <- delta.inits 
				}
	
	#****** indicator matrices
	datindw <- list(1:maxK)
	for (kk in 1:maxK){ 
		datindw[[kk]] <- (resp == kk - 1 ) * resp.ind * pweights
						}
	
	#*************
	# gammaslope constraints
	if ( ! is.null(gammaslope.constr.V) ){
			V <- gammaslope.constr.V
			e2 <- matrix( gammaslope.constr.c , nrow=ncol(V) , ncol=1 )
			V1 <- solve( crossprod(V) )
										}	
	
  	 gammaslope <- .mml.3pl.gammaslope.center( gammaslope , gammaslope.center.index  ,
	         			gammaslope.center.value  )

	
	#******
	# prior distribution guessing parameter
	if ( ! is.null(guess.prior) ){
		guess.mean <- guess.prior[,1] / rowSums( guess.prior )
		i1 <- which( guess.prior[,1] > 0 )
		guess[ i1 ] <- guess.mean[i1]
		guess.prior[ guess.prior == 0 ] <- .001
								}
	#******
	# prior distribution slope parameter
	if ( ( ! is.null(gammaslope.prior) ) & ( ! init.gammaslope) ){		
		i1 <- which( gammaslope.prior[,2] < 10 )
		gammaslope[ i1 ] <- gammaslope.prior[i1,1]	
										}
			
	#******
	# prior distribution slope parameter
	if ( ! is.null(xsi.prior) ){		
		i1 <- which( xsi.prior[,2] < 10 )
		xsi[ i1 ] <- xsi.prior[i1,1]
										}
										
	  #******
	  # compute F design matrix for loadings

	  Fdes <- .mml.3pl.computeFdes( E , gammaslope , theta )	  
      # use simplified design for F
	  dimFdes <- dim(Fdes)
 
	  res <- .Call("mml3_calc_Fdes" , as.vector(Fdes) , dimFdes=dimFdes ,
				       PACKAGE="TAM" )	  
	  FdesM <- res$FdesM[ 1:res$NFdesM , ]
	
	#*****
	#@@@@ 2015-06-26
#	Avector <- as.vector(A)
#	Avector[ is.na(Avector) ] <- 0
	#@@@@
	
    ##**SE
    se.xsi <- 0*xsi
    se.B <- 0*B 
    
    YSD <- max( apply( Y , 2 , stats::sd ) )
    if (YSD > 10^(-15) ){ YSD <- TRUE } else { YSD <- FALSE }

	#*****
	#@@@@ speed gains, further auxiliary objects, 2015-06-26
	Avector <- as.vector(A)
	Avector[ is.na(Avector) ] <- 0
	unidim_simplify <- TRUE
	if (G > 1){ unidim_simplify <- FALSE }
	if ( YSD){ unidim_simplify <- FALSE }	
	if (  is.null(beta.fixed) ){ unidim_simplify <- FALSE }
	#@@@@

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
	gammaslope_acceleration <- list( "acceleration" = acceleration , "w" = .35 ,
							"w_max" = .95 , 
							parm_history = cbind( gammaslope, gammaslope , gammaslope) ,
							"beta_new" = 0 ,
							"beta_old" = 0 
						)
	d1 <- as.vector(delta)	
	#delta_acceleration <- list( "acceleration" = "none" , "w" = .35 ,
	delta_acceleration <- list( "acceleration" = acceleration , "w" = .35 ,
							"w_max" = .95 , 
							parm_history = cbind( d1, d1 , d1) ,
							"beta_new" = 0 ,
							"beta_old" = 0 
						)
    #v1 <- qlogis( guess + 1E-4 )
	v1 <- guess
    ind.guess <- which( est.guess > 0 )	
	guess_acceleration <- list( "acceleration" = acceleration , "w" = .35 ,
							"w_max" = .95 , 
							parm_history = cbind( v1, v1 , v1 ) ,
							"beta_new" = 0 ,
							"beta_old" = 0 ,
							"ind_guess" = ind.guess
						)							
	#@@@@AAAA@@@@@		
	
	
	
    
    hwt.min <- 0
    rprobs.min <- 0
    AXsi.min <- 0
    B.min <- 0
    deviance.min <- 1E100
    itemwt.min <- 0
    se.xsi.min <- se.xsi
    se.B.min <- se.B
	gammaslope.change <- guess.change <- 0
	a4 <- 0
    a44 <- 1000
	old.increment.guess <- .6
	se.gammaslope <- NULL
	
	se.guess <- 0*guess
	
	
    # display
    disp <- "....................................................\n"
    # define progress bar for M step
    mpr <- round( seq( 1 , np , len = 10 ) )

	if ( ! is.null(gammaslope.fixed ) ){
		gammaslope.fixed <- as.matrix( gammaslope.fixed )
								}
	# inits delta parameters
    if ( ! is.null( delta.inits) ){
		for ( gg in 1:G){
			pi.k[,gg] <- exp( delta.designmatrix %*% delta.inits[,gg] )
						}
					}

	#*****
	#@@@@ speed gains, further auxiliary objects, 2015-06-26
	# Avector <- as.vector(A)
	# Avector[ is.na(Avector) ] <- 0
	unidim_simplify <- TRUE
	if (G > 1){ unidim_simplify <- FALSE }
	if ( YSD){ unidim_simplify <- FALSE }
	#@@@@
							

							
	##############################################################
	##############################################################
    ##############################################################   
    # EM loop starts here
    while ( ( (!betaConv | !varConv)  | ((a1 > conv) | (a44 > conv) | (a02 > convD)) )  & (iter < maxiter) ) { 
a0 <- Sys.time()
	  delta0 <- delta
	  
      iter <- iter + 1
      if (progress){ 
        cat(disp)	
        cat("Iteration" , iter , "   " , paste( Sys.time() ) )
        cat("\nE Step\n") ; flush.console()
      }
      # calculate nodes for Monte Carlo integration	
      if ( snodes > 0){
        #      theta <- beta[ rep(1,snodes) , ] +  
        #				t ( t(chol(variance)) %*% t(theta0.samp) )
        theta <- beta[ rep(1,snodes) , ] + theta0.samp %*% chol(variance) 
        
        # calculate density for all nodes
        thetasamp.density <- mvtnorm::dmvnorm( theta , mean = as.vector(beta[1,]) , sigma = variance )
        # recalculate theta^2
        #      theta2 <- matrix( theta.sq(theta) , nrow=nrow(theta) , ncol=ncol(theta)^2 )   
        theta2 <- matrix( theta.sq2(theta) , nrow=nrow(theta) , ncol=ncol(theta)^2 )   
      }			
      olddeviance <- deviance
	  
      # calculation of probabilities
      res <- .mml.3pl.calc_prob.v5(iIndex=1:nitems , A=A , AXsi=AXsi , B=B , xsi=xsi , theta=theta , 
                          nnodes=nnodes , maxK=maxK , recalc=TRUE ,
						  guess=guess )	
      rprobs <- res[["rprobs"]]
	  rprobs0 <- res$rprobs0
      AXsi <- res[["AXsi"]]
 
 
# cat("calc_prob") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1

    

	  #***********************************
	  # student's prior distribution
	  if (skillspace == "normal" ){	  	 	# normal distribution
		  gwt <- stud_prior.v2(theta=theta , Y=Y , beta=beta , 
					variance=variance , nstud=nstud , 
					nnodes=nnodes , ndim=ndim,YSD=YSD , unidim_simplify=unidim_simplify , 
					snodes = snodes )
						   							   
		 if ( snodes == 0 ){
			gwt <- gwt / rowSums( gwt ) 
							}
							
									}

      if ( skillspace != "normal" ){		# non-normal distribution			
		  for (gg in 1:G){
			ind.gg <- group1.list[[gg]]
#			calcpik <- ( is.null(delta.designmatrix) ) | ( iter <= 3) 
			calcpik <- FALSE
#			calcpik <- TRUE
			if ( calcpik ){
				pi.k[,gg] <- colSums( ( pweights*hwt )[ind.gg,] )
				pi.k[,gg] <- pi.k[,gg] / sum( pi.k[,gg] )
							}
			gwt[ind.gg,] <- matrix( pi.k[,gg] , nrow=length(ind.gg) , 
					ncol=ntheta , byrow=TRUE )						
                					  }
							}
									
#	 cat("stud prior") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1				
	  
	  
	  #******************************************
	  # likelihood	and posterior  
				
      # calculate student's likelihood	
	   gwt1 <- gwt
	   res.hwt <- calc_posterior.v2(rprobs=rprobs , gwt=gwt1 , resp=resp , nitems=nitems , 
	                               resp.ind.list=resp.ind.list , normalization=FALSE , 
                                   thetasamp.density=thetasamp.density , snodes=snodes ,
                                   resp.ind=resp.ind , logprobs=TRUE , avoid.zerosum=TRUE )
#	   hwt0 <- hwt <- res.hwt$hwt * gwt   
	   hwt0 <- hwt <- res.hwt[["hwt"]]  	   
       hwt <- hwt / rowSums( hwt )	   	   
	   
# cat("posterior v2") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1
      
      if (progress){ cat("M Step Intercepts   |"); flush.console() }
      # collect old values for convergence indication

      oldxsi <- xsi
      oldbeta <- beta
      oldvariance <- variance 
     
	    #******************************************
        # M step: distribution parameter estimation of beta and variance
	    if ( skillspace == "normal" ){
		  resr <- mstep.regression( resp=resp , hwt=hwt , resp.ind=resp.ind , pweights=pweights , 
									pweightsM=pweightsM , Y=Y , theta=theta , theta2=theta2 , YYinv=YYinv , 
									ndim=ndim , nstud=nstud , beta.fixed=beta.fixed , variance=variance , 
									Variance.fixed=variance.fixed , group=group ,  G=G , snodes = snodes ,
									thetasamp.density=thetasamp.density,nomiss=nomiss)		  
		  beta <- resr$beta
		  # cat("m step regression") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1		  
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
		  if (max(abs(variance-oldvariance)) < conv) varConv <- TRUE	  
			}  # end skillspace == "normal"
		#******************************************
		
		
		#******************************************
		# skill space estimation non-normal distribution	
		if ( skillspace != "normal" ){
			itemwt <- crossprod( hwt , resp.ind * pweightsM  )			
			
		  for (gg in 1:G){
			ind.gg <- group1.list[[gg]]
			pi.k[,gg] <- colSums( ( pweights*hwt )[ind.gg,] )
			pi.k[,gg] <- pi.k[,gg] / sum( pi.k[,gg] )
							}
						
			# log-linear smoothing of skill space
			res <- .mml.3pl.skillspace( Ngroup, pi.k , 
						delta.designmatrix , G , delta , delta.fixed )
			pi.k <- res$pi.k
			delta <- res$delta	
			covdelta <- res$covdelta
			
        #@@@@AAAA@@@@@
		# acceleration
		if ( delta_acceleration$acceleration != "none" ){		
			delta_acceleration <- accelerate_parameters( xsi_acceleration=delta_acceleration , 
							xsi=delta , iter=iter , itermin=3)
			delta <- matrix( delta_acceleration$parm , nrow=nrow(delta) , ncol=ncol(delta) )
								}
	    #@@@@AAAA@@@@@	  			
			
			varConv <- TRUE
			betaConv <- TRUE			
									}
      
		if (skillspace == "normal"){
		   if (G==1){  diag(variance) <- diag(variance)+10^(-14)	 }
			if ( ! est.variance ){ 
			  if ( G == 1 ){ variance <- cov2cor(variance)  } # fix variance at 1  
			  if ( G > 1 ){ variance[ group == 1 ] <- 1 }     
			# fix variance of first group to 1
							}
						}
#      }
      
      #---end 2PL---
      # stop("er")
# cat("sufficient statistics 2PL") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1	

	#******
	# generate input for fixed parameters
	xsi.fixed.estimated <- generate.xsi.fixed.estimated( xsi , A )
	B.fixed.estimated <- generate.B.fixed.estimated(B)
      	  
      
	  ######################################
	  # calculation of expected counts
	  
		res <- .mml.3pl.expected.counts( datindw , nitems , maxK , ntheta , hwt)
		n.ik <- res$n.ik
		N.ik <- res$N.ik
		
# cat("expected counts") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1	
		
      ######################################
      # M-step item intercepts
	  res <- .mml.3pl.est.intercepts( max.increment , np , est.xsi.index0 , 
			Msteps , nitems , A , AXsi , B , xsi , guess , theta , nnodes , maxK ,
			progress , itemwt , indexIP.no , indexIP.list2 ,
			ItemScore , fac.oldxsi , rprobs , xsi.fixed , convM , rprobs0 ,
			n.ik , N.ik , xsi.prior )
	  xsi <- res$xsi
	  se.xsi <- res$se.xsi
      
        #@@@@AAAA@@@@@
		# acceleration
		if ( xsi_acceleration$acceleration != "none" ){		
			xsi_acceleration <- accelerate_parameters( xsi_acceleration=xsi_acceleration , 
							xsi=xsi , iter=iter , itermin=3)
			xsi <- xsi_acceleration$parm
								}
	    #@@@@AAAA@@@@@	  
	  
# cat("M steps intercepts") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1		

      ###############################################
	  # M-step item slopes	  
      if ( est.some.slopes){
          oldgamma <- gammaslope		  
		  res <- .mml.3pl.est.slopes( max.increment , np , 
				Msteps , nitems , A , AXsi , B , xsi , guess , theta , nnodes , maxK ,
				progress ,ItemScore , fac.oldxsi , rprobs , xsi.fixed , convM , rprobs0 ,
				n.ik , N.ik , gammaslope , E , FdesM , dimFdes ,
				gammaslope.fixed , gammaslope.prior , maxgamma = maxgamma , Edes )	
				
		  gammaslope <- res$gammaslope	
		  se.gammaslope <- res$se.gammaslope
		  gammaslope.change <- res$gammachange	
		  
		 if ( ! is.null(gammaslope.constr.V) ){
            e1 <- matrix( gammaslope , ncol=1 )
			gammaslope <- ( e1 + V %*% V1 %*% ( e2 - t(V) %*% e1 ) )[,1]		
									}

		gammaslope <- .mml.3pl.gammaslope.center( gammaslope , gammaslope.center.index  ,
				gammaslope.center.value  )
									
		  # userfunction gammaslope
		  if ( ! is.null( userfct.gammaslope ) ){
				gammaslope <- do.call( userfct.gammaslope , list(gammaslope) )
#				B <- .mml.3pl.computeB( E , gammaslope )		
					}	  		  
		  gammaslope <- fac.oldxsi	* oldgamma + ( 1 - fac.oldxsi)*gammaslope		

        #@@@@AAAA@@@@@
		# acceleration
		if ( gammaslope_acceleration$acceleration != "none" ){		
			gammaslope_acceleration <- accelerate_parameters( xsi_acceleration=gammaslope_acceleration , 
							xsi=gammaslope , iter=iter , itermin=3)
			gammaslope <- gammaslope_acceleration$parm
								}
	    #@@@@AAAA@@@@@			  
		  
		  # B <- .mml.3pl.computeB( E , gammaslope )	  
		  B <- .mml.3pl.computeB.v2( Edes , gammaslope , E )		  
		  
				}
					
	  #*********************
	  # 3PL estimation
	  if ( est.some.guess ){
	      oldguess <- guess
		  res <- .mml.3pl.est.guessing( guess , Msteps , convM , 
						nitems , A , AXsi , B, xsi , theta , nnodes , maxK ,
						n.ik , N.ik , est.guess ,  old.increment.guess ,
						guess.prior  , progress	)	  
		  guess <- res$guess
		  guess.change <- res$guess.change
		  se.guess <- res$se.guess
		  
        #@@@@AAAA@@@@@
		# acceleration
		if ( guess_acceleration$acceleration != "none" ){		
			# g1 <- qlogis( guess + 1E-5 )
			g1 <- guess
			guess_acceleration <- accelerate_parameters( xsi_acceleration=guess_acceleration , 
							xsi= g1 , iter=iter , itermin=3 , 
							ind = guess_acceleration$ind_guess)
			#guess <- plogis( guess_acceleration$parm )
			guess <- guess_acceleration$parm
			guess[ guess < 0 ] <- 1E-5
			guess.change <- max( abs(guess - oldguess))
								}
	    #@@@@AAAA@@@@@	
		
		  
#		  old.increment.guess <- guess.change
					}
# cat("M steps slopes") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1		
      
      #***
      # decrease increments in every iteration
      if( increment.factor > 1){max.increment <-  1 / increment.factor^iter }

	    
      # calculate deviance
      if ( snodes == 0 ){ 
#        deviance <- - 2 * sum( pweights * log( res.hwt$rfx * thetawidth ) )
		 # rfx <- rowSums( hwt * gwt )
		 rfx <- rowSums( hwt0 )
		 deviance <- - 2 * sum( pweights * log( rfx ) )
      } else {
        #      deviance <- - 2 * sum( pweights * log( res.hwt$rfx ) )
#        deviance <- - 2 * sum( pweights * log( rowMeans( res.hwt$swt ) ) )
# Revalpr("sum(res.hwt$rfx)")		
       #	   deviance <- - 2 * sum( pweights * log( rowMeans( res.hwt$hwt ) ) )
        deviance <- - 2 * sum( pweights * log( res.hwt$rfx   ) )		
		
      }
      deviance.history[iter,2] <- deviance
      a01 <- abs( ( deviance - olddeviance ) / deviance  )
      a02 <- abs( ( deviance - olddeviance )  )	
	  if (con$dev_crit == "relative" ){ a02 <- a01 }
 
      if( ( deviance < deviance.min ) | ( iter == 1)  ){ 
        xsi.min.deviance <- xsi 
        beta.min.deviance <- beta
        variance.min.deviance <- variance	
        hwt.min <- hwt	
        AXsi.min <- AXsi	
        B.min <- B
		gammaslope.min <- gammaslope
		se.gammaslope.min <- se.gammaslope
        deviance.min <- deviance
		delta.min <- delta
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
        if ( est.some.slopes ){				
			cat( "\n  Maximum slope parameter change:" , round( gammaslope.change , 6 ) )
							}										
        if ( est.some.guess ){				
			cat( "\n  Maximum guessing parameter change:" , round( guess.change , 6 ) )
							}
		a44 <- max( a4 , gammaslope.change , guess.change )
		
		#***********************		
		#**** skill space == "normal"		
		if ( skillspace == "normal"){		
			cat( "\n  Maximum regression parameter change:" , round( a2 , 6 ) )  
			if ( G == 1 ){ 
			  cat( "\n  Variance: " , round( variance[ ! lower.tri(variance)] , 4 ) , " | Maximum change:" , round( a3 , 6 ) )  
			} else {
			  cat( "\n  Variance: " , round( variance[var.indices] , 4 ) ,
				   " | Maximum change:" , round( a3 , 6 ) )  		
			}					
			cat( "\n  beta ",round(beta,4)  )
					}
		if ( skillspace != "normal" ){		
              a31 <- max( abs( delta - delta0 ))		
			  cat( "\n  Maximum delta parameter change: " , round( a31 , 6 ) )  						
					}
		cat( "\n" )					
        flush.console()
		
				
      }
    } # end of EM algorithm loop
    ############################################################################
    ############################################################################
    ############################################################################
# stop("check here")  

    xsi.min.deviance -> xsi 
    beta.min.deviance -> beta
    variance.min.deviance -> variance	
    hwt.min -> hwt	
    AXsi.min -> AXsi	
    B.min -> B
	gammaslope.min -> gammaslope
	se.gammaslope.min -> se.gammaslope
	delta.min -> delta
    deviance.min -> deviance
    itemwt.min -> itemwt
    se.xsi.min -> se.xsi	
    se.B.min -> se.B		
	
# Revalpr("deviance")	
    #******
    ##**SE  
    # standard errors of AXsi parameters
    # check for missing entries in A
    se.AXsi <- 0*AXsi
    A1 <- A
    A1[ is.na(A) ] <- 0
    #	se.xsiD <- diag( se.xsi^2 )
    if ( length( se.xsi) > 1){
      se.xsiD <- diag( se.xsi^2 )
    } else {
      se.xsiD <- matrix( se.xsi^2,1,1)
    }	
    
    for (kk in 1:maxK){  # kk <- 1
      # se.AXsi[,kk] <- sqrt( diag( A1[,kk,] %*% se.xsiD %*% t( A1[,kk,]) ) )
      #**** bugfix
	  dim_A1 <- dim(A1)
      A1_kk <- A1[,kk,]
      if ( is.vector(A1_kk) ){
        A1_kk <- matrix( A1_kk , nrow=dim_A1[1] , ncol=dim_A1[3] )
      }
      se.AXsi[,kk] <- sqrt( diag( A1_kk %*% se.xsiD %*% t( A1_kk ) ) )	
      #****		
    } 
    ##*** Information criteria
    ic <- .mml.3pl.TAM.ic( nstud=nstud100 , deviance , xsi , xsi.fixed ,
                   beta , beta.fixed , ndim , variance.fixed , G ,
                   irtmodel ,B_orig=B_orig ,  B.fixed , E , est.variance , resp ,
                   est.slopegroups , skillspace , delta , delta.fixed , est.guess ,
				   fulldesign , est.some.slopes , gammaslope ,
				   gammaslope.fixed , gammaslope.constr.V , gammaslope.constr.Npars ,
				   gammaslope.center.index  , gammaslope.prior , 
				   numdiff.parm=5*1E-4 )
# cat("q200\n")
    #***
    # calculate counts
    res <- .tam.calc.counts( resp, theta , resp.ind , 
                             group , maxK , pweights , hwt )
    if ( skillspace == "normal"){							 
		n.ik <- res$n.ik
		pi.k <- res$pi.k   
				}
    #****
    # collect item parameters
    item1 <- .mml.3pl.TAM.itempartable( resp , maxK , AXsi , B , ndim ,
                                resp.ind , rprobs,n.ik,pi.k , guess , est.guess )
# cat("q300\n")
								
	# distribution moments
	if ( skillspace != "normal" ){
		D <- ncol(theta.k)
		moments <- .mml.3pl.distributionmoments( D =D , G =G , pi.k=pi.k , theta.k=theta.k )
								} else { moments <- NULL }
								
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
      #		hwtE <- hwt / snodes 
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
	  #*****
	  # skillspace == normal
	  if (skillspace=="normal"){	
		  cat("...................................\n")
		  cat("Regression Coefficients\n")
		  print( beta , 4  )
		  cat("\nVariance:\n" ) # , round( varianceM , 4 ))
		  if ( ( G==1 ) ){ 
			varianceM <- matrix( variance , nrow=ndim , ncol=ndim ) 
			print( varianceM , 4 )	
		  } else { 
			print( variance[ var.indices] , 4 )	}
		  if ( ndim > 1){ 
			cat("\nCorrelation Matrix:\n" ) # , round( varianceM , 4 ))	
			print( cov2cor(varianceM) , 4 )	
				}
		  }
	  #****	  
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
    
	# labels gammaslope parameters
	names(gammaslope) <- dimnames(E)[[4]]
	
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
				 "moments" = moments ,
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
                 "hwt" = hwt , "like"= res.like , "ndim" = ndim ,
                 "xsi.fixed" = xsi.fixed , 
				 "xsi.fixed.estimated" = xsi.fixed.estimated , 				 
				 "beta.fixed" = beta.fixed , "Q" = Q, 
                 "B.fixed" = B.fixed , 
				 "B.fixed.estimated" = B.fixed.estimated , 				 
				 "est.slopegroups" = est.slopegroups , "E" = E , "basispar" = basispar,
                 "variance.fixed" = variance.fixed ,
                 "nnodes" = nnodes , "deviance" = deviance ,
                 "ic" = ic , 
                 "deviance.history" = deviance.history ,
                 "control" = con1a , "irtmodel" = irtmodel ,
                 "iter" = iter ,
                 "printxsi"=printxsi 	, "YSD"=YSD		,
				 "skillspace"= skillspace ,
				 "delta" = delta , "delta.designmatrix" = delta.designmatrix , 
				 "gammaslope" = gammaslope , "se.gammaslope" = se.gammaslope ,
				 "guess" = guess ,  "se.guess" = se.guess , 
				 "E"= E , "Edes" = Edes , CALL = CALL 
                 #			   "design"=design				
                 #			   "xsi.min.deviance" = xsi.min.deviance ,
                 #			   "beta.min.deviance" = beta.min.deviance , 
                 # "variance.min.deviance" = variance.min.deviance 
    )
    class(res) <- "tam.mml.3pl"
    return(res)
  }

