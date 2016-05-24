tam.mml.mfr <-
  function( resp , Y=NULL , group = NULL ,  irtmodel ="1PL" ,
            formulaY = NULL , dataY = NULL , 
            ndim = 1 , pid = NULL ,
            xsi.fixed=NULL ,  xsi.setnull = NULL , 
			xsi.inits = NULL , 			
            beta.fixed = NULL , beta.inits = NULL , 
            variance.fixed = NULL , variance.inits = NULL , 
            est.variance = FALSE , formulaA=~item+item:step, constraint="cases",
            A=NULL , B=NULL , B.fixed = NULL , 
            Q=NULL , facets=NULL, est.slopegroups=NULL , E = NULL , 
            pweights = NULL , control = list() ,
            delete.red.items=TRUE
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
    
	CALL <- match.call()
    a0 <- Sys.time()    
    s1 <- Sys.time()
    # display
    disp <- "....................................................\n"
    increment.factor <- progress <- nodes <- snodes <- ridge <- xsi.start0 <- QMC <- NULL
    maxiter <- conv <- convD <- min.variance <- max.increment <- Msteps <- convM <- NULL 
    resp_orig <- resp
    B00 <- B 
	B <- NULL

    # attach control elements
    e1 <- environment()
    con <- list( nodes = seq(-6,6,len=21) , snodes = 0 , QMC=TRUE , 
                 convD = .001 ,conv = .0001 , convM = .0001 , Msteps = 4 ,            
                 maxiter = 1000 , max.increment = 1 , 
                 min.variance = .001 , progress = TRUE , ridge=0,seed= NULL ,
                 xsi.start0= 0 , increment.factor=1 , fac.oldxsi=0 , acceleration="none" ,
				 dev_crit = "absolute"  )  	
    con[ names(control) ] <- control  
    Lcon <- length(con)
    con1a <- con1 <- con ; 
    names(con1) <- NULL

    for (cc in 1:Lcon ){
      assign( names(con)[cc] , con1[[cc]] , envir = e1 ) 
    }
    if ( !is.null(con$seed)){ set.seed( con$seed )	 }
	
	
	acceleration <- con$acceleration

    #***
    fac.oldxsi <- max( 0 , min( c( fac.oldxsi , .95 ) ) )
    
    if ( constraint=="items" ){ beta.fixed <- FALSE }
    
    pid0 <- pid <- unname(c(unlist(pid)))
    if (progress){ 
      cat(disp)	
      cat("Processing Data     ", paste(Sys.time()) , "\n") ; utils::flush.console()
    }    
    if ( ! is.null(group) ){ 
      con1a$QMC <- QMC <- FALSE
      con1a$snodes <- snodes <- 0
    }
    
	resp <- as.matrix(resp)
	resp <- add.colnames.resp(resp)		
    itemnames <- colnames(resp)
	
    nullY <- is.null(Y)
    
    if ( ! is.null(facets) ){ facets <- as.data.frame(facets) }
    
# cat("read data" ) ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1
    # create design matrices
    if(ncol(resp)>1){
      maxKi <- apply( resp , 2 , max , na.rm=TRUE )
    } else {
		if(ncol(resp)==1){
		  item.ind <- grep("Item", names(facets), ignore.case=TRUE)
		  if(!is.null(item.ind)){
				if ( length(item.ind) == 0 ){
					item.ind <- NULL 
								}
							}		   
		  if(!is.null(item.ind)){ 
			maxKi <- stats::aggregate( resp , facets[,item.ind,drop=FALSE] , 
								max, na.rm=TRUE )[,2]
		  }else{
			maxKi <- stats::aggregate( resp , facets[,1,drop=FALSE] , 
								max, na.rm=TRUE )[,2]
			
		  }
		}		
	}
# cat("w200\n")
# Revalpr("maxKi")
# Revalpr("head(resp)")
	
	#*****************
	# handle formula and facets
    resp00 <- resp

	res <- mfr.dataprep( formulaA , xsi.setnull , B , Q ,
				resp, pid, facets , beta.fixed  )				
	formulaA <- res$formula_update
	xsi.setnull <- res$xsi.setnull	
	beta.fixed <- res$beta.fixed
	facets <- res$facets
	PSF <- res$PSF
	pid <- res$pid

# Revalpr("head(resp00)")
# cat("s100\n")

# maxKi <- 10


	
	diffKi <- FALSE
# Revalpr("maxKi")	



	var_ki <- var( maxKi )
	if ( is.na( var_ki) ){ var_ki <- 0 }

# Revalpr("var_ki")	
	
	if ( var_ki > .001 ){ 
	     diffKi <- TRUE
# cat("heere1\n")		 
		 design <- designMatrices.mfr2(resp, formulaA=formulaA, facets=facets,  
                                 constraint=constraint, ndim=ndim,
                                 Q=Q, A=A, B=B , progress=progress)
		 xsi.elim <- design$xsi.elim
		if ( ! is.null(xsi.elim) ){	
			 if ( nrow(xsi.elim) > 0 ){
				 xsi.elim2 <- cbind( xsi.elim[,2] , 99 )		 
				 xsi.fixed <- rbind( xsi.fixed , xsi.elim2 )
										}
									}
			# set first beta coefficient to zero
			if ( is.null( beta.fixed ) ){
				dimB <- dim(design$B$B.3d.0	)	
			    beta.fixed <- cbind( 1 , 1:dimB[3] , 0)
							}
					} else {				
         design <- designMatrices.mfr(resp, formulaA=formulaA, facets=facets,  
                                 constraint=constraint, ndim=ndim,
                                 Q=Q, A=A, B=B , progress=progress)	
								 
							}	 
																										
													
    A <- design$A$A.3d.0	
    cA <- design$A$A.flat.0	
    B <- design$B$B.3d.0
    Q <- design$Q$Q.flat.0
    X <- design$X$X
    X.red <- design$X$X.noStep
    gresp <- design$gresp$gresp
    gresp.noStep <- design$gresp$gresp.noStep
    xsi.constr <- design$xsi.constr
	
	#****************************
	items00 <- colnames(resp00)
	I00 <- length(items00)
	D <- dim(B00)[3]
	if ( ! is.null(B00) ){		
			rownames_A <- dimnames(A)[[1]]
			for (ii in 1:I00){
				#		ii <- 1
				ind <- grep( paste0( items00[ii] , "-" ) , rownames_A )
				if ( length(ind) > 0 ){
					I2 <- length(ind)
					for (vv in 1:I2){
						B[ ind[vv] , , 1:D] <- B00[ii , , 1:D ,drop=FALSE] * ( B[ ind[vv] , , 1:D] != 0  )
								}						
							}		
						}									
					}   # end is.null()
			#*****************

	
    if ( is.null( pid ) ){ pid <- 1:(nrow(gresp) ) }
    
	
#Revalpr(" head(gresp.noStep)")	
	
    design <- NULL
    
    if (progress){ 
      cat("    * Created Design Matrices   (", 
          paste(Sys.time()) , ")\n") ; utils::flush.console()	  
    }    
    
#    cat(" --- design matrix ready" ) ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1
    
    #***
    # preprocess data if multiple person IDs do exist
    tp <- max( table( pid ))

    if ( tp > 1){
      persons <- sort( unique( pid ) )
      NP <- length( persons )
      # cat("*** multiple persons start" ) ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1		
      #****
      # ARb 2013-08-23: added simplify=TRUE
      person.ids <- sapply( persons , FUN = function( pp){ which( pid == pp ) } ,
                            simplify=FALSE)
      #print(person.ids[[5]] )					
      #cat("*** multiple persons sapply function" ) ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1				
      PP <- matrix( NA , nrow=NP , ncol=tp)
      for (pos in 1:tp){
        #pos <- 1
        PP[,pos] <- unlist( lapply( person.ids , FUN = function( vv){ vv[pos] } ) )
      }

     
#     cat("*** multiple persons lapply function" ) ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1
      gresp0 <- matrix( NA , nrow=NP , ncol= ncol(gresp) )
      colnames(gresp0) <- colnames(gresp)
      gresp0.noStep <- matrix( NA , nrow=NP , ncol= ncol(gresp.noStep) )
      colnames(gresp0.noStep) <- colnames(gresp.noStep)
      grespNA <- ( ! is.na( gresp ) )
      grespnoStepNA <- ( ! is.na( gresp.noStep ) )		
      #cat("*** multiple persons start pos" ) ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1			
	  
	  #***
	  # check multiple rows
	  m1 <- rowsum( 1-is.na(gresp.noStep) , pid )
	  	  
# Revalpr("table( rowsum( m1 > 1 )[,1] )")	  
	  
	  h1 <- sum(m1>1)
	  if (h1>0){
		cat("* Combinations of person identifiers and facets are not unique.\n")
		cat("* Use an extended 'formulaA' to include all \n")
		cat("  relevant facets and the argument 'xsi.setnull'.\n")
		cat("  See the help page of 'tam.mml' (?tam.mml) Example 10a.\n") 
		stop()			
				}
	  	  
      for (pos in 1:tp){
        ind.pos <- which( ! is.na( PP[,pos]  ) )
        PP.pos <- PP[ind.pos,pos]
        g1 <- gresp[ PP.pos , ]
        g0 <- gresp0[ ind.pos , ]
        #			ig1 <- ( ! is.na(g1) )
        ig1 <- grespNA[ PP.pos , ]
        # * this check is time-consuming! release it to rcpp
        g0[ ig1 ] <- g1[ ig1 ]
        gresp0[ ind.pos , ] <- g0
        g1 <- gresp.noStep[ PP.pos , ]
        g0 <- gresp0.noStep[ ind.pos , ]
        #			ig1 <- ( ! is.na(g1) )
        ig1 <- grespnoStepNA[ PP.pos , ]
        g0[ ig1 ] <- g1[ ig1 ]
        gresp0.noStep[ ind.pos , ] <- g0

	
      }
      #cat("*** multiple persons loop over pos" ) ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1				
      gresp0 -> gresp
      gresp0.noStep -> gresp.noStep
      pid <- persons	
      if (progress){ 
        cat("    * Arranged Response Data with Multiple Person Rows   (", 
            paste(Sys.time()) , ")\n") ; utils::flush.console()	  
      }  		
    }
    ###################################################
    
	# set some xsi effects to null
	if ( ! is.null(xsi.setnull) ){	
		xsi.labels <- dimnames(A)[[3]]
		xsi0 <- NULL	
		N1 <- length(xsi.setnull)
		for (nn in 1:N1){
			ind.nn <- grep( xsi.setnull[nn] , xsi.labels )	
			l1 <- cbind( ind.nn , 0 )
			xsi0 <- rbind( xsi0 , l1 )	
			colnames(xsi0) <- NULL
						}
		xsi.fixed <- rbind( xsi.fixed , xsi0 )
		i2 <- duplicated(xsi.fixed[,1])
		if ( sum(i2) > 0 ){
			xsi.fixed <- xsi.fixed[ - i2  , ]
								}							
						}

#cat("process data in case of multiple persons" ) ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1


    nitems <- nrow( X.red )
    nstud <- nrow(gresp)        # number of students
    if ( is.null( pweights) ){
      pweights <- rep(1,nstud) # weights of response pattern
    }
    
    if (progress){ 
      cat("    * Response Data:" , nstud , "Persons and " , 
          ncol(gresp.noStep) , "Generalized Items (" , paste(Sys.time()) ,")\n" )  ;
      utils::flush.console()	  
    }  	
    #!! check dim of person ID pid
    if ( is.null(pid) ){ pid <- seq(1,nstud) }
    
    # normalize person weights to sum up to nstud
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
      #   cat( paste0( l1 , nnodes , " nodes\n") )
      if (nnodes > 8000){
        cat("      @ Are you sure that you want so many nodes?\n")
        cat("      @ Maybe you want to use Quasi Monte Carlo integration with fewer nodes.\n")		
      }
    }
    #********* 
    # maximum no. of categories per item. Assuming dichotomous
    maxK <- max( resp , na.rm=TRUE ) + 1 
    
    ################################
    # number of parameters
    np <- dim(A)[[3]]
    
    # xsi inits
    if ( ! is.null(xsi.inits) ){
      #      xsi <- xsi.inits 
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
    
    #    if ( ! is.null(Y) ){ 
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
    # redefine Y in case of multiple persons
    if (tp>1){
      if ( nrow(gresp) != nrow(Y) ){
        Ypid <- rowsum( Y , pid0 )
        Y <- Ypid / Ypid[,1]
      }
    }
       
    #    W <- t(Y * pweights) %*% Y
    W <- crossprod(Y * pweights ,  Y )
    
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
    # ARb 2013-08-20: Handling of no beta constraints	
    # ARb 2013-08-24: correction	
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
       
    #***
    #*** ARb 2013-03-26   o Change setup of calculation
    
    # define response indicator matrix for missings
    #    resp.ind <- 1 - is.na(resp)
    #    resp.col.ind <- as.numeric(X.red[,ifelse("item" %in% colnames(X.red), "item", 1)]) 
    resp.ind.list <- list( 1:nitems )
    gresp.ind <- 1 - is.na( gresp ) 
    gresp.noStep.ind <- 1 - is.na( gresp.noStep )	 
    #	nomiss <- sum( is.na(gresp.noStep) == 0 )  	#*** included nomiss in M step regression
    resp.ind <- gresp.noStep.ind
    nomiss <- mean( gresp.noStep.ind ) == 1
    #***
    miss.items <- rep(0,nitems)
    for (i in 1:nitems){ 
      #      resp.ind.list[[i]] <- which( resp.ind[,resp.col.ind[i]] == 1)  
      resp.ind.list[[i]] <- which( gresp.noStep.ind[,i] == 1)  
      miss.items[i] <- i * ( length(resp.ind.list[[i]]) == 0 )
    }
    #    resp[ is.na(resp) ] <- 0 	# set all missings to zero
    gresp0.noStep <- gresp.noStep
    gresp[ is.na( gresp) ] <- 0
    gresp.noStep[ is.na( gresp.noStep) ] <- 0
    #    gresp.noStep.ind <- resp.ind[ ,resp.col.ind]
    #    gresp.ind <- resp.ind[ ,as.numeric(X[,ifelse("item" %in% colnames(X), "item", 1)])]
        
    #*****
    # ARb 2013-09-09: deletion of items
    miss.items <- miss.items[ miss.items > 0 ]	
    if ( length(miss.items) == 0 ){ delete.red.items <- FALSE }
    if (delete.red.items){					
      miss.itemsK <- NULL
      for (kk in 1:maxK ){
        miss.itemsK <- c( miss.itemsK , ( miss.items - 1 )* maxK + kk )
      }
      
      miss.itemsK <- sort(miss.itemsK)
      gresp <- gresp[ , - miss.itemsK ]
      gresp.noStep <- gresp.noStep[ , - miss.items ]
      gresp.noStep.ind <- gresp.noStep.ind[ , - miss.items ]		
      A <- A[ - miss.items , , , drop=FALSE]
      B <- B[ - miss.items , , ,drop=FALSE]	
      resp.ind.list <- resp.ind.list[ - miss.items ]
      resp.ind <- resp.ind[ , - miss.items ]
      nitems <- ncol(gresp.noStep)
      pweightsM <- outer( pweights , rep(1,nitems) )
      
      if (progress){ 
        cat("    * Reduced Response Data:" , nstud , "Persons and " , 
            ncol(gresp.noStep) , "Generalized Items (" , paste(Sys.time()) ,")\n" )  ;
        utils::flush.console()	  
      } 		
      
    }
    
    #****
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
        
    #***********************	
    #...TK - 24.08.2012:  cA returned from designMatrices
    # sufficient statistics
    #    ItemScore <- (cResp %*% cA) %t*% pweights
    col.index <- rep( 1:nitems , each = maxK )
    
    cResp <- (gresp.noStep+1)*gresp.noStep.ind	
    cResp <- cResp[ , col.index  ]
    
    
    cResp <- 1 * ( cResp == matrix( rep(1:(maxK), nitems) , nrow(cResp) , 
                                    ncol(cResp) , byrow=TRUE ) )
									
    cA <- t( matrix( aperm( A , c(2,1,3) ) , nrow = dim(A)[3] , byrow = TRUE ) )
    cA[is.na(cA)] <- 0		
    if ( stats::sd(pweights) > 0 ){ 
      ItemScore <- as.vector( t( colSums( cResp * pweights ) ) %*% cA )
    } else { 
      ItemScore <- as.vector( t( colSums( cResp) ) %*% cA )			
    }
    
    # cat("calc ItemScore" ) ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1
       
    if (progress){ 
      cat("    * Calculated Sufficient Statistics   (", 
          paste(Sys.time()) , ")\n") ; utils::flush.console()	  
    }   			
    # starting values for xsi
    gresp.ind.tmp <- gresp.noStep.ind[ , col.index  ]
    #    gresp.ind.tmp[,- grep(paste0("-step",(maxK-1)),colnames(gresp))] <- 0
    # ItemMax <- (gresp.ind.tmp %*% cA) %t*% pweights
	ItemMax <- crossprod(gresp.ind.tmp %*% cA , pweights )
    ItemMax <- as.vector( t( colSums( gresp.ind.tmp * pweights ) ) %*% cA )    
    
    
    # cat("calc ItemMax" ) ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1
    
    
    #    xsi[est.xsi.index] <- - 
    #  log(abs(ItemScore[est.xsi.index]/(ItemMax[est.xsi.index]-ItemScore[est.xsi.index])))  
    xsi[est.xsi.index] <- - log(abs(( ItemScore[est.xsi.index]+.5)/
                                      (ItemMax[est.xsi.index]-ItemScore[est.xsi.index]+.5) ) )							  
    # starting values of zero
    if( xsi.start0 == 1){ 
			xsi <- 0*xsi 
					}
    if( xsi.start0 == 2){ 
		ind1 <- which( dimnames(A)[[3]] %in% colnames(resp) )
		ind2 <- which( dimnames(A)[[3]] %in% paste0( "step" ,1:9) )
		ind3 <- setdiff( seq(1,length(xsi) ) , union(ind1,ind2) )
		xsi[ind3] <- 0
					}
					
    
    
    #log of odds ratio of raw scores  
    xsi[ is.na(xsi) ] <- 0
    if ( ! is.null(xsi.inits) ){  
      #			xsi <- xsi.inits  
      xsi[ xsi.inits[,1] ] <- xsi.inits[,2]			
    }
    if ( ! is.null(xsi.fixed) ){   xsi[ xsi.fixed[,1] ] <- xsi.fixed[,2] }
    
    xsi.min.deviance <- xsi
    beta.min.deviance <- beta
    variance.min.deviance <- variance
    
    # nodes
    if ( snodes == 0 ){ 
      theta <- as.matrix( expand.grid( as.data.frame( matrix( rep(nodes, ndim) , ncol = ndim ) ) ) )
      #we need this to compute sumsig2 for the variance
      #      theta2 <- matrix(theta.sq(theta), nrow=nrow(theta),ncol=ncol(theta)^2)            
      theta2 <- matrix(theta.sq2(theta), nrow=nrow(theta),ncol=ncol(theta)^2)            
      # grid width for calculating the deviance
      thetawidth <- diff(theta[,1] )
      thetawidth <- ( ( thetawidth[ thetawidth > 0 ])[1] )^ndim 
      thetasamp.density <- NULL
    } else {
      # sampled theta values
      if (QMC){			
        r1 <- sfsmisc::QUnif (n=snodes, min = 0, max = 1, n.min = 1, p=ndim, leap = 409)
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
    
    hwt.min <- 0
    rprobs.min <- 0
    AXsi.min <- 0
    B.min <- 0
    deviance.min <- 1E100
    itemwt.min <- 0
    
	#*****
	#@@@@ 2015-06-26
	Avector <- as.vector(A)
	Avector[ is.na(Avector) ] <- 0
	unidim_simplify <- TRUE
    YSD <- max( apply( Y , 2 , stats::sd ) )
    if (YSD > 10^(-15) ){ YSD <- TRUE } else { YSD <- FALSE }		
	if (G > 1){ unidim_simplify <- FALSE }
	if ( YSD){ unidim_simplify <- FALSE }	
	if (  is.null(beta.fixed) ){ unidim_simplify <- FALSE }
	#@@@@	
	
	#@@@@AAAA@@@@@
	# xsi acceleration
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
	
    ##**SE
    se.xsi <- 0*xsi
    se.B <- 0*B
    se.xsi.min <- se.xsi
    se.B.min <- se.B
    

    
    devch <- 0
    
    # display
    disp <- "....................................................\n"
    # define progress bar for M step        
 # cat("rest  " ) ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1								 

   
    ##############################################################   
    ##############################################################   
    ##############################################################   
    #Start EM loop here
    while ( ( (!betaConv | !varConv)  | ((a1 > conv) | (a4 > conv) | (a02 > convD)) )  & 
              (iter < maxiter) ) { 

      # a0 <- Sys.time()	
      iter <- iter + 1
      if (progress){ 
        cat(disp)	
        cat("Iteration" , iter , "   " , paste( Sys.time() ) )
        cat("\nE Step\n") ; utils::flush.console()
      }
      # calculate nodes for Monte Carlo integration	
      if ( snodes > 0){
        #        theta <- beta[ rep(1,snodes) , ] +  t ( t(chol(variance)) %*% t(theta0.samp) )
        theta <- beta[ rep(1,snodes) , ] + theta0.samp %*% chol(variance) 
        # calculate density for all nodes
        thetasamp.density <- mvtnorm::dmvnorm( theta , mean = as.vector(beta[1,]) , sigma = variance )
        # recalculate theta^2
        #        theta2 <- matrix( theta.sq(theta) , nrow=nrow(theta) , ncol=ncol(theta)^2 )   
        theta2 <- matrix( theta.sq2(theta) , nrow=nrow(theta) , ncol=ncol(theta)^2 )   
      }			
      olddeviance <- deviance
      # calculation of probabilities
      res <- calc_prob.v5(iIndex=1:nitems , A=A , AXsi=AXsi , B=B , xsi=xsi , theta=theta , 
                          nnodes=nnodes , maxK=maxK , recalc=TRUE )	
      # cat("calc prob") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1							  

      rprobs <- res[["rprobs"]]
      AXsi <- res[["AXsi"]]
      # calculate student's prior distribution
      gwt <- stud_prior.v2(theta=theta , Y=Y , beta=beta , variance=variance , nstud=nstud , 
                           nnodes=nnodes , ndim=ndim,YSD=YSD , unidim_simplify=unidim_simplify ,
							snodes = snodes )
      #print( head(gwt))						   
      
      # calculate student's likelihood
      res.hwt <- calc_posterior.v2(rprobs=rprobs , gwt=gwt , resp=gresp.noStep , 
					nitems=nitems , resp.ind.list= resp.ind.list , normalization=TRUE , 
                    thetasamp.density=thetasamp.density , snodes=snodes ,
                    resp.ind=resp.ind	, avoid.zerosum=TRUE)	
      hwt <- res.hwt[["hwt"]] 

# cat("calc posterior") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1		  
      if (progress){ cat("M Step Intercepts   |"); utils::flush.console() }
      # collect old values for convergence indication
      oldxsi <- xsi
      oldbeta <- beta
      oldvariance <- variance 
      # M step: estimation of beta and variance
      resr <- mstep.regression( resp=gresp.noStep , hwt=hwt , resp.ind=gresp.noStep.ind , 
                                pweights=pweights ,  pweightsM=pweightsM , Y=Y , theta=theta , 
                                theta2=theta2 , YYinv=YYinv , 
                                ndim=ndim , nstud=nstud , beta.fixed=beta.fixed , variance=variance , 
                                Variance.fixed=variance.fixed , group=group ,  G=G , snodes = snodes ,
                                nomiss=nomiss)
														
	if ( ( iter < 2 ) & is.na(resr$variance) ){
		stop("Choose argument control=list( xsi.start0=TRUE, ...) ")
						}												
														
      beta <- resr$beta
      #cat("m step regression") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1		        
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
        #		if( det(variance) < 10^(-20) ){
        #		  stop("\n variance is close to singular or zero. Estimation cannot proceed")
        #@@ARb: I would not prefer to stop the program but adding a small
        #       constant in the diagonal.
        #				} 
      }
	  
        #@@@@AAAA@@@@@
		# variance acceleration
		if ( variance_acceleration$acceleration != "none" ){		
			variance_acceleration <- accelerate_parameters( xsi_acceleration=variance_acceleration , 
							xsi=as.vector(variance) , iter=iter , itermin=3)
			variance <- matrix( variance_acceleration$parm , nrow= nrow(variance) , ncol=ncol(variance) )
								}
	    #@@@@AAAA@@@@@		  	  
	  
      if (max(abs(variance-oldvariance)) < conv){ varConv <- TRUE      }
      ######################################
      #M-step item parameters
      converge <- FALSE
      Miter <- 1	  
      old_increment <- rep( max.increment , np )
      #	  if (TRUE & (iter>1) ){
      #	     old_increment <- xsi.change
      #				}

     
      est.xsi.index <- est.xsi.index0	  
      while (!converge & ( Miter <= Msteps ) ) {	  
        # Only compute probabilities for items contributing to param p
        if (Miter > 1){ 
          res.p <- calc_prob.v5( iIndex=1:nitems , A=A , AXsi=AXsi , B=B , 
                                 xsi=xsi , theta=theta , nnodes=nnodes, maxK=maxK)					
          rprobs <- res.p[["rprobs"]]            
        }
        #	indexIP.list2 <- unlist(indexIP.list)
        # ==> Vector with item indices for parameter estimation
        #	indexIP.no <- cbind( c(1 , lipl[-length(lipl)]+1 ) , lipl )
        # ==> Vector with start and end indices for item parameter estimation
        
        res <- calc_exp_TK3( rprobs , A , np , est.xsi.index , itemwt ,
                             indexIP.no , indexIP.list2 , Avector )
        xbar <- res$xbar
        xbar2 <- res$xbar2
        xxf <- res$xxf
        
        
        # Compute the difference between sufficient statistic and expectation
        diff <- as.vector(ItemScore) - xbar
        #Compute the Newton-Raphson derivative for the equation to be solved
        deriv <- xbar2 - xxf 			
        increment <- diff*abs(1/( deriv + 10^(-20) ) )
        if ( !is.null( xsi.fixed) ){ increment[ xsi.fixed[,1] ] <- 0 } 
        #!!!	  necessary to include statement to control increment?
        ci <- ceiling( abs(increment) / ( abs( old_increment) + 10^(-10) ) )
        a1 <- abs(old_increment)
#		a1 <- max.increment
#***
choice1 <- TRUE
#choice1 <- FALSE
if (choice1){		
        increment <- ifelse( abs( increment) > a1  , 
                             increment/(2*ci) , increment )					 
        old_increment <- abs(increment)							 
			}
  	    #****
        # inclusion
if (!choice1){		
                increment <- ifelse( abs( increment) > old_increment  , 
                                     sign(increment)*old_increment / 1.5 , increment )
        		max.increment <- max( abs(increment) )
#				old_increment <- max( abs(increment) )				
			}	
        #****
               
#        w <- 1		
#		old_increment <- w * abs(increment) + (1-w)*abs(old_increment)

		
        xsi <- xsi+increment   # update parameter p
        
        # stabilizing the algorithm | ARb 2013-09-10
        if (fac.oldxsi > 0 ){
		  fac.oldxsi1 <- (devch>0)*fac.oldxsi	  
		  # fac.oldxsi1 <- fac.oldxsi
          xsi <-  (1-fac.oldxsi1) * xsi + fac.oldxsi1 *oldxsi
        }
        
        
        #		est.xsi.index <- which( abs(increment) > convM )		
        if ( max(abs(increment)) < convM ) { 
          converge <- TRUE 
        }
        Miter <- Miter + 1						
        
        ##**SE
        se.xsi <- sqrt( 1 / abs(deriv) )
        if ( ! is.null( xsi.fixed) ){ se.xsi[ xsi.fixed[,1] ] <- 0 } 
        ##**        
        # progress bar
        if (progress){ 
          #    cat( paste( rep("-" , sum( mpr == p ) ) , collapse="" ) )
          cat("-")
          utils::flush.console()
        }
      } # end of all parameters loop
      
      
        #@@@@AAAA@@@@@
		# acceleration
		if ( xsi_acceleration$acceleration != "none" ){		
			xsi_acceleration <- accelerate_parameters( xsi_acceleration=xsi_acceleration , 
							xsi=xsi , iter=iter , itermin=3)
			xsi <- xsi_acceleration$parm
								}
	    #@@@@AAAA@@@@@	  
	  
      #***
      # decrease increments in every iteration
      if( increment.factor > 1){   max.increment <-  max.increment / increment.factor }	  
      

      
# cat("m step item parameters") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1		        
      # calculate deviance
      if ( snodes == 0 ){ 
        deviance <- - 2 * sum( pweights * log( res.hwt$rfx * thetawidth ) )
      } else {
        #        deviance <- - 2 * sum( pweights * log( res.hwt$rfx ) )
        # deviance <- - 2 * sum( pweights * log( rowMeans( res.hwt$swt ) ) )
		deviance <- - 2 * sum( pweights * log( res.hwt$rfx   ) )
      }
      deviance.history[iter,2] <- deviance
      a01 <- abs( ( deviance - olddeviance ) / deviance  )
      a02 <- abs( ( deviance - olddeviance )  )	

	if (con$dev_crit == "relative" ){ a02 <- a01 }
      
      if( deviance > deviance.min ){ 
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
      # xsi.change <- abs( xsi - oldxsi )
      
      a2 <- max( abs( beta - oldbeta ))	
      a3 <- max( abs( variance - oldvariance ))
	  devch <- -( deviance - olddeviance )
      if (progress){ 
        cat( paste( "\n  Deviance =" , round( deviance , 4 ) ))        
        cat( " | Deviance change:", round( devch  , 4 ) )		
	cat( " | Relative deviance change:", round( a01  , 8 ) )
        if ( devch < 0 & iter > 1 ){ 
          cat ("\n!!! Deviance increases!                                        !!!!") 
          cat ("\n!!! Choose maybe fac.oldxsi > 0 and/or increment.factor > 1    !!!!") 			
        }
        cat( "\n  Maximum intercept parameter change:" , round( a1 , 6 ) )
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
      

      
    } # end of EM loop
    #############################################################
    #############################################################
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
    #***
    resp <- gresp0.noStep
    resp.ind <- gresp.noStep.ind
    
	#****
	# look for non-estimable xsi parameters
#    xsi[ xsi == 99 ] <- NA	

	#******
	# generate input for fixed parameters
	xsi.fixed.estimated <- generate.xsi.fixed.estimated( xsi , A )
	B.fixed.estimated <- generate.B.fixed.estimated(B)
	
	
    ##**SE  
    # standard errors of AXsi parameters
    # check for missing entries in A
    se.AXsi <- 0*AXsi
    A1 <- A
    A1[ is.na(A) ] <- 0
    se.xsiD <- diag( se.xsi^2 )
    for (kk in 1:maxK){  # kk <- 1
      #	se.AXsi[,kk] <- sqrt( diag( A1[,kk,] %*% se.xsiD %*% t( A1[,kk,]) ) )
	  dim_A1 <- dim(A1)
      A1_kk <- A1[,kk,]
      if ( is.vector(A1_kk) ){
        A1_kk <- matrix( A1_kk , nrow=dim_A1[1] , ncol=dim_A1[3] )
      }
      se.AXsi[,kk] <- sqrt( diag( A1_kk %*% se.xsiD %*% t( A1_kk ) ) )	
      #****		
    }
    
    ##*** Information criteria
    ic <- .TAM.ic( nstud , deviance , xsi , xsi.fixed ,
                   beta , beta.fixed , ndim , variance.fixed , G ,
                   irtmodel , B_orig=NULL , B.fixed , E , est.variance =TRUE ,
                   resp , est.slopegroups=NULL , 
				   variance.Npars= NULL , group )
    
    #***
    # calculate counts
    res <- .tam.calc.counts( resp = gresp.noStep, theta , 
				resp.ind=gresp.noStep.ind , 
                group , maxK , pweights , hwt )
    n.ik <- res$n.ik
    pi.k <- res$pi.k 
    
    
    #****
    # collect item parameters
    
    #	item1 <- .TAM.itempartable( resp , maxK , AXsi , B , ndim ,
    #				 resp.ind , rprobs,n.ik,pi.k)
    item1 <- .TAM.itempartable( resp=gresp.noStep , maxK , AXsi , B , ndim ,
                                resp.ind=gresp.noStep.ind , rprobs,n.ik,pi.k)
    
  
    #####################################################
    # post ... posterior distribution	
    # create a data frame person	
    person <- data.frame( "pid"=pid , "case" = 1:nstud , "pweight" = pweights )
    
    #    person$score <- rowSums( resp * resp.ind )
    #    resp2 <- resp
    resp2 <- gresp.noStep
    resp2[ is.na(resp2) ] <- 0
    #    person$score <- rowSums( resp * resp.ind , na.rm=TRUE)
    person$score <- rowSums( gresp.noStep * gresp.noStep.ind , na.rm=TRUE)
    
    # use maxKi here; from "design object"
    nstudl <- rep(1,nstud)
    
    #    person$max <- rowSums( outer( nstudl , apply( resp2 ,2 , max , na.rm=T) ) * resp.ind )
    person$max <- rowSums( outer( nstudl , apply( resp2 ,2 , max , na.rm=T) ) * gresp.noStep.ind )
    
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
      EAP.variance <- stats::weighted.mean( person$EAP^2 , pweights ) - 
				( stats::weighted.mean( person$EAP , pweights ) )^2
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
        EAP.variance <- stats::weighted.mean( person$EAP^2 , pweights ) - 
					( stats::weighted.mean( person$EAP , pweights ) )^2
        EAP.error <- stats::weighted.mean( person$SD.EAP^2 , pweights )
        EAP.rel[dd] <- EAP.variance / ( EAP.variance + EAP.error )	
        colnames(person)[ which( colnames(person) == "EAP" ) ] <- paste("EAP.Dim" , dd , sep="")
        colnames(person)[ which( colnames(person) == "SD.EAP" ) ] <- paste("SD.EAP.Dim" , dd , sep="")				
      }
#      person <- data.frame( "pid" = pid , person )
    }
    ############################################################
    s2 <- Sys.time()
    
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
        print( stats::cov2cor(varianceM) , 4 )	
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
    
    ####################################
    # collect xsi parameters
    xsiFacet <- as.data.frame( (xsi.constr$xsi.table)[,1:2]	)
    obji <- data.frame( "parameter" = dimnames(A)[[3]] , 
                        "xsi"=xsi , "se.xsi"=se.xsi ) 		
    rownames(obji) <- paste(obji$parameter)
    rownames(xsiFacet) <- paste( xsi.constr$xsi.table[,1] )

    xsi1 <- merge( x = xsiFacet , y= obji , by="parameter" , all=TRUE )
    A1 <- xsi.constr$xsi.constraint %*% xsi
	
    xsi1[ match( rownames(xsi.constr$xsi.constraint) , xsi1$parameter) , "xsi" ] <- A1

    xsi1 <- xsi1[ match( xsiFacet$parameter , xsi1$parameter) , ]
    xsi.facets <- xsi1
    rownames(xsi.facets) <- NULL
    i1 <- grep( "Intercept" , xsi.facets$parameter)

    if ( length(i1) > 0 ){
      xsi.facets <-  xsi.facets[ - i1 , ] 
    }

	#@@@@@@@@@@@@@@@@@@@@@@@@ control xsi.facets
	if( xsi.constr$intercept_included ){	
		ind <- which( paste(xsi.facets$facet) == "item" )
		n1 <- length(ind)
			if ( n1 > 0 ){
				itemc <- itemnames[n1]
				itemo <- paste0("item" , n1 )		
				g1 <- which( paste(xsi.facets$parameter) == itemo )
				if ( length(g1) > 0 ){
					xsi.facets$parameter[g1] <- itemc
									}

				g1 <- grep( paste0(itemo , ":") , paste(xsi.facets$parameter)  )
				if ( length(g1) > 0 ){
					xsi.facets$parameter <- gsub( paste0(itemo , ":") , paste0(itemc , ":")  ,
										paste(xsi.facets$parameter) )
									}

				g1 <- grep( paste0(itemo , "-") , dimnames(A)[[1]]  )
				if ( length(g1) > 0 ){
					dimnames(A)[[1]] <- gsub( paste0(itemo , "-") , paste0(itemc , "-")  ,
										dimnames(A)[[1]] )
									}

									
						}
				}	
	#@@@@@@@@@@@@@@@@@@@@@@@@@
	
    xsi <- obji[,-1]
    rownames(xsi) <- dimnames(A)[[3]]
    
    if(delete.red.items) resp <- resp[,-miss.items]
    colnames(resp) <- dimnames(A)[[1]]
    
	 res.hwt <- calc_posterior.v2(rprobs=rprobs , gwt=1+0*gwt , resp=resp , nitems=nitems , 
                                   resp.ind.list=resp.ind.list , normalization=FALSE , 
                                   thetasamp.density=thetasamp.density , snodes=snodes ,
                                   resp.ind=resp.ind	)	
      res.like <- res.hwt[["hwt"]] 		

	
    # Output list
    deviance.history <- deviance.history[ 1:iter , ]
    res <- list( "xsi" = xsi , "xsi.facets" = xsi.facets , 
                 "beta" = beta , "variance" = variance ,
                 "item" = item1 , 
                 "person" = person , pid = pid , "EAP.rel" = EAP.rel , 
                 "post" = hwt ,  "rprobs" = rprobs , "itemweight" = itemwt ,
                 "theta" = theta , 
                 "n.ik" = n.ik , "pi.k" = pi.k ,
                 "Y" = Y , "resp" = resp , 
                 "resp.ind" = resp.ind , "group" = group , 
                 "G" = if ( is.null(group)){1} else { length(unique( group ) )} , 
                 "groups" = if ( is.null(group)){1} else { groups } , 			   			   
                 "formulaY" = formulaY , "dataY" = dataY , 
                 "pweights" = pweights , 
                 "time" = c(s1,s2,s2-s1) , "A" = A , "B" = B  ,
                 "se.B" = se.B , 
                 "nitems" = nitems , "maxK" = maxK , "AXsi" = AXsi ,
                 "AXsi_" = - AXsi ,			   
                 "se.AXsi" = se.AXsi , 
                 "nstud" = nstud , "resp.ind.list" = resp.ind.list ,
                 "hwt" = hwt , "like" = res.like , "ndim" = ndim ,
                 "xsi.fixed" = xsi.fixed , 
				 "xsi.fixed.estimated" = xsi.fixed.estimated , 
				 "B.fixed.estimated" = B.fixed.estimated , 
				 "beta.fixed" = beta.fixed , "Q"=Q,
                 "formulaA"=formulaA , "facets"=facets ,
				 "xsi.constr" = xsi.constr , 
                 "variance.fixed" = variance.fixed ,
                 "nnodes" = nnodes , "deviance" = deviance ,
                 "ic" = ic , 
                 "deviance.history" = deviance.history ,
                 "control" = con1a , "irtmodel" = irtmodel ,
                 "iter" = iter , "resp_orig" = resp_orig ,
                 "printxsi"=TRUE , "YSD"=YSD , "PSF" = PSF ,
				 CALL = CALL 
                 #			   "design"=design
                 #			   "xsi.min.deviance" = xsi.min.deviance ,
                 #			   "beta.min.deviance" = beta.min.deviance , 
                 # "variance.min.deviance" = variance.min.deviance 
    )
    class(res) <- "tam.mml"
    return(res)
  }
