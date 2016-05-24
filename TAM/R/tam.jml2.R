
##################################################################
##################################################################
##################################################################

tam.jml2 <-
  function( resp , group = NULL , adj=.3 , disattenuate = FALSE ,
            bias = TRUE, xsi.fixed=NULL ,  xsi.inits = NULL ,  
            A=NULL , B=NULL , Q=NULL , ndim=1 ,
            pweights = NULL , control = list() 
            # control can be specified by the user 
  ){
    
    #------------------------------------
    # INPUT:
    # control:
    #    control = list( nodes = seq(-6,6,len=15) , 
    #		          			convD = .001 ,conv = .0001 , convM = .0001 , Msteps = 30 ,            
    #                   maxiter = 1000 , progress = TRUE) 
    # progress ... if TRUE, then display progress
    #-------------------------------------

	maxiter <- conv <- progress <- tamobj <- convM <- Msteps <- NULL 
    R <- NULL	
	
#    adj <- 0.3   # adjustment for perfect and zero scores
    s11 <- Sys.time()
    # attach control elements
    e1 <- environment()
    con <- list( nodes = seq(-6,6,len=21) , snodes = 0 ,
                 convD = .001 ,conv = .0001 , convM = .0001 , Msteps = 10 ,            
                 maxiter = 1000 , progress = TRUE )  	
    con[ names(control) ] <- control  
    Lcon <- length(con)
    con1a <- con1 <- con ; 
    names(con1) <- NULL
    for (cc in 1:Lcon ){
      assign( names(con)[cc] , con1[[cc]] , envir = e1 ) 
    }
	
	resp <- add.colnames.resp(resp)
	
    # maximum no. of categories per item.
    maxK <- max( resp , na.rm=TRUE ) + 1
    
    resp <- as.matrix(resp)
    nitems <- ncol(resp)       # number of items
    nstud <- nrow(resp)        # number of students
    
    nitems <- ncol(resp)       # number of items
    nstud <- nrow(resp)        # number of students
    ################################
    # create design matrices
    
    design <- designMatrices( modeltype = "PCM" , maxKi = NULL , resp = resp , 
                              A = A , B = B , Q = Q , R = R, ndim = ndim )
    A <- design$A
    A.0 <- A
    A.0[ is.na(A.0) ] <- 0  
    
    
    B <- design$B
    B.0 <- B
    B.0[ is.na(B.0) ] <- 0    
    cA <- design$flatA
    cA[is.na(cA)] <- 0
    
    # number of parameters
    np <- dim(A)[[3]]
    errorP <- rep(0,np)
    
    AXsi <- matrix(0, nrow=nitems, ncol=maxK) 
    
    if ( is.null( pweights) ){
      pweights <- rep(1,nstud) # weights of response pattern
    }
    
    # normalize person weights to sum up to nstud
    pweights <- nstud * pweights / sum(pweights)
    # a matrix version of person weights
    pweightsM <- outer( pweights , rep(1,nitems) )
    
    # xsi inits
    if ( ! is.null(xsi.inits) ){
      xsi <- xsi.inits 
    } else { xsi <- rep(0,np)   } 
    if ( ! is.null( xsi.fixed ) ){
      xsi[ xsi.fixed[,1] ] <- xsi.fixed[,2]
      est.xsi.index <- setdiff( 1:np , xsi.fixed[,1] )
    } else { est.xsi.index <- 1:np }
    
    
    # group indicators for variance matrix
    if ( ! is.null(group) ){ 
      groups <- sort(unique(group))
      G <- length(groups)
      # user must label groups from 1, ... , G
      if ( length( setdiff( 1:G , groups)  ) > 0 ){
        stop("Label groups from 1, ...,G\n")
      }
    } else { G <- 1 }
    
    # define response indicator matrix for missings
    resp.ind <- 1 - is.na(resp)
    resp.ind.list <- list( 1:nitems )
    for (i in 1:nitems){ resp.ind.list[[i]] <- which( resp.ind[,i] == 1)  }
    resp[ is.na(resp) ] <- 0 	# set all missings to zero
    
    # Create an index linking items and parameters
    indexIP <- colSums(aperm(A, c(2,1,3)) != 0, na.rm = TRUE)
    # define list of elements for item parameters
    indexIP.list <- list( 1:np )
    for ( kk in 1:np ){
      indexIP.list[[kk]] <- which( indexIP[,kk] > 0 )
    }
    
    
    # These sufficient statistics must be changed
    # to make it more general
    # First extension:  pweights and dependent on A; needs to be further extended (e.g., different number of categories)
    # Second extension: multiple category option       -> resp \in 0:maxKi (see method definition calc_posterior_TK)
    #                                                  -> length(ItemScore) = np (see diff computation in M Step)
    #                   multiple category option Bugfix
    #                                                  -> dim(cResp) = (nstud, nitems*maxK)
    #                                                  -> adapt dim(A) to dim(cResp) for sufficient statistic (cf. print.designMatrices)
    
    col.index <- rep( 1:nitems , each = maxK )
    cResp <- resp[ , col.index  ]*resp.ind[ , col.index ]
    # This line does not take missings into account
    cResp <- 1 * t( t(cResp) == rep(0:(maxK-1), nitems) )
    ##@@@##
    # ARb: I added this line
    cResp <- cResp * resp.ind[ , col.index ] 
    cB <- t( matrix( aperm( B , c(2,1,3) ) , nrow = dim(B)[3] , byrow = TRUE ) )
    cB[is.na(cB)] <- 0
    
    # Item sufficient statistics
    # ItemScore <- (cResp %*% cA) %t*% pweights
	ItemScore <- crossprod(cResp %*% cA , pweights )
    
    # Computer possible maximum parameter score for each person
    maxAi <-  - (apply(-(A) , 3 , rowMaxs , na.rm=TRUE))  
    personMaxA <- resp.ind %*% maxAi
    # ItemMax <- personMaxA %t*% pweights
	ItemMax <- crossprod( personMaxA , pweights )
    
    #Adjust perfect and zero scores for the parameters
    ItemScore[ItemScore==ItemMax] <- ItemScore[ItemScore==ItemMax] + adj #..."+" sign, because ItemScore is -ve)
    ItemScore[ItemScore==0] <- ItemScore[ItemScore==0] - adj  
    
    #Initialise xsi
    xsi[est.xsi.index] <- - log(abs(ItemScore[est.xsi.index]/(ItemMax[est.xsi.index]-ItemScore[est.xsi.index])))  #log of odds ratio of raw scores

    
    #Compute person sufficient statistics (total score on each dimension)
    PersonScores <- cResp %*% cB
    # define response pattern
    rp3 <- as.data.frame( resp.pattern3( resp.ind ) )
    rp3$caseid <- 1:nstud
    #rp3 <- as.data.frame( rp3 )
    rp3$PersonScores <- PersonScores
    rp3$mp.scores <- paste( rp3$mp.index , rp3$PersonScores, sep="-")
    rp3$theta.index <- match( rp3$mp.scores , unique(rp3$mp.scores ) )
    rp3 <- rp3[ order(rp3$theta.index) , ]
    rp3.sel <- rp3[ c(TRUE , diff(rp3$theta.index) == 1 ) , ]
    rp3 <- rp3[ order( rp3$caseid) , ]
    
    rp3.pweightsM  <- rowsum( pweightsM ,  rp3$theta.index )
    
    #Compute possible maximum score for each item on each dimension
    maxBi <- apply(B , 3 , rowMaxs , na.rm = TRUE)
    
    #Compute possible maximum score for each person on each dimension
    PersonMaxB <- resp.ind %*% maxBi
    
    #Adjust perfect scores for each person on each dimension
    PersonScores[PersonScores==PersonMaxB] <- PersonScores[PersonScores==PersonMaxB] - adj
    
    #Adjust zero scores for each person on each dimension
    PersonScores[PersonScores==0] <- PersonScores[PersonScores==0] + adj
    
    #Initialise theta (WLE) values for all students
    theta <- log(PersonScores/(PersonMaxB-PersonScores)) #log of odds ratio of raw score
    
    deviance <- 0  
    deviance.history <- matrix( 0 , nrow=maxiter , ncol = 2)
    colnames(deviance.history) <- c("iter" , "deviance")
    deviance.history[,1] <- 1:maxiter
    
    
    iter <- 0 
    maxthetachange <- meanChangeWLE <- maxChangeWLE <- maxChangeP <- 999	# item parameter change
    # display
    disp <- "....................................................\n"
    
    ##############################################################   
    #Start convergence loop here
    while ( (( maxthetachange > conv) | (maxChangeP > conv))  & (iter < maxiter) ) { 
      
      iter <- iter + 1
      if (progress){ 
        cat(disp)	
        cat("Iteration" , iter , "   " , paste( Sys.time() ) )
        #      cat( "\n" )
        utils::flush.console()
      }
      olddeviance <- deviance
      
      
      #update theta, ability estimates
      #    jmlAbility <- tam.jml.WLE ( resp , resp.ind, A, B, nstud, nitems, maxK, convM, 
      #                                PersonScores, theta, xsi, Msteps, WLE=FALSE)
      theta_old <- theta
      jmlAbility <- tam.jml.WLE ( resp=resp , resp.ind=resp.ind[ rp3.sel$caseid,] , 
                                  A=A, B=B, 
                                  nstud=nrow(rp3.sel) , 
                                  nitems=nitems, maxK=maxK, convM=convM, 
                                  PersonScores=PersonScores[ rp3.sel$caseid ] , 
                                  theta=theta[ rp3.sel$caseid , , drop=FALSE] , 
                                  xsi=xsi, Msteps=Msteps, WLE=FALSE)		
      theta <- jmlAbility$theta								
      theta <- theta[ rp3$theta.index , , drop=FALSE]
      
      if (is.null( xsi.fixed))  theta <- theta - mean(theta)
      meanChangeWLE <- jmlAbility$meanChangeWLE
      maxthetachange <- max( abs( theta - theta_old ) )
      errorMLE <- jmlAbility$errorWLE

      
      #update xsi, item parameters
      #    jmlxsi0 <- tam.jml.xsi ( resp , resp.ind, A, B, nstud, nitems, maxK, convM, 
      #                            ItemScore, theta, xsi, Msteps, pweightsM,
      #                            est.xsi.index)   
      jmlxsi <- tam.jml.xsi2( resp , resp.ind, A=A,A.0=A.0 ,  B=B, nstud, nitems, maxK, convM, 
                              ItemScore, theta, xsi, Msteps, pweightsM,
                              est.xsi.index , rp3 , rp3.sel , rp3.pweightsM	)
      xsi[est.xsi.index] <- jmlxsi$xsi[est.xsi.index]
      maxChangeP <- jmlxsi$maxChangeP
      errorP[est.xsi.index] <- jmlxsi$errorP[est.xsi.index]
      #Deviance
      #Calculate Axsi. Only need to do this once for ability estimates.
      for (k in 1:maxK){ 
        AXsi[ , k ] <- A[ ,k , ] %*% xsi
      }
      theta.unique <- unique( theta )
      nstud1 <- length( theta.unique )
      res <- calc_prob.v5(iIndex = 1:nitems , A , AXsi , 
                          B , xsi , theta.unique , nstud1, maxK , recalc=FALSE )      	
      rprobs <- res[["rprobs"]] 
      crprobs <- t( matrix( aperm( rprobs , c(2,1,3) ) , nrow = dim(rprobs)[3] , byrow = TRUE ) )
      crprobs <- crprobs[ , match( theta , theta.unique) ]
      #    cr <- log(crprobs + 10^(-80)) * t(cResp)
      cr <- log(crprobs ) * t(cResp)
      deviance <- -2 * sum(cr, na.rm=TRUE)
      deviance.history[iter,2] <- deviance
      # progress bar
      if (progress){ 
        cat( paste( "\n  Deviance =" , round( deviance , 4 ) ))
        if (iter>1){ cat( " | Deviance change:", -round( deviance-olddeviance , 4 ) ) }
        #      cat( "\n  Mean MLE/WLE change:" , round( meanChangeWLE , 6 ) )
        cat( "\n  Maximum MLE/WLE change:" , round( maxthetachange , 6 ) )	  
        cat( "\n  Maximum item parameter change:" , round( maxChangeP , 6 ) )  
        cat( "\n" )
        utils::flush.console()
      }
      #@ARb 2012-08-27
      # stop loop (break) if there is no change in deviation
      if ( abs( deviance-olddeviance) < 10^(-10) ){ break }
    }# end of all convergence 
    
    s1 <- Sys.time()  
    #After convergence, compute final WLE (WLE set to TRUE)
    jmlWLE <- tam.jml.WLE ( tamobj , resp , resp.ind[ rp3.sel$caseid,], A, B, 
                            nrow(rp3.sel) , nitems, maxK, convM, 
                            PersonScores[ rp3.sel$caseid ] , 
                            theta[ rp3.sel$caseid , , drop=FALSE] , xsi, Msteps, WLE=TRUE)					  
    thetaWLE <- jmlWLE$theta[rp3$theta.index,1]
    # cat("\n WLE \n"); s2 <- Sys.time(); print(s2-s1) ; s1 <- s2    
    
    meanChangeWLE <- jmlWLE$meanChangeWLE
    errorWLE <- jmlWLE$errorWLE[rp3$theta.index]
    
    #WLE person separation reliability	
	WLEreliability <- WLErel(thetaWLE , errorWLE, pweights)
	
    # varWLE <- stats::var(thetaWLE)
    # WLEreliability <- (varWLE - mean(errorWLE^2)) / varWLE
    
    if (progress){ cat("\n Item fit calculation \n") }  
#     #Compute fit statistics
#     fit <- tam.jml.fit ( tamobj , resp , resp.ind, A, B, nstud, nitems, maxK, 
#                          ItemScore, theta, xsi, Msteps, pweightsM,
#                          est.xsi.index)
#     # cat("\n fit \n"); s2 <- Sys.time(); print(s2-s1) ; s1 <- s2    					   
#     outfitPerson <- fit$outfitPerson
#     outfitItem <- fit$outfitItem
#     infitPerson <- fit$infitPerson
#     infitItem <- fit$infitItem
#     outfitPerson_t <- fit$outfitPerson_t
#     outfitItem_t <- fit$outfitItem_t
#     infitPerson_t <- fit$infitPerson_t
#     infitItem_t <- fit$infitItem_t 
    
    #disattenuate
    if (disattenuate == TRUE) {
      thetaWLE <- sqrt(WLEreliability) * thetaWLE
      theta <- sqrt(WLEreliability) * theta
    }
    
    #bias
    if (bias == TRUE) {
      xsi[est.xsi.index] <- (nitems - 1)/nitems * xsi[est.xsi.index]     #Check this for more complex models
    }
    
  # collect item statistics
  item <- data.frame( "xsi.label" = dimnames(A)[[3]] ,
		"xsi.index" = 1:( length(xsi) ) , "xsi" = xsi ,
		"se.xsi" = errorP 
#     , "outfit" = outfitItem ,		"infit"=infitItem 
    )
	
    ############################################################
    s2 <- Sys.time()
    if (progress){
      cat(disp)
      cat( "\nStart: " , paste(s11))
      cat( "\nEnd: " , paste(s2),"\n")
      print(s2-s11)
      cat( "\n" )
    }
    
    # Output list
    deviance.history <- deviance.history[ 1:iter , ]
    res <- list( "item"=item , "xsi" = xsi ,  "errorP" = errorP , 
                 "theta" = theta[,1] ,"errorWLE" = errorWLE , "WLE" = thetaWLE , 
                 "WLEreliability" = WLEreliability ,
                 "PersonScores" = PersonScores , "ItemScore" = ItemScore ,             
                 "PersonMax" = PersonMaxB , "ItemMax" = ItemMax , 
#                  "outfitPerson" = outfitPerson , "outfitItem" = outfitItem, 
#                  "infitPerson" = infitPerson , "infitItem" = infitItem, 
#                  "outfitPerson_t" = outfitPerson_t , "outfitItem_t" = outfitItem_t, 
#                  "infitPerson_t" = infitPerson_t , "infitItem_t" = infitItem_t,
                 "deviance" = deviance, "deviance.history" = deviance.history, 
                 "resp" = resp , "resp.ind" = resp.ind , "group" = group ,
                 "pweights" = pweights , "A" = A , "B" = B  ,               
                 "nitems" = nitems , "maxK" = maxK , 
                 "nstud" = nstud , "resp.ind.list" = resp.ind.list ,
                 "xsi.fixed" = xsi.fixed , "deviance" = deviance ,
                 "deviance.history" = deviance.history ,
                 "control" = con1a , "iter"=iter)
    res$time <-  c(s11,s2,s2-s11)
    class(res) <- "tam.jml"
    return(res)
  }
