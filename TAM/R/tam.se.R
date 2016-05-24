tam.se <-
  function( tamobj , ...){
    SE.quick <- TRUE
    if(SE.quick){
      if(class(tamobj) == "tam.mml"){
        res <- tam.mml.se( tamobj, ...)
      }
      if(class(tamobj) == "tam.jml"){
        # res <- tam.jml.se( tamobj, ...)
      }
      
    }else{
      
    }
    
    return( res )
  }

tam.mml.se <-
  function( tamobj , numdiff.parm = .001){
    h <- numdiff.parm
    B <- tamobj$B
    A <- tamobj$A
    Y <- tamobj$Y
    YSD <- tamobj$YSD
    nitems <- tamobj$nitems
    xsi <- ( tamobj$xsi )[,1]	
    beta <- tamobj$beta
    variance <- tamobj$variance
    nstud <- tamobj$nstud
    AXsi <- tamobj$AXsi
    resp <- tamobj$resp
    ndim <- tamobj$ndim	
    theta <- tamobj$theta
    maxK <- tamobj$maxK
    thetawidth <- diff(theta[,1] )
    thetawidth <- ( ( thetawidth[ thetawidth > 0 ])[1] )^ndim
    hwt <- tamobj$hwt
    resp <- tamobj$resp
    resp.ind <- tamobj$resp.ind
    resp.ind.list <- tamobj$resp.ind.list 
    nnodes <- tamobj$nnodes	
    ntheta <- length(theta)
    irtmodel <- tamobj$irtmodel
    est.slopegroups <- tamobj$est.slopegroups
    # multiplication parameters for numerical differentiation
    mult <- c(0,1,-1)
    ##############################################
    # Item parameters xsi
    ##############################################	
    ip <- length(xsi)
    # create result object for item parameters
    se.xsi <- rep( 0 , ip )
    cat("Item parameters\n|")
    VP <- min( ip, 10 )
    cat(paste( rep("*" , VP) , collapse=""))
    if (VP<10){ disp_progress <- 1:ip } else {
      disp_progress <- 100* ( 1:ip ) / (ip+1)
      disp_progress <- sapply( seq(5,95,10) , FUN = function(pp){ # pp <- 5
        which.min( abs( disp_progress - pp ) )[1] }
      )
    }
    cat("|\n|")    
    # compute likelihood
    # prior distribution for each student (normal density)
    res0a <- stud_prior.v2( theta=theta , Y=Y , beta=beta , variance=variance , 
                            nstud=nstud , nnodes=nnodes , ndim=ndim , YSD=YSD ,
							unidim_simplify=FALSE )
    
    ll <- matrix( 0 , nrow=nstud , ncol=3 )	
    vv <- 1
    for (pp in 1:ip){
      # pp <- 10  # parameter pp
      if (pp == 1){ vec <- 1:3 } else { vec <- 2:3 }
      for (mm in vec){
        xsi0 <- xsi	
        xsi0[ pp ] <- xsi0[ pp] + mult[mm] * numdiff.parm
        ll0 <- 0
        # calculate probabilities
        res0 <- calc_prob.v5( iIndex=1:nitems , A=A , AXsi=AXsi , B=B , 
                              xsi=xsi0 , theta=theta , nnodes=nnodes , maxK=maxK )
        rprobs <- res0[["rprobs"]]
        # posterior distribution
        # calculate student's likelihood
        res0b <- calc_posterior.v2(rprobs=rprobs , gwt=res0a , resp=resp , 
                                   nitems=nitems , resp.ind.list=resp.ind.list , 
                                   normalization = FALSE , thetasamp.density = NULL , 
                                   snodes = 0 )$hwt
        # calculate individual log likelihood																						
        ll[,mm] <- log( rowSums( res0b * thetawidth ) )
      }
      se.xsi[pp] <- sqrt( - 1 / sum( tamobj$pweights * ( ll[,2] + ll[,3] - 2*ll[,1] )/h^2 ) )
      if ( ( pp==disp_progress[vv] ) & ( vv <=10) ){
        cat("-") ; 
        flush.console() 
        vv <- vv+1
      }
    }
    cat("|\n")	
    
    ##############################################
    # Item parameters B
    ##############################################
    se.B <- 0*B
    if( get.se.B <- irtmodel %in% c("2PL", "GPCM", "GPCM.design", "2PL.groups") ){
      cat("Loading parameters\n|")
      
      if(irtmodel=="GPCM.design"){
        basispar <- tamobj$basispar
        E <- tamobj$E
        
        pair.ind <- which(lower.tri(diag(length(basispar)), TRUE), arr.ind=TRUE)
        se.basispar <- rep(0,length(basispar))
        
        warning("SE for discrimination in irtmodel='GPCM.design' is experimental.")
      } 
      
      ip <- switch(irtmodel,
                   "2PL.groups"=length(unique(est.slopegroups)),
                   "GPCM.design"=length(basispar),
                   length(B[,-1,]) # default
      )
      
      # create result object for item parameters
      VP <- min( ip, 10 )
      cat(paste( rep("*" , VP) , collapse=""))
      if (VP<10){ disp_progress <- 1:ip } else {
        disp_progress <- 100* ( 1:ip ) / (ip+1)
        disp_progress <- sapply( seq(5,95,10) , FUN = function(pp){ # pp <- 5
          which.min( abs( disp_progress - pp ) )[1] }
        )
      }
      cat("|\n|")    
      
      ll <- matrix( 0 , nrow=nstud , ncol=3 )  
      vv <- 1
      for (pp in 1:ip){
        # pp <- 10  # parameter pp
        if (pp == 1){ vec <- 1:3 } else { vec <- 2:3 }
        pp.ind <- switch(irtmodel,
                         "2PL.groups"= which(est.slopegroups == sort(unique(est.slopegroups))[pp]),
                         #"GPCM.design"= pair.ind[pp,,drop=FALSE],
                         pp # default)
        )
        for (mm in vec){
          
          B0 <- B  
          if(irtmodel %in% c("2PL", "GPCM", "2PL.groups")){
            B0[,-1,][ pp.ind ] <- B0[,-1,][ pp.ind ] + mult[mm] * numdiff.parm  
          }
          if(irtmodel=="GPCM.design"){
            basispar0 <- basispar
            basispar0[ pp.ind ] <- basispar0[ pp.ind ] + mult[mm] * numdiff.parm  
            B0[,-1,] <- E%*%basispar0
          }
          
          ll0 <- 0
          # calculate probabilities
          res0 <- calc_prob.v5( iIndex=1:nitems , A=A , AXsi=AXsi , B=B0 , 
                                xsi=xsi , theta=theta , nnodes=nnodes , maxK=maxK )
          rprobs <- res0[["rprobs"]]
          # posterior distribution
          # calculate student's likelihood
          res0b <- calc_posterior.v2(rprobs=rprobs , gwt=res0a , resp=resp , 
                                     nitems=nitems , resp.ind.list=resp.ind.list , 
                                     normalization = FALSE , thetasamp.density = NULL , 
                                     snodes = 0 )$hwt
          # calculate individual log likelihood  																					
          ll[,mm] <- log( rowSums( res0b * thetawidth ) )
        }
        
        if(irtmodel %in% c("2PL", "GPCM", "2PL.groups")){
          se.B[,-1,][pp.ind] <- sqrt( - 1 / sum( tamobj$pweights * ( ll[,2] + ll[,3] - 2*ll[,1] )/h^2 ) )
        }
        if(irtmodel=="GPCM.design"){
          var.basispar[ pp.ind ] <- ( - 1 / sum( tamobj$pweights * ( ll[,2] + ll[,3] - 2*ll[,1] )/h^2 ) )
        } 
        
        if ( ( pp==disp_progress[vv] ) & ( vv <=10) ){
          cat("-") ; 
          flush.console() 
          vv <- vv+1
        }
      }
      
      if(irtmodel=="GPCM.design"){
        # se.B[,-1,] <- sqrt( diag((E%*%diag(var.basispar)%*t%E)) )
		se.B[,-1,] <- sqrt( diag(( E%*% diag(var.basispar) %*% t(E) )) )
      }
      
      cat("|\n")	
    }
    
    ##############################################
    # Regression parameters
    ##############################################	
    # create result object for item parameters
    se.beta <- 0*beta
    nreg <- nrow(beta)
    cat("Regression parameters\n|")
    #	cat(paste( rep("*" , nreg*ndim) , collapse=""))
    ip <- nreg*ndim
    VP <- min( ip, 10 )
    cat(paste( rep("*" , VP) , collapse=""))
    if (VP<10){ disp_progress <- 1:ip } else {
      disp_progress <- 100* ( 1:ip ) / (ip+1)
      disp_progress <- sapply( seq(5,95,10) , FUN = function(pp){ # pp <- 5
        which.min( abs( disp_progress - pp ) )[1] }
      )
    }
    
    cat("|\n|")    
    vv <- 1
    pp1 <- 1
    # compute response probabilities
    for (pp in 1:nreg){
      for (dd in 1:ndim){	
        #		ll <- matrix( 0 , nrow=nstud , ncol=3 )
        for (mm in 2:3){
          beta0 <- beta 
          beta0[ pp ,dd] <- beta0[pp,dd] + mult[mm] * numdiff.parm
          # compute likelihood
          # prior distribution for each student (normal density)
          res0a <- stud_prior.v2( theta=theta , Y=Y , beta=beta0 , variance=variance , 
                                  nstud=nstud , nnodes=nnodes , ndim=ndim, YSD=YSD ,
								  unidim_simplify=FALSE  )
          # calculate probabilities
          res0 <- calc_prob.v5( iIndex=1:nitems , A=A , AXsi=AXsi , B=B , 
                                xsi=xsi , theta=theta , nnodes=nnodes, maxK=maxK )
          rprobs <- res0[["rprobs"]]
          # posterior distribution
          res0b <- calc_posterior.v2( rprobs=rprobs , gwt=res0a , resp=resp , 
                                      nitems=nitems , resp.ind.list=resp.ind.list ,
                                      normalization=FALSE , thetasamp.density=NULL , 
                                      snodes=0)$hwt
          ll[,mm] <- log( rowSums( res0b *thetawidth ) )
        }
        se.beta[pp,dd] <- sqrt( - 1 / sum( tamobj$pweights * ( ll[,2] + ll[,3] - 2*ll[,1] )/h^2 ) )
        if ( ( pp1==disp_progress[vv] ) & ( vv <=10) ){
          cat("-") ; 
          flush.console() 
          vv <- vv+1
        }
        pp1 <- pp1+1 ;
      } }
    
    #-----------------------------------------------------------
    # handle fixed parameters
    if ( ! is.null( tamobj$xsi.fixed ) ){
      se.xsi[ tamobj$xsi.fixed[,1] ] <- 0
    }
    if ( ! is.null( tamobj$beta.fixed ) ){
      se.beta[ tamobj$beta.fixed[,1:2] ] <- 0
    }
    
    if ( ! is.null( tamobj$B.fixed ) ){
      se.B[ tamobj$B.fixed[,1:3] ] <- 0
    }
    if ( ! is.null( tamobj$Q ) & tamobj$ndim>1 ){
      Q.ind <- which(tamobj$Q==0, arr.ind=TRUE)
      Q.ind <- cbind(Q.ind[rep(1:nrow(Q.ind), maxK),], rep(1:maxK, each=nrow(Q.ind)))
      se.B[ Q.ind[,c(1,3,2)] ] <- 0
    }
    
    #-----------------------------------------------------------
    cat("|\n")	
    #***
    # ARb 2013-08-19: big fix for faceted models
    N1 <- nrow(tamobj$item)
    if (N1 != length(xsi) ){
      xsi <- data.frame( "item" = rownames(tamobj$xsi) , "N"=NA ,
                         "est" = xsi , "se" = se.xsi )
    } else {
      #		xsi <- data.frame( tamobj$item[,1:2] , 
      #					"est" = xsi , "se" = se.xsi )		
      xsi <- data.frame( "item" = rownames(tamobj$xsi)  , 
                         "est" = xsi , "se" = se.xsi )		
      
    }
    
    beta <- data.frame( "beta" = beta , "se" = se.beta )
    colnames(beta) <- c( paste("est.Dim" , 1:ndim , sep="")	, paste("se.Dim" , 1:ndim , sep="")	)
    
    B.out <- data.frame( "item" = dimnames(B)[[1]] )
    for (kk in 1:(maxK-1)){ # kk <- 1
      for (dd in 1:ndim){
        B.out[ , paste0("B.Cat" , kk,".Dim",dd) ] <- B[,kk+1,dd]
        B.out[ , paste0("se.B.Cat" , kk,".Dim",dd) ] <- se.B[,kk+1,dd]
      }
    }				
    
    
    flush.console()
    res <- list( "xsi" = xsi , "beta" = beta, "B"=B.out )
    if(irtmodel=="GPCM.design"){
      basispar.res <- data.frame("basispar"=1:length(basispar),
                                 "gamma"=basispar, "se"=sqrt(var.basispar) )
      res <- c(res, list("basispar"=basispar.res))
    } 
    
    return(res)
  }
