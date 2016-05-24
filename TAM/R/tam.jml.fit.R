


tam.jml.fit <-
  function ( tamobj ){
    #####################################################
    # INPUT:
    # tamobj ... result from tam.jml analysis
    ####################################################
    
    resp<-tamobj$resp
    resp.ind<-tamobj$resp.ind
    
    A<-tamobj$A
    B<-tamobj$B
    nstud<-tamobj$nstud
    nitems<-tamobj$nitems
    maxK<-tamobj$maxK
    ItemScore<-tamobj$ItemScore
    theta<-tamobj$theta
    xsi<-tamobj$xsi
    #     Msteps=tamobj
    pweightsM<-tamobj$pweights
    #     est.xsi.index=tamobj
    
    AXsi <- matrix(0, nrow=nitems, ncol=maxK) 
    BB <- array (0, dim=c(nitems,maxK))
    B_Sq <- array(0,dim=c(nstud, nitems))
    C4 <- array(0,dim=c(nitems,nstud))
    for (k in 1:maxK){ 
      AXsi[,k] <- A[,k,] %*% xsi
    } 
    B.0 <- B
    B.0[ is.na(B.0) ] <- 0					
    B1 <- B.0[,,1]					
    BB <- B1^2
    theta.unique <- unique( theta )
    NU <- length(theta.unique)
    B_bari <- array(0,dim=c(NU, nitems))
    BB_bari <- array(0, dim=c(NU, nitems))  
    res <- calc_prob.v5(iIndex = 1:nitems , A , AXsi , 
                        B , xsi , theta= matrix( theta.unique , nrow=NU , ncol=1) , 
                        NU, maxK , recalc=FALSE )        
    rprobs <- res[["rprobs"]] 
    #  rprobs <- rprobs[ , , match( theta[,1] , theta.unique)  ]
    rprobs[ is.na( rprobs) ] <- 0 
    # rprobs [ nitems , maxK , nstud ]
    # B1	[ nitems , maxK ]
    # BB	[ nitems , maxK ]
    # B_bari, BB_bari 	[nstud , nitems ]
    
    for (kk in 1:maxK){ 
      B_bari <- B_bari + t( B1[,kk]*rprobs[,kk,] )
      BB_bari <- BB_bari + t( BB[,kk] * rprobs[ , kk , ] )
    }
    ind.theta <- match( theta , theta.unique)				 
    rprobs <- rprobs[ , ,  ind.theta ]				 
    B_bari <- B_bari[ ind.theta , ] 
    BB_bari <- BB_bari[ ind.theta , ] 				 
    B_bari <- B_bari * resp.ind 		
    BB_bari <- BB_bari  * resp.ind
    #  B_bari <- sapply(1:nitems, function(i) colSums(B1[i,] * rprobs[i,,] , na.rm = TRUE)) * resp.ind
    #  BB_bari <- sapply(1:nitems, function(i) colSums(BB[i,] * rprobs[i,,] , na.rm = TRUE))  * resp.ind
    B_var <- BB_bari - (B_bari^2)
    z_sq <- (resp - B_bari)^2/B_var
    zsd <- sd(as.numeric(z_sq),na.rm=TRUE)
    z_sq[z_sq > 10*zsd] <-  10*zsd   #Trim extreme values
    B_bariM <- aperm(outer(B_bari,rep(1,maxK)),c(2,3,1))
    B1M <- outer(B1,rep(1,nstud))
    # B1M	[ nitems , maxK , nstud ]
    # B_bariM [ nitems , maxK , nstud ]
    #  C4 <- sapply(1:nitems, function(i) colSums((B1M[i,,]-B_bariM[i,,])^4 * rprobs[i,,] , na.rm = TRUE)) * resp.ind  
    # rprobs[i,,]  [ maxK , nstud ]
    # C4	[ nstud , nitems ]
    # B1M	[ nitems , maxK , nstud ]
    for (kk in 1:maxK){
      C4 <- C4 + ( B1M[,kk,] - B_bariM[,kk,] )^4 * rprobs[,kk,]
    }
    
    C4 <- t(C4) * resp.ind			 
    #  outfitPerson <- apply(z_sq, 1, mean, na.rm = TRUE)
    outfitPerson <- rowMeans( z_sq , na.rm=TRUE )
    #  outfitItem <- apply(z_sq * pweightsM, 2, mean, na.rm = TRUE)
    outfitItem <- colMeans(z_sq * pweightsM, na.rm = TRUE)
    
    infitPerson <- rowSums((resp - B_bari)^2, na.rm=TRUE)/rowSums(B_var, na.rm=TRUE)
    infitItem <- colSums((resp - B_bari)^2 * pweightsM, na.rm=TRUE)/colSums(B_var * pweightsM, na.rm=TRUE) 
    
    z_sq.ind <- !is.na(z_sq)
    #  var_outfit <- rowSums(C4/(B_var^2), na.rm=TRUE)/(nitems^2) - 1/nitems
    var_outfit <- rowSums(C4/(B_var^2), na.rm=TRUE)/(rowSums(z_sq.ind)^2) - 1/rowSums(z_sq.ind)
    outfitPerson_t <- (outfitPerson^(1/3) - 1) * (3/sqrt(var_outfit)) + sqrt(var_outfit)/3
    
    #  var_outfit <- colSums(C4/(B_var^2), na.rm=TRUE)/(nstud^2) - 1/nstud
    var_outfit <- colSums(C4/(B_var^2), na.rm=TRUE)/(colSums(z_sq.ind)^2) - 1/colSums(z_sq.ind)
    outfitItem_t <- (outfitItem^(1/3) - 1) * (3/sqrt(var_outfit)) + sqrt(var_outfit)/3
    
    var_infit <- rowSums(C4-B_var^2, na.rm=TRUE)/((rowSums(B_var, na.rm=TRUE))^2)
    infitPerson_t <- (infitPerson^(1/3) - 1) * (3/sqrt(var_infit)) + sqrt(var_infit)/3  
    
    var_infit <- colSums(C4-B_var^2, na.rm=TRUE)/((colSums(B_var, na.rm=TRUE))^2)
    infitItem_t <- (infitItem^(1/3) - 1) * (3/sqrt(var_infit)) + sqrt(var_infit)/3
    
    # collect item statistics
    fit.item <- data.frame("item" = colnames(tamobj$resp) , "outfitItem" = outfitItem , 
                           "outfitItem_t" = outfitItem_t, "infitItem" = infitItem , 
                           "infitItem_t" = infitItem_t)
    
    fit.person <- data.frame("outfitPerson" = outfitPerson , 
                             "outfitPerson_t" = outfitPerson_t, "infitPerson" = infitPerson , 
                             "infitPerson_t" = infitPerson_t)
    
    return( list("fit.item"=fit.item, "fit.person"=fit.person) )  
  }
