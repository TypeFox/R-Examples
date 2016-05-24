tam.fit <- function( tamobj, ... ){
  if(class(tamobj) == "tam.mml"){
    res <- tam.mml.fit( tamobj, ...)
  }
  if(class(tamobj) == "tam.jml"){
    res <- tam.jml.fit( tamobj, ...)
  }
  class(res) <- "tam.fit"
  return(res)
}

tam.mml.fit <-
  function( tamobj, FitMatrix=NULL , Nsimul=NULL , progress = TRUE ,
     useRcpp=TRUE , seed=123 , fit.facets=TRUE ){
    #####################################################
    # INPUT:
    # tamobj ... result from tam analysis
    # FitMatrix is the fit design matrix. If it's NULL, then we will use the A matrix.
    # MW: We need to check whether this works when there are missing responses.
    # progress ... fit progress
    ####################################################
    
	if ( ! is.null(seed) ){
	   set.seed(seed)  }
	
	if ( is.null(Nsimul) ){
			nstud <- tamobj$nstud
			Nsimul <- 5
			if ( nstud < 3000 ){ Nsimul <- 15 }			
			if ( nstud < 1000 ){ Nsimul <- 40 }
			if ( nstud < 400 ){ Nsimul <- 100 }						
						}
	 if (progress){
	    cat( paste0("Item fit calculation based on " , Nsimul , " simulations\n") )
					}
		
    resp <- tamobj$resp  
    rprobs <- tamobj$rprobs    	
	
    indexIP.list <- tamobj$indexIP.list
    hwt <- tamobj$hwt
    resp.ind <- tamobj$resp.ind
    nnodes <- nrow(tamobj$theta)
    pweights <- tamobj$pweights
    nstud <- tamobj$nstud
    nitems <- tamobj$nitems
    maxK <- tamobj$maxK
    
#    if( "formulaA" %in% names(tamobj) ) warning("tam.fit is experimental for objects obtained by tam.mml.mfr.")
    if ( is.null(FitMatrix) ) {
      FitMatrix <- tamobj$A 
               }
	if( ( "formulaA" %in% names(tamobj) ) & ( fit.facets )){
		FitMatrix2 <- FitMatrix.facets(tamobj)
		F1 <- dim(FitMatrix)[[3]]
		F2 <- dim(FitMatrix2)[[3]]
		FitMatrix3 <- array( 0 , dim = c( dim(FitMatrix)[[1]] ,
				dim(FitMatrix)[[2]] , F1+F2) )
		dimnames(FitMatrix3)[[1]] <- dimnames(FitMatrix)[[1]]
		dimnames(FitMatrix3)[[2]] <- dimnames(FitMatrix)[[2]]
		dimnames(FitMatrix3)[[3]] <- c( dimnames(FitMatrix)[[3]] ,
				dimnames(FitMatrix2)[[3]] )
		FitMatrix3[ ,,1:F1] <- FitMatrix
		FitMatrix3[ ,,(F1+1):(F1+F2)] <- FitMatrix2
		FitMatrix3 <- FitMatrix3[ ,, paste(tamobj$xsi.facets$parameter) ]
		FitMatrix <- FitMatrix3
					}
	
					
    col.index <- rep( 1:nitems , each = maxK )
    cResp <- resp[ , col.index  ]*resp.ind[ , col.index ]
    cResp <- 1 * t( t(cResp) == rep(0:(maxK-1), nitems) )
    cF <- t( matrix( aperm( FitMatrix , c(2,1,3) ) , nrow = dim(FitMatrix)[3] , byrow = TRUE ) )
    cF[is.na(cF)] <- 0
    cResp[ is.na(cResp) ] <- 0
    # sufficient statistics by person by parameter, for Fit design matrix.
    ParamScore <- (cResp %*% cF)
    
    np <- dim(FitMatrix)[3]  #The parameter dimension
    indexIP <- colSums( aperm( FitMatrix, c(2,1,3) ) != 0, na.rm = TRUE )
    
    # define list of elements for item parameters
    indexIP.list <- list( 1:np )
    for ( kk in 1:np ){ 
      indexIP.list[[kk]] <- which( indexIP[,kk] > 0 )
    }
    
    Outfit <- rep(0,np)
    Infit <- rep(0,np)
    Outfit_t <- rep(0,np)
    Infit_t <- rep(0,np)
    if (progress){
      cat(paste( "|" , paste(rep("*" , 10 ), collapse="") , "|\n|" ,sep="") )
	  prbar <- round( seq( 1 , np , len = 10 ) )
    }

	N <- nrow(hwt)
    rn1M <- matrix( stats::runif(N*Nsimul) , nrow=N , ncol= Nsimul )        

    
    for (p in 1:np) {
      ip <- indexIP.list[[p]]
      xbari <- sapply( ip, function(i) colSums(FitMatrix[i,,p] * rprobs[i,,] , na.rm = TRUE ))
      #... TK: multiple category option -> na.rm = TRUE
      xxfi <- sapply( ip, function(i) colSums(FitMatrix[i,,p]^2 * rprobs[i,,] , na.rm = TRUE ))
      vari <- xxfi - xbari^2
      xxxfi <- sapply( ip, function(i) colSums((FitMatrix[i,,p])^3 * rprobs[i,,] , na.rm = TRUE ))
      xxxxfi <- sapply( ip, function(i) colSums(FitMatrix[i,,p]^4 * rprobs[i,,] , na.rm = TRUE ))
      C4i <- xxxxfi - 4*xbari*xxxfi + 6*(xbari^2)*xxfi - 3*(xbari^4)
      Vz2i <- C4i - vari^2
      Uz2i <- C4i/(vari^2) - 1
      
      # xbar <- resp.ind[,ip]%*%t( xbari )      
	  xbar <- tcrossprod( resp.ind[,ip] ,  xbari )
	  # var1 <- resp.ind[,ip]%*%t( vari )
	  var1 <- tcrossprod( resp.ind[,ip] ,  vari )	  	  
      # Vz2 <- resp.ind[,ip]%*%t( Vz2i )
	  Vz2 <- tcrossprod( resp.ind[,ip] ,  Vz2i )	  	  
      # Uz2 <- resp.ind[,ip]%*%t( Uz2i )
	  Uz2 <- tcrossprod( resp.ind[,ip] ,  Uz2i )	
      
      Ax <- matrix(rep(ParamScore[,p],nnodes),nrow=nstud, ncol=nnodes)
      
      
      #	c_hwt0 <- aperm(apply (hwt,1,cumsum),c(2,1))
      c_hwt <- rowCumsums.TAM(hwt)
      
	  Outfit_SIM <- Infit_SIM <- rep(NA, Nsimul )
	  Infit_t_SIM <- Outfit_t_SIM <- rep(NA,Nsimul )

      #calculate number of students per item parameter	
#      nstud.ip <- sum( rowMeans( resp.ind[ , ip , drop=FALSE],na.rm = TRUE ), na.rm=TRUE )
      if (TRUE){
		nstud.ip <- rowSums( resp.ind[ , ip , drop=FALSE],na.rm = TRUE )
		nstud.ip <- sum( 1*(nstud.ip > 0))
				}
	  
	  
  
	  if ( ! useRcpp ){
		  
		  #*******
		  # start simulation
		  for ( hh in 1:Nsimul ){
				#hh <- 1
			  rn1 <- rn1M[,hh]
			  nthetal <- rep(1,ncol(c_hwt))      
			  j <- rowSums( c_hwt < rn1)
			  j[ j == 0 ] <- 1
			  NW <- ncol(c_hwt)
			  j <- j + 1
			  j[ j > NW ] <- NW	
			  s <- cbind(seq(1:N),j)
			  wt_numer <- ( Ax[s] - xbar[s] )^2
			  wt_denom <- var1[s]
			  z2 <- wt_numer/wt_denom
			  varz2 <- Uz2[s]
			  wt_var <- Vz2[s]	  	 
				  #**********************************
				  # simulated quantities:
				  # wt_numer
				  # wt_denom
				  # z2
				  # varz2
				  # wt_var
				  #**********************************	  
			  
			  # end simulation
			  #****************************************
			  
			  #Outfit MNSQ (unweighted fit)
			  
			  #	z2[z2 > 10*sd(z2)] <- 10*sd(z2)  #Trim extreme values
			  
			  Outfit[p] <- sum( z2*pweights, na.rm = TRUE  ) / nstud.ip
			  Outfit_SIM[hh] <- Outfit[p]
			  
			  #Infit MNSQ (weighted fit)
			  Infit[p] <- sum( wt_numer*pweights,na.rm = TRUE )/sum(wt_denom*pweights,na.rm = TRUE  )
			  Infit_SIM[hh] <- Infit[p]
			  
			  #Infit t
			  vf <- sum(wt_var*pweights,na.rm = TRUE )/(sum(wt_denom*pweights,na.rm = TRUE)^2 ) 
			  Infit_t[p] <- (Infit[p]^(1/3)-1) * 3/sqrt(vf) + sqrt(vf)/3  
			  Infit_t_SIM[hh] <- Infit_t[p]
			  
			  #Outfit t
			  vf2 <- sum(varz2*pweights,na.rm = TRUE )/(nstud.ip^2)
			  Outfit_t[p] <- (Outfit[p]^(1/3)-1) * 3/sqrt(vf2) + sqrt(vf2)/3
			  Outfit_t_SIM[hh] <- Outfit_t[p]
		  }
				}
	
      # calculate fit in Rcpp	
	  if ( useRcpp ){
		res0 <- .Call( "tam_fit_simul" , 
					rn1M , c_hwt , Ax , xbar , var1 , Uz2 , Vz2 , nstud.ip ,
					pweights , 
					PACKAGE="TAM" )
		Outfit_SIM <- res0$Outfit_SIM
		Infit_SIM <- res0$Infit_SIM
		Infit_t_SIM <- res0$Infit_t_SIM
		Outfit_t_SIM <- res0$Outfit_t_SIM
					}
				
	  Outfit[p] <- mean( Outfit_SIM)
	  Infit[p] <- mean( Infit_SIM)
	  Infit_t[p] <- mean( Infit_t_SIM )
	  Outfit_t[p] <- mean( Outfit_t_SIM )	  
	  
	  
	if (progress){
		if ( p %in% prbar ){
		       cat("-") ; utils::flush.console()
							}
						}
	  
    }
    if (progress){ cat("|\n") ; utils::flush.console() }
    res <- data.frame(
	    "parameter" = dimnames(FitMatrix)[[3]] ,
					  "Outfit" = Outfit , 
                      "Outfit_t" = Outfit_t, 
					  "Outfit_p" = stats::pnorm(-abs(Outfit_t)) /2 ,
					  "Outfit_pholm" =NA , 
					  "Infit" = Infit , 
                      "Infit_t" = Infit_t,  
					  "Infit_p" = stats::pnorm(-abs(Infit_t)) /2 ,
					  "Infit_pholm" =NA 
							)
	res$Outfit_pholm <- stats::p.adjust( res$Outfit_p , method="holm")
	res$Infit_pholm <- stats::p.adjust( res$Infit_p , method="holm")
    #data.frame( "Outfit" = round(Outfit,2) , "Outfit_t" = round(Outfit_t,1), "Infit" = round(Infit,2), Infit_t = round(Infit_t,1) )    
    res <- list( "itemfit" = res )
    class(res) <- "tam.fit"	
	return(res)
  }


########################################
# create FitMatrix for facets

FitMatrix.facets <- function(tamobj){
    xsi.constraints <- tamobj$xsi.constr$xsi.constraints
    A <- tamobj$A
    XX <- nrow(xsi.constraints)
    A.ext <- array( 0 , dim=c( dim(A)[1],dim(A)[2] , XX ) )
    dimnames(A.ext)[[1]] <- dimnames(A)[[1]]
    dimnames(A.ext)[[2]] <- dimnames(A)[[2]]
    dimnames(A.ext)[[3]] <- rownames(xsi.constraints)    
    for (pp in 1:XX){
        # pp <- 2
        xx <- xsi.constraints[pp,]
        VV <- dim(A)[3]
        Axx <- 0
        for (vv in 1:VV){
            Axx <- Axx + xx[vv] * A[,,vv]
                        }
        A.ext[,,pp] <- Axx
                }
    return(A.ext)
        }  
 
 
