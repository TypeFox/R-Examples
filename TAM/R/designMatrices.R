designMatrices <-
  function( modeltype = c( "PCM" , "RSM" ) , 
            maxKi = NULL , resp = resp , ndim = 1 ,
            A = NULL , B = NULL , Q = NULL , R = NULL , 
			constraint="cases" , ... ){
    
    modeltype <- match.arg(modeltype)
#a0 <- Sys.time();
    I <- ncol(resp)
	if ( ! is.null(A) ){
		constraint <- "cases" 
					}
	
    A.draft <- A	# if ! is.null(A), it is necessary
	
    if( is.null(maxKi) ){
      if( !is.null(resp) ){
        resp[is.na(resp)] <- 0
        maxKi <- apply( resp , 2 , max , na.rm=TRUE )
      } else 
        #... TK: 24.07.2012 -- check     
        if( !is.null(A) ){   
          np <- ncol(A)
          maxKi <- -colSums(A)
          maxKi <- maxKi[ - (which( (maxKi - 1) > 0) + maxKi[ which( (maxKi - 1) > 0) ]-1) ]
        } else return( warning("Not enough information to generate design matrices") )
    }
	# stop processing if there are items with a maximum score of 0
	i11 <- names(maxKi)[ maxKi == 0 ]
    if ( length(i11) > 0 ){
#		stop( cat( "Items with maximum score of 0:" , paste(i11 , collapse=" " ) ) )
 		cat( "Items with maximum score of 0:" , paste(i11 , collapse=" " ) ) 
					}
	
    nI <- length(maxKi)
    maxK <- max(maxKi)
    item <- rep( 1:nI , maxKi+1 )
#cat("g100"); a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1      
    if(modeltype %in%  c("PCM","RSM") ){
      
      cat <- unlist( lapply ( maxKi, seq , from=0 ) )
      np <- sum(maxKi)
      repnP <- cat[ cat != 0 ]
      revCat <- unlist( lapply ( maxKi, seq , to=1 ) )
      
      # Q Matrix
      if( is.null(Q) ){
#        if( ndim > 1 ) warning("random q matrix")
        Q.draft <- matrix( 0 , nrow = nI , ncol = ndim )
        Q.draft[ cbind( 1:nI , sample(1:ndim, nI, replace=TRUE) ) ] <- 1 
      }else{
        Q.draft <- Q
      } 
      ndim <- dim(Q.draft)[2]
#cat("g150"); a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1            
      # B Matrix
      if( is.null(B) ){      
        B.draft <- array( 0 , dim = c(nI, maxK+1 , ndim) , 
                          dimnames = list(
                            #							paste( "Item", ifelse( 1:nI<10 , "0" , "" ), 1:nI , sep = "" ) , 
                            colnames(resp)  , 
                            paste( "Cat", 0:maxK , sep = "" ) , 
                            paste( "Dim" , ifelse( (1:ndim) < 10 , "0" , "" ), 1:ndim , sep = "")))
        
        for(dd in ndim){
          ind <- cbind( rep( 1:nI , maxKi+1 ) ,
                        cat + 1 , 
                        rep( dd , sum(maxKi+1) ) )
          B.draft[ind] <- cat*rep(Q.draft[,dd], maxKi+1)
        }
        
      } else { 
        B.draft <- B 
      }
      
      if ( ! is.null(Q) ){
        for (dd in 1:dim(Q)[2] ){
          for (zz in 1:dim(B.draft)[2] ){ 
            B.draft[ , zz , dd ] <- (zz-1)*Q[,dd]
          } 
        }
      }
#cat("g200"); a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1         
      # A Matrix
      if( is.null(A) || length( dim(A) ) < 3 ){
        A.draft <- array(, dim = c( nI, maxK+1 , np ) ,
                         dimnames = list( paste( "Item", ifelse( 1:nI<10 , "0" , "" ), 1:nI , sep = "" ) , 
                                          paste( "Category", 0:maxK , sep = "" ) , 
                                          paste( "Xsi", ifelse( (1:np) < 10 , "0" , "" ), 1:np , sep = "" ) ))
        
        ind.0 <- cbind( rep( item , np ) ,
                        rep( cat + 1 , np ) ,
                        rep( 1 : np , each = nI+np ) )
        
        i_pars <- cbind( pars_start <- rep(c(0, cumsum(maxKi))+1, c(maxKi, 0)),
                         pars_end <- pars_start+repnP-1 )
        
        ind.1 <- cbind( "item" = rep( item[ cat != 0 ] , repnP ), 
                        "category" = rep( repnP+1 , repnP) ,
                        "xsi" = unlist( apply( i_pars , 1 , 
								function(i_par) seq(from = i_par[1], to = i_par[2]) ) )
					)
        
        A.draft[ ind.0 ] <- 0
        A.draft[ ind.1 ] <- -1
        
        
        # item labels for Partial Credit Model
        l0 <- unlist(sapply( maxKi , FUN = function(cc){ seq( 1, cc)  } ))
        l1 <- paste( rep( colnames(resp) , maxKi) , "_Cat" , l0 , sep="" )
        dimnames(A.draft)[[3]] <- l1
        
        # item labels for Rasch model
        if (maxK == 1 ){
          dimnames(A.draft)[[3]] <- colnames(resp)
        }
        
      }    
      
    }

	#*****************************
	# constraint = "items"

    if ( constraint == "items"){ 
		unidim <- is.null(Q)
	    if ( ! is.null(Q) ){
		      unidim <- ncol(Q) == 1
							}	
		#***** dichotomous items unidimensionality
		if ( ( maxK==1 ) & (unidim ) ){	
			I <- dim(A.draft)[3]
			A.draft[ I , 2 , ] <- +1
			A.draft <- A.draft[ ,, seq(1,I-1) ]
					}
		#***** dichotomous items multidimensionality
		if ( ( maxK==1 ) & (! unidim ) ){	
			rem.pars <- NULL
			D <- ncol(Q)
			for ( dd in 1:D){
#				dd <- 1
				ind.dd <- which(Q[,dd] != 0 )
				ndd <- length(ind.dd)
				rem.pars <- c(rem.pars , ind.dd[ndd] )
				A.draft[ ind.dd[ndd] , 2 , ind.dd ] <- 1
							}
			A.draft <- A.draft[ , , - rem.pars ]
							}
				}
	
#cat("g300"); a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1       
    if(modeltype == "RSM"){
      if( is.null(A) ) 
#        return( warning("Not enough information to generate design matrices") )
#          nP <- sum(maxKi) + length(rater)
	Nxsi <- I + maxK - 1
	Kitem <- maxKi+1
	A <- array( 0 , dim=c( I , maxK+1 , Nxsi ) )
	vv <- 1
	for (ii in 1:I){
		A[ ii , 2:Kitem[ii] , vv ] <- - ( 2:Kitem[ii] - 1 )
		if ( Kitem[ii] <= maxK ){
			A[ ii , ( Kitem[ii] + 1 ):(maxK+1) ,  ] <- NA                
							}
		vv <- vv+1
					}
	Kitem2 <- maxK+1 + 0*Kitem					
    for (ii in 1:I){					
	  if ( Kitem2[ii] > 2 ){
		for (kk in 1:(Kitem2[ii] - 2) ){
			A[ ii , 1 + ( kk:(Kitem2[ii]-2) ) , I+kk ] <- - 1
						}
					}
				}
				
	dimnames(A)[[1]] <- colnames(resp)
	vars <- colnames(resp)
	vars <- c(vars , paste0( "Cat"  , 1:(maxK-1) ) )
	dimnames(A)[[3]] <- vars			
	A.draft <- A


	
    }
    
    flatA <- t( matrix( aperm( A.draft , c(2,1,3) ) , nrow = dim(A.draft)[3] , byrow = TRUE ) )
    colnames(flatA) <- dimnames(A.draft)[[3]]
    
    flatB <- t( matrix( aperm( B.draft , c(2,1,3) ) , nrow = dim(B.draft)[3] , byrow = TRUE ) )
    colnames(flatB) <- dimnames(B.draft)[[3]]
    rownames(flatB) <- rownames(flatA) <- t(outer(dimnames(B.draft)[[1]] , 
					dimnames(B.draft)[[2]] , paste , sep ="."))
    
    out <- list( "item" = item , "maxKi" = maxKi , "cat" = cat , 
                 "A" = A.draft , "flatA" = flatA , "B" = B.draft , 
                 "flatB" = flatB , "Q" = Q , "R" = R )
    class(out) <- "designMatrices"
    return(out)
  }
