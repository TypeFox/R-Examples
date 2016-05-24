####################################################
# create E matrix
.mml.3pl.create.E <- function( resp , E , Q , gammaslope.des ,
		Q.fixed = NULL ){  
  Qdes <- NULL		
  gammaslope.fixed <- NULL
  if ( is.null(E) ){
	maxKi <- apply( resp , 2 , max , na.rm=TRUE )
	I <- ncol(resp)
	if ( is.null(Q) ){ Q <- matrix( 1 , nrow=I , ncol=1 ) }
	D <- ncol(Q)
	maxK <- max( maxKi ) + 1
    if ( gammaslope.des == "2PL" ){
		Ngam <- sum( abs(Q) > 0 )
								}
    if ( gammaslope.des == "2PLcat" ){
		Ngam <- sum( rowSums(( abs(Q) > 0 )  )*maxKi )
								}																					
    ng <- 1
	kk <- 1
    vv <- 1	
	Qdes <- matrix( 0 , nrow=maxK*I*D , ncol=5 )	
	colnames(Qdes) <- c("gammapar" , "item" , "dim" , "category" , "Qval")
	for (ii in 1:I){
		for (dd in 1:D){
		 if ( Q[ii,dd] != 0 ){
		  for (kk in 1:maxKi[ii]){
		    Qdes[vv,1] <- ng
			Qdes[vv,2:3] <- c(ii,dd)
			Qdes[vv,4] <- kk
			if (  gammaslope.des=="2PL" ){
					Qdes[vv,5] <- Q[ii,dd]*kk
								}
			if (  gammaslope.des=="2PLcat" ){
					Qdes[vv,5] <- Q[ii,dd]
								}								
			vv <- vv + 1
			if ( ( kk==maxKi[ii] ) & ( gammaslope.des=="2PL") ){
						ng <- ng + 1
							}
			if (  gammaslope.des=="2PLcat" ){
						ng <- ng + 1
							}							
							}
							}
					   }  # end dd   					   					   
				}  # end ii
	  Qdes <- as.data.frame( Qdes[ 1:(vv-1) , ] )
	  Ngam <- max( Qdes$gammapar )
	  gammaslope.fixed <- NULL
	  # fixed gammaslope parameters
	  Qdes$gamma.fixed <- NA	  	  	  
	  if ( ! is.null(Q.fixed) ){
	      for (dd in 1:D){
			# dd <- 1
			Q1 <- Q.fixed[ , dd ]
			ind.dd <- which( ! is.na( Q1) )
			if ( length(ind.dd) > 0 ){
			    I1 <- length(ind.dd)
				for (ii in 1:I1){
						i2 <- which( ( Qdes$item == ind.dd[ii] ) & 
								  (	Qdes$dim == dd )	)
						Qdes[i2,"gamma.fixed"] <- Q1[ ind.dd[ii] ]	
								  }	
								}  # end if len(ind.dd) > 0	  
					} # end dd
		  gam1 <- stats::aggregate( Qdes$gamma.fixed , list(Qdes$gammapar) , mean )
		  gam1 <- stats::na.omit(gam1)
		  gammaslope.fixed <- gam1[ , c(1,2) ]
		  colnames(gammaslope.fixed) <- NULL
				}  # end ! is.null(Q.fixed)


	#****  
	# c("gammapar" , "item" , "dim" , "category" , "Qval")
	E <- array( 0 , dim=c(I,maxK , D , Ngam ) )
	for (ee in 1:(nrow(Qdes)) ){
		# ee <- 1
		Qdes.ii <- Qdes[ee,]		
		E[ Qdes.ii$item , Qdes.ii$category + 1, Qdes.ii$dim , 
				Qdes.ii$gammapar ] <- Qdes.ii$Qval
					}
			}
	res <- list(E=E , Qdes=Qdes , gammaslope.fixed=gammaslope.fixed )
	return(res)
	}
####################################################################	


#####################################################################
# compute F design matrix
.mml.3pl.computeFdes <- function( E , gammaslope , theta ){	  
	   #*********
	   #  SETUP: 3PL Structured latent class analysis
       #  design matrix E ... E[ii, cc , dd , pp ] with loading parameters \gamma
       #  B[ ii , cc , dd ] = E[ii,cc,dd,] %*% gamma[]
       #  f[ii,cc,tt,kk] = \sum_dd E[ii,cc,dd,kk] * theta[jj,dd]
       #*********	   
	   Ngam <- length(gammaslope)
	   dimE <- dim(E)
	   I <- dimE[1]
	   maxK <- dimE[2]
	   D <- dimE[3]
	   TP <- nrow(theta)	   	   
	   Fdes <- array( 0 , dim=c(I , maxK , TP , Ngam ) )
	   ttheta <- t(theta)

       for (pp in 1:Ngam){
	   for (cc in 1:maxK){
         if (D==1){		
			Fdes[,cc,,pp] <- E[ , cc , , pp , drop=FALSE] %*% ttheta
							}
         if (D >1){
			E.cc.pp <- E[ , cc , , pp  ]
		  #***** Fdes <- array( 0 , dim=c(I , maxK , TP , Ngam ) )
			Fdes[,cc,,pp] <- E.cc.pp %*% ttheta
			# Fdes[,cc,,pp] <- E.cc.pp %*% ttheta
							}							
						}
					}
					
      return(Fdes)
		}


####################################################
# create a dummy A matrix
.mml.3pl.create.notA <- function( E , notA ){  
    dimE <- dim(E)
	A <- array( 0 , dim=c(dimE[1] , dimE[2] , 2 ) )
    xsi.fixed <- cbind( c(1,2) , 0 )
	res <- list("A"=A , "xsi.fixed" = xsi.fixed )
	return(res)
	}
####################################################################	
		