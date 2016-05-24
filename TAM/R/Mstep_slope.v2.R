


############################################################################
############################################################################
Mstep_slope.v2 <-
  function (B_orig, B, B_obs, B.fixed , max.increment, nitems, A, 
            AXsi, xsi, theta, nnodes, maxK, itemwt, Msteps, ndim, convM ,
            irtmodel , progress , est.slopegroups , E , basispar , se.B ,
			equal.categ 
  ){     
    # begin loop of slope parameters (items and categories, where B coeff is not zero)
    #        old_increment <- rep(max.increment, nitems )

    old_increment <- array( max.increment , dim= c(nitems , maxK , ndim) )

    xbar2 <- xxf <- xbar <- array(0,dim=c( nitems , maxK ) )
    converge <- FALSE
    Biter <- 1
    mK <- 1:maxK		
    items.temp <- 1:nitems
    items.conv <- NULL
		
    if (irtmodel == "GPCM"){
      old_increment.temp <- matrix( .3 , nrow=nitems , ncol= ndim )	
    }
    if (irtmodel == "GPCM.design" ){
      Nlambda <- ncol(E)	# number of lambda parameters							
      old_increment.temp <- matrix( .3 , nrow=Nlambda , ncol= ndim )				
    }	
    if (irtmodel == "2PL.groups"){
      ES <- length( unique( est.slopegroups) )
      old_increment.temp <- array( .3 , dim=c(ES , maxK ,ndim )	)
    }							
    while (!converge & ( Biter <= Msteps ) ) {  
      #compute expectation
      res.p <- calc_prob.v5( iIndex= items.temp , A=A , AXsi=AXsi , B=B, 
                             xsi=xsi , theta=theta , nnodes=nnodes, maxK=maxK)    			
      rprobs <- res.p[["rprobs"]] 
      #########################################
      ######     D I M E N S I O N S     ######
      for (dd in 1:ndim){
        #***
        if ( irtmodel %in% c("GPCM","GPCM.design") ){ 
          xtemp <- matrix(0 , length(items.temp) , nrow(theta) )
        }
        if ( irtmodel == "2PL.groups"){ 
          xtemp <- array(0 , dim=c(length(items.temp) , nrow(theta) , maxK) )
        }									
        ###################################
        ##### C A T E G O R I E S #########	
		if ( ! equal.categ ){ rprobs[ is.na(rprobs ) ] <- 0 }
        for ( k in 1:maxK ){	
          #			k <- 1	
		  # SIMPLIFY HERE!!!!
		  #    transpose and theta[,dd]!!
		  # t(A) * t(B) = t( B * A )
#          xbar[items.temp,k] <- t(theta[,dd]) %*% t( rprobs[,k,]*t(itemwt[,items.temp]) )
#		  rpr.it <- rprobs[,k,] * t( itemwt[,items.temp] )
#		  rpr.it <- t( rprobs[,k,] ) * itemwt[,items.temp,drop=FALSE]
		  rprobs.k <- matrix( rprobs[,k,] , nrow=length(items.temp) , ncol=nrow(theta) )
#		  rpr.it <- t( rprobs[,k,] ) * itemwt[,items.temp]
 		  rpr.it <- t( rprobs.k ) * itemwt[,items.temp]
		  ttheta.dd2 <- t( theta[,dd,drop=FALSE]^2)
		  ttheta.dd <- t( theta[,dd,drop=FALSE] )
#		  xbar[items.temp,k] <- t( theta[,dd] ) %*% rpr.it
		  xbar[items.temp,k] <- ttheta.dd %*% rpr.it

#          xxf[items.temp,k] <- t( theta[,dd]^2) %*% t( rprobs[,k,]*t(itemwt[,items.temp]) )
#          xxf[items.temp,k] <- t( theta[,dd]^2) %*% rpr.it
          xxf[items.temp,k] <- ttheta.dd2 %*% rpr.it
          if ( irtmodel == "2PL" ){ 	# begin usual 2PL
#			xbar2[items.temp,k] <- 
#              t( theta[,dd]^2 ) %*% t( rprobs[,k,]^2*t(itemwt[,items.temp]) )
#			xbar2[items.temp,k] <- 
#              ttheta.dd2 %*% t( rprobs[,k,]^2*t(itemwt[,items.temp]) )
#			xbar2[items.temp,k] <- 
#              ttheta.dd2 %*% ( t( rprobs[,k,]^2 ) * itemwt[,items.temp] )
			xbar2[items.temp,k] <- 
              ttheta.dd2 %*% ( t( rprobs.k^2 ) * itemwt[,items.temp] )

					}
          if ( irtmodel %in% c("GPCM","GPCM.design") ){ 	# begin GPCM
            xxf[items.temp,k] <- xxf[items.temp,k] * (k-1)^2	
            xtemp <- xtemp + matrix(theta[,dd],length(items.temp),nrow(theta),byrow=TRUE) * 
              rprobs[,k,] * ( k - 1 )  									
			}
          if ( irtmodel == "2PL.groups"){ 	# begin 2PL groups
            #	xtemp[items.temp,,k] <- matrix(theta[,dd],length(items.temp),
			#           nrow(theta),byrow=TRUE) * rprobs[,k,]
            xbar2[items.temp,k] <- 
              t( theta[,dd]^2 ) %*% t( rprobs[,k,]^2*t(itemwt[,items.temp]) )
			  # t(A) %*% t( B * t(C) ) -> simplify this formula!!!
          }
          
          if ( ! is.null( items.conv) ){
            xbar[ items.conv , ] <- 0
            xxf[ items.conv , ] <- 0
            xbar2[ items.conv , ] <- 0
          }		  
        }	# end categories k
 

		
        ###################################
        if ( irtmodel %in% c("GPCM","GPCM.design")){  # begin GPCM / GPCM.design
          B_obs.temp <- B_obs
          B_obs.temp[ items.temp,,dd] <- 
            matrix( mK-1 , length(items.temp) , maxK , byrow=TRUE ) * B_obs[ items.temp ,,dd]
		  #***
		  xbar[ is.na(xbar) ] <- 0		
          xbar.temp <- matrix( mK-1 , length(items.temp) , maxK , byrow=TRUE ) *xbar					
#          diff.temp <- rowSums(B_obs.temp) - rowSums(xbar.temp)	
          diff.temp <- rowSums(B_obs.temp[,,dd]) - rowSums(xbar.temp)	
          xbar2.temp <- diag( xtemp^2 %*% itemwt[ , items.temp  ] )
          xxf.temp <- rowSums( xxf )		
          if ( ! is.null( items.conv) ){ diff.temp[ items.conv, ] <- 0   }
  
          deriv.temp <- xbar2.temp - xxf.temp
 		  
        }
        if (irtmodel == "GPCM"){	  # begin GPCM
		  #***
		  deriv.temp[ is.na(deriv.temp)] <- 10^(45)	
          increment.temp <- diff.temp*abs(1/( deriv.temp + 10^(-20) ) )
          ci <- ceiling( abs(increment.temp) / ( abs( old_increment.temp[,dd]) + 10^(-10) ) )
          increment.temp <- ifelse( abs( increment.temp) > abs(old_increment.temp[,dd])  , 
                                    increment.temp /(2*ci) , increment.temp )
          old_increment.temp[,dd] <- increment.temp	
          increment <- outer( increment.temp , mK - 1) 
#		  dt1 <- outer( deriv.temp , mK - 1)
		  if (Biter==1){ 
					se.B[,,dd]  <- outer( sqrt( 1 / abs( deriv.temp )) , mK-1 ) 
					}
        }  # end GPCM
        if ( irtmodel == "GPCM.design"){						
          diff.temp <- ( t(E) %*% diff.temp )[,1]
          deriv.temp <- ( t(E^2) %*% deriv.temp )[,1]
          increment.temp <- diff.temp*abs(1/( deriv.temp + 10^(-20) ) )  
          ci <- ceiling( abs(increment.temp) / ( abs( old_increment.temp[,dd]) + 10^(-10) ) )
          increment.temp <- ifelse( abs( increment.temp) > abs(old_increment.temp[,dd])  , 
                                    increment.temp /(2*ci) , increment.temp )
          old_increment.temp[,dd] <- increment.temp	
          basispar[,dd] <- basispar[,dd] + increment.temp
          increment.temp <- ( E %*% increment.temp	)[,1]
          increment <- outer( increment.temp , mK - 1) 
		  d1 <- outer(  1 / abs( deriv.temp ) , mK-1 ) 

          LL <- ncol(d1)		
		  for (ll in 1:LL){
	#		  ll <- 2
			  m1 <- sqrt( diag( E %*% d1[,ll] %*% t( d1[,ll] ) %*% t(E) ) )
			  if (Biter==1){ se.B[,ll,dd]  <- m1 }	  
				}
		#**** Bug fix ARb 2015-12-16
		nB <- dim(B)	
		# B_ind <- 1 * ( B_obs != 0 )
		B_ind <- 1 * ( B_orig != 0 )
		for (dd in 1:nB[3]){
			# dd <- 1
			EB <- E %*% basispar[,dd] 
			for (cc in 1:(nB[2]-1)){ 
			#	cc <- 1
				B[,cc+1,dd] <- cc * EB * B_ind[,cc,dd]
								}
						}
		B00 <- B

        } # end GPCM.design																			
        
        #**********
        if (irtmodel == "2PL.groups"){  # begin 2PL slopegroups
          a1 <- stats::aggregate( B_obs[,,dd] - xbar  , list(est.slopegroups) , sum )		
          a2 <- stats::aggregate( xbar2 - xxf  , list(est.slopegroups) , sum )				
          deriv.temp <- as.matrix(a2[,-1])
          diff.temp <- as.matrix(a1[,-1])
          increment.temp <- diff.temp*abs(1/( deriv.temp + 10^(-20) ) )  
          ci <- ceiling( abs(increment.temp) / ( abs( old_increment.temp[,,dd]) + 10^(-10) ) )
          increment.temp <- ifelse( abs( increment.temp) > abs(old_increment.temp[,,dd])  , 
                                    increment.temp /(2*ci) , increment.temp )
          old_increment.temp[,,dd] <- increment.temp	
		  ind.match <- match( est.slopegroups , a1[,1] )
          increment <- increment.temp[ ind.match , ]
	  
		  if (Biter==1){ se.B[ ,,dd]  <- sqrt( 1 / abs( deriv.temp[ind.match , ] )) }	  
        } # end 2PL slope groups
        #****
        if (irtmodel == "2PL"){		# begin 2PL					
          diff <- B_obs[,,dd] - xbar	
          if ( ! is.null( items.conv) ){  diff[ items.conv, ] <- 0	}
          deriv <- xbar2 - xxf
	  
          increment <- diff*abs(1/( deriv + 10^(-20) ) )        
		  		  
          ci <- ceiling( abs(increment) / ( abs( old_increment[,,dd]) + 10^(-10) ) )
          increment <- ifelse( abs( increment) > abs(old_increment[,,dd])  , 
                               increment/(2*ci) , 
                               increment )	
#          increment <- ifelse( abs( increment) > abs(old_increment[,,dd])  , 
#                               sign(increment)*abs(old_increment[,,dd] ) , increment )	


		  if (Biter==1){ se.B[,,dd] <- sqrt( 1 / abs( deriv )) }
        }   # end 2PL
        ###########################################
        increment[B_orig[,,dd]==0] <- 0  # B[i,k,] could be zero for some dimensions
        old_increment[,,dd] <- increment        	
        B[,,dd] <- B[,,dd] + increment   # update B parameter
		
if (irtmodel=="GPCM.design"){
     B <- B00
				}		
		
      }  # end dimensions
      ###############################################
      if ( irtmodel == "2PL" ){
        items.temp <- which( apply( old_increment , 1 , 
                                    FUN = function(ll){ ! ( max(abs(ll)) < convM ) } ) )
        items.conv <- setdiff( 1:nitems , items.temp )	
      }
      #          if ( max(abs(increment)) < convM ) { converge <- TRUE }
      if ( max(abs(old_increment)) < convM ) { converge <- TRUE }
      Biter <- Biter + 1
      if (progress){ cat( "-" ) ; utils::flush.console()	}		  
    }		# end while loop
	se.B[ B_orig == 0 ] <- 0
    res <- list( "B" = B , "basispar" = basispar , "se.B" = se.B )
    return(res) 
  }
