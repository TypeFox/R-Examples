
#----------------------------------------------------------------------------------------------
.ll.rasch.copula320 <- function( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
		CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
		bdat2.li , bdat2.li.resp , pattern , GG , copula.type , Ncat.ld ,
        Qmatrix=Qmatrix , ... ){
		#-------------------------------------------------------
		# use this function for log likelihood calculation
		# calculation of terms of independent itemclusters?
		calc.ind <- length(itemcluster0) > 0
		G <- GG			# number of groups
		eps1 <- 10^(-14)
		ndat2 <- nrow(dat2.ld)
		M1 <- rep(1,ndat2)
		ntheta <- nrow(theta.k)
		M2 <- rep( 1, ntheta)
 		pjk.theta.k <- .prob.raschtype.genlogis( theta.k , b  , alpha1 , alpha2 , a ,
				Qmatrix=Qmatrix )
		pjk.theta.k01 <- cbind( pjk.theta.k , 1 - pjk.theta.k , 1  )
		
		################################################
		# E step
		################################################
		#.............................................		
		# probabilities for indepedent items
		if ( calc.ind ){
			pjk.theta.k0 <- pjk.theta.k01[ , c( itemcluster0 , itemcluster0 + I ) ]
								} else  {
			pjk.theta.k0 <- NULL
								}
		# probabilities for dependent items
		pjk.theta.kCC <- as.list( 1:CC )
	
		for (cc in 1:CC){
			# cc <- 2	# itemcluster cc
			dp.ld.cc <- dp.ld[[cc]]
			m1.cc <- pjk.theta.k01[ , dp.ld.cc$independent$items ]		
			v1.cc <- dp.ld.cc$independent$N.Index1
			#--------------------------------------------
			# Boundary Mixture Copula (Braeken, 2011)
			if (copula.type[cc] == "bound.mixt" ){
				# likelihood under independence				
				F0pjk.cc <- .rowProds2.bundle( m1 = m1.cc , v1 = v1.cc)
				# likelihood under dependence
				m1.cc <- pjk.theta.k01[ , dp.ld.cc$dependent$items ]		
				v1.cc <- dp.ld.cc$dependent$N.Index1
				pjk.cc <- .rowMins2cpp.bundle( m1 = m1.cc , v1 = v1.cc)
#				F1pjk.cc <- t( dp.ld.cc$calc %*% t( pjk.cc ) )
				F1pjk.cc <- tcrossprod(  pjk.cc ,  dp.ld.cc$calc  )
				pjk.theta.kCC[[cc]] <- ( 1 - delta[cc] ) * F0pjk.cc + delta[cc] * F1pjk.cc
							}
			#-----------------------------------------------
			# include other Copula models here
			# Cook-Johnson copula (Braeken et al., 2007)
			if (copula.type[cc] == "cook.johnson" ){			
				F.Xr <-  1- pjk.theta.k[ , dp.ld.cc$items ]
				R <- ncol(F.Xr)
				delta.cc <- delta[cc]
				patt.cc <- dp.ld.cc$patt
				pjk.cc  <- matrix( 0 , nrow=ntheta , ncol= nrow(patt.cc) )
				for (pp in 1:( nrow(patt.cc) ) ){
					ppcc <- 1 - patt.cc[pp,]
					pjk.cc[ , pp ] <- ( rowSums( ( F.Xr^(-delta.cc))^( outer( rep(1,ntheta) , ppcc )) ) 
													- R + 1 )^(-1/delta.cc)											
								}
#					pjk.theta.kCC[[cc]] <- t( dp.ld.cc$calc %*% t( pjk.cc ) )
#					temp1 <- t( dp.ld.cc$calc %*% t( pjk.cc ) )												
					temp1 <- tcrossprod( pjk.cc , dp.ld.cc$calc )
					temp1[ temp1 < 0 ] <- eps1					
					pjk.theta.kCC[[cc]]	<- temp1					
							}
			#******************************************
			# Frank copula
			if (copula.type[cc] == "frank" ){			
				F.Xr <-  1- pjk.theta.k[ , dp.ld.cc$items ]
				R <- ncol(F.Xr)
				delta.cc <- delta[cc]
				prod.delta <- ( 1 - exp( - delta.cc ) )^(R-1)
				patt.cc <- dp.ld.cc$patt
				pjk.cc  <- matrix( 0 , nrow=ntheta , ncol= nrow(patt.cc) )
				for (pp in 1:( nrow(patt.cc) ) ){
					ppcc <- 1 - patt.cc[pp,]
					g1 <- rowProds2( ( 1 - exp( - delta.cc * F.Xr^( outer( rep(1,ntheta) , ppcc ))  ) ) )
					pjk.cc[,pp] <- - log( 1 - g1 / prod.delta ) / delta.cc										
								}										
					temp1 <- tcrossprod( pjk.cc , dp.ld.cc$calc )
					temp1[ temp1 < 0 ] <- eps1					
					pjk.theta.kCC[[cc]]	<- temp1					
							}  # end Frank copula
						}
								
#print( str(pjk.theta.kCC))				
# cat( "\n probabilities \n" ) ; aa1 <- Sys.time() ; print(aa1-aa0) ; aa0 <- aa1	
	
		#****
		post0 <- matrix( 1 , nrow=ndat2 , ncol=ntheta )
		I0 <- length(itemcluster0)
		if ( calc.ind ){
			pjkL <- array( NA , dim=c(2 , nrow(theta.k) , I0 ) )
			# pjkL: [ #categories , #thetagrid , #items ]
			pjkL[1,,] <- 1 - pjk.theta.k0[,1:I0]
			pjkL[2,,] <- pjk.theta.k0[,1:I0]
			for (ii in 1:I0 ){
				ind.ii <- which( bdat2.li.resp[,ii] == 1 )
				post0[ind.ii,] <- post0[ind.ii,] * pjkL[ bdat2.li[ind.ii,ii]+1 , ,ii]
						}
						}	
									
		#****
		post2 <- matrix( 1 , nrow=ndat2 , ncol=ntheta )
		pjkL <- array( NA , dim=c(Ncat.ld , nrow(theta.k) , CC ) )
	    for (cc in 1:CC){
			#	cc <- 1
			p1.cc <- t( pjk.theta.kCC[[cc]] )	
			pjkL[ seq( 1 , nrow(p1.cc) )  ,, cc ] <- p1.cc
						}
		for (cc in 1:CC ){
			ind.ii <- which( dat2.ld.resp[,cc] == 1 )
			post2[ind.ii,] <- post2[ind.ii,] * pjkL[ dat2.ld[ind.ii,cc] , ,cc]
					}						
						
				
		post <- post0 * post2		# product of independent and dependent parts
		post.unnorm <-  post 		
		post <- post * outer( M1 , wgt.theta )		
		post <- post / rowSums( post)	# standardization of posterior distribution

# cat( "\n posterior\n" ) ; aa1 <- Sys.time() ; print(aa1-aa0) ; aa0 <- aa1			
	
		#....................................................
		# Calculate expected counts
		# expected counts independent item responses
		njk0 <- rjk0.0 <- rjk0.1 <- NULL
		# pattern[,gg+1] is the frequency weight of a response pattern in group gg
		gg <- 1
		if ( calc.ind ){
			njk0 <- crossprod( post , pattern[,gg+1] * dat2.li.resp )			
			rjk0.1 <- crossprod( post ,  pattern[,gg+1] * dat2.li.resp * dat2.li )
			rjk0.0 <- njk0 - rjk0.1
									}
					
		#*****************
		# expected counts dependent items
		rjkCC <- as.list( 1:CC )
		gg <- 1		
		for ( cc in 1:CC){	
			rjkCC[[cc]] <- crossprod( post ,  pattern[,gg+1] * dat3.ld[[cc]] * dat2.ld.resp[,cc] )
					}
					
		# total counts
		gg <- 1
		tc <- outer( pattern[,gg+1] , rep( 1 , ntheta ) )
		nk <- colSums( tc * post)
		nk <- nk + 1e-30
		pik <- nk / sum(nk)

# cat( "\n expected counts \n" ) ; aa1 <- Sys.time() ; print(aa1-aa0) ; aa0 <- aa1			
		
		#---------------------------------------------------------------
		# calculate log likelihood
		#.......................
		# one group 
		
		if (G == 1){ 
			 ll0 <- sum( nk * log(pik) )					
#			 ll0 <- sum( nk * log(wgt.theta) )					
			 # likelihood part from independent items		 
			 if ( calc.ind ){	
				if ( nrow( rjk0.1) == 1 ){ 
					rjk.temp <- cbind( t( rjk0.1) , t(rjk0.0) )
							} else {
					rjk.temp <- cbind( rjk0.1 , rjk0.0 )
								}
					ll0 <- ll0 + sum( log(pjk.theta.k0 ) * rjk.temp )
							}	
			 # likelihood part from dependent items
			for (cc in 1:CC){
				ll0 <- ll0 + sum( rjkCC[[cc]] * log( pjk.theta.kCC[[cc]] + 10^(-15) ) )
								}
								
					}
# cat( "\n Likelihood \n" ) ; aa1 <- Sys.time() ; print(aa1-aa0) ; aa0 <- aa1	
				
		#.......................
		# multiple groups 
		#	... to do	...



			lli <- ll0
		################################
		# arrange output
		res <- list( "ll"=ll0 ,  "b" = b , "a" = a , "delta"=delta ,
						"alpha1" = alpha1 , "alpha2" = alpha2 , "lli" = lli ,
						"post" = post , "post.unnorm" = post.unnorm ,
						"rjk0.1" = rjk0.1 ,
						"rjk0.0" = rjk0.0 , "rjkCC" = rjkCC ,
						"pjk.theta.kCC" = pjk.theta.kCC , 
						"pjk.theta.k0" = pjk.theta.k0 ,
						"nk" = nk , "pik" = pik ,"calc.ind" = calc.ind
													)
		return(res)
				}
#----------------------------------------------------------------------------------------------				




#----------------------------------------------------------------------------------------------
.update.ll.rasch.copula321 <- function( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
		CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
		bdat2.li , bdat2.li.resp , rescopula , itemcluster , pattern , GG , copula.type , 
		Qmatrix=Qmatrix , ... ){
		#-------------------------------------------------------
		# use this function for log likelihood calculation
		# look for items which change parameters for necessary update
		G <- GG
		# calculation of terms of independent itemclusters?
		calc.ind <- length(itemcluster0) > 0	
		eps1 <- 10^(-14)		
		# calculate necessary updates
		ind.b <- which( b != rescopula$b )
		ind.a <- which( a != rescopula$a )
		ind.delta <- which( delta != rescopula$delta )
		ind.alpha1 <- ( alpha1 != rescopula$alpha1 )	+ ( alpha2 != rescopula$alpha2 )
		if (ind.alpha1 > 0){ ind.alpha <- seq(1 , ncol(dat2.ld) ) } else { ind.alpha <- NULL }
		itemset <- union( ind.b , ind.a )
		itemset <- union( itemset , ind.alpha )	
		# update term local independence
		li.update <- 1 * ( sum( itemcluster0 %in% itemset ) > 0 )
		# update terms item dependence parameters
		ld.update <- sapply( 1:CC , FUN = function(cc){ 
				g1 <- intersect( which( itemcluster == cc )  , itemset )
				if ( length(g1)){ v1 <- cc } else { v1 <- NULL }
				v1
					} )
		ld.update <- unique( union( ind.delta , unlist( ld.update) ) )
		###########################################################################
		ndat2 <- nrow(dat2.ld)
		M1 <- rep(1,ndat2)
		ntheta <- nrow(theta.k)
		M2 <- rep( 1, ntheta)
		pjk.theta.k <- .prob.raschtype.genlogis( theta.k , b  , alpha1 , alpha2 , a ,
			Qmatrix=Qmatrix)
		pjk.theta.k01 <- cbind( pjk.theta.k , 1 - pjk.theta.k , 1  )
		#.............................................		
		# probabilities for indepedent items
		if ( calc.ind ){
			pjk.theta.k0 <- pjk.theta.k01[ , c( itemcluster0 , itemcluster0 + I ) ]
								} else  {
			pjk.theta.k0 <- NULL
								}
		# probabilities for dependent items
		pjk.theta.kCC <- rescopula$pjk.theta.kCC
		for (cc in ld.update){
			# cc <- 2	# itemcluster cc
			dp.ld.cc <- dp.ld[[cc]]
			m1.cc <- pjk.theta.k01[ , dp.ld.cc$independent$items ]		
			v1.cc <- dp.ld.cc$independent$N.Index1
			#--------------------------------------------
			# Boundary Mixture Copula (Braeken, 2011)
			if (copula.type[cc] == "bound.mixt" ){			
				# likelihood under independence				
				F0pjk.cc <- .rowProds2.bundle( m1 = m1.cc , v1 = v1.cc)
				# likelihood under dependence
				m1.cc <- pjk.theta.k01[ , dp.ld.cc$dependent$items ]		
				v1.cc <- dp.ld.cc$dependent$N.Index1
				pjk.cc <- .rowMins2cpp.bundle( m1 = m1.cc , v1 = v1.cc)
#				F1pjk.cc <- t( dp.ld.cc$calc %*% t( pjk.cc ) )
                # t( A * t(B) ) = B * t(A) 
				F1pjk.cc <- tcrossprod(  pjk.cc ,  dp.ld.cc$calc  ) 
				pjk.theta.kCC[[cc]] <- ( 1 - delta[cc] ) * F0pjk.cc + delta[cc] * F1pjk.cc
										}
			#-----------------------------------------------
			# Cook-Johnson Copula
			if (copula.type[cc] == "cook.johnson" ){			
				F.Xr <-  1- pjk.theta.k[ , dp.ld.cc$items ]
				R <- ncol(F.Xr)
				delta.cc <- delta[cc]
				patt.cc <- dp.ld.cc$patt
				pjk.cc  <- matrix( 0 , nrow=ntheta , ncol= nrow(patt.cc) )
				for (pp in 1:( nrow(patt.cc) ) ){
					ppcc <- 1 - patt.cc[pp,]
					pjk.cc[ , pp ] <- ( rowSums( ( F.Xr^(-delta.cc))^( outer( rep(1,ntheta) , ppcc )) ) 
													- R + 1 )^(-1/delta.cc)
		#			temp1 <- t( dp.ld.cc$calc %*% t( pjk.cc ) )		
					temp1 <- tcrossprod( pjk.cc , dp.ld.cc$calc )
					temp1[ temp1 < 0 ] <- eps1					
					pjk.theta.kCC[[cc]]	<- temp1
					
								}
							}
			#******************************************
			# Frank copula
			if (copula.type[cc] == "frank" ){			
				F.Xr <-  1- pjk.theta.k[ , dp.ld.cc$items ]
				R <- ncol(F.Xr)
				delta.cc <- delta[cc]
				prod.delta <- ( 1 - exp( - delta.cc ) )^(R-1)
				patt.cc <- dp.ld.cc$patt
				pjk.cc  <- matrix( 0 , nrow=ntheta , ncol= nrow(patt.cc) )
				for (pp in 1:( nrow(patt.cc) ) ){
					ppcc <- 1 - patt.cc[pp,]
					g1 <- rowProds2( ( 1 - exp( - delta.cc * F.Xr^( outer( rep(1,ntheta) , ppcc ))  ) ) )
					pjk.cc[,pp] <- - log( 1 - g1 / prod.delta ) / delta.cc										
								}
#					temp1 <- t( dp.ld.cc$calc %*% t( pjk.cc ) )
					temp1 <- tcrossprod( pjk.cc , dp.ld.cc$calc )
					temp1[ temp1 < 0 ] <- eps1					
					pjk.theta.kCC[[cc]]	<- temp1
							} # end Frank copula
						}							
		#---------------------------------------------------------------

		
		#####################################
		# rearrange output
				res <- list( "pjk.theta.kCC"=pjk.theta.kCC , "pjk.theta.k0" = pjk.theta.k0  )
				return(res)
				}
#----------------------------------------------------------------------------------------------				




#----------------------------------------------------------------------------------------------
.update.ll.rasch.copula320 <- function( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
		CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
		bdat2.li , bdat2.li.resp , rescopula , itemcluster , pattern , GG , copula.type , 
		Qmatrix=Qmatrix , ... ){
		#-------------------------------------------------------
		# use this function for log likelihood calculation
		# look for items which change parameters for necessary update
		G <- GG
		# calculation of terms of independent itemclusters?
		calc.ind <- length(itemcluster0) > 0	
		eps1 <- 10^(-14)		
		# calculate necessary updates
		ind.b <- which( b != rescopula$b )
		ind.a <- which( a != rescopula$a )
		ind.delta <- which( delta != rescopula$delta )
		ind.alpha1 <- ( alpha1 != rescopula$alpha1 )	+ ( alpha2 != rescopula$alpha2 )
		if (ind.alpha1 > 0){ ind.alpha <- seq(1 , ncol(dat2.ld) ) } else { ind.alpha <- NULL }
		itemset <- union( ind.b , ind.a )
		itemset <- union( itemset , ind.alpha )	
		# update term local independence
		li.update <- 1 * ( sum( itemcluster0 %in% itemset ) > 0 )
		# update terms item dependence parameters
		ld.update <- sapply( 1:CC , FUN = function(cc){ 
				g1 <- intersect( which( itemcluster == cc )  , itemset )
				if ( length(g1)){ v1 <- cc } else { v1 <- NULL }
				v1
					} )
		ld.update <- unique( union( ind.delta , unlist( ld.update) ) )
		###########################################################################
		ndat2 <- nrow(dat2.ld)
		M1 <- rep(1,ndat2)
		ntheta <- nrow(theta.k)
		M2 <- rep( 1, ntheta)
		pjk.theta.k <- .prob.raschtype.genlogis( theta.k , b  , alpha1 , alpha2 , a , Qmatrix=Qmatrix)
		pjk.theta.k01 <- cbind( pjk.theta.k , 1 - pjk.theta.k , 1  )
		#.............................................		
		# probabilities for indepedent items
		if ( calc.ind ){
			pjk.theta.k0 <- pjk.theta.k01[ , c( itemcluster0 , itemcluster0 + I ) ]
								} else  {
			pjk.theta.k0 <- NULL
								}
		# probabilities for dependent items
		pjk.theta.kCC <- rescopula$pjk.theta.kCC
		for (cc in ld.update){
			# cc <- 2	# itemcluster cc
			dp.ld.cc <- dp.ld[[cc]]
			m1.cc <- pjk.theta.k01[ , dp.ld.cc$independent$items ]		
			v1.cc <- dp.ld.cc$independent$N.Index1
			#--------------------------------------------
			# Boundary Mixture Copula (Braeken, 2011)
			if (copula.type[cc] == "bound.mixt" ){			
				# likelihood under independence				
				F0pjk.cc <- .rowProds2.bundle( m1 = m1.cc , v1 = v1.cc)
				# likelihood under dependence
				m1.cc <- pjk.theta.k01[ , dp.ld.cc$dependent$items ]		
				v1.cc <- dp.ld.cc$dependent$N.Index1
				pjk.cc <- .rowMins2cpp.bundle( m1 = m1.cc , v1 = v1.cc)
#				F1pjk.cc <- t( dp.ld.cc$calc %*% t( pjk.cc ) )
				F1pjk.cc <- tcrossprod(  pjk.cc ,  dp.ld.cc$calc  )
				pjk.theta.kCC[[cc]] <- ( 1 - delta[cc] ) * F0pjk.cc + delta[cc] * F1pjk.cc
										}
			#-----------------------------------------------
			# Cook-Johnson Copula
			if (copula.type[cc] == "cook.johnson" ){			
				F.Xr <-  1- pjk.theta.k[ , dp.ld.cc$items ]
				R <- ncol(F.Xr)
				delta.cc <- delta[cc]
				patt.cc <- dp.ld.cc$patt
				pjk.cc  <- matrix( 0 , nrow=ntheta , ncol= nrow(patt.cc) )
				for (pp in 1:( nrow(patt.cc) ) ){
					ppcc <- 1 - patt.cc[pp,]
					pjk.cc[ , pp ] <- ( rowSums( ( F.Xr^(-delta.cc))^( outer( rep(1,ntheta) , ppcc )) ) 
													- R + 1 )^(-1/delta.cc)
#					temp1 <- t( dp.ld.cc$calc %*% t( pjk.cc ) )												
					temp1 <- tcrossprod( pjk.cc , dp.ld.cc$calc )
					temp1[ temp1 < 0 ] <- eps1					
					pjk.theta.kCC[[cc]]	<- temp1
					
								}
							}
			#******************************************
			# Frank copula
			if (copula.type[cc] == "frank" ){			
				F.Xr <-  1- pjk.theta.k[ , dp.ld.cc$items ]
				R <- ncol(F.Xr)
				delta.cc <- delta[cc]
				prod.delta <- ( 1 - exp( - delta.cc ) )^(R-1)
				patt.cc <- dp.ld.cc$patt
				pjk.cc  <- matrix( 0 , nrow=ntheta , ncol= nrow(patt.cc) )
				for (pp in 1:( nrow(patt.cc) ) ){
					ppcc <- 1 - patt.cc[pp,]
					g1 <- rowProds2( ( 1 - exp( - delta.cc * F.Xr^( outer( rep(1,ntheta) , ppcc ))  ) ) )
					pjk.cc[,pp] <- - log( 1 - g1 / prod.delta ) / delta.cc										
								}
#					temp1 <- t( dp.ld.cc$calc %*% t( pjk.cc ) )
					temp1 <- tcrossprod( pjk.cc , dp.ld.cc$calc )
					temp1[ temp1 < 0 ] <- eps1					
					pjk.theta.kCC[[cc]]	<- temp1
							} # end Frank copula
						}							
		#---------------------------------------------------------------
		# calculate log likelihood
		rjk0.1 <- rescopula$rjk0.1
		rjk0.0 <- rescopula$rjk0.0
		rjkCC <- rescopula$rjkCC
		#.......................
		# one group
		if (G == 1){ 
			ll0 <- 	sum( rescopula$nk * log(rescopula$pik) )			
			 # likelihood part from independent items
			 if ( calc.ind ){	
				if ( nrow( rjk0.1) == 1 ){ 
					rjk.temp <- cbind( t( rjk0.1) , t(rjk0.0) )
							} else {
					rjk.temp <- cbind( rjk0.1 , rjk0.0 )
								}			 
					ll0 <- ll0 + sum( log(pjk.theta.k0) * rjk.temp )
							}
			 # likelihood part from dependent items
			for (cc in 1:CC){
				ll0 <- ll0 + sum( rjkCC[[cc]] * log( pjk.theta.kCC[[cc]] + 10^(-15) ) )
							}

								}
#					}
		#.......................
		# multiple groups 
		#	... to do	...		
			lli <- ll0		
		
		#####################################
		# rearrange output
				res <- list( "ll" = ll0 , "lli" = lli )
				return(res)
				}
#----------------------------------------------------------------------------------------------				




