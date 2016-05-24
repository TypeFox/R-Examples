



#--------------------------------------------------------------------------------------
# Function calculates necessary patterns for copula IRT models (Braeken, 2011)
.calc.copula.itemcluster <- function(D){
#aa0 <- Sys.time()
    res <- gtools::permutations(n=2, r=D, v=0:1, repeats.allowed=TRUE)
    rownames(res) <- apply( res , 1 , FUN = function(ll){ paste("P" , paste( ll ,collapse="") ,sep="") } )
    RR <- nrow(res) 
    matr <- matrix( 0 , RR , RR )
    rownames(matr) <- colnames(matr) <- rownames(res) 
    colnames(matr) <- gsub( "P" , "F" , colnames(matr) )
#cat("   ***  permutations") ; aa1 <- Sys.time(); print(aa1-aa0) ; aa0 <- aa1	
    vec <- 1:RR
    # calculation of formulas
	# This loop must be simplified!!!
    for (rr in vec){
        # rr <- 2
        res.rr <- outer( rep(1,nrow(res)) , res[rr,] ) - res
        a1.rr <- apply( res.rr , 1 , FUN = function(ll){ paste("F" , paste( ll ,collapse="") ,sep="") } )
        g1.rr <- ( (-1)^rowSums( res ))
        ind.rr <- which( apply( res.rr , 1 , min ) > -1 )
        a1.rr <- a1.rr[ind.rr]
        g1.rr <- g1.rr[ind.rr]
        matr[ rr , a1.rr ]  <- g1.rr
            }
# cat("   ***  after outer") ; aa1 <- Sys.time(); print(aa1-aa0) ; aa0 <- aa1			
# print(matr)
    res1 <- list( "patt" = res , "calc" = matr )
    return(res1)
    }
#--------------------------------------------------------------------------------------



#--------------------------------------------------------------------------------------
# Function calculates necessary patterns for copula IRT models (Braeken, 2011)
.calc.copula.itemcluster2 <- function(D){
#aa0 <- Sys.time()
    res <- gtools::permutations(n=2, r=D, v=0:1, repeats.allowed=TRUE)
    rownames(res) <- apply( res , 1 , FUN = function(ll){ paste("P" , paste( ll ,collapse="") ,sep="") } )
    RR <- nrow(res) 
    matr <- matrix( 0 , RR , RR )
    rownames(matr) <- colnames(matr) <- rownames(res) 
    colnames(matr) <- gsub( "P" , "F" , colnames(matr) )
#cat("   ***  permutations") ; aa1 <- Sys.time(); print(aa1-aa0) ; aa0 <- aa1	
    matr <- .Call( "calc_copula_itemcluster_C" , D , res , package="sirt" )$matr	
# cat("   ***  after outer") ; aa1 <- Sys.time(); print(aa1-aa0) ; aa0 <- aa1			
# print(matr)
    res1 <- list( "patt" = res , "calc" = matr )
    return(res1)
    }
#--------------------------------------------------------------------------------------




#----------------------------------------------------------------------------------------------
.update.ll.rasch.copula21 <- function( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
		CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
		bdat2.li , bdat2.li.resp , rescopula , itemcluster , pattern , GG , copula.type , ... ){
		#-------------------------------------------------------
		# use this function for log likelihood calculation
		# look for items which change parameters for necessary update
		G <- GG
# vv0 <- Sys.time()		
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
# cat("    ---- ld.update") ; vv1 <- Sys.time(); print(vv1-vv0) ; vv0 <- vv1		
		###########################################################################
		ndat2 <- nrow(dat2.ld)
		M1 <- rep(1,ndat2)
		ntheta <- length(theta.k)
		M2 <- rep( 1, ntheta)
		pjk.theta.k <- .prob.raschtype.genlogis( theta.k , b  , alpha1 , alpha2 , a)
		pjk.theta.k01 <- cbind( pjk.theta.k , 1 - pjk.theta.k , 1  )
		#.............................................		
		# probabilities for indepedent items
		if ( calc.ind ){
			pjk.theta.k0 <- pjk.theta.k01[ , c( itemcluster0 , itemcluster0 + I ) ]
								} else  {
			pjk.theta.k0 <- NULL
								}
# cat("    ---- probs independence") ; vv1 <- Sys.time(); print(vv1-vv0) ; vv0 <- vv1										
		# probabilities for dependent items
		pjk.theta.kCC <- rescopula$pjk.theta.kCC
		for (cc in ld.update){
# cat("     ***" , 'cc=' , cc )
# hh0 <- Sys.time()		
			# cc <- 2	# itemcluster cc
			dp.ld.cc <- dp.ld[[cc]]
			m1.cc <- pjk.theta.k01[ , dp.ld.cc$independent$items ]		
			v1.cc <- dp.ld.cc$independent$N.Index1
# cat("    ---- cc begin probs") ; hh1 <- Sys.time(); print(hh1-hh0) ; hh0 <- hh1													
			#--------------------------------------------
			# Boundary Mixture Copula (Braeken, 2011)
			if (copula.type[cc] == "bound.mixt" ){			
				# likelihood under independence				
				F0pjk.cc <- .rowProds2.bundle( m1 = m1.cc , v1 = v1.cc)
				# likelihood under dependence
				m1.cc <- pjk.theta.k01[ , dp.ld.cc$dependent$items ]		
				v1.cc <- dp.ld.cc$dependent$N.Index1

#cat("    ---- cc before bundle") ; hh1 <- Sys.time(); print(hh1-hh0) ; hh0 <- hh1				
											
				pjk.cc <- .rowMins2cpp.bundle( m1 = m1.cc , v1 = v1.cc)
#cat("    ---- cc after bundle") ; hh1 <- Sys.time(); print(hh1-hh0) ; hh0 <- hh1				

# This function is most time consuming!!!!
#				F1pjk.cc <- t( dp.ld.cc$calc %*% t( pjk.cc ) )
                # t( A * t(B) ) = B * t(A) 
				F1pjk.cc <- tcrossprod(  pjk.cc ,  dp.ld.cc$calc  ) 
				pjk.theta.kCC[[cc]] <- ( 1 - delta[cc] ) * F0pjk.cc + delta[cc] * F1pjk.cc
# cat("    ---- cc clac pjk.theta.kCC") ; hh1 <- Sys.time(); print(hh1-hh0) ; hh0 <- hh1																											
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
					pjk.cc[ , pp ] <- ( rowSums( 
						( F.Xr^(-delta.cc))^( outer( rep(1,ntheta) , ppcc )) ) 	- R + 1 )^(-1/delta.cc)
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
# cat("    ---- probs dependence") ; vv1 <- Sys.time(); print(vv1-vv0) ; vv0 <- vv1										
		
		#####################################
		# rearrange output
				res <- list( "pjk.theta.kCC"=pjk.theta.kCC , "pjk.theta.k0" = pjk.theta.k0  )
				return(res)
				}
#----------------------------------------------------------------------------------------------				



#----------------------------------------------------------------------------------------------
.update.ll.rasch.copula20 <- function( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
		CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
		bdat2.li , bdat2.li.resp , rescopula , itemcluster , pattern , GG , copula.type , ... ){
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
		ntheta <- length(theta.k)
		M2 <- rep( 1, ntheta)
		pjk.theta.k <- .prob.raschtype.genlogis( theta.k , b  , alpha1 , alpha2 , a)
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




#----------------------------------------------------------------------------------------------
.ll.rasch.copula20 <- function( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
		CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
		bdat2.li , bdat2.li.resp , pattern , GG , copula.type , Ncat.ld , ... ){
		#-------------------------------------------------------
		# use this function for log likelihood calculation
		# calculation of terms of independent itemclusters?
		calc.ind <- length(itemcluster0) > 0
		G <- GG			# number of groups
		eps1 <- 10^(-14)
		ndat2 <- nrow(dat2.ld)
		M1 <- rep(1,ndat2)
		ntheta <- length(theta.k)
		M2 <- rep( 1, ntheta)
 		pjk.theta.k <- .prob.raschtype.genlogis( theta.k , b  , alpha1 , alpha2 , a)
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
#					pjk.theta.kCC[[cc]] <- t( dp.ld.cc$calc %*% t( pjk.cc ) )	
#					temp1 <- t( dp.ld.cc$calc %*% t( pjk.cc ) )												
					temp1 <- tcrossprod( pjk.cc , dp.ld.cc$calc )
					temp1[ temp1 < 0 ] <- eps1					
					pjk.theta.kCC[[cc]]	<- temp1					
							}  # end Frank copula
						}
								
#print( str(pjk.theta.kCC))				
# cat( "\n probabilities \n" ) ; aa1 <- Sys.time() ; print(aa1-aa0) ; aa0 <- aa1	
						
		#.............................................
		# Calculate posterior distribution	
#		post0 <- matrix( 1 , nrow=ndat2 , ncol=ntheta )
		# posterior distribution independent items
#		if ( length(itemcluster0) > 0){ 
#			post0 <- sapply( 1:ntheta , FUN = function(tt){ 
#					# tt <- 11
#					g1 <- outer( M1 , pjk.theta.k0[ tt , ] )
#					rowProds2( g1^( bdat2.li * bdat2.li.resp ) )
#						} )
#							}
							
		#****
		post0 <- matrix( 1 , nrow=ndat2 , ncol=ntheta )
		I0 <- length(itemcluster0)
		if ( calc.ind ){
			pjkL <- array( NA , dim=c(2 , length(theta.k) , I0 ) )
			# pjkL: [ #categories , #thetagrid , #items ]
			pjkL[1,,] <- 1 - pjk.theta.k0[,1:I0]
			pjkL[2,,] <- pjk.theta.k0[,1:I0]
			for (ii in 1:I0 ){
				ind.ii <- which( bdat2.li.resp[,ii] == 1 )
				post0[ind.ii,] <- post0[ind.ii,] * pjkL[ bdat2.li[ind.ii,ii]+1 , ,ii]
						}
						}	
					
		# posterior distribution dependent items
#		post2 <- sapply( 1:ntheta , FUN = function(tt){ 		
				# tt <- 11
#				h1 <- 1
#				for (cc in 1:CC){
					# cc <- 1
#					pcc <- pjk.theta.kCC[[cc]]
#					h1 <- h1 * ( (pcc[tt,])[ dat2.ld[,cc] ] )^( dat2.ld.resp[,cc] )
#								}
#					h1
#							} )								
		#****
		post2 <- matrix( 1 , nrow=ndat2 , ncol=ntheta )
		pjkL <- array( NA , dim=c(Ncat.ld , length(theta.k) , CC ) )
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
#			FW1 <- outer( pattern[,gg+1] , rep( 1 , ncol(dat2.li.resp) ) )
			# students at items
#			njk0 <- t( sapply( 1:ntheta , FUN = function(tt){
#						colSums( FW1 * dat2.li.resp * outer( post[,tt] , rep( 1 , ncol(dat2.li.resp) ) ) )
#							} )	)
#			njk0 <- t( post ) %*% ( pattern[,gg+1] * dat2.li.resp )
			njk0 <- crossprod( post , pattern[,gg+1] * dat2.li.resp )			
			
			# how many students solved correctly items
#			rjk0.1 <- t( sapply( 1:ntheta , FUN = function(tt){
#						colSums( dat2.li * FW1 * dat2.li.resp * outer( post[,tt] , rep( 1 , ncol(dat2.li.resp) ) ) )
#							} )	)
#			rjk0.1 <- t( post ) %*% ( pattern[,gg+1] * dat2.li.resp * dat2.li )
			rjk0.1 <- crossprod( post ,  pattern[,gg+1] * dat2.li.resp * dat2.li )
			rjk0.0 <- njk0 - rjk0.1
									}
		# use loops in case of multiple groups !!!
				
		
		#*****************
		# expected counts dependent items
		rjkCC <- as.list( 1:CC )
		gg <- 1		
		for ( cc in 1:CC){	
			#	cc <- 1		
#			rjkCC[[cc]] <- t(post) %*% ( pattern[,gg+1] * dat3.ld[[cc]] * dat2.ld.resp[,cc] )
			rjkCC[[cc]] <- crossprod( post ,  pattern[,gg+1] * dat3.ld[[cc]] * dat2.ld.resp[,cc] )
					}
					
		# total counts
		gg <- 1
		tc <- outer( pattern[,gg+1] , rep( 1 , ntheta ) )
		nk <- colSums( tc * post)
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
					ll0 <- ll0 + sum( log(pjk.theta.k0) * rjk.temp )
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
		#~~~~~~~~~~~~~~
		# Output Version < 1.0
		#				res <- list( "ll"=ll1 , "g1"=res , "b" = b , "a" = a , "delta"=delta ,
		#									"alpha1" = alpha1 , "alpha2" = alpha2 , "lli" = lli ,
		#									"ll.theta.post" = ll.theta.post ,
		#									"exp.theta.post" = exp.theta.post
		#											)
				}
#----------------------------------------------------------------------------------------------				



###############################################################################################
# auxiliary function: person parameter estimation
.likelihood.rasch.copula <- function( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
        CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , I , 
        bdat2.li , bdat2.li.resp , eps=10^(-20) , ... ){
        pjk.theta.k.tt <- .prob.raschtype.genlogis( theta.k , b  , alpha1 , alpha2   , a)
        g1 <- matrix( 0 , nrow(pjk.theta.k.tt) , CC + 1 )
		ndat2 <- nrow(dat2.ld)
        M1 <- rep(1,ndat2)
        pqjk.theta.k.tt <- cbind( pjk.theta.k.tt , 1 - pjk.theta.k.tt , 1 )
                    if ( length(itemcluster0) > 0 ){
                pqjk.theta.k.tt0 <- pqjk.theta.k.tt[ , c( itemcluster0 , itemcluster0+I) ]
                # likelihood for independent items at theta tt
                ll.tt <- ( pqjk.theta.k.tt0 ^ bdat2.li )^bdat2.li.resp
                g1[,1] <- rowProds2( ll.tt )        } else { g1[,1] <- 1 }
 
            g1i <- g1    
            # likelihood for dependent items
            for (cc in 1:CC){
                #                cc <- 2
                dat3.ld.cc <- dat3.ld[[cc]] 
                dp.ld.cc <- dp.ld[[cc]]
                m1.cc <- pqjk.theta.k.tt[ , dp.ld.cc$independent$items ]
                v1.cc <- dp.ld.cc$independent$N.Index1
                # product under independence                
                Fpjk.cc <- .rowProds2.bundle( m1 = m1.cc , v1 = v1.cc)
                # evaluate likelihood
                g1i[,cc+1] <- g1.tt <- ( rowSums(Fpjk.cc * dat3.ld.cc ) )^dat2.ld.resp[  ,cc]
                # product under dependence  
                m1.cc <-  pqjk.theta.k.tt[ , dp.ld.cc$dependent$items ] 
                v1.cc <- dp.ld.cc$dependent$N.Index1
                F0pjk.cc <- .rowMins2cpp.bundle( m1 = m1.cc , v1 = v1.cc)
#                F1pjk.cc <- F0pjk.cc %*% t(dp.ld.cc$calc)
                F1pjk.cc <- tcrossprod( F0pjk.cc , dp.ld.cc$calc )
                g2.tt <- ( rowSums(F1pjk.cc * dat3.ld.cc) )^dat2.ld.resp[ ,cc]
                g3.tt <- ( 1 - delta[cc] ) * g1.tt + delta[cc] * g2.tt               
                g1[,cc+1] <- g3.tt
                            }
                res <- g1   
                g1 <- rowProds2( g1 )
                g1i <- rowProds2( g1i ) 
				# calculate log likelihood
				g1[ g1 < eps] <- eps
				g1i[ g1i < eps] <- eps				
				g1 <- log( g1 )
				g1i <- log( g1i )
                res <- list( "loglike.dep" = g1 , "loglike.ind" = g1i )
                return(res)
                    }                
###############################################################################################



##################################################################
# function for person parameter estimation in rasch copula models
person.parameter.rasch.copula <- function( raschcopula.object , numdiff.parm = .001 , 
					conv.parm = .001 , maxiter = 20 , stepwidth = 1 , 
					print.summary = TRUE , 					... ){
        dat2.li <- NULL					
        dat2 <- raschcopula.object$datalist$dat2
        dat2.resp <- raschcopula.object$datalist$dat2.resp
        dat2.ld <- raschcopula.object$datalist$dat2.ld
        dat2.ld.resp <- raschcopula.object$datalist$dat2.ld.resp
        dat2.li.resp <- raschcopula.object$datalist$dat2.li.resp
        dat3.ld <- raschcopula.object$datalist$dat3.ld
        bdat2.li <- raschcopula.object$datalist$bdat2.li
        bdat2.li.resp <- raschcopula.object$datalist$bdat2.li.resp
        CC <- raschcopula.object$datalist$CC
        dp.ld <- raschcopula.object$datalist$dp.ld
        itemcluster0 <- raschcopula.object$datalist$itemcluster0
        b <- raschcopula.object$b
        a <- raschcopula.object$a
        alpha1 <- raschcopula.object$alpha1
        alpha2 <- raschcopula.object$alpha2
        delta <- raschcopula.object$delta
		pattern <- raschcopula.object$pattern
        ndat2 <- nrow(dat2)
        ntheta <- nrow(dat2)
        I <- ncol(dat2)
		######################
		# missing response pattern
		mp <- paste("R" , dat2.resp[,1] , sep="")
		for (vv in seq( 2 , ncol(dat2) ) ){
			mp <- paste( mp , dat2.resp[,vv] , sep="")
					}
        # initial estimate of theta
#        theta0 <- rowMeans( dat2 * dat2.resp  ) 
        dat20 <- dat2 * dat2.resp
		dat20[ dat2.resp == 0 ] <- NA
		theta0 <- rowMeans( dat20 , na.rm=T )
        theta0 <- stats::qlogis( theta0 )
        # theta0[1:2] <- c(0,1)
        theta0[ theta0 == - Inf] <- -9999
        theta0[ theta0 ==  Inf] <- 9999
        theta.init <- theta0i <- theta0
        ii <- 0
        h <- numdiff.parm
        a1m <- 990
        while( a1m > conv.parm & ii < maxiter ){
            # evaluation of likelihood at theta0
            rescop0 <- .likelihood.rasch.copula( theta.k = theta0 , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
                        CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , I , 
                        bdat2.li , bdat2.li.resp  )
            rescop1 <- .likelihood.rasch.copula( theta.k = theta0 + h, b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
                        CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , I , 
                        bdat2.li , bdat2.li.resp  )
            rescop2 <- .likelihood.rasch.copula( theta.k = theta0 - h, b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
                        CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , I , 
                        bdat2.li , bdat2.li.resp  )
            rescop0i <- .likelihood.rasch.copula( theta.k = theta0i , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
                        CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , I , 
                        bdat2.li , bdat2.li.resp  )
            rescop1i <- .likelihood.rasch.copula( theta.k = theta0i + h, b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
                        CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , I , 
                        bdat2.li , bdat2.li.resp  )
            rescop2i <- .likelihood.rasch.copula( theta.k = theta0i - h, b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
                        CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , I , 
                        bdat2.li , bdat2.li.resp  )                        
            #******
            # estimation assuming dependence
            ll1 <- rescop1$loglike.dep
            ll2 <- rescop2$loglike.dep
            ll0 <- rescop0$loglike.dep            
            d1 <- ( ll1 - ll2  ) / ( 2 * h )    
            # second order derivative
            # f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
            d2d <- d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2       
            theta.change <- - d1 / d2
            theta.change[ abs( theta.init ) == 9999 ] <- 0
			theta.change[ is.na( theta.change ) ] <- 0
            a1t1 <- theta.change <- ifelse( abs( theta.change ) > stepwidth , stepwidth*sign(theta.change) , theta.change )                        
            theta0 <- theta0 + theta.change
			ind1 <- ( abs( theta.change ) < conv.parm )

            #******
            # estimation assuming independence
            ll1 <- rescop1i$loglike.ind
            ll2 <- rescop2i$loglike.ind
            ll0 <- rescop0i$loglike.ind            
            d1 <- ( ll1 - ll2  ) / ( 2 * h )    
            # second order derivative
            # f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
            d2i <- d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2       
            theta.change <- - d1 / d2
            theta.change[ abs( theta.init ) == 9999 ] <- 0
			theta.change[ is.na( theta.change ) ] <- 0
            a1t2 <- theta.change <- ifelse( abs( theta.change ) > stepwidth , stepwidth*sign(theta.change) , theta.change )                        
            theta0i <- theta0i + theta.change
            a1m <- max( abs(a1t1) , abs(a1t2) )
            ii <- ii+1
            cat("Iteration" , ii , ":   max. parm. change" , round( a1m , 5))
			cat("   |" , sum(ind1) , "out of" , ndat2 , "cases converged")
			cat(" (", round(100*sum(ind1)/length(ind1),1) , "%)\n")
			utils::flush.console()
                    }
		theta0[ abs( theta.init  ) == 9999 ] <- NA
		theta0i[ abs( theta.init  ) == 9999 ] <- NA
        res <- data.frame( "pattern" = pattern[,1] , 
					"missing.pattern" = match( mp , unique(mp) ) , 
					"freqwgt" = pattern$freqwgt , 
					"converged" = 1*(ind1 == 1) , 
				"score" = rowSums(dat2) , "max" = rowSums( dat2.resp) , 
                 "theta.dep" = theta0 , "theta.ind" = theta0i )		
		res$setheta.dep <- sqrt( - 1 / d2d )
		res$setheta.ind <- sqrt( - 1 / d2i )
		res$setheta.dep[ is.na(theta0i) ] <- NA
		res$setheta.ind[ is.na(theta0i) ] <- NA
		res$seinflat <- res$setheta.dep / res$setheta.ind
		res[ is.na(res$theta.dep) , "converged" ] <- NA
		x1 <- seq( grep( "theta" , colnames(res) )[1] , ncol(res) ) 
		for (vv in x1){ 
				res[ paste(res$converged) == 0, vv] <- NA 
						}	
		res0 <- res <- res[ order( paste( 10000+ res$missing.pattern , 10000+ res$score)) , ]
		# calculate a summary
		index <- rep( seq(1,nrow(res0)) , res$freqwgt )
		res <- res[ index , ]
		a1 <- stats::aggregate( res[ , c("theta.dep" , "theta.ind")] , list( res$missing.pattern ,  res$score , res$max) , mean , na.rm=T )
		colnames(a1) <- c("missing.pattern" , "score" , "max" , "M.theta.dep" , "M.theta.ind" )
		a1$N <- stats::aggregate( 1+0*res[ , c("theta.dep")] , list( res$missing.pattern ,res$score, res$max) , sum  , na.rm=T )[,4]	
		a1 <- a1[ , c(1:3,6,4,5) ]
		a1$SD.theta.dep <- stats::aggregate( res[ , c("theta.dep")] , 
				list( res$missing.pattern ,res$score, res$max) , stats::sd  , na.rm=T)[,4]		
		a1$Min.theta.dep <- stats::aggregate( res[ , c("theta.dep")] , 
				list( res$missing.pattern ,res$score, res$max) , min  , na.rm=T)[,4]
		a1$Max.theta.dep <- stats::aggregate( res[ , c("theta.dep" )] , 
				list( res$missing.pattern ,res$score, res$max) , max , na.rm=T )[,4]
		if ( abs(alpha1) + abs( alpha2) > 0 ){
			a1$SD.theta.ind <- stats::aggregate( res[ , c("theta.ind")] , 
					list( res$missing.pattern ,res$score, res$max) , stats::sd  , na.rm=T)[,4]		
			a1$Min.theta.ind <- stats::aggregate( res[ , c("theta.ind")] , 
						list( res$missing.pattern ,res$score, res$max) , min  , na.rm=T )[,4]
			a1$Max.theta.ind <- stats::aggregate( res[ , c("theta.ind" )] , 
						list( res$missing.pattern ,res$score, res$max) , max  , na.rm=T)[,4]
					}		
		a1$M.seinflat <- stats::aggregate( res[ , c("seinflat")] , 
				list( res$missing.pattern ,res$score, res$max) , mean  , na.rm=T)[,4]
		a1$M.setheta.dep <- stats::aggregate( res[ , c("setheta.dep" )] , 
					list( res$missing.pattern ,res$score, res$max) , mean , na.rm=T )[,4]		
		a1$M.setheta.ind <- stats::aggregate( res[ , c("setheta.ind" )] , 
					list( res$missing.pattern ,res$score, res$max) , mean  , na.rm=T)[,4]		
		if (print.summary){
			cat("\n..................................................................\n")
			cat("Mean percentage standard error inflation\n\n")
			a4 <- stats::aggregate( res[,"seinflat"] , list( res$missing.pattern) , mean , na.rm=T)
			a4[,2] <- round( 100*(a4[,2] - 1) , 2 )
			colnames(a4) <- c("missing.pattern" , "Mperc.seinflat")
			print(a4)	
			cat("\n..................................................................\n")
			cat("Summary theta estimation\n\n")
			a1b <- a1
			a1b[,seq(5,ncol(a1))] <- round( a1[ , seq(5,ncol(a1)) ] , 4 )
			print(a1b)
				}
		# person parameter estimates		
				
		ind <- match( raschcopula.object$datalist$pattern.in.data , res0$pattern )		
 
		# results		
		res <- list( "person" = res0[ ind , ] , "se.inflat" = a4 , 
						"theta.table" = res0 , "pattern.in.data" = raschcopula.object$datalist$pattern.in.data ,
						"summary.theta.table" = a1)
        return(res)
           }
#########################################################################################



								
##################################################################################
# product of rows in a matrix m1 | bundlewise calculated by a vector v1
.rowProds.bundle <- function( m1 , v1 ){
	L1 <- length(v1)
	m1prod <- matrix( 0 , nrow=nrow(m1) , ncol= L1 )
	v1min <- c(1 , cumsum(v1)[ - L1 ]+1 )
	v1max <- cumsum(v1)
	for (ll in 1:L1){
		if (v1[ll] > 1 ){ 
				m1prod[,ll] <- rowProds( m1[ , v1min[ll]:v1max[ll] ] )
					} else { m1prod[ll] <- m1[ , v1min[ll] ] }
						}
	m1prod
		}
##################################################################################
.rowProds2.bundle <- function( m1 , v1 ){
	L1 <- length(v1)
	m1prod <- matrix( 0 , nrow=nrow(m1) , ncol= L1 )
	v1min <- c(1 , cumsum(v1)[ - L1 ]+1 )
	v1max <- cumsum(v1)
	for (ll in 1:L1){
		if (v1[ll] > 1 ){ 
				m1prod[,ll] <- rowProds2( m1[ , v1min[ll]:v1max[ll] ] )
					} else { m1prod[ll] <- m1[ , v1min[ll] ] }
						}
	m1prod
		}
#*********************************************************************************
##################################################################################
# product of rows in a matrix m1 | bundlewise calculated by a vector v1
.rowMins2.bundle <- function( m1 , v1 ){
	L1 <- length(v1)
	m1min <- matrix( 0 , nrow=nrow(m1) , ncol= L1 )
	v1min <- c(1 , cumsum(v1)[ - L1 ]+1 )
	v1max <- cumsum(v1)
	m1min[ , which(v1==1)] <- m1[ , v1min[ v1 == 1 ] ]
	for (ll in (1:L1)[ v1 > 1] ){
				m1min[,ll] <- rowMins2( m1[ , v1min[ll]:v1max[ll] ] )
						}
	m1min
		}
##################################################################################

.rowMins2cpp.bundle <- function( m1 , v1 ){
	m1min <- .Call( "rowmins2_bundle_C" , m1 , v1 , PACKAGE="sirt" )
	return(m1min)
		}



#************************************************************************
# Function rowProds2 
rowProds2 <- function(matr){
	y <- matr[,1]
	for (ii in 2:dim(matr)[2]){
		y <- y * matr[,ii] }
	return(y)
		}
#...................................................................
