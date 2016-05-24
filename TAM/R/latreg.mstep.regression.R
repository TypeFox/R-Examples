

###############################################
# latent regression mstep
latreg.mstep.regression <- function( hwt , 
	  pweights , pweightsM , Y , theta , theta2 , YYinv , ndim , 
	  nstud , beta.fixed , variance , Variance.fixed , group , G , 
	  snodes = 0 , thetasamp.density=NULL , nomiss=FALSE){
	    	  
	  # calculate item weights
	  variance.fixed <- Variance.fixed
#  a0 <- Sys.time()	

	#*****
	# numerical integration	
	if ( snodes == 0){	
		# hwt ... N x q matrix				 			
		thetabar <- hwt %*% theta
			## -----------
			## Functions from tensor package
			## %*t% 	Matrix product A %*% t(B)  -> tcrossprod	
			## %t*% 	Matrix product t(A) %*% B  -> crossprod	
			## %t*t% 	Matrix product t(A) %*% t(B)		
			##-----------
			## Functions base::crossprod and base::tcrossprod
			## t(x) %*% y (crossprod)
			## x %*% t(y) (tcrossprod). 
		# sumbeta <- Y %t*% ( thetabar*pweights )
		sumbeta <- base::crossprod( Y , ( thetabar*pweights ) )
        # sumsig2 <- as.vector( t( colSums( pweights * hwt ) ) %*% theta2 )		
		sumsig2 <- as.vector( base::crossprod( colSums( pweights * hwt ) , theta2 ) )		
 #cat("- sum sig2") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1		
					}

	
	#****
	# Monte Carlo integration

	if ( snodes > 0 ){
# cat("- start monte carlo ") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1								
		hwtS <- hwt
#		hwtS <- hwt / outer( rep(1,nrow(hwt) ) , thetasamp.density )
		hwtS <- hwtS / rowSums( hwtS )   # maybe this can be fastened
# cat("- pure matrix R ") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1				
		#*** included		
		# hwt ... N x q matrix
# cat("- calc itemwt ") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1				
		hwtS <- hwt		
		thetabar <- hwtS %*% theta
# cat("- thetabar ") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1						
		# sumbeta <- Y %t*% ( thetabar*pweights )
		sumbeta <- base::crossprod( Y ,  thetabar*pweights )
# cat("- tensor operation ") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1								
		sumsig2 <- as.vector( crossprod( colSums( pweights * hwtS ) , theta2 ) )
# cat("- sum sig2 modified ") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1	
					}				
	# calculation of variance and regression coefficients					
    beta <- YYinv%*%sumbeta                     #new beta
    sumsig2 <- matrix(sumsig2,ndim,ndim)
    if (G==1){ 
		variance <- (sumsig2- crossprod( sumbeta , beta) )/nstud  #new variance
			}
			
	# fixed beta coefficients
	if ( ! is.null( beta.fixed )){ 
		beta[ beta.fixed[,1:2,drop=FALSE] ] <- beta.fixed[,3] 
		beta <- as.matrix( beta , ncol=ndim )	
					}
	# fixed covariance matrix entries
	if ( ! is.null(variance.fixed) ){ 
		variance[ variance.fixed[,1:2,drop=FALSE] ] <- variance.fixed[,3]  	
		variance[ variance.fixed[,c(2,1),drop=FALSE] ] <- variance.fixed[,3]  		
				}		
	if ( G> 1){	# begin multiple groups
		if ( snodes > 0 ){ 
				hwt <- hwt / snodes 
#				hwt <- hwt / outer( rep(1,nrow(hwt) ) , thetasamp.density )	
				hwt <- hwt / rowSums( hwt )				
						}
		for (gg in 1:G){
			#gg <- 1
			ind.gg <- which( group == gg )
			thetabar <- hwt[ind.gg,]%*%theta
			# sumbeta <- Y[ind.gg,]%t*%( thetabar*pweights[ind.gg] )
			sumbeta <- base::crossprod( Y[ind.gg,] , thetabar*pweights[ind.gg] )
			# -- sumsig2 <- sum( (pweights*hwt) %*% theta2 )
			sumsig2 <- colSums((pweights[ind.gg]*hwt[ind.gg,]) %*% theta2)   
			sumsig2 <- matrix(sumsig2,ndim,ndim)
			variance[ind.gg] <- (sumsig2- base::crossprod( sumbeta , beta) ) / 
										sum(pweights[ind.gg]) #new variance
				}
			}		# end multiple groups
    res <- list( "beta" = beta , "variance" = variance  )
	return(res)
			}

############################################################################
############################################################################