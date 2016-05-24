###############################################
# mstep.regression
mstep.regression <-
function( resp , hwt ,  resp.ind , 
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
		
		if (!nomiss){
		#	itemwt <- t( hwt ) %*% ( resp.ind * pweightsM )
			itemwt <- crossprod( hwt , resp.ind * pweightsM  )			
					}
		if ( nomiss ){
			itemwt <- matrix( colSums(hwt*pweights) , nrow=ncol(hwt) , ncol=ncol(resp.ind) )
					}						 			
 		# original implementation without missings
		#  -- itemwt0 <- matrix(rep(colSums(hwt), nitems), nrow=nnodes, ncol=nitems)
		thetabar <- hwt %*% theta
		# sumbeta <- Y %t*% ( thetabar*pweights )
		sumbeta <- crossprod( Y , thetabar*pweights )
		# -- sumsig2 <- sum( (pweights*hwt) %*% theta2 )
		# sumsig2 <- colSums((pweights*hwt) %*% theta2)
        sumsig2 <- as.vector( crossprod( colSums( pweights * hwt ) , theta2 ) )		
 #cat("- sum sig2") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1		
					}

	#****
	# Monte Carlo integration

	if ( snodes > 0 ){
# cat("- start monte carlo ") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1								
#		hwtS <- hwt / snodes
		hwtS <- hwt
		
#****		
#		hwtS <- hwt / outer( rep(1,nrow(hwt) ) , thetasamp.density )

		hwtS <- hwtS / rowSums( hwtS )   # maybe this can be fastened
# Revalpr("rowSums(hwtS)")		
		
# cat("- pure matrix R ") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1				
		#*** included		
		# hwt ... N x q matrix
#		itemwt <- t( hwtS ) %*% ( resp.ind * pweightsM )
		itemwt <- crossprod( hwtS , resp.ind * pweightsM  )
		# make this formula easier and make some speed checks!!
# cat("- calc itemwt ") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1				

#		hwtS <- hwt			
		
		#  -- itemwt0 <- matrix(rep(colSums(hwt), nitems), nrow=nnodes, ncol=nitems)
		# This formula is not faster!
		thetabar <- hwtS %*% theta
# cat("- thetabar ") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1						
		#sumbeta <- Y %t*% ( thetabar*pweights )
		sumbeta <- crossprod( Y ,  thetabar*pweights )
# cat("- tensor operation ") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1								
		# -- sumsig2 <- sum( (pweights*hwt) %*% theta2 )
		# sumsig2 <- colSums((pweights*hwtS) %*% theta2)      		
		# sumsig2 <- as.vector( t( colSums( pweights * hwtS ) ) %*% theta2 )
		sumsig2 <- as.vector( crossprod( colSums( pweights * hwtS ) , theta2 ) )
# cat("- sum sig2 modified ") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1	
					}				
	# calculation of variance and regression coefficients					
    beta <- YYinv%*%sumbeta                     #new beta
    sumsig2 <- matrix(sumsig2,ndim,ndim)
    if (G==1){ 
		variance <- (sumsig2- crossprod( sumbeta , beta ) )/nstud  #new variance
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
	#		sumbeta <- Y[ind.gg,]%t*%( thetabar*pweights[ind.gg] )
			sumbeta <- crossprod( Y[ind.gg,] ,  thetabar*pweights[ind.gg] )
			# -- sumsig2 <- sum( (pweights*hwt) %*% theta2 )
			sumsig2 <- colSums((pweights[ind.gg]*hwt[ind.gg,]) %*% theta2)   
			sumsig2 <- matrix(sumsig2,ndim,ndim)
			# variance[ind.gg] <- (sumsig2-sumbeta%t*%beta)/sum(pweights[ind.gg]) #new variance
			variance[ind.gg] <- (sumsig2- crossprod( sumbeta , beta) )/sum(pweights[ind.gg]) #new variance
				}
			}		# end multiple groups
    res <- list( "beta" = beta , "variance" = variance , "itemwt" = itemwt )
			}

############################################################################
############################################################################