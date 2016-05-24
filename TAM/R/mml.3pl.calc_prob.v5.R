


#####################################################################
# calc_prob
# Calculation of probabilities
.mml.3pl.calc_prob.v5 <-
  function(iIndex, A, AXsi, B, xsi, theta, 
           nnodes, maxK, recalc=TRUE , guess){
		   
    if(recalc){
      AXsi.tmp <- array( tensor::tensor( A[iIndex,,, drop = FALSE], xsi, 3, 1 ) , 
                         dim = c( length(iIndex) , maxK , nnodes ) )
      AXsi[iIndex,] = AXsi.tmp[,,1]
    } else {
      AXsi.tmp <- array( AXsi, dim = c( length(iIndex) , maxK , nnodes ) )
    }

	
    Btheta <- array(0, dim = c(length(iIndex) , maxK , nnodes) )
    for( dd in 1:ncol(theta) ){  
      Btheta <- Btheta + array(B[iIndex,,dd ,drop = FALSE] %o% theta[,dd] , dim = dim(Btheta))
							}
							
    rprobs <- ( rr <- exp(Btheta+AXsi.tmp) )/aperm( array( rep( colSums( aperm( rr , c(2,1,3) ) ,
					dims=1 , na.rm = TRUE) ,    maxK ), dim=dim(rr)[c(1,3,2)] ) , c(1,3,2) )
					
	# include guessing	
	rprobs0 <- rprobs
	ind <- which(guess > 1E-6 )
	if ( length(ind) > 0 ){
		rprobs[ ind , 2 , ] <- guess[ind] + ( 1-guess[ind] ) * rprobs0[ind,2,]	
		# include guessing here
		rprobs[ ind , 1 , ] <- 1 - rprobs[ ind , 2 , ]
							}						
    return(list("rprobs" = rprobs, "AXsi" = AXsi , "rprobs0"=rprobs0 ))
  }
##############################################################################