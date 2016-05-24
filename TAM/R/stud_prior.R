


#############################################################
#############################################################
stud_prior.v2 <-
  function(theta , Y , beta , variance , nstud , 
           nnodes , ndim , YSD , unidim_simplify , snodes=0 ){	      
# unidim_simplify <- FALSE		   
    if(ndim == 1) {
      ##################################
      # SINGLE DIMENSION
	if ( ! unidim_simplify  ){	  
		gwt <- matrix( stats::dnorm(rep(theta, each = nstud), mean= Y%*%beta, sd = sqrt(variance)),
                    nrow = nstud)
						} else {
	  TP <- nrow(theta)				
      gwt <- matrix(   stats::dnorm(theta[,1] , mean= 0, sd = sqrt(variance[1,1])),  
							nrow = nstud , ncol=TP , byrow=TRUE)					
								}
	  
    } else {
      ###################################
      # MULTIPLE DIMENSIONS
      mu <- Y%*%beta     #mean vector for each student: dimensions nstud by ndim 
      eps <- 10^(-7)
	  # eps <- 1E-3
	  
	    
	  
	# @@@ ARb 2014-10-19: Stabilization of the covariance matrix
		svd_var <- svd(variance)
		d0 <- d <- svd_var$d		
		eps2 <- .05
		# eps2 <- 1E-8
		ind <- which( d < eps2)
		if (length(ind)>0){
			d[ ind ] <- eps2
			d <- d / sum(d) * sum(d0)
			variance <- svd_var$u %*% diag(d) %*% t( svd_var$v  )
							}


      #	variance[ diag(variance) ] <- diag(variance) + eps
#      diag(variance) <- diag(variance) + eps
      varInverse <- solve(variance)
	  detvar <- det(variance)

      coeff <- 1/sqrt( (2*pi)^ndim * detvar ) 
	  # coeff <- 1 
      gwt <- matrix( 0 , nrow=nstud , ncol=nnodes )  
	  

					
						
#  cat(" * prior Ysd") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1	
    if ( YSD ){
	  if (snodes >= 0 ){	
		gwt <- prior.normal.density.R( theta_=theta , mu_=mu , 
			     varInverse_=varInverse ,  coeff_=coeff) 
						}
						
		if ( snodes < 0 ){ 		
# Revalpr("cbind( theta[ rep(1:snodes , nstud ) , ] , mu[ rep(1:nstud , each=snodes)   , ])")		
			gwt <- dmvnorm_TAM( x = theta[ rep(1:snodes , nstud ) , ] , 
							mean = mu[ rep(1:nstud , each=snodes)   , ] , sigma = variance )	
			gwt <- matrix( gwt , nrow=nstud , byrow=TRUE)
# gwt <- matrix( gwt , nrow=nstud , ncol=nnodes , byrow=TRUE )			
# gwt <- gwt / matrix( rowSums(gwt) , nrow=nstud , ncol=nnodes , byrow=TRUE )
				}
			 
			   }	
#  cat(" * prior nnodes2") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1						
    if ( ! YSD){
	  if (snodes >= 0 ){
		gwt <- prior.normal.densityALL.R( theta_=theta , mu_=mu , 
			           varInverse_=varInverse ,  coeff_=coeff) 			 
					}				
	  if (snodes < 0 ){
#       gwt1 <- mvtnorm::dmvnorm( theta , mean=as.vector(mu[1,1:2]) , sigma=variance )	  
		gwt <- dmvnorm_TAM( x=theta , mean= as.vector(mu[1,1:2]) , sigma=variance )
						}	
		# gtw <- gwt / sum(gwt) 
		gwt <- matrix( gwt , nrow=nstud , ncol=nnodes , byrow=TRUE )
				}

# gwt <- gwt / matrix( rowSums(gwt) , nrow=nstud , ncol=nnodes , byrow=TRUE )			
				
#  cat(" * prior nnodes3") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1						
    }
	
    return(gwt)
  }
#####################################################################  
  
#*************************
# auxiliary functions for calculation of prior functions  
prior.normal.density.R <- 
function( theta_ , mu_ , varInverse_ , coeff_){
	.Call("prior_normal_density_C",  theta_ , mu_ , 
			varInverse_ , coeff_ , PACKAGE = "TAM")
} 

prior.normal.densityALL.R <- 
 function( theta_ , mu_ , varInverse_ , coeff_){
 	.Call("prior_normal_densityALL_C",  theta_ , mu_ , 
		varInverse_  , coeff_ , PACKAGE = "TAM")
 } 
