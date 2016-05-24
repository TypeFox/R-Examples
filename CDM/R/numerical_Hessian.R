
##############################################################################
# numerical computation of the Hessian matrix
numerical_Hessian <- function(par , FUN , h = 1E-5, gradient=FALSE){

	NP <- length(par)
	
	f0 <- FUN( x = par )
	fh <- rep(NA,NP)    # f(x+h)
	f2h <- rep(NA,NP)   # f(x+2*h)
	hess <- matrix(NA,nrow=NP , ncol=NP)
	fhh <- hess
	
#	if ( length(h) == 1 ){
#	   h <- rep(h,NP)
#		}
		
	#--- loop for computing f(x+h)
	for (ii in 1:NP){
		# ii <- 1
		par1 <- par
		par1[ii] <- par[ii] + h
		fh[ii] <- FUN( x = par1 )
		}
    
	#--- computation of the gradient
	if (gradient){
		res <- ( fh - f0 ) / h
				}
	
	#------
	# second partial derivatives
	# d F / dx dy
	# dF/dx = g(x,y) = ( F(x+h,y) - F(x,y) ) / h
	# (dF/dx)/dy = ( g(x,y+h) - g(x,y) ) / h
	#            = ( F(x+h,y+h) - F(x,y+h) - F(x+h,y) + F(x,y+h) ) / h^2 

	#---- hessian
	if (!gradient){
		
		fh1 <- matrix( fh , nrow=NP , ncol=NP, byrow=TRUE)
		fh2 <- matrix( fh , nrow=NP , ncol=NP, byrow=FALSE)
		
		#--- computation f(x+2*h)
		for (ii in 1:NP){
			par1 <- par
			par1[ii] <- par[ii] + 2*h
			f2h[ii] <- FUN( x = par1 )
			}	
		#--- computation f(x+h,y+h) 	
		for (ii in 1:NP){
		for (jj in 1:NP){
		if (ii < jj){
			par1 <- par
			par1[ii] <- par[ii] + h
			par1[jj] <- par[jj] + h
			fhh[ii,jj] <- fhh[jj,ii] <- FUN( x = par1 )
				}
			}	
			}
	    

		hess <- ( fhh - fh1 - fh2 + f0 ) / h^2
		diag(hess) <- ( f2h - 2*fh + f0)/ h^2		
		res <- hess
		}
	
	return(res)

		}
##############################################################################		