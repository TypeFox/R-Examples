
################################################
# estimating polychoric correlation using the Olsson method
# of maximum likelihood estimation
polychoric2 <- function( dat , maxiter = 100 , cor.smooth = TRUE ){
	dat1 <- as.matrix(dat)
	# maximum number of categories
	maxK <- max(dat1 , na.rm=TRUE )
	# perform polychoric correlation
	res0 <- .Call("polychoric2_aux_rcpp" ,
				dat1 , maxK , maxiter , PACKAGE="sirt")
	# include cpp Call here
	#***
	# output cleaning
	# thresholds
	tau <- res0$thresh[,c(2:(maxK+1))]
	tau[ tau == 99 ] <- Inf
	rho <- res0$rho
	rownames(rho) <- colnames(rho) <- colnames(dat1)
	if ( maxK > 1 ){ 
		rownames(tau) <- rownames(rho)
		colnames(tau) <- paste0("Cat" , 1:maxK)
				}
	if ( maxK == 1 ){
        names(tau) <- rownames(rho) 
					}				
	# handling missing entries in rho
	rho[ res0$Nobs == 0 ] <- NA 
	diag(rho) <- 1
	if ( sum( is.na(rho) > 0 ) ){ cor.smooth <- FALSE }
	if (cor.smooth){ 
			rho <- psych::cor.smooth(rho) 
					}
	# output list
	res <- list("tau"=tau , "rho"= rho , "Nobs"=res0$Nobs , "maxcat" = res0$maxcat )
	return(res)
	}
#########################################################################
