
#*********************************************
# Clustering for continuous fuzzy data
fuzcluster <- function(dat_m , dat_s , K=2 , nstarts= 7 , 
        seed = NULL , maxiter=100 , parmconv=.001 , 
		fac.oldxsi=0.75 , progress=TRUE ){
	s1 <- Sys.time()
	dev0 <- 10^(200)
    dat_resp <- 1 - is.na(dat_m)
	if ( ! is.null(seed) ){ nstarts <- 1 }
	for ( rr in 1:nstarts ){
		res1 <- fuzcluster_estimate(K , dat_m , dat_s , dat_resp ,
			maxiter=maxiter , parmconv=parmconv , progress=progress ,
			seed= seed , fac.oldxsi=fac.oldxsi)
		if ( res1$deviance < dev0 ){ 
			res <- res1 
			dev0 <- res1$deviance
			}
			}
	s2 <- Sys.time()
	### end random starts
	### computation of information criteria
	dev <- res$deviance
    ic <- list( "deviance" = dev , "n" = nrow(dat_m) )
	I <- ncol(dat_m)
	ic$np <- (K-1) + 2*K*I
	# AIC
        ic$AIC <- dev + 2*ic$np
        # BIC
        ic$BIC <- dev + ( log(ic$n) )*ic$np
        # CAIC (conistent AIC)
        ic$CAIC <- dev + ( log(ic$n) + 1 )*ic$np
		# corrected AIC
        ic$AICc <- ic$AIC + 2*ic$np * ( ic$np + 1 ) / ( ic$n - ic$np - 1 )		
	res$ic <- ic
	res$K <- K
	res$s1 <- s1
	res$s2 <- s2
	res$nstarts <- nstarts
	class(res) <- "fuzcluster"
	###
	return(res)	
	}