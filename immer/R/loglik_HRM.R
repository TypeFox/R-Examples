
##############################################################
# log likelihood in the HRM
loglik_HRM <- function( dat , dat_ind , est_pars , theta_like ,
			rater , pid , maxK ){	
    useRcpp <- FALSE 			
	TP <- length(theta_like)
	b <- est_pars$b
	a <- est_pars$a
	phi <- est_pars$phi
	psi <- est_pars$psi
	K <- ncol(b)
	I <- nrow(b)
	N <- length( unique(pid))
	ND <- nrow(dat)
	dat <- as.matrix(dat)	
	#****
	# log likelihood (based on dyads)
	ll1 <- matrix( 1 , nrow=ND , ncol= TP )	
	for (tt in 1:TP){
		# tt <- 24 # theta	
		theta_tt <- theta_like[tt]
		theta0 <- rep(theta_tt,ND)		
		for (ii in 1:I){	
			# ii <- 1     # item ii
			phi_ii <- as.numeric( phi[ii,] )
			psi_ii <- as.numeric( psi[ii,] )
			K_ii <- maxK[ii]
			pr_tot <- rep(0,ND)			
			for (hh in seq(0,K_ii) ){
				# hh <- 1     # category hh for item ii and theta tt
				x0 <- rep(hh,ND)		
				pr1_hh <- probs_gpcm( x=x0 , theta= theta0 , b = as.numeric(b[ii,]) , a = a[ii] , 
							K = K_ii , useRcpp=useRcpp)		
				pr2_x <- probs_hrm( x = dat[,ii] , xi = x0 , 
							phi= phi_ii[ rater ]  , psi = psi_ii[ rater ] ,
								K = K_ii , useRcpp=useRcpp )
				pr_tot <- pr_tot + pr1_hh * pr2_x
						}
			ll1[,tt] <- ifelse( dat_ind[,ii] == 1 , ll1[,tt] * pr_tot , ll1[,tt] )
						}
				}
	#*****
	# aggregation for persons
	# individual likelihood
	eps <- 1E-200
	ll1 <- log( ll1 + eps )
	ll2 <- rowsum( ll1 , pid )	
	f.yi.qk <- exp( ll2 )

	#****
	# individual posterior
	pi.k <- stats::dnorm( theta_like , mean=est_pars$mu , sd = est_pars$sigma )
	pi.k <- pi.k / sum( pi.k )
	piM <- matrix( pi.k , nrow=N , ncol=TP , byrow=TRUE )
	f.qk.yi <- f.yi.qk * piM
	f.qk.yi <- f.qk.yi / rowSums(f.qk.yi)
	
	#***
	# computation of log-likelihood
	ll <- sum( log( rowSums( piM * f.yi.qk ) ) )
	
	#****
	# output
	res <- list( f.yi.qk = f.yi.qk , f.qk.yi=f.qk.yi , 
					theta_like = theta_like , pi.k = pi.k , ll = ll )
	return(res)	
				}