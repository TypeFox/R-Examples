
fuzcluster_estimate <- function(K , dat_m , dat_s , dat_resp ,
	maxiter=1000 , parmconv=.0001 , progress=TRUE , seed=NULL ,
	fac.oldxsi=0){

	dat.resp <- dat_resp
	# inits
	m1 <- colMeans( dat_m , na.rm=TRUE)
	s1 <- apply( dat_m , 2 , stats::sd , na.rm=TRUE )
    
	# items
	I <- ncol(dat_m)
	N <- nrow(dat_m)
    if ( is.null(seed) ){ seed <- round( stats::runif(1 ,1000, 9999 ) )}
	set.seed( seed )
	# initial mu estimate
	mu_est <- mvtnorm::rmvnorm( K , mean=m1 , sigma = 1.0*diag(s1) )
# mu_est[1,] <- 4
# mu_est[2,] <- 2

	# initial SD estimate
	sd_est <- t(sapply( 1:K , FUN = function(kk){ 
			stats::runif( I , 1.05*s1 , 2*s1 ) } ))
	# initial pi estimate
	pi_est <- stats::runif(K)
#	pi_est <- rep(1/K , K )
	pi_est <- pi_est / sum( pi_est )
    dev <- 0
	iter <- 0
	eps <- 10^(-10)
	parmchange <- 1
	
	mu_estmin <- mu_est
	sd_estmin <- sd_est
	pi_estmin <- pi_est
	dev_min <- 10^90
	itermin <- iter
		if (progress){
		   cat("****************************************************************\n")
				}
	#*************************************************
	# begin algorithm
	while ( ( iter < maxiter ) & ( parmchange > parmconv) ){
		mu_est0 <- mu_est
		sd_est0 <- sd_est	
		pi_est0 <- pi_est
		dev0 <- dev
		
		#*** E-step

		# matrix inits
		t_ik <- matrix( 0 , nrow=N , ncol=K)
		phi_ik <- xi_ik <- sig2_ik <- mu_ik <- array( 0 , dim=c(N,I,K) )

		for (kk in 1:K){
			# mu_k
			mu_k <- matrix( mu_est[kk,] , nrow=N , ncol=I , byrow=TRUE )
			# sd_k
			sd_k <- matrix( sd_est[kk,] , nrow=N , ncol=I , byrow=TRUE )
			# numerator
			num1 <- dat.resp*(dat_s^2 + sd_k^2 + eps)
			# mu estimate
			mu_ik[,,kk]  <- dat.resp*( mu_k * dat_s^2 + dat_m * sd_k^2  ) / num1
			# sigma2 estimate
			sig2_ik[,,kk] <- dat.resp*(dat_s^2 * sd_k^2 / num1)
			# xi estimate
			xi_ik[,,kk] <- mu_ik[,,kk]^2 + sig2_ik[,,kk]
			# t_ik: posterior distribution
			phi_ik[,,kk] <- stats::dnorm( dat.resp*( dat_m - mu_k ) / sqrt( sd_k^2 + dat_s^2 + eps) )
			t_ik[,kk] <- pi_est[kk] * rowProds2( phi_ik[,,kk] )
					}
		t_ik <- t_ik / rowSums( t_ik )

		#*** M-step

		# sum t_ik
		sum_tik <- colSums(t_ik)+eps

		# pi estimate
		pi_est <- colMeans(t_ik )+eps
#		pi_est[ pi_est < .001 ] <- .001
		pi_est <- pi_est / sum( pi_est )
		if ( fac.oldxsi > 0 ){
			# fac.oldxsi <- .75
			pi_est <- (1-fac.oldxsi) * pi_est + fac.oldxsi *pi_est0
			mu_est <- (1-fac.oldxsi) * mu_est + fac.oldxsi *mu_est0
			sd_est <- (1-fac.oldxsi) * sd_est + fac.oldxsi *sd_est0
					}
		# mu estimate
		for (kk in 1:K){
			# kk <- 1
			mu_est[kk,] <- colSums( mu_ik[,,kk] * t_ik[,kk] ) / sum_tik[kk]
						}

		# sd estimate
		for (kk in 1:K){
			# kk <- 1
			mu_k <- matrix( mu_est[kk,] , nrow=N , ncol=I , byrow=TRUE )
			h1 <- colSums( t_ik[,kk] * ( xi_ik[,,kk] - 2*mu_k*mu_ik[,,kk] + mu_k^2 ) ) / sum_tik[kk]
			h1 <- ifelse( h1 < eps , eps , sqrt(h1) )	
			sd_est[kk,] <- h1
						}
        sd_est[ sd_est < .01 ] <- .01
		# calculate deviance
		dev <- sum( colSums( t_ik ) * log( pi_est +eps) )
		dev <- dev - N*I / 2 * log(2*pi) 
		for (kk in 1:K){
			dev <- dev - sum( t_ik[,kk] * sum( log( sd_est[kk,] + eps) ) )
				}
		for (kk in 1:K){
			sd2_k <- matrix( sd_est[kk,]^2 , nrow=N , ncol=I , byrow=TRUE )
			mu_k <- matrix( mu_est[kk,] , nrow=N , ncol=I , byrow=TRUE )
			dev <- dev - .5 * sum( t_ik[,kk] * rowSums( ( xi_ik[,,kk] - 2*mu_k*mu_ik[,,kk] + mu_k^2 ) / (sd2_k+eps) ) )
						}
		dev <- -2*dev
		iter <- iter+1
		devchange <- abs( dev - dev0 )
		devchange1 <- - ( dev - dev0 )
		parmchange <- max( abs( mu_est - mu_est0 ) , abs( sd_est - sd_est0) )
		if (progress){
		   cat("............ seed =" , seed , " | Iteration",iter , ".............\n")
		   cat(paste("Deviance =", round( dev,3) ))
		   cat(paste(" | Deviance Change =" , round(devchange1,3) , "\n"))
		   cat("Class Probabilities =", paste(round( pi_est,4) ))		   
		   cat(paste(" | Maximum Parameter Change =" , round(parmchange,3) ,"\n"))
		   utils::flush.console()	
				}
		if ( dev < dev_min){
			t_ikmin <- t_ik
			mu_estmin <- mu_est
			sd_estmin <- sd_est
			pi_estmin <- pi_est
			itermin <- iter
			dev_min <- dev
			}			
	}
	#*************** end algorithm
	res <- list( "deviance" = dev_min , "iter" = itermin , "pi_est" = pi_estmin ,
		"mu_est"=mu_estmin , "sd_est" = sd_estmin , "posterior" = t_ikmin ,
		"seed"=seed 
		 )
	return(res)
		}