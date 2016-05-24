BFGS_special <-
function(init, knobj, fun_like, verbose = FALSE){
  ## Same as BFGS, designed fo our setting
  ## Special stoping criterion and tracking of
  ## solver ability to solve or infinite or na values in gradient
  
  theta <- init
	names(theta) <- knobj$global_parameters$param_names

	k <- 0
	lk <- 1
	#theta <- init
	Bk <- -diag(length(init))
	Hk <- Bk

	temp_valf <- fun_like(theta, knobj, fail_incoming = TRUE)
	valf <- temp_valf$res
	fail <- temp_valf$fail
	next_it <- TRUE
	
	if(!fail){
		gradf_new <- compute_gradient(theta, fun_like, knobj = knobj)
		pk <- gradf_new
		fail <- fail | (sum(is.infinite(gradf_new)) != 0)
		fail <- fail | (sum(is.na(gradf_new)) != 0)
		next_it <- (sqrt(sum(gradf_new^2)) > knobj$global_parameters$tol)
	}

	next_it <- next_it & !fail

	while(next_it){
		k <- k+1
		theta_prev <- theta
		gradf <- gradf_new
		temp_valf <- fun_like(theta, knobj, fail_incoming = TRUE)
		valf <- temp_valf$res
		fail <- temp_valf$fail
	
		pk <- as.numeric(-Hk %*% gradf)
		names(pk) <- names(gradf)
	
		lk <- armijo(theta, fun_like, pk, gradf, valf, l = lk, beta = knobj$global_parameters$beta, c = knobj$global_parameters$c,  knobj = knobj)
		theta <- theta + lk * pk
		#theta[theta > 100] <- runif(0,100,n=1)
		#theta[theta < 0] <- runif(0,100,n=1)
		names(theta) <- names(init)
		sk <- lk * pk
		gradf_new <- compute_gradient(theta, fun_like, knobj = knobj)		
		#gradf_new[is.infinite(gradf_new)] <- 0
		yk <- gradf_new - gradf
	
		temp_denom <- sum(sk * yk)
		if ( temp_denom < 0 & !is.na(temp_denom)){
			temp1 <- Bk %*% sk
			Bk <- Bk - temp1 %*% t(temp1) / sum(sk * temp1) + yk %*% t(yk) / temp_denom
			temp2 <- Hk %*% yk
			Hk <- Hk + as.numeric(sum(sk * yk) + t(yk) %*% temp2) * (sk %*% t(sk)) / temp_denom^2 - (temp2 %*% t(sk) + sk %*% t(temp2)) / temp_denom
		}
	
		fail <- fail | (sum(is.infinite(gradf_new)) != 0)
		fail <- fail | (sum(is.na(gradf_new)) != 0)
		next_it <- ((sqrt(sum(gradf_new^2)) > knobj$global_parameters$tol) & (k <= knobj$global_parameters$max_it) & (lk > 1e-12) & !fail)
	
		temp_valf <- fun_like(theta, knobj, fail_incoming = TRUE)
		valf <- temp_valf$res
		fail <- temp_valf$fail | fail
	
		if(verbose){
		  #write.table(k, paste(home, "follow_up", sep=""))
			print(paste(sqrt(sum(gradf_new^2)),sqrt(sum(sk^2))))
			print(lk)
			print(paste("Iteration",k," : Objective",valf))
		}
	
		if (lk < 1e-7){
			Bk <- -diag(length(init))
			Hk <- Bk
		}
	}
	res <- c()
	res$theta <- theta
	res$fail <- fail
	res
}
