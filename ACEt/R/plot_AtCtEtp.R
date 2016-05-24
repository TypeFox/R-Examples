plot_AtCtEtp <- function(AtCtEtp_mcmc)
{
	if(class(AtCtEtp_mcmc)!='AtCtEtp_mc_model')
	{
		stop('The first parameter must be an object obtained from the acetp_mcmc function.')
	}
	
	model_cur <- AtCtEtp_mcmc

	#pheno_m <- c(t(data_m[,1:2]))
	#pheno_d <- c(t(data_d[,1:2]))
	#T_m <- rep(data_m[,3], each=2)
	#T_d <- rep(data_d[,3], each=2)
	p_n <- 500

	order <- 3
	x <- seq(from=model_cur$min_t, to=model_cur$max_t, length.out=p_n)
	t_int <- model_cur$max_t-model_cur$min_t
	l_m_1 <- (model_cur$max_t-x)/t_int
	l_m_2 <- (x-model_cur$min_t)/t_int
	
	n_a <- length(model_cur$beta_a_mc)
	n_c <- length(model_cur$beta_c_mc)
	n_e <- length(model_cur$beta_e_mc)

	if(n_a>2)
	{
		bb_a <- splineDesign(model_cur$knots_a, x = x, ord=order, outer.ok = TRUE)
	}else{
		if(n_a==2)
		{
			bb_a <- matrix(NA, p_n, 2)
			bb_a[,1] <- l_m_1
			bb_a[,2] <- l_m_2
		}else{
			bb_a <- matrix(1, p_n, 1)
		}
	}
	points_a <- exp(bb_a%*%model_cur$beta_a_mc)

	if(n_c>2)
	{
		bb_c <- splineDesign(model_cur$knots_c, x = x, ord=order, outer.ok = TRUE)
	}else{
		if(n_c==2)
		{
			bb_c <- matrix(NA, p_n, 2)
			bb_c[,1] <- l_m_1
			bb_c[,2] <- l_m_2
		}else{
			bb_c <- matrix(1, p_n, 1)
		}
	}
	points_c <- exp(bb_c%*%model_cur$beta_c_mc)

	if(n_e>2)
	{
		bb_e <- splineDesign(model_cur$knots_e, x = x, ord=order, outer.ok = TRUE)
	}else{
		if(n_e==2)
		{
			bb_e <- matrix(NA, p_n, 2)
			bb_e[,1] <- l_m_1
			bb_e[,2] <- l_m_2
		}else{
			bb_e <- matrix(1, p_n, 1)
		}
	}
	points_e <- exp(bb_e%*%model_cur$beta_e_mc)

	plot(range(x), c(0,max(points_c, points_a, points_e)+1), type = "n", xlab = "Age", ylab = "Variance",main =  "Variance curves of the A, C and E components")
	
	#bb <- splineDesign(model_cur$knots_a, x = x, ord=order, outer.ok = TRUE)
	lines(x, points_a, col = "red", lwd = 2)

	fisher_a <- model_cur$cov_mc[1:n_a,1:n_a]
	lower <- rep(NA, length(x))
	upper <- rep(NA, length(x))
	sd <- rep(NA, length(x))

	for(i in 1:length(x))
	{
		sd[i] <- sqrt((t(bb_a[i,])%*%fisher_a%*%bb_a[i,]))
		lower[i] <- sum(bb_a[i,]*model_cur$beta_a_mc) - 1.96*sd[i]
		upper[i] <- sum(bb_a[i,]*model_cur$beta_a_mc) + 1.96*sd[i]
	}

	lines(x, exp(lower), col = "orange" ,lty = 2 , lwd = 0.6)
	lines(x, exp(upper), col = "orange" ,lty = 2 , lwd = 0.6)
	polygon(c(x, rev(x)),c(exp(upper), rev(exp(lower))),col='grey',border = NA, lty=3, density=20)

	
	lines(x, points_c, col = "blue", lwd = 2)
	fisher_c <- model_cur$cov_mc[(n_a+1):(n_a+n_c),(n_a+1):(n_a+n_c)]
	lower <- rep(NA, length(x))
	upper <- rep(NA, length(x))
	sd <- rep(NA, length(x))

	for(i in 1:length(x))
	{
		sd[i] <- sqrt((t(bb_c[i,])%*%fisher_c%*%bb_c[i,]))
		lower[i] <- sum(bb_c[i,]*model_cur$beta_c_mc) - 1.96*sd[i]
		upper[i] <- sum(bb_c[i,]*model_cur$beta_c_mc) + 1.96*sd[i]
	}

	lines(x, exp(lower), col = "green" ,lty = 2 , lwd = 0.6)
	lines(x, exp(upper), col = "green" ,lty = 2 , lwd = 0.6)
	polygon(c(x, rev(x)),c(exp(upper), rev(exp(lower))),col='grey',border = NA, lty=3, density=20)

	lines(x, points_e, col = "pink", lwd = 2)

	fisher_e <- model_cur$cov_mc[(n_a+n_c+1):(n_a+n_c+n_e),(n_a+n_c+1):(n_a+n_c+n_e)]
	lower <- rep(NA, length(x))
	upper <- rep(NA, length(x))
	sd <- rep(NA, length(x))

	for(i in 1:length(x))
	{
		sd[i] <- sqrt((t(bb_e[i,])%*%fisher_e%*%bb_e[i,]))
		lower[i] <- sum(bb_e[i,]*model_cur$beta_e_mc) - 1.96*sd[i]
		upper[i] <- sum(bb_e[i,]*model_cur$beta_e_mc) + 1.96*sd[i]
	}

	lines(x, exp(lower), col = "yellow" ,lty = 2 , lwd = 0.6)
	lines(x, exp(upper), col = "yellow" ,lty = 2 , lwd = 0.6)
	polygon(c(x, rev(x)),c(exp(upper), rev(exp(lower))),col='grey',border = NA, lty=3, density=20)
	

	legend(x[1], max(points_c, points_a, points_e)+1, c('Additive genetic component','Common environmental component', 'Unique environmental component'), col = c('red','blue','pink'), lty=c(1,1,1), lwd=c(2,2,2))

}