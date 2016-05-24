test_AtCtEtp <-
function(acetp, comp, sim = 200)
{
	if(!(class(acetp) %in% c('AtCtEtp_model')))
	{
		stop('The first parameter must be an acetp object.')
	}

	if(!(comp %in% c('a','c','e')))
	{
		stop('The variable \'comp\' must be \'a\',\'c\' or \'e\' to specify which component to test linearity.')
	}

	if(acetp$mod[match(comp,c('a','c','e'))]!='d')
	{
		stop('The component to test is linear or constant.')
	}

	mod_a <- acetp$mod
	mod_n <- mod_a
	mod_n[match(comp,c('a','c','e'))] <- 'l' 
	
	k_a <- length(acetp$beta_a)-1
	k_c <- length(acetp$beta_c)-1
	k_e <- length(acetp$beta_e)-1
	data_m <- cbind(acetp$pheno_m[seq(from=1,to=nrow(acetp$pheno_m),by=2)],acetp$pheno_m[seq(from=2,to=nrow(acetp$pheno_m),by=2)])
	data_m <- cbind(data_m,acetp$T_m[seq(from=1,to=length(acetp$T_m),by=2)])
	data_d <- cbind(acetp$pheno_d[seq(from=1,to=nrow(acetp$pheno_d),by=2)],acetp$pheno_d[seq(from=2,to=nrow(acetp$pheno_d),by=2)])
	data_d <- cbind(data_d,acetp$T_d[seq(from=1,to=length(acetp$T_d),by=2)])	
	

	m_a <- AtCtEtp(data_m, data_d, knot_a=k_a, knot_c=k_c, knot_e=k_e, mod=mod_a)

	m_n <- AtCtEtp(data_m, data_d, knot_a=k_a, knot_c=k_c, knot_e=k_e, mod=mod_n)

	llr <- m_a$lik - m_n$lik

	beta_a_n <- m_n$beta_a
	beta_c_n <- m_n$beta_c
	beta_e_n <- m_n$beta_e
	var_b_a_n <- m_n$var_b_a
	var_b_c_n <- m_n$var_b_c
	var_b_e_n <- m_n$var_b_e

	order <- 3
	max_t <- max(m_n$T_m,m_n$T_d)
	min_t <- min(m_n$T_m,m_n$T_d)
	t_int <- max_t-min_t

	l_m_1 <- (max_t-m_n$T_m)/t_int
	l_m_2 <- (m_n$T_m-min_t)/t_int
	l_d_1 <- (max_t-m_n$T_d)/t_int
	l_d_2 <- (m_n$T_d-min_t)/t_int

	num_m <- length(m_n$T_m)
	num_d <- length(m_n$T_d)

	if(m_n$mod[1]=='d')
	{
		bb_a_m <- splineDesign(m_n$knot_a, x = m_n$T_m, ord=order, outer.ok = TRUE)
		a_m <- exp(bb_a_m%*%beta_a_n)
		bb_a_d <- splineDesign(m_n$knot_a, x = m_n$T_d, ord=order, outer.ok = TRUE)
		a_d <- exp(bb_a_d%*%beta_a_n)
	}else{
		bb_a_m <- matrix(NA, num_m, 2)
		bb_a_m[,1] <- l_m_1
		bb_a_m[,2] <- l_m_2
		bb_a_d <- matrix(NA, num_d, 2)
		bb_a_d[,1] <- l_d_1
		bb_a_d[,2] <- l_d_2
		a_m <- exp(beta_a_n[1]+(beta_a_n[2]-beta_a_n[1])/(max_t-min_t)*(m_n$T_m-min_t))
		a_d <- exp(beta_a_n[1]+(beta_a_n[2]-beta_a_n[1])/(max_t-min_t)*(m_n$T_d-min_t))
	}
	if(m_n$mod[2]=='d')
	{
		bb_c_m <- splineDesign(m_n$knot_c, x = m_n$T_m, ord=order, outer.ok = TRUE)
		c_m <- exp(bb_c_m%*%beta_c_n)
		bb_c_d <- splineDesign(m_n$knot_c, x = m_n$T_d, ord=order, outer.ok = TRUE)
		c_d <- exp(bb_c_d%*%beta_c_n)
	}else{
		bb_c_m <- matrix(NA, num_m, 2)
		bb_c_m[,1] <- l_m_1
		bb_c_m[,2] <- l_m_2
		bb_c_d <- matrix(NA, num_d, 2)
		bb_c_d[,1] <- l_d_1
		bb_c_d[,2] <- l_d_2
		c_m <- exp(beta_c_n[1]+(beta_c_n[2]-beta_c_n[1])/(max_t-min_t)*(m_n$T_m-min_t))
		c_d <- exp(beta_c_n[1]+(beta_c_n[2]-beta_c_n[1])/(max_t-min_t)*(m_n$T_d-min_t))
	}
	if(m_n$mod[3]=='d')
	{
		bb_e_m <- splineDesign(m_n$knot_e, x = m_n$T_m, ord=order, outer.ok = TRUE)
		e_m <- exp(bb_e_m%*%beta_e_n)
		bb_e_d <- splineDesign(m_n$knot_e, x = m_n$T_d, ord=order, outer.ok = TRUE)
		e_d <- exp(bb_e_d%*%beta_e_n)
	}else{
		bb_e_m <- matrix(NA, num_m, 2)
		bb_e_m[,1] <- l_m_1
		bb_e_m[,2] <- l_m_2
		bb_e_d <- matrix(NA, num_d, 2)
		bb_e_d[,1] <- l_d_1
		bb_e_d[,2] <- l_d_2
		e_m <- exp(beta_e_n[1]+(beta_e_n[2]-beta_e_n[1])/(max_t-min_t)*(m_n$T_m-min_t))
		e_d <- exp(beta_e_n[1]+(beta_e_n[2]-beta_e_n[1])/(max_t-min_t)*(m_n$T_d-min_t))
	}

	if(comp=='a')
	{
		bb_a_m_a <- splineDesign(m_a$knot_a, x = m_a$T_m, ord=order, outer.ok = TRUE)
		bb_a_d_a <- splineDesign(m_a$knot_a, x = m_a$T_d, ord=order, outer.ok = TRUE)
		bb_c_m_a <- bb_c_m
		bb_c_d_a <- bb_c_d
		bb_e_m_a <- bb_e_m
		bb_e_d_a <- bb_e_d
	}
	if(comp=='c')
	{
		bb_c_m_a <- splineDesign(m_a$knot_c, x = m_a$T_m, ord=order, outer.ok = TRUE)
		bb_c_d_a <- splineDesign(m_a$knot_c, x = m_a$T_d, ord=order, outer.ok = TRUE)
		bb_a_m_a <- bb_a_m
		bb_a_d_a <- bb_a_d
		bb_e_m_a <- bb_e_m
		bb_e_d_a <- bb_e_d
	}
	if(comp=='e')
	{
		bb_e_m_a <- splineDesign(m_a$knot_e, x = m_a$T_m, ord=order, outer.ok = TRUE)
		bb_e_d_a <- splineDesign(m_a$knot_e, x = m_a$T_d, ord=order, outer.ok = TRUE)
		bb_c_m_a <- bb_c_m
		bb_c_d_a <- bb_c_d
		bb_a_m_a <- bb_a_m
		bb_a_d_a <- bb_a_d
	}

	sim_m <- data_m
	sim_d <- data_d

	low_a <- -20
	upp_a <- 20
	low_c <- -20
	upp_c <- 20
	low_e <- -10
	upp_e <- 20
 
	llr_sim <- rep(NA, sim)
 
	for(i in 1:sim)
	{
		for(j in 1:(num_m/2))
		{
			sigma <- matrix(c(a_m[j]+c_m[j]+e_m[j],a_m[j]+c_m[j],a_m[j]+c_m[j],a_m[j]+c_m[j]+e_m[j]),2,2)
			sim_m[j,1:2] <- mvrnorm(1, rep(0,2), sigma)
		}
		for(j in 1:(num_d/2))
		{
			sigma <- matrix(c(a_d[j]+c_d[j]+e_d[j],0.5*a_d[j]+c_d[j],0.5*a_d[j]+c_d[j],a_d[j]+c_d[j]+e_d[j]),2,2)
			sim_d[j,1:2] <- mvrnorm(1, rep(0,2), sigma)
		}

		n_a <- ncol(bb_a_m)
		n_c <- ncol(bb_c_m)
		n_e <- ncol(bb_e_m)

		sim_m_l <- rep(0,nrow(bb_a_m))
		sim_m_l[seq(from=1, to=nrow(bb_a_m),by=2)] <- sim_m[,1]
		sim_m_l[seq(from=2, to=nrow(bb_a_m),by=2)] <- sim_m[,2]
		sim_d_l <- rep(0,nrow(bb_a_d))
		sim_d_l[seq(from=1, to=nrow(bb_a_d),by=2)] <- sim_d[,1]
		sim_d_l[seq(from=2, to=nrow(bb_a_d),by=2)] <- sim_d[,2]

		beta_a_t <- beta_a_n
		beta_c_t <- beta_c_n
		beta_e_t <- beta_e_n
		result <- optim(c(beta_a_t,beta_c_t,beta_e_t), loglik_AtCtEt_epsp_g, gr_AtCtEt_epsp_g, pheno_m = sim_m_l, pheno_d = sim_d_l, B_des_a_m=bb_a_m, B_des_a_d=bb_a_d, B_des_c_m=bb_c_m, B_des_c_d=bb_c_d, B_des_e_m=bb_e_m, B_des_e_d=bb_e_d, var_b_a=var_b_a_n, var_b_c=var_b_c_n, var_b_e=var_b_e_n, D_a=m_n$D_a, D_c=m_n$D_c, D_e=m_n$D_e, lower = c(rep(low_a,n_a),rep(low_c,n_c),rep(low_e,n_e)), upper = c(rep(upp_a,n_a),rep(upp_c,n_c),rep(upp_e,n_e)), method = "L-BFGS-B", control=list(maxit = 3000))
		for(j in 1:2)
		{
			init <- c(beta_a_t,beta_c_t,beta_e_t) + runif(n_a+n_c+n_e,min=-0.2,max=0.2)
			result_r <- optim(init, loglik_AtCtEt_epsp_g, gr_AtCtEt_epsp_g, pheno_m = sim_m_l, pheno_d = sim_d_l, B_des_a_m=bb_a_m, B_des_a_d=bb_a_d, B_des_c_m=bb_c_m, B_des_c_d=bb_c_d, B_des_e_m=bb_e_m, B_des_e_d=bb_e_d, var_b_a=var_b_a_n, var_b_c=var_b_c_n, var_b_e=var_b_e_n, D_a=m_n$D_a, D_c=m_n$D_c, D_e=m_n$D_e, lower = c(rep(low_a,n_a),rep(low_c,n_c),rep(low_e,n_e)), upper = c(rep(upp_a,n_a),rep(upp_c,n_c),rep(upp_e,n_e)), method = "L-BFGS-B", control=list(maxit = 3000))
			if(result_r$value < result$value)
			{
				result <- result_r
			}
		}
		result <- loglik_AtCtEt_epsp(c(var_b_a_n,var_b_c_n,var_b_e_n), pheno_m=sim_m_l, pheno_d=sim_d_l, B_des_a_m=bb_a_m, B_des_a_d=bb_a_d, beta_a=result$par[1:n_a], D_a=m_n$D_a, B_des_c_m=bb_c_m, B_des_c_d=bb_c_d, beta_c=result$par[(n_a+1):(n_a+n_c)], D_c=m_n$D_c, B_des_e_m=bb_e_m, B_des_e_d=bb_e_d, beta_e=result$par[(n_a+n_c+1):(n_a+n_c+n_e)], D_e=m_n$D_e)

		var_b_a_1 <- m_a$var_b_a
		var_b_c_1 <- m_a$var_b_c
		var_b_e_1 <- m_a$var_b_e
		beta_a_t <- m_a$beta_a
		beta_c_t <- m_a$beta_c
		beta_e_t <- m_a$beta_e
		n_a <- ncol(bb_a_m_a)
		n_c <- ncol(bb_c_m_a)
		n_e <- ncol(bb_e_m_a)

		result_2 <- optim(c(beta_a_t,beta_c_t,beta_e_t), loglik_AtCtEt_epsp_g, gr_AtCtEt_epsp_g, pheno_m = sim_m_l, pheno_d = sim_d_l, B_des_a_m=bb_a_m_a, B_des_a_d=bb_a_d_a, B_des_c_m=bb_c_m_a, B_des_c_d=bb_c_d_a, B_des_e_m=bb_e_m_a, B_des_e_d=bb_e_d_a, var_b_a=var_b_a_1, var_b_c=var_b_c_1, var_b_e=var_b_e_1, D_a=m_a$D_a, D_c=m_a$D_c, D_e=m_a$D_e, lower = c(rep(low_a,n_a),rep(low_c,n_c),rep(low_e,n_e)), upper = c(rep(upp_a,n_a),rep(upp_c,n_c),rep(upp_e,n_e)), method = "L-BFGS-B", control=list(maxit = 3000))
		for(j in 1:2)
		{
			init <- c(beta_a_t,beta_c_t,beta_e_t) + runif(n_a+n_c+n_e,min=-0.2,max=0.2)
			result_2r <- optim(init, loglik_AtCtEt_epsp_g, gr_AtCtEt_epsp_g, pheno_m = sim_m_l, pheno_d = sim_d_l, B_des_a_m=bb_a_m_a, B_des_a_d=bb_a_d_a, B_des_c_m=bb_c_m_a, B_des_c_d=bb_c_d_a, B_des_e_m=bb_e_m_a, B_des_e_d=bb_e_d_a, var_b_a=var_b_a_1, var_b_c=var_b_c_1, var_b_e=var_b_e_1, D_a=m_a$D_a, D_c=m_a$D_c, D_e=m_a$D_e, lower = c(rep(low_a,n_a),rep(low_c,n_c),rep(low_e,n_e)), upper = c(rep(upp_a,n_a),rep(upp_c,n_c),rep(upp_e,n_e)), method = "L-BFGS-B", control=list(maxit = 3000))
			if(result_2r$value < result_2$value)
			{
				result_2 <- result_2r
			}
		}
		result_2 <- loglik_AtCtEt_epsp(c(var_b_a_1,var_b_c_1,var_b_e_1), pheno_m=sim_m_l, pheno_d=sim_d_l, B_des_a_m=bb_a_m_a, B_des_a_d=bb_a_d_a, beta_a=result_2$par[1:n_a], D_a=m_a$D_a, B_des_c_m=bb_c_m_a, B_des_c_d=bb_c_d_a, beta_c=result_2$par[(n_a+1):(n_a+n_c)], D_c=m_a$D_c, B_des_e_m=bb_e_m_a, B_des_e_d=bb_e_d_a, beta_e=result_2$par[(n_a+n_c+1):(n_a+n_c+n_e)], D_e=m_a$D_e)
		llr_sim[i] <- result_2/2-result/2
 
	}

	p <- sum(llr_sim<llr)/sim
	return(p)
}