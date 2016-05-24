acetp_mcmc <- function(acetp, iter_num = 10000, sd = 0.1, burnin =1000)
{
	if(!(class(acetp) %in% c('AtCtEp_model','AtEtp_model','AtCtEtp_model')))
	{
		stop('The first parameter must be an acetp object.')
	}

	
	if(burnin >= iter_num)
	{
		stop('The number of burnins must be smaller than the number of MCMC iterations.')
	}


	if(class(acetp)=='AtCtEp_model')
	{
		res <- AtCtEp_mcmc(acetp, iter_num, sd, burnin)
		return(res)
	}

	if(class(acetp)=='AtCtEtp_model')
	{
		res <- AtCtEtp_mcmc(acetp, iter_num, sd, burnin)
		return(res)
	}

	if(class(acetp)=='AtEtp_model')
	{
		res <- AtEtp_mcmc(acetp, iter_num, sd, burnin)
		return(res)
	}
	
}



AtCtEp_mcmc <-
function(AtCtEp, iter_num = 10000, sd = 0.1, burnin =1000)
{

if(class(AtCtEp)!='AtCtEp_model')
{
	stop('The first parameter must be an object obtained from the AtCtEp function.')
}

order <- 3

B_des_a_m <- splineDesign(AtCtEp$knot_a, x=AtCtEp$T_m, ord=order)
B_des_a_d <- splineDesign(AtCtEp$knot_a, x=AtCtEp$T_d, ord=order)
B_des_c_m <- splineDesign(AtCtEp$knot_c, x=AtCtEp$T_m, ord=order)
B_des_c_d <- splineDesign(AtCtEp$knot_c, x=AtCtEp$T_d, ord=order)

result <- mcmc_epsp(AtCtEp$pheno_m, AtCtEp$pheno_d, B_des_a_m, B_des_a_d, B_des_c_m, B_des_c_d, AtCtEp$var, AtCtEp$var_b_a, AtCtEp$var_b_c, AtCtEp$D_a, AtCtEp$D_c, iter_num, burnin, sd)

AtCtEp_mc_mod <- list(beta_a_mc=result$beta_a_mc, beta_c_mc=result$beta_c_mc, cov_a_mc = result$cov_a, cov_c_mc = result$cov_c, knots_a=AtCtEp$knot_a, knots_c=AtCtEp$knot_c, min_t = min(AtCtEp$T_m, AtCtEp$T_d), max_t = max(AtCtEp$T_m, AtCtEp$T_d))

class(AtCtEp_mc_mod) <- 'AtCtEp_mc_model'

return(AtCtEp_mc_mod)
}

AtCtEtp_mcmc <-
function(AtCtEtp, iter_num = 10000, sd = 0.1, burnin =1000)
{

if(class(AtCtEtp)!='AtCtEtp_model')
{
	stop('The first parameter must be an object obtained from the AtCtEtp function.')
}

T_m <- AtCtEtp$T_m
num_m <- length(T_m)
T_d <- AtCtEtp$T_d
num_d <- length(T_d)

t_int <- max(c(T_m,T_d))-min(c(T_m,T_d))
l_m_1 <- (max(c(T_m,T_d))-T_m)/t_int
l_m_2 <- (T_m-min(c(T_m,T_d)))/t_int
l_d_1 <- (max(c(T_m,T_d))-T_d)/t_int
l_d_2 <- (T_d-min(c(T_m,T_d)))/t_int

order <- 3
if(length(AtCtEtp$beta_a)>2)
{
	B_des_a_m <- splineDesign(AtCtEtp$knot_a, x=AtCtEtp$T_m, ord=order)
	B_des_a_d <- splineDesign(AtCtEtp$knot_a, x=AtCtEtp$T_d, ord=order)
}else{
	if(length(AtCtEtp$beta_a)==2)
	{
		B_des_a_m <- matrix(NA, num_m, 2)
		B_des_a_m[,1] <- l_m_1
		B_des_a_m[,2] <- l_m_2
		B_des_a_d <- matrix(NA, num_d, 2)
		B_des_a_d[,1] <- l_d_1
		B_des_a_d[,2] <- l_d_2
	}else{
		B_des_a_m <- matrix(1, num_m, 1)
		B_des_a_d <- matrix(1, num_d, 1)
	}
}
if(length(AtCtEtp$beta_c)>2)
{
	B_des_c_m <- splineDesign(AtCtEtp$knot_c, x=AtCtEtp$T_m, ord=order)
	B_des_c_d <- splineDesign(AtCtEtp$knot_c, x=AtCtEtp$T_d, ord=order)
}else{
	if(length(AtCtEtp$beta_c)==2)
	{
		B_des_c_m <- matrix(NA, num_m, 2)
		B_des_c_m[,1] <- l_m_1
		B_des_c_m[,2] <- l_m_2
		B_des_c_d <- matrix(NA, num_d, 2)
		B_des_c_d[,1] <- l_d_1
		B_des_c_d[,2] <- l_d_2
	}else{
		B_des_c_m <- matrix(1, num_m, 1)
		B_des_c_d <- matrix(1, num_d, 1)
	}
}
if(length(AtCtEtp$beta_e)>2)
{
	B_des_e_m <- splineDesign(AtCtEtp$knot_e, x=AtCtEtp$T_m, ord=order)
	B_des_e_d <- splineDesign(AtCtEtp$knot_e, x=AtCtEtp$T_d, ord=order)
}else{
	if(length(AtCtEtp$beta_e)==2)
	{
		B_des_e_m <- matrix(NA, num_m, 2)
		B_des_e_m[,1] <- l_m_1
		B_des_e_m[,2] <- l_m_2
		B_des_e_d <- matrix(NA, num_d, 2)
		B_des_e_d[,1] <- l_d_1
		B_des_e_d[,2] <- l_d_2
	}else{
		B_des_e_m <- matrix(1, num_m, 1)
		B_des_e_d <- matrix(1, num_d, 1)
	}
}

result <- mcmc_epsp_AtCtEt(AtCtEtp$pheno_m, AtCtEtp$pheno_d, B_des_a_m, B_des_a_d, B_des_c_m, B_des_c_d, B_des_e_m, B_des_e_d, AtCtEtp$var_b_a, AtCtEtp$var_b_c, AtCtEtp$var_b_e, AtCtEtp$D_a, AtCtEtp$D_c, AtCtEtp$D_e, iter_num, burnin, sd)

AtCtEtp_mc_mod <- list(beta_a_mc=result$beta_a_mc, beta_c_mc=result$beta_c_mc, beta_e_mc=result$beta_e_mc, cov_mc = result$cov, knots_a=AtCtEtp$knot_a, knots_c=AtCtEtp$knot_c, knots_e=AtCtEtp$knot_e, min_t = min(AtCtEtp$T_m, AtCtEtp$T_d), max_t = max(AtCtEtp$T_m, AtCtEtp$T_d))

class(AtCtEtp_mc_mod) <- 'AtCtEtp_mc_model'

return(AtCtEtp_mc_mod)
}

AtEtp_mcmc <-
function(AtEtp, iter_num = 10000, sd = 0.1, burnin =1000)
{

if(class(AtEtp)!='AtEtp_model')
{
	stop('The first parameter must be an object obtained from the AtEtp function.')
}

model_cur <- AtEtp

order <- 3

B_des_a_m <- splineDesign(model_cur$knot_a, x=model_cur$T_m, ord=order)
B_des_a_d <- splineDesign(model_cur$knot_a, x=model_cur$T_d, ord=order)
B_des_e_m <- splineDesign(model_cur$knot_e, x=model_cur$T_m, ord=order)
B_des_e_d <- splineDesign(model_cur$knot_e, x=model_cur$T_d, ord=order)

result <- mcmc_epsp_AtEt(model_cur$pheno_m, model_cur$pheno_d, B_des_a_m, B_des_a_d, B_des_e_m, B_des_e_d, model_cur$var_b_a, model_cur$var_b_e, model_cur$D_a, model_cur$D_e, iter_num, burnin, sd)

AtEtp_mc_mod <- list(beta_a_mc=result$beta_a_mc, beta_e_mc=result$beta_e_mc, cov_a_mc = result$cov_a, cov_e_mc = result$cov_e, knots_a=model_cur$knot_a, knots_e=model_cur$knot_e, min_t = min(AtEtp$T_m, AtEtp$T_d), max_t = max(AtEtp$T_m, AtEtp$T_d))

class(AtEtp_mc_mod) <- 'AtEtp_mc_model'

return(AtEtp_mc_mod)
}

