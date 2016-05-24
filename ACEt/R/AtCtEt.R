AtCtEt <-
function(data_m, data_d, mod = c('d','d','d'), knot_a=5, knot_c=5, knot_e=5, boot=FALSE, num_b = 100, init = rep(0,3), robust = 0)
{

pheno_m <- c(t(data_m[,1:2]))
pheno_d <- c(t(data_d[,1:2]))
T_m <- rep(data_m[,3], each=2)
T_d <- rep(data_d[,3], each=2)

mag <- var(pheno_m)
init_max <- log(mag)
init_min <- log(mag) - abs(log(mag))*1.3

if((is.vector(mod)==FALSE) | (length(mod)!=3) )
{stop('The model parameter must be a vector of length 3.')}

if(!(mod[1] %in% c('d','c','n')))
{stop('The \'mod\' parameter for the A component must be \'d\'(dynamic), \'c\'(constant) or \'n\'(NA).')}

if(!(mod[2] %in% c('d','c','n')))
{stop('The \'mod\' parameter for the C component must be \'d\'(dynamic), \'c\'(constant) or \'n\'(NA).')}

if(!(mod[3] %in% c('d','c')))
{stop('The \'mod\' parameter for the E component must be \'d\'(dynamic), \'c\'(constant).')}

order <- 3
if(mod[1]=='d')
{
	order <- 3

	if(knot_a < 3)
	{stop('The number of interior knots must be no less than 3.')}

}else
{
	order <- 1
}
#knot <- 8
if(mod[1]=='d')
{
	knots_a <- seq(from=min(T_m, T_d), to=max(T_m, T_d), length.out=knot_a)
	interval_a <- knots_a[2] - knots_a[1]
	knots_a <- c(c(min(T_m, T_d)-interval_a*2,min(T_m, T_d)-interval_a), knots_a)
	knots_a <- c(knots_a, c(max(T_m, T_d)+interval_a,max(T_m, T_d)+interval_a*2))
	B_des_a_m <- splineDesign(knots_a, x=T_m, ord=order)
	B_des_a_d <- splineDesign(knots_a, x=T_d, ord=order)
}else{
	knots_a <- c(min(T_m, T_d),max(T_m, T_d))
	B_des_a_m <- splineDesign(knots_a, x=T_m, ord=order)
	B_des_a_d <- splineDesign(knots_a, x=T_d, ord=order)
}

if(mod[2]=='d')
{
	order <- 3

	if(knot_c < 3)
	{stop('The number of interior knots must be no less than 3.')}

}else
{
	order <- 1
}
if(mod[2]=='d')
{
	knots_c <- seq(from=min(T_m, T_d), to=max(T_m, T_d), length.out=knot_c)
	interval_c <- knots_c[2] - knots_c[1]
	knots_c <- c(c(min(T_m, T_d)-interval_c*2,min(T_m, T_d)-interval_c), knots_c)
	knots_c <- c(knots_c, c(max(T_m, T_d)+interval_c,max(T_m, T_d)+interval_c*2))
	B_des_c_m <- splineDesign(knots_c, x=T_m, ord=order)
	B_des_c_d <- splineDesign(knots_c, x=T_d, ord=order)
}else{
	knots_c <- c(min(T_m, T_d),max(T_m, T_d))
	B_des_c_m <- splineDesign(knots_c, x=T_m, ord=order)
	B_des_c_d <- splineDesign(knots_c, x=T_d, ord=order)
}

if(mod[3]=='d')
{
	order <- 3

	if(knot_e < 3)
	{stop('The number of interior knots must be no less than 3.')}

}else
{
	order <- 1
}
if(mod[3]=='d')
{
knots_e <- seq(from=min(T_m, T_d), to=max(T_m, T_d), length.out=knot_e)
interval_e <- knots_e[2] - knots_e[1]
knots_e <- c(c(min(T_m, T_d)-interval_e*2,min(T_m, T_d)-interval_e), knots_e)
knots_e <- c(knots_e, c(max(T_m, T_d)+interval_e,max(T_m, T_d)+interval_e*2))
B_des_e_m <- splineDesign(knots_e, x=T_m, ord=order)
B_des_e_d <- splineDesign(knots_e, x=T_d, ord=order)
}else{
knots_e <- c(min(T_m, T_d),max(T_m, T_d))
B_des_e_m <- splineDesign(knots_e, x=T_m, ord=order)
#knots_e <- c(min(T_d),max(T_d))
B_des_e_d <- splineDesign(knots_e, x=T_d, ord=order)

}

n_c <- ncol(B_des_c_m)
n_a <- ncol(B_des_a_m)
n_e <- ncol(B_des_e_m)

init_a <- rep(init[1],n_a)
init_c <- rep(init[2],n_c)
init_e <- rep(init[3],n_e)

up_a <- up_c <- up_e <- 10
lo_a <- lo_c <- -50
lo_e <- -15

if(mod[1]=='n')
{
up_a <- lo_a <- -50
init_a <- -50
}

if(mod[1]=='c')
{
up_a <- 20
lo_a <- -50
#init_a <- 1
}

if(mod[2]=='n')
{
up_c <- lo_c <- -50
init_c <- -50
}

if(mod[2]=='c')
{
up_c <- 20
lo_c <- -50
#init_c <- 1
}

if(mod[3]=='c')
{
up_e <- 20
lo_e <- -10
#init_e <- 1
}

result <- optim(c(init_a,init_c,init_e), loglik_AtCtEt_esp, gr_AtCtEt_esp, pheno_m = matrix(pheno_m), pheno_d = matrix(pheno_d), B_des_a_m = B_des_a_m, B_des_a_d = B_des_a_d, B_des_c_m = B_des_c_m, B_des_c_d = B_des_c_d, B_des_e_m = B_des_e_m, B_des_e_d = B_des_e_d,lower = c(rep(lo_a, n_a),rep(lo_c, n_c),rep(lo_e, n_e)), upper = c(rep(up_a, n_a),rep(up_c, n_c),rep(up_e, n_e)), method = "L-BFGS-B", hessian = TRUE, control=list(maxit = 3000))
#while(result$con>0)
#{
#	init <- runif(n_a+n_c+n_e,min=init_min,max=init_max)
#	if(mod[1]!='n')
#	{init_a <- init[1:n_a]}
#	if(mod[2]!='n')
#	{init_c <- init[(n_a+1):(n_a+n_c)]}
#	init_e <- init[(n_a+n_c+1):(n_a+n_c+n_e)]
#	result <- optim(c(init_a,init_c,init_e), loglik_AtCtEt_esp, gr_AtCtEt_esp, pheno_m = matrix(pheno_m), pheno_d = matrix(pheno_d), B_des_a_m = B_des_a_m, B_des_a_d = B_des_a_d, B_des_c_m = B_des_c_m, B_des_c_d = B_des_c_d, B_des_e_m = B_des_e_m, B_des_e_d = B_des_e_d,lower = c(rep(lo_a, n_a),rep(lo_c, n_c),rep(lo_e, n_e)), upper = c(rep(up_a, n_a),rep(up_c, n_c),rep(up_e, n_e)), method = "L-BFGS-B", hessian = TRUE, control=list(maxit = 3000))
#}

if(robust>0)
{
	for(i in 1:ceiling(robust))
	{
		init <- runif(n_a+n_c+n_e,min=init_min,max=init_max)
		if(mod[1]!='n')
		{init_a <- init[1:n_a]}
		if(mod[2]!='n')
		{init_c <- init[(n_a+1):(n_a+n_c)]}
		init_e <- init[(n_a+n_c+1):(n_a+n_c+n_e)]
		result_r <- optim(c(init_a,init_c,init_e), loglik_AtCtEt_esp, gr_AtCtEt_esp, pheno_m = matrix(pheno_m), pheno_d = matrix(pheno_d), B_des_a_m = B_des_a_m, B_des_a_d = B_des_a_d, B_des_c_m = B_des_c_m, B_des_c_d = B_des_c_d, B_des_e_m = B_des_e_m, B_des_e_d = B_des_e_d,lower = c(rep(lo_a, n_a),rep(lo_c, n_c),rep(lo_e, n_e)), upper = c(rep(up_a, n_a),rep(up_c, n_c),rep(up_e, n_e)), method = "L-BFGS-B", hessian = TRUE, control=list(maxit = 3000))
		if(result_r$value < result$value)
		{
			result <- result_r
		}
	}
}

res_a <- result$par[1:n_a]
res_c <- result$par[(1+n_a):(n_c+n_a)]
res_e <- result$par[(1+n_a+n_c):(n_e+n_c+n_a)]

hes <- .Call('hessian_AtCtEt_esp_c', res_a, res_c, res_e, matrix(pheno_m), matrix(pheno_d), B_des_a_m, B_des_a_d, B_des_c_m, B_des_c_d, B_des_e_m, B_des_e_d)

n_t <- n_a+n_c+n_e
hes_m <- matrix(0, n_t, n_t)
k <- 1
for(i in 1:n_t)
{
  for(j in i:n_t)
  {
    hes_m[i,j] <- hes[k]
    k <- k + 1
  }
}

hes_m_t <- t(hes_m)
diag(hes_m_t) <- 0
hes_m <- hes_m_t + hes_m

if(mod[1]=='n')
{res_a <- -Inf}
if(mod[2]=='n')
{res_c <- -Inf}

AtCtEt_model <- list(n_beta_a=n_a, n_beta_c=n_c, n_beta_e=n_e, beta_a=res_a, beta_c=res_c, beta_e=res_e, hessian_ap=result$hessian, hessian=hes_m, con=result$convergence, lik=result$value, knots_a =knots_a, knots_c = knots_c, knots_e = knots_e, min_t = min(T_m, T_d), max_t = max(T_m, T_d), boot = NULL )

#if(AtCtEt_model$con!=0)
#{
#	warning('The optimization algorithm might not converge. Try different initial values.')
#}

class(AtCtEt_model) <- 'AtCtEt_model'

if(boot==TRUE)
{
	boot_res <- AtCtEt_boot(res = AtCtEt_model, mod, data_m, data_d, knot_a, knot_c, knot_e, B=num_b,alpha=0.05,m=500)
	AtCtEt_model$boot <- boot_res
}

return(invisible(AtCtEt_model))

}
