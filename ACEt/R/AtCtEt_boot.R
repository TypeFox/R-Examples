AtCtEt_boot <- function(res, model, data_m, data_d, knot_a, knot_c, knot_e, B=100,alpha=0.05,m=500)
{
spline.main <- sp.spline.estimator(data_m, data_d, model, knot_a, knot_c, knot_e, m)
spline.boots_a <- matrix(NA,m,B)
spline.boots_c <- matrix(NA,m,B)
spline.boots_e <- matrix(NA,m,B)
spline.boots_h <- matrix(NA,m,B)
b_c <- res$beta_c
b_a <- res$beta_a
b_e <- res$beta_e
order <- 3
if(res$n_beta_a==1)
{order <- 1}
B_des_a_m <- splineDesign(res$knots_a, x=data_m[,3], ord=order)	
B_des_a_d <- splineDesign(res$knots_a, x=data_d[,3], ord=order)
order <- 3
if(res$n_beta_c==1)
{order <- 1}
B_des_c_m <- splineDesign(res$knots_c, x=data_m[,3], ord=order)	
B_des_c_d <- splineDesign(res$knots_c, x=data_d[,3], ord=order)
order <- 3
if(res$n_beta_e==1)
{order <- 1}
B_des_e_m <- splineDesign(res$knots_e, x=data_m[,3], ord=order)	
B_des_e_d <- splineDesign(res$knots_e, x=data_d[,3], ord=order)
num_m <- nrow(data_m)
num_d <- nrow(data_d)

for(i in 1:B)
{
pheno_mr <- matrix(0, num_m,2)
for(j in 1:num_m)
{
	var_ct <- sum(B_des_c_m[j,]*b_c)
	var_at <- sum(B_des_a_m[j,]*b_a)
	var_et <- sum(B_des_e_m[j,]*b_e)
	#if(res$n_beta_a>1)
	#{var_at <- exp(var_at)}
	var_at <- exp(var_at)
	#if(res$n_beta_c>1)
	#{var_ct <- exp(var_ct)}
	var_ct <- exp(var_ct)
	#if(res$n_beta_e>1)
	#{var_et <- exp(var_et)}
	var_et <- exp(var_et)
	# var_ct <- var_c
	sigma <- matrix(c(var_at+var_ct+var_et,var_at+var_ct,var_at+var_ct,var_at+var_ct+var_et),2,2)

	pheno_mr[j,1:2] <- mvrnorm(1, rep(0,2), sigma)
}

pheno_dr <- matrix(0, num_d,2)	
		
for(j in 1:num_d)
{
	var_ct <- sum(B_des_c_d[j,]*b_c)
	var_at <- sum(B_des_a_d[j,]*b_a)
	var_et <- sum(B_des_e_d[j,]*b_e)
	#if(res$n_beta_a>1)
	#{var_at <- exp(var_at)}
	var_at <- exp(var_at)
	#if(res$n_beta_c>1)
	#{var_ct <- exp(var_ct)}
	var_ct <- exp(var_ct)
	#if(res$n_beta_e>1)
	#{var_et <- exp(var_et)}
	var_et <- exp(var_et)
	sigma <- matrix(c(var_at+var_ct+var_et,0.5*var_at+var_ct,0.5*var_at+var_ct,var_at+var_ct+var_et),2,2)

	pheno_dr[j,1:2] <- mvrnorm(1, rep(0,2), sigma)
}


spline.boots <- sp.spline.estimator(cbind(pheno_mr,data_m[,3]),cbind(pheno_dr, data_d[,3]),model=model, knot_a, knot_c, knot_e, m=m)
spline.boots_a[,i] <- spline.boots$est_a
spline.boots_c[,i] <- spline.boots$est_c
spline.boots_e[,i] <- spline.boots$est_e
spline.boots_h[,i] <- spline.boots$est_a/(spline.boots$est_a+spline.boots$est_c+spline.boots$est_e)
}
# Result has m rows and B columns
#cis.lower_a <- 2*spline.main$est_a - apply(spline.boots_a,1,quantile,probs=1-alpha/2)
#cis.upper_a <- 2*spline.main$est_a - apply(spline.boots_a,1,quantile,probs=alpha/2)
#cis.lower_c <- 2*spline.main$est_c - apply(spline.boots_c,1,quantile,probs=1-alpha/2)
#cis.upper_c <- 2*spline.main$est_c - apply(spline.boots_c,1,quantile,probs=alpha/2)
#cis.lower_e <- 2*spline.main$est_e - apply(spline.boots_e,1,quantile,probs=1-alpha/2)
#cis.upper_e <- 2*spline.main$est_e - apply(spline.boots_e,1,quantile,probs=alpha/2)

# percentile method
cis.upper_a <- apply(spline.boots_a,1,quantile,probs=1-alpha/2)
cis.lower_a <- apply(spline.boots_a,1,quantile,probs=alpha/2)
cis.upper_c <- apply(spline.boots_c,1,quantile,probs=1-alpha/2)
cis.lower_c <- apply(spline.boots_c,1,quantile,probs=alpha/2)
cis.upper_e <- apply(spline.boots_e,1,quantile,probs=1-alpha/2)
cis.lower_e <- apply(spline.boots_e,1,quantile,probs=alpha/2)
cis.upper_h <- apply(spline.boots_h,1,quantile,probs=1-alpha/2)
cis.lower_h <- apply(spline.boots_h,1,quantile,probs=alpha/2)

# return(list(lower.ci_a=cis.lower_a,upper.ci_a=cis.upper_a,lower.ci_c=cis.lower_c,upper.ci_c=cis.upper_c, lower.ci_e=cis.lower_e,upper.ci_e=cis.upper_e, x=seq(from=res$min_t, to=res$max_t, length.out=m),boots_a=spline.boots_a, boots_c=spline.boots_c, boots_e=spline.boots_e))
return(list(lower.ci_a=cis.lower_a,upper.ci_a=cis.upper_a,lower.ci_c=cis.lower_c,upper.ci_c=cis.upper_c, lower.ci_e=cis.lower_e,upper.ci_e=cis.upper_e, lower.ci_h=cis.lower_h, upper.ci_h=cis.upper_h, x=seq(from=res$min_t, to=res$max_t, length.out=m)))

}


sp.spline.estimator <- function(data_m, data_d, model, knot_a, knot_c, knot_e, m) {
# Fit spline to data, with cross-validation to pick lambda
fit <- AtCtEt(data_m, data_d, model, knot_a, knot_c, knot_e)
T_m <- rep(data_m[,3], each=2)
T_d <- rep(data_d[,3], each=2)
eval.grid <- seq(from=min(T_m, T_d), to=max(T_m, T_d), length.out=m)
order <- 3
if(fit$n_beta_a==1)
{order <- 1}
bb_a <- splineDesign(fit$knots_a, x = eval.grid, ord=order, outer.ok = TRUE)
order <- 3
if(fit$n_beta_c==1)
{order <- 1}
bb_c <- splineDesign(fit$knots_c, x = eval.grid, ord=order, outer.ok = TRUE)
order <- 3
if(fit$n_beta_e==1)
{order <- 1}
bb_e <- splineDesign(fit$knots_e, x = eval.grid, ord=order, outer.ok = TRUE)

est_a <- bb_a%*%fit$beta_a
#if(fit$n_beta_a>1)
#{est_a <- exp(est_a)}
est_a <- exp(est_a)

est_c <- bb_c%*%fit$beta_c
#if(fit$n_beta_c>1)
#{est_c <- exp(est_c)}
est_c <- exp(est_c)

est_e <- bb_e%*%fit$beta_e
#if(fit$n_beta_e>1)
#{est_e <- exp(est_e)}
est_e <- exp(est_e)

return(list(e = fit$beta_e, c = fit$beta_c, a = fit$beta_a, est_c = est_c,est_a = est_a, est_e = est_e))

}
