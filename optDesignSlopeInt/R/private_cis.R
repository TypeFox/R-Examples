normal_approx_homo_ci = function(alpha, b0, b1, XtXinv, s_e){
	z = qnorm(1 - alpha / 2)
	est = b1 / b0
	del_g_beta_hat = t(as.matrix(c(-est, 1)))
	moe = z * s_e^2 / b0^2 * as.numeric(del_g_beta_hat %*% XtXinv %*% t(del_g_beta_hat))
	c(est - moe, est + moe)	
}

normal_approx_hetero_ci = function(alpha, b0, b1, X, XtXinv, n, xs, es){
	z = qnorm(1 - alpha / 2)
	est = b1 / b0
	del_g_beta_hat = t(as.matrix(c(-est, 1)))
	
	omega_diag = array(NA, n)
	list_es_to_xs = match_uniques_to_ests(xs, es)
	for (i in 1 : n){
		es_x = list_es_to_xs[[as.character(xs[i])]]
		omega_diag[i] = mean(es_x^2)
	}	
	Omega = diag(omega_diag)
	#note: we don't multiply by s_e^2 since it is baked into the Omega matrix
	moe = z / b0^2 * as.numeric(del_g_beta_hat %*% XtXinv %*% t(X) %*% Omega %*% X %*% XtXinv %*% t(del_g_beta_hat))
	c(est - moe, est + moe)
}

ratio_density_ci = function(alpha, b0, b1, XtXinv, s_e){
	c(0,0)
}

bayesian_weighting_ci = function(alpha, B, n, xs, ys){
	ws = rdirichlet(B, rep(1, n)) #stick breaking in n dimensions
	thetahat_bs = array(NA, B)
	for (b in 1 : B){
		mod = lm(ys ~ xs, weights = ws[b, ])
		b0_b = coef(mod)[1]
		b1_b = coef(mod)[2]
		thetahat_bs[b] = b1_b / b0_b 
	}
	confidence_endpoints(B, thetahat_bs, alpha)	
}

param_homo_ci = function(alpha, B, b0, b1, s_e, n, xs){
	thetahat_bs = array(NA, B)
	for (b in 1 : B){
		ys_b = b0 + b1 * xs + rnorm(n = n, mean = 0, sd = s_e)
		mod = lm(ys_b ~ xs)
		b0_b = coef(mod)[1]
		b1_b = coef(mod)[2]
		thetahat_bs[b] = b1_b / b0_b 
	}
	confidence_endpoints(B, thetahat_bs, alpha)	
}

nonparam_homo_ci = function(alpha, B, xs, ys, es){
	thetahat_bs = array(NA, B)
	for (b in 1 : B){
		ys_b = ys + sample(es, replace = T)
		mod = lm(ys_b ~ xs)
		b0_b = coef(mod)[1]
		b1_b = coef(mod)[2]
		thetahat_bs[b] = b1_b / b0_b 
	}
	confidence_endpoints(B, thetahat_bs, alpha)
}

nonparam_hetero_ci = function(alpha, B, n, xs, ys, es){
	list_es_to_xs = match_uniques_to_ests(xs, es)
	thetahat_bs = array(NA, B)
	for (b in 1 : B){
		ys_b = array(NA, n)
		for (i in 1 : n){
			es_x = list_es_to_xs[[as.character(xs[i])]]
			e = sample(es_x, 1)
			ys_b[i] = ys[i] + e * sample(c(-1, 1), 1)
		}
		mod = lm(ys_b ~ xs)
		b0_b = coef(mod)[1]
		b1_b = coef(mod)[2]
		thetahat_bs[b] = b1_b / b0_b 
	}
	confidence_endpoints(B, thetahat_bs, alpha)	
}

match_uniques_to_ests = function(keys, ests){
	l = list()
	for (i in 1 : length(keys)){
		k = as.character(keys[i])
		if (is.null(l[[k]])){
			l[[k]] = c()
		}
		l[[k]] = c(l[[k]], ests[i])
	}
	l
}

confidence_endpoints = function(B, ests, alpha){
	ests = sort(ests)
	l = ests[round(B * alpha / 2)]	
	r = ests[round(B * (1 - alpha / 2))]
	c(l, r)
}