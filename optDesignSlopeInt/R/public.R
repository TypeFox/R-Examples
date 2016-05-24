#' Create an optimal design for measuring the slope divided by the intercept
#' 
#' @param n 				The number of experimental runs.
#' @param xmin 				The minimum value of the independent variable.
#' @param xmax 				The maximum value of the independent variable.
#' @param theta0 			The guess of the true value of the slope / intercept.
#' @param f_hetero			Specification of heteroskedasticity: the h(x) which relates the value of the 
#' 							independent variable to the variance in the response around the line at that place
#' 							or the proportional variance at that point. If \code{NULL}, homoskedasticity is
#' 							assumed (this is the default behavior).
#' @param MaxIter			For the heteroskedastic design, a Nelder-Mead search is used (via the function \code{fminbnd}). 
#' 							This is the \code{MaxIter} value for the search. Default is \code{6000}. Lower if \code{n} is high.
#' @param MaxFunEvals		For the heteroskedastic design, a Nelder-Mead search is used (via the function \code{fminbnd}). 
#' 							This is the \code{MaxFunEvals} value for the search. Default is \code{6000}. Lower if \code{n} is high.
#' @param TolFun			For the heteroskedastic design, a Nelder-Mead search is used (via the function \code{fminbnd}). 
#' 							This is the \code{TolFun} value for the search. Default is \code{1e-6}. Increase for faster execution.
#' @param NUM_RAND_STARTS	For the heteroskedastic design, a Nelder-Mead search is used (via the function \code{fminbnd}). 
#' 							The Nelder-Mead search must be given a starting location. Our implementation uses many
#' 							starting locations. This parameter controls the number of additional random starting
#' 							locations in the space \code{[xmin, xmax]}. Default is \code{50}.
#' 
#' @return 					An n-vector of x-values which specifies the optimal design
#' 
#' @author 					Adam Kapelner
#' @examples 
#' xmin = 5 / 15
#' xmax = 19 / 1
#' n = 10
#' theta0 = 0.053
#' opt_homo_design = oed_for_slope_over_intercept(n, xmin, xmax, theta0)
#' table(opt_homo_design)
#' @export
oed_for_slope_over_intercept = function(n, xmin, xmax, theta0, f_hetero = NULL, MaxIter = 6000, MaxFunEvals = 6000, TolFun = 1e-6, NUM_RAND_STARTS = 50){
	#throw errors here
	if (is.null(f_hetero)){
		oed_for_slope_over_intercept_homo(n, xmin, xmax, theta0)
	} else {
		oed_for_slope_over_intercept_hetero(n, xmin, xmax, theta0, f_hetero, MaxIter, MaxFunEvals, TolFun, NUM_RAND_STARTS)
	}
}

#' Plots a standard error estimate of thetahat (slope over intercept) over a range of possible theta0 values
#' in order to investigate robustness of the the initial theta0 guess.
#' 
#' @param n 			The number of experimental runs.
#' @param xmin 			The minimum value of the independent variable.
#' @param xmax 			The maximum value of the independent variable.
#' @param theta			The putative true value. This is used to see how much efficiency given up by designing it for \code{theta0}.
#' @param theta0 		The guess used to construct the experimental design. Specify only if you wish to see this 
#' 						value plotted. Default is \code{NULL}.
#' @param theta0_min 	Simulating over different guesses of theta0, this is the minimum guess.
#' @param theta0_max 	Simulating over different guesses of theta0, this is the maximum guess.
#' @param beta0	 	A guess to be used for the intercept. Defaults to \code{1}.
#' @param sigma 		A guess to be used for the homoskedastic variance of the measurement errors. If known accurately,
#' 						then the standard errors (i.e. the y-axis on the plot) will be accurate. Otherwise, the standard
#' 						errors are useful only when compared to each other in a relative sense. Defaults to \code{1}.
#' @param error_est 	The error metric for the estimates. The sample standard deviation (i.e. \code{sd}) 
#' 						is unstable at low sample sizes. The default is the 90 percentile minus the 10 percentile.
#' @param RES 			The number of points on the x-axis to simulate. Higher numbers will give smoother results. Default is \code{20}.
#' @param Nsim 			The number of models to be simulated for estimating the standard error at each value on the x-axis. Default is \code{1000}.
#' @param theta_logged	Should the values of theta be logged? Default is \code{TRUE}.
#' @param error_pct		Plot error as a percentage increase from minimum. Default is \code{TRUE}.
#' @param plot_rhos		Plot an additional graph of rho by theta0. Default is \code{FALSE}.
#' @param ...			Additional arguments passed to the \code{plot} function.
#' 
#' @return 				A list with original parameters as well as data from the simulation
#' 
#' @author 				Adam Kapelner
#' @export
err_vs_theta0_plot_for_homo_design = function(n, xmin, xmax, theta, theta0_min, theta0_max, 
		theta0 = NULL, beta0 = 1, sigma = 1, RES = 500, Nsim = 5000, 
		error_est = function(est){quantile(est, 0.99) - quantile(est, 0.01)}, #99%ile - 1%ile
		theta_logged = TRUE, error_pct = TRUE, plot_rhos = FALSE, ...){
	if (theta_logged){
		theta0s = seq(log(theta0_min), log(theta0_max), length.out = RES) / log(10)
	} else {
		theta0s = seq(theta0_min, theta0_max, length.out = RES)
	}
	

	homo_designs = matrix(NA, nrow = 0, ncol = n)
	for (i in 1 : RES){
		homo_design = oed_for_slope_over_intercept_homo(n, xmin, xmax, ifelse(theta_logged, 10^(theta0s[i]), theta0s[i]))
		homo_designs = rbind(homo_designs, homo_design)
	}
	
	unique_homo_designs = unique(homo_designs)
	unique_homo_designs
	unique_rhos = apply(unique_homo_designs, 1, function(row){sum(row == xmin) / n})
	
	unique_error_ests = array(NA, nrow(unique_homo_designs))
	for (d in 1 : nrow(unique_homo_designs)){
		ests_opt = array(NA, Nsim)
		xs_opt = unique_homo_designs[d, ]
		for (nsim in 1 : Nsim){						
			ys_opt = beta0 + theta  * beta0 * xs_opt + rnorm(n, 0, sigma)			
			mod_opt = lm(ys_opt ~ xs_opt)			
			ests_opt[nsim] = coef(mod_opt)[2] / coef(mod_opt)[1]
		}
		unique_error_ests[d] = error_est(ests_opt)
	}
	
	#now we need to realign the uniques
	error_ests = array(NA, RES)
	rhos = array(NA, RES)
	for (i in 1 : RES){
		homo_design = homo_designs[i, ]
		#get index in the unique set
		ind = apply(unique_homo_designs, 1, function(row){all(row == homo_design)}) 
		error_ests[i] = unique_error_ests[ind]
		rhos[i] = unique_rhos[ind]
	}
	
	if (plot_rhos){
		par(mfrow = c(1, 2))
	}
	
	if (error_pct){
		error_ests_pct = (error_ests - min(error_ests)) / min(error_ests) * 100
		plot(theta0s, error_ests_pct, 
				type = "l", 
				ylab = "% increase in optimal error of theta-hat", 
				xlab = ifelse(theta_logged, "log_10(theta0) hypothesis", "theta0 hypothesis"),
				...)		
	} else {
		plot(theta0s, error_ests, 
				type = "l", 
				ylab = "error of theta-hat", 
				xlab = ifelse(theta_logged, "log_10(theta0) hypothesis", "theta0 hypothesis"),
				...)		
	}
		
		

	
	if (!is.null(theta0)){
		abline(v = ifelse(theta_logged, log(theta0) / log(10), theta0), col = "red")
	}
	abline(v = ifelse(theta_logged, log(theta) / log(10), theta), col = "green")
	
	if (plot_rhos){
		plot(theta0s, rhos, 
				type = "l", 
				ylab = "rho", 
				xlab = ifelse(theta_logged, "log_10(theta0) hypothesis", "theta0 hypothesis"),
				...)
		
		if (!is.null(theta0)){
			abline(v = ifelse(theta_logged, log(theta0) / log(10), theta0), col = "red")
		}
		abline(v = ifelse(theta_logged, log(theta) / log(10), theta), col = "green")
	}	
	
	#return the simulated results only if user cares
	invisible(list(
		n = n,
		xmin = xmin,
		xmax = xmax,
		theta = theta, 
		theta_logged = theta_logged,
		theta_vs_design = cbind(ifelse(rep(theta_logged, RES), 10^(theta0s), theta0s), homo_designs),
		theta_vs_error_ests = cbind(ifelse(rep(theta_logged, RES), 10^(theta0s), theta0s), unique_error_ests)
	))		
}





#' A visualiation for comparing slope-divided-by-intercept estimates
#' for a number of designs
#' 
#' @param xmin 						The minimum value of the independent variable.
#' @param xmax 						The maximum value of the independent variable.
#' @param designs 					A d x n matrix where each of the d rows is a design (the x values
#' 									used to run the experiment).
#' @param gen_resp 					A model for the response which takes the design as its parameter.
#' @param Nsim 						The number of estimates per design. Default is \code{1000}.
#' @param l_quantile_display 		The lowest quantile of the simulation estimates displayed. Default is \code{0.025}.
#' @param u_quantile_display 		The highest quantile of the simulation estimates displayed. Default is \code{0.975}.
#' @param error_est 				The error metric for the estimates. The sample standard deviation (i.e. \code{sd}) 
#' 									is unstable at low sample sizes. The default is the 90 percentile minus the 10 percentile.
#' @param num_digits_round			The number of digits to round the error results. Default is 2.
#' @param draw_theta_at 			If the user wishes to draw a horizontal line marking theta (to checked biasedness)
#' 									it is specified here. The default is \code{NULL} with no line being drawn.
#' @param ...						Additional arguments passed to the \code{boxplot} function.
#' 
#' @return 							A list with the simulated estimates and error estimates for each design.
#' 
#' @author 							Adam Kapelner
#' @export
design_bakeoff = function(xmin, xmax, designs, 
		gen_resp = function(xs){1 + 2 * xs + rnorm(length(xs), 0, 1)}, 
		Nsim = 1000, l_quantile_display = 0.025, u_quantile_display = 0.975,
		error_est = function(est){quantile(est, 0.99) - quantile(est, 0.01)}, #99%ile - 1%ile i.e the interdecile range
		num_digits_round = 3,
		draw_theta_at = NULL, ...){
	num_designs = nrow(designs)
	n = ncol(designs)
	
	ests = matrix(NA, nrow = Nsim, ncol = num_designs)
	
	for (nsim in 1 : Nsim){
		for (d in 1 : num_designs){
			xs = designs[d, ]
			ys = gen_resp(xs)
			mod = lm(ys ~ xs)
			ests[nsim, d] = coef(mod)[2] / coef(mod)[1]
		}
	}
	error_ests = apply(ests, 2, error_est)
	
	l = list()
	for (d in 1 : num_designs){
		l[[d]] = ests[, d]
	}
	boxplot(l, ylim = quantile(c(ests), c(l_quantile_display, u_quantile_display)), 
			main = "Bakeoff: Error Estimates for Many Designs",  
			ylab = "theta-hat", 
			xlab = paste("Error Results:", paste(round(error_ests, num_digits_round), collapse = ", ")), 
			names = 1 : num_designs,
			...)
	if (!is.null(draw_theta_at)){
		abline(h = draw_theta_at, col = "blue")
	}
	invisible(list(ests = l, error_ests = error_ests))	
}

#' Report the results of the experiment as well as confidence intervals.
#' 
#' @param xs 		The design
#' @param ys 		The measurements of the response
#' @param alpha		\code{1 - alpha} is the confidence of the computed intervals. Default is \code{0.05}.
#' @param B			For the confidence interval methods with an embedded bootstrap (or resampling), the number
#' 					of resamples (defaults to \code{1000}).
#' 
#' @return			A list object containing the estimate as well as confidence intervals and parameters.
#' 
#' @author			Adam Kapelner
#' @export
experimental_results = function(xs, ys, alpha = 0.05, B = 1000){
	mod = lm(ys ~ xs)
	b0 = coef(mod)[1]
	b1 = coef(mod)[2]
	thetahat = b1 / b0 
	es = mod$residuals
	s_e = summary(mod)$sigma
	n = length(xs)
	
	X = cbind(rep(1, n), xs)
	XtXinv = solve(t(X) %*% X)
	
	cis = data.frame(matrix(NA, nrow = 0, ncol = 2))
	colnames(cis) = c("lower", "upper")
	#now generate them outside
	cis = rbind(
		normal_approx_homo_ci(alpha, b0, b1, XtXinv, s_e),
		normal_approx_hetero_ci(alpha, b0, b1, X, XtXinv, n, xs, es),
		bayesian_weighting_ci(alpha, B, n, xs, ys),
		param_homo_ci(alpha, B, b0, b1, s_e, n, xs),
		nonparam_homo_ci(alpha, B, xs, ys, es),
		nonparam_hetero_ci(alpha, B, n, xs, ys, es)
	)
	rownames(cis) = c(
		"Normal Approx Homo (poor coverage)",
		"Normal Approx Hetero (poor coverage)",
		"Bayesian Weighting",
		"Parametric Homo",
		"Nonparametric Homo",
		"Nonparametric Hetero"
	)	
	list(
		thetahat = thetahat,
		cis = cis,
		n = n,
		alpha = alpha,
		B = B
	)
}