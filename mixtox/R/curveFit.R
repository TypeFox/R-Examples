curveFit <- function(x, expr, eq , param, effv, fig = TRUE, ylimit, xlabel = "log[concentration, mol/L]", 
							ylabel = "Inhibition [%]", sigLev = 0.05, noec = TRUE, algo = "default"){
	# NLS curve fitting for monotonic and non-monotonic equations
	# x is a vector 
	# expr is a vector or matrix
	# for non-monotonic curve fitting, Brain_Consens, BCV, and Biphasic are highly recommended.
	## Jacobian matrix calculation
	jacobian <- function(eq, x, paraHat){
		n <- length(x)
		mpara <- length(paraHat)
		Alpha <- paraHat[1]
		Beta <- paraHat[2]
		if (mpara == 3) Gamma <- paraHat[3]
		if (mpara == 4) {Gamma <- paraHat[3]; Delta <- paraHat[4]}
		if (mpara == 5) {Gamma <- paraHat[3]; Delta <- paraHat[4]; Epsilon <- paraHat[5]}
		
		jac <- matrix(rep(0, n * mpara), n, mpara)
		
		jacFun <- switch(eq,
			Hill = c('-1 / (1 + (Alpha / x)^Beta)^2 * (Alpha / x)^Beta*Beta / Alpha', '-1 / (1 + (Alpha / x)^Beta)^2*(Alpha / x)^Beta * log(Alpha / x)'),
			Hill_two = c('x / (Alpha + x)', '-Beta * x / (Alpha + x)^2'),
			Hill_three = c('-Gamma / (1 + (Alpha / x)^Beta)^2 * (Alpha / x)^Beta * Beta / Alpha', '-Gamma / (1 + (Alpha / x)^Beta)^2 * (Alpha / x)^Beta * log(Alpha / x)', '1 / (1 + (Alpha / x)^Beta)'),
			Hill_four = c('-(Gamma - Delta) / (1 + (Alpha / x)^Beta)^2 * (Alpha / x)^Beta * Beta / Alpha', '-(Gamma - Delta) / (1 + (Alpha / x)^Beta)^2 * (Alpha / x)^Beta * log(Alpha / x)', 
							'1 /(1 + (Alpha / x)^Beta)', '1 - 1 / (1 + (Alpha / x)^Beta)'),
			Hill_six = c('-Gamma / (1 + (Alpha / x)^Beta)^2 * Gamma_one / (1 + (Alpha_one / x)^Beta_one) * (Alpha / x)^Beta * Beta / Alpha', '-Gamma / (1 + (Alpha / x)^Beta)^2 * Gamma_one / (1 + (Alpha_one / x)^Beta_one) * (Alpha / x)^Beta * log(Alpha / x)',
							'1 / (1 + (Alpha / x)^Beta) * Gamma_one / (1 + (Alpha_one / x)^Beta_one)', '-Gamma / (1 + (Alpha / x)^Beta) * Gamma_one / (1 + (Alpha_one / x)^Beta_one)^2 * (Alpha_one / x)^Beta_one * Beta_one / Alpha_one',
							'-Gamma / (1 + (Alpha / x)^Beta) * Gamma_one / (1 + (Alpha_one / x)^Beta_one)^2 * (Alpha_one / x)^Beta_one * log(Alpha_one / x)', 'Gamma / (1 + (Alpha / x)^Beta) / (1 + (Alpha_one / x)^Beta_one)'),
			Weibull = c('exp(Alpha + Beta * log(x) / log(10)) * exp(-exp(Alpha + Beta * log(x) / log(10)))', 
						'log(x) / log(10) * exp(Alpha + Beta * log(x) / log(10)) * exp(-exp(Alpha + Beta * log(x) / log(10)))'),
			
			Weibull_three = c('Gamma * exp(Alpha + Beta * log(x) / log(10)) * exp( -exp(Alpha + Beta * log(x) / log(10)))',
								'Gamma * log(x) / log(10) * exp(Alpha + Beta * log(x) / log(10)) * exp( -exp(Alpha + Beta * log(x) / log(10)))', '1 - exp( -exp(Alpha + Beta * log(x) / log(10)))'),
			
			Weibull_four = c(' -(Delta - Gamma) * exp(Alpha + Beta * log(x) / log(10)) * exp( -exp(Alpha + Beta * log(x) / log(10)))',
								' -(Delta - Gamma) * log(x) / log(10) * exp(Alpha + Beta * log(x) / log(10)) * exp( -exp(Alpha + Beta * log(x) / log(10)))', 
								'1 -exp( -exp(Alpha + Beta * log(x) / log(10)))', 'exp( -exp(Alpha + Beta * log(x) / log(10)))'),
			
			Logit = c('1 / (1 + exp(-Alpha - Beta * log(x) / log(10)))^2 * exp(-Alpha - Beta * log(x) / log(10))', 
						'1 / (1 + exp(-Alpha - Beta * log(x) / log(10)))^2 * log(x) / log(10) * exp(-Alpha - Beta * log(x) / log(10))'),
			
			Logit_three = c('Gamma / (1 + exp( -Alpha - Beta * log(x) / log(10)))^2 * exp( -Alpha - Beta * log(x) / log(10))',
								'Gamma / (1 + exp( -Alpha - Beta * log(x) / log(10)))^2 * log(x) / log(10) * exp( -Alpha - Beta * log(x) / log(10))',
								'1 / (1 + exp( -Alpha - Beta * log(x) / log(10)))'), 
 
			Logit_four = c('( -Delta + Gamma) / (1 + exp( -Alpha - Beta * log(x) / log(10)))^2 * exp( -Alpha - Beta * log(x) / log(10))',
							'( -Delta + Gamma) / (1 + exp( -Alpha - Beta * log(x) / log(10)))^2 * log(x) / log(10) * exp( -Alpha - Beta * log(x) / log(10))',
							'1 / (1 + exp( -Alpha - Beta * log(x) / log(10)))', '1 - 1 / (1 + exp( -Alpha - Beta * log(x) / log(10)))'), 
			BCW = c('exp(Alpha + Beta * (x^Gamma - 1) / Gamma) * exp(-exp(Alpha + Beta * (x^Gamma - 1) / Gamma))', 
					'(x^Gamma - 1) / Gamma * exp(Alpha + Beta * (x^Gamma - 1) / Gamma) * exp(-exp(Alpha + Beta * (x^Gamma - 1) / Gamma))', 
					'(Beta * x^Gamma * log(x) / Gamma - Beta * (x^Gamma - 1) / Gamma^2) * exp(Alpha + Beta * (x^Gamma - 1) / Gamma) * exp(-exp(Alpha + Beta * (x^Gamma - 1) / Gamma))'),
			BCL = c('1 /(1 + exp(-Alpha - Beta * (x^Gamma - 1) / Gamma))^2 * exp(-Alpha - Beta * (x^Gamma - 1) / Gamma)', 
					'1 / (1 + exp(-Alpha - Beta * (x^Gamma - 1) / Gamma))^2 * (x^Gamma - 1) / Gamma * exp(-Alpha - Beta * (x^Gamma - 1) / Gamma)', 
					'-1 / (1 + exp(-Alpha - Beta * (x^Gamma - 1) / Gamma)) ^ 2 * (-Beta * x^Gamma * log(x) / Gamma + Beta * (x^Gamma - 1) / Gamma^2) * exp(-Alpha - Beta * (x^Gamma - 1) / Gamma)'),
			GL = c('1 / ((1 + exp(-Alpha - Beta * log(x) / log(10)))^Gamma) * Gamma * exp(-Alpha - Beta * log(x) / log(10)) / (1 + exp(-Alpha - Beta * log(x) / log(10)))', 
					'1 / ((1 + exp(-Alpha - Beta * log(x) / log(10)))^Gamma) * Gamma * log(x) / log(10) * exp(-Alpha - Beta * log(x) / log(10)) / (1 + exp(-Alpha - Beta * log(x) / log(10)))', 
					'-1 / ((1 + exp(-Alpha - Beta * log(x) / log(10)))^Gamma) * log(1 + exp(-Alpha - Beta * log(x) / log(10)))'),
		
			Brain_Consens = c('-x / (1 + exp(Beta * Gamma) * x^Beta)', 
								'(1 + Alpha * x) / (1 + exp(Beta * Gamma) * x^Beta)^2 * (Gamma * exp(Beta * Gamma) * x^Beta + exp(Beta * Gamma) * x^Beta * log(x))', 
								'(1 + Alpha * x) / (1 + exp(Beta * Gamma) * x^Beta)^2 * Beta * exp(Beta * Gamma) * x^Beta'),
			
			BCV = c('-(1  +  Beta * x) / (1 + (1 + 2 * Beta * Gamma) * (x / Gamma)^Delta)', 
					'-Alpha * x / (1 + (1 + 2 * Beta * Gamma) * (x / Gamma)^Delta) - 2 * Alpha * (1 + Beta * x) / (1 + (1 + 2 * Beta * Gamma) * (x / Gamma)^Delta)^2 * Gamma * (x / Gamma)^Delta',
					'Alpha * (1 + Beta * x) / (1 + (1 + 2 * Beta * Gamma) * (x / Gamma)^Delta)^2 * (2 * Beta * (x / Gamma)^Delta - (1 + 2 * Beta * Gamma) * (x / Gamma)^Delta * Delta / Gamma)', 
					'Alpha * (1 + Beta * x) / (1 + (1 + 2 * Beta * Gamma) * (x / Gamma)^Delta)^2 * (1 + 2 * Beta * Gamma) * (x / Gamma)^Delta * log(x / Gamma)'),
			  
			Cedegreen = c('-exp(-1 / (x^Beta)) / (1 + exp(Gamma * (log(x) - log(Delta))))',
							'-Alpha / (x^Beta) * log(x) * exp(-1 / (x^Beta)) / (1 + exp(Gamma * (log(x) - log(Delta))))',
							'(1 + Alpha * exp(-1 / (x^Beta))) / (1 + exp(Gamma * (log(x) - log(Delta))))^2 * (log(x) - log(Delta)) * exp(Gamma * (log(x) - log(Delta)))', 
							'-(1 + Alpha * exp(-1 / (x^Beta))) / (1 + exp(Gamma * (log(x) - log(Delta))))^2 * Gamma / Delta * exp(Gamma * (log(x) - log(Delta)))'),

			Beckon = c('(1 - 1 / (1 + (Beta / x)^Gamma)) / (1 + (x / Delta)^Epsilon)', 
						'Alpha / (1 + (Beta / x)^Gamma)^2 * (Beta / x)^Gamma * Gamma / Beta / (1 + (x / Delta)^Epsilon)',
						'Alpha / (1 + (Beta / x)^Gamma)^2 * (Beta / x)^Gamma * log(Beta / x) / (1 + (x / Delta)^Epsilon)', 
						'(Alpha + 1 - Alpha / (1 + (Beta / x)^Gamma)) / (1 + (x / Delta)^Epsilon)^2 * (x / Delta)^Epsilon * Epsilon / Delta',
						'-(Alpha + 1 - Alpha / (1 + (Beta / x)^Gamma)) / (1 + (x / Delta)^Epsilon)^2 * (x / Delta)^Epsilon * log(x / Delta)'),
	
			Biphasic = c('1 - 1 / (1 + 10^((x-Beta) * Gamma)) - 1 / (1 + 10^((Delta - x) * Epsilon))', 
							'-Alpha / (1 + 10^((x-Beta) * Gamma))^2 * 10^((x - Beta) * Gamma) * Gamma * log(10)',
							'Alpha / (1 + 10^((x - Beta) * Gamma))^2 * 10^((x - Beta) * Gamma) * (x - Beta) * log(10)',
							'-(1 - Alpha) / (1 + 10^((Delta - x) * Epsilon))^2 * 10^((Delta - x) * Epsilon) * Epsilon * log(10)',
							'-(1 - Alpha) / (1 + 10^((Delta - x) * Epsilon))^2 * 10^((Delta - x) * Epsilon) * (Delta - x) * log(10)')
		)
		
		for (i in seq(mpara)) jac[, i] <- eval(parse(text = jacFun[i]))
		return(jac)
	}
	
	#############################################################
	ecxCI <- function(ciInfo, effv){
		# effect concentration and associated confidence intervals calculation
		oci.up <- spline(ciInfo[, 3], ciInfo[, 1], method = 'fmm', xout = effv)$y # effect concentration upper bound
		oci.low <- spline(ciInfo[, 4], ciInfo[, 1], method = 'fmm', xout = effv)$y # effect concentration lower bound
		oci.low[which(oci.low < 0)] = 0
		
		fci.up <- spline(ciInfo[, 5], ciInfo[, 1], method = 'fmm', xout = effv)$y # effect concentration upper bound
		fci.low <- spline(ciInfo[, 6], ciInfo[, 1], method = 'fmm', xout = effv)$y # effect concentration lower bound
		fci.low[which(fci.low < 0)] = 0
		
		ec.CI <- cbind(oci.low, oci.up, fci.low, fci.up)
		colnames(ec.CI) <- c('OCI.low', 'OCI.up', 'FCI.low', 'FCI.up')
		return(ec.CI)
	}
	
	#############################################################
	## confidence intervals for effect
	effvCI <- function(ciInfo, effv, ecx){
		# confidence interval for effect based on spline interpolation
		eoci.low <- spline(ciInfo[, 1], ciInfo[, 3], method = 'fmm', xout = ecx)$y
		eoci.up <- spline(ciInfo[, 1], ciInfo[, 4], method = 'fmm', xout = ecx)$y
		eoci.low[which(eoci.low < 0)] = 0
		
		efci.up <- spline(ciInfo[, 1], ciInfo[, 6], method = 'fmm', xout = ecx)$y # effect concentration upper bound
		efci.low <- spline(ciInfo[, 1], ciInfo[, 5], method = 'fmm', xout = ecx)$y # effect concentration lower bound
		efci.low[which(efci.low < 0)] = 0
		
		effv.CI <- cbind(eoci.low, eoci.up, efci.low, efci.up)
		colnames(effv.CI) <- c('eOCI.low', 'eOCI.up', 'eFCI.low', 'eFCI.up')
		return(effv.CI)
	}

	#############################################################
	figPlot <- function(crcInfo, ylimit, xlabel, ylabel){
		# plot the concentration-response curves
		#tiff(file = paste(root_name, "_04-12.tiff", sep = ""), res = 100)
		size <- dim(crcInfo)
		x <- crcInfo[, 1]
		yhat <- crcInfo[, 2]
		expr <- crcInfo[, 3 : (size[2] - 4)]
		if(is.vector(expr)) expr <- as.matrix(expr)
		oci <- crcInfo[, (size[2] - 3) : (size[2] - 2)]
		fci <- crcInfo[,  (size[2] - 1) : size[2]]
		
		if(missing(ylimit)) ylimit <- c((min(expr) * 100 -20), (max(expr) * 100 + 20))
		par(mar=c(5,5,1,1))
		plot(rep(log10(x), ncol(expr)), expr * 100, ylim = ylimit, pch = 16, xlab = xlabel, ylab = ylabel, cex = 1.8, cex.lab = 1.8, cex.axis = 1.8)
		lines(log10(x), yhat * 100, col = 1, lwd = 1.9)
		lines(log10(x), oci[, 1] * 100, col = 'blue', lwd = 1.9)
		lines(log10(x), oci[, 2] * 100, col = 'blue', lwd = 1.9)
		lines(log10(x), fci[, 1] * 100, col = 'red', lwd = 1.9)
		lines(log10(x), fci[, 2] * 100, col = 'red', lwd = 1.9)
		
		#legend("topleft", inset = 0.01, root_name, box.col = 'white', cex = 1.9) 
		#dev.off()
	}
	#############################################################	
	## source('ECx.R')
	
	## checking experimental data (expr)
	## main
	
	if (missing(x) || missing(expr) || missing(eq) || missing(param)) stop('argument missing')
	
	n <- length(x) # the number of concentrations
	mode(param) <- "numeric"
	m <- length(param) # the number of parameters
	
	if (is.vector(expr)){
		if (n != length(expr)) stop("x and y should be in the same length")
		y <- expr
		expr <- as.matrix(expr)
		nrep <- 1		
	}else if (is.matrix(expr)){
		size <- dim(expr)
		nrep <- size[2]
		y <- rowMeans(expr)
		if(n != size[1]) stop("x and dim(y)[1] should be in the same length")
	}
	
	## deploying the equation
	# if(eq == 'Hill'){
		# if(m == 3) eq = "Hill_three" else if (m == 4) eq = "Hill_four"
	# }
	
	# nonmonotonic or monotonic
	if(eq == 'Brain_Consens' || eq == 'BCV' || eq == 'Cedergreen' || eq == 'Beckon' || eq == 'Biphasic' || eq == 'Hill_six') Hormesis <- TRUE else Hormesis <- FALSE
	
	# define equation expression
	fun <- switch(eq,
		# For Hill equation: Alpha = EC50; Beta = m(Hill coefficient); Gamma = Top; Delta = Bottom
		Hill = 'y ~ 1 / (1 + (Alpha / x)^Beta)',
		# Howard GJ, Webster TF. 2009. Generalized concentration addition: A method for examining mixtures containing partial agonists. J. Theor. Biol. 259:469~477
		# Hill function with slope parameter 1. Alpha is EC50 here.
		Hill_two = 'y ~ Beta * x / (Alpha + x)',
		Hill_three = 'y ~ Gamma /(1 + (Alpha / x)^Beta)',
		Hill_four = 'y ~ Delta + (Gamma - Delta) / (1 + (Alpha / x)^Beta)',
		Hill_six = 'y ~ (Gamma / (1 + (Alpha / x)^Beta)) * (Gamma_one / (1 + (Alpha_one / x)^Beta_one))',
		# Hill_nine = 'y ~ (Gamma / (1 + (Alpha / x)^Beta)) * (Gamma_one / (1 + (Alpha_one / x)^Beta_one)) * (Gamma_two / (1 + (Alpha_two / x)^Beta_two))',
		Weibull = 'y ~ 1 - exp(-exp(Alpha + Beta * log10(x)))',
		Weibull_three = 'y ~ Gamma * (1 - exp(-exp(Alpha + Beta * log10(x))))',
		Weibull_four = 'y ~ Gamma + (Delta - Gamma) * exp(-exp(Alpha + Beta * log10(x)))',		
		Logit = 'y ~ 1/(1 + exp((-Alpha) - Beta * log10(x)))',
		Logit_three = 'y ~ Gamma / (1 + exp((-Alpha) - Beta * log10(x)))',
		Logit_four = 'y ~ Delta + (Gamma - Delta) / (1 + exp((-Alpha) - Beta * log10(x)))',
		BCW = 'y ~ 1 - exp(-exp(Alpha + Beta * ((x^Gamma - 1) / Gamma)))',
		BCL = 'y ~ (1 + exp(-Alpha - Beta *((x^Gamma - 1) / Gamma)))^(-1)',
		GL = 'y ~ 1 / (1 + exp(-Alpha - Beta * log10(x)))^Gamma',
		# An equation to describe dose responses where there isstimulatin of growth at low doses. 1989. Weed Research.
		Brain_Consens = 'y ~ 1 - (1 + Alpha * x) / (1 + exp(Beta * Gamma) * x^Beta)',
		# Vanewijk, P. H. and Hoekstra, J.A. Calculation of the EC50 and its confidence interval when subtoxic stimulus is present. 1993, Ecotoxicol. Environ. Saf.
		BCV = 'y ~ 1 - Alpha * (1 + Beta * x) / (1 + (1 + 2 * Beta * Gamma) * (x / Gamma)^Delta)',
		# Cedergreen, N., Ritz, C., Streibig, J.C., 2005. Improved empirical models describing hormesis. Environ. Toxicol. Chem. 24, 3166~3172
		Cedergreen = 'y ~ 1 - (1 + Alpha * exp(-1 / (x^Beta))) / (1 + exp(Gamma * (log(x) - log(Delta))))',
		# Beckon, W. et.al. 2008. A general approach to modeling biphasic relationships. Environ. Sci. Technol. 42, 1308~1314.
		Beckon = 'y ~ (Alpha + (1 - (Alpha) / (1 + (Beta / x)^Gamma))) / (1 + (x / Delta)^Epsilon)',
		# Zhu X-W, et.al . 2013. Modeling non-monotonic dose-response relationships: Model evaluation and hormetic quantities exploration. Ecotoxicol. Environ. Saf. 89:130~136;
		Biphasic = 'y ~ Alpha - Alpha / (1 + 10^((x - Beta) * Gamma)) + (1 - Alpha) / (1 + 10^((Delta - x) * Epsilon))'
	)
	
	## checking nls2 package, use the nls2 or built-in nls for curve fitting
	#if(require(nls2)){
	dframe <- data.frame(x, y)
	
	if(requireNamespace("nls2", quietly = TRUE)){
		print("use the nls2 package")
		
		if(eq == "Weibull" || eq == "Logit" || eq == "Hill" || eq == "Hill_two"){
			fit <- nls2::nls2(fun, data = dframe, start = list(Alpha = param[1], Beta = param[2]), control = nls.control(maxiter = 1000), algorithm = algo)
		
		}else if(eq == "BCW" || eq == "BCL" || eq == "GL" || eq == 'Brain_Consens' || eq == "Hill_three" || eq == 'Weibull_three' || eq == 'Logit_three'){
			fit <- nls2::nls2(fun, data = dframe, start = list(Alpha = param[1], Beta = param[2], Gamma = param[3]), control = nls.control(maxiter = 1000), algorithm = algo)
		
		}else if(eq == 'BCV'|| eq == 'Cedergreen' || eq == "Hill_four" || eq == 'Weibull_four' || eq == 'Logit_four'){
			#fit <- nls2::nls2(fun, data = dframe, start = list(Alpha = param[1], Beta = param[2], Gamma = param[3], Delta = param[4]), control = nls.control(maxiter = 1000), algorithm = "grid-search")
			fit <- nls2::nls2(fun, data = dframe, start = list(Alpha = param[1], Beta = param[2], Gamma = param[3], Delta = param[4]), control = nls.control(maxiter = 1000), algorithm = algo)
		}else if(eq == 'Beckon' || eq == 'Biphasic'){
			fit <- nls2::nls2(fun, data = dframe, start = list(Alpha = param[1], Beta = param[2], Gamma = param[3], Delta = param[4], Epsilon = param[5]), control = nls.control(maxiter = 1000), algorithm = algo)
			#fit <- nls2::nls2(fun, data = dframe, start = list(Alpha = param[1], Beta = param[2], Gamma = param[3], Delta = param[4], Epsilon = param[5]))
		}else if(eq == 'Hill_six'){
			fit <- nls2::nls2(fun, data = dframe, start = list(Alpha = param[1], Beta = param[2], Gamma = param[3], Alpha_one = param[4], Beta_one = param[5], Gamma_one = param[6]), control = nls.control(maxiter = 1000), algorithm = algo)
		} else{
			stop('input the right equation name')
		}
		
		#detach(package: nls2)
	}else {
		stop('please install package nls2')
	}
	
	fitInfo <- summary(fit) # fitting information
	
	yhat <- predict(fit, x) # y prediction
	res <- y - yhat
	sst <- sum((y - mean(y))^2) # total sum of squares
	sse <- sum((res)^2) # sum of squared errors
	r2 <- 1 - sse / sst # coefficient of determination
	adjr2 <- 1 - sse * (n - 1) / (sst * (n - m)) # adjusted coefficient of determination
	rmse <- sqrt(sse / (n - m)) # root-mean-square error
	mae <- sum(abs(res)) / n # mean absolute error
	#Spiess A-N, Neumeyer N. 2010. An evaluation of R2 as an inadequate measure for nonlinear models in pharmacological and biochemical research: A Monte Carlo approach. BMC Pharmacol. 10: 11.
	lnL <- 0.5 * (-n * (log(2 * pi) + 1 - log(n) + log(sse)))
	aic <- 2 * m - 2 * lnL # Akaike information criterion 
	aicc <- aic + 2 * m * (m + 1) / (n - m - 1)
	bic <- m * log(n) - 2 * lnL # Bayesian information criterion
	sta <- t(c(r2, adjr2, mae, rmse, aic, aicc, bic))
	colnames(sta) <- c('r2', 'adjr2', 'MAE', 'RMSE', 'AIC', 'AICc', 'BIC')
	
	paramHat <- t(as.matrix(summary(fit)$parameters[, 1]))	
	jac <- jacobian(eq, x, paramHat) # jacobian matrix calculation
	probT <- qt(1 - sigLev / 2, n - m) # the student t distribution
	mse <- rmse^2  # squared residual standard error
	covPara <- mse * solve(t(jac) %*% jac)  # covariance matrix of the parameter estimates
	
	gap.OCI <- sqrt(mse + diag(jac %*% covPara %*% t(jac))) # observation based confidence intervals
	gap.FCI <- sqrt(diag(jac %*% covPara %*% t(jac))) # function based confidence intervals
	
	OCI.up <- yhat + probT * gap.OCI # OCI upper bound
	OCI.low <- yhat - probT * gap.OCI # OCI lower bound
	FCI.up <- yhat + probT * gap.FCI # FCI upper bound
	FCI.low <- yhat - probT * gap.FCI # FCI lower bound
	
	crcInfo <- cbind(x, yhat, expr, OCI.low, OCI.up, FCI.low, FCI.up)
	ciInfo <- cbind(x, yhat, OCI.low, OCI.up, FCI.low, FCI.up)
	
	# compute highest stimulation (minimum effect) of the J-shaped curve and associated concentration. 
	# Brain_Consens, BCV, Cedergreen, Beckon, Biphasic
	if(Hormesis == TRUE){
		if(eq == 'Brain_Consens') Alpha = paramHat[1]; Beta = paramHat[2]; Gamma = paramHat[3]
		if(eq == 'BCV' || eq == 'Cedergreen') Alpha <- paramHat[1]; Beta <- paramHat[2]; Gamma <- paramHat[3]; Delta <- paramHat[4]
		if(eq == 'Beckon' || eq == 'Biphasic') Alpha <- paramHat[1]; Beta <- paramHat[2]; Gamma <- paramHat[3]; Delta <- paramHat[4]; Epsilon <- paramHat[5]
		if(eq == 'Hill_six') Alpha <- paramHat[1]; Beta <- paramHat[2]; Gamma <- paramHat[3]; Alpha_one <- paramHat[4]; Beta_one <- paramHat[5]; Gamma_one <- paramHat[6]
		
		if (eq == 'Brain_Consens') f <- function(x) 1 - (1 + Alpha * x) / (1 + exp(Beta * Gamma) * x^Beta)
		if(eq == 'BCV') f <- function(x) 1 - Alpha * (1 + Beta * x) / (1 + (1 + 2 * Beta * Gamma) * (x / Gamma)^Delta)
		if(eq == 'Cedergreen') f <- function(x) 1 - (1 + Alpha * exp(-1 / (x^Beta))) / (1 + exp(Gamma * (log(x) - log(Delta))))
		if(eq == 'Beckon') f <- function(x) (Alpha + (1 - (Alpha) / (1 + (Beta / x)^Gamma))) / (1 + (x / Delta)^Epsilon)
		if(eq == 'Biphasic') f <- function(x) Alpha - Alpha / (1 + 10^((x - Beta) * Gamma)) + (1 - Alpha) / (1 + 10^((Delta - x) * Epsilon))
		if(eq == 'Hill_six') f <- function(x) (Gamma / (1 + (Alpha / x)^Beta)) * (Gamma_one / (1 + (Alpha_one / x)^Beta_one))
		
		# intervals for finding the minimum
		intv <- c(x[2], x[length(x) - 1])
		
		minxy <- tryCatch({
			minxy <- optimize(f, intv)
		}, warning = function(w){
			message("Input an optimal intv")
		}, finally = {
			minxy <- list(minimum = NULL, objective = NULL)
		})
		minx <- minxy$minimum
		miny <- minxy$objective
	}

	## confidence intervals for effect concentration and effect
	## checking argument
	if(!missing(effv) && Hormesis == FALSE) {
		## effect concentration and confidence intervals 
		ecx <- ECx(eq, paramHat, effv)
		ecx.ci <- ecxCI(ciInfo, effv)
		ecx.ci <- cbind(t(ecx), ecx.ci)
		
		## effect confidence intervals 
		effv.ci <- effvCI(ciInfo, effv, ecx)
		effv.vec <- t(t(effv))
		rownames(effv.vec) <- paste('E', effv * 100, sep = '')
		colnames(effv.vec) <- 'Effect'
		effv.ci <- cbind(effv.vec, effv.ci)
	}else{
		ecx.ci = NULL
		effv.ci = NULL
	}
	
	if(!missing(effv) && Hormesis == TRUE){
		effv <- sort(effv)
		effv_pos <- effv[effv > 0]
		ecx <- nmECx(eq, paramHat, effv, minx)
		if(length(effv_pos) > 0){
			ecx.ci <- ecxCI(ciInfo, effv_pos)
			rownames(ecx.ci) <- paste('EC', effv_pos * 100, sep = '')
		}
	}
	
	## non-observed effect concentration
	## least observed effect concentration
	## at least 3 repetition
	if(nrep >= 3 && noec == TRUE){
		## source('NOEC.R')
		noecInfo <- NOEC(x, expr, sigLev)	
	}else {
		noecInfo <- NULL
	}
		
	## show concentration-response curve
	if(fig == TRUE){
		## source('figPlot.R')
		figPlot(crcInfo, ylimit, xlabel, ylabel)
	}
	
	if(Hormesis == FALSE){
		list(fitInfo = fitInfo, p = paramHat, res = res, sta = sta, crcInfo = crcInfo, eci = ecx.ci, effvci = effv.ci, noecInfo = noecInfo)
	}else{
		list(fitInfo = fitInfo, p = paramHat, res = res, sta = sta, minx = minx, miny = miny, crcInfo = crcInfo, ecx = ecx, eci = ecx.ci, noecInfo = noecInfo)
	}
}
