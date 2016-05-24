tuneFit <- function(conc, rspn, eq = 'Weibull', effv, highBar = 5000, bar = 1000, sav = FALSE){

# File with the first line as header. assay name in the first column, casrn the second, and compounds' name the third column.
# the following columns should be conc an rspn in a same number.
	#library(nls2)
	#source('fitKenel.R')
	#source('ECx.R')
	#source('nmECx.R')
	#load('staval.rda')
	##########################################
	fitKenel <- function(x, expr, eq , param, effv, algo = "default"){
		# NLS curve fitting for monotonic and non-monotonic equations
		# x is a vector 
		# expr is a vector or matrix
		# for non-monotonic curve fitting, Brain_Consens, BCV, and Biphasic are highly recommended.

		fun <- switch(eq,
			# For Hill equation: Alpha = EC50; Beta = m(Hill coefficient); Gamma = Top; Delta = Bottom 
			Hill = 'y ~ 1 / (1 + (Alpha / x)^Beta)',
			# Howard GJ, Webster TF. 2009. Generalized concentration addition: A method for examining mixtures containing partial agonists. J. Theor. Biol. 259:469-477 
			#Hill function with slope parameter 1. Alpha is EC50 here.
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
			# Cedergreen, N., Ritz, C., Streibig, J.C., 2005. Improved empirical models describing hormesis. Environ. Toxicol. Chem. 24, 3166-3172
			Cedergreen = 'y ~ 1 - (1 + Alpha * exp(-1 / (x^Beta))) / (1 + exp(Gamma * (log(x) - log(Delta))))',
			# Beckon, W. et.al. 2008. A general approach to modeling biphasic relationships. Environ. Sci. Technol. 42, 1308-1314.
			Beckon = 'y ~ (Alpha + (1 - (Alpha) / (1 + (Beta / x)^Gamma))) / (1 + (x / Delta)^Epsilon)',
			# Zhu X-W, et.al . 2013. Modeling non-monotonic dose-response relationships: Model evaluation and hormetic quantities exploration. Ecotoxicol. Environ. Saf. 89:130-136;
			Biphasic = 'y ~ Alpha - Alpha / (1 + 10^((x - Beta) * Gamma)) + (1 - Alpha) / (1 + 10^((Delta - x) * Epsilon))'
		)
		
		## checking nls2 package, use the nls2 or built-in nls for curve fitting
		#if(require(nls2)){
		if (missing(x) || missing(expr) || missing(eq) || missing(param)) stop('argument missing')
		n <- length(x) # the number of concentrations
		mode(param) <- "numeric"
		m <- length(param) # the number of parameters
		n.a. <- 1000001
		
		if (is.vector(expr)){
			if (n != length(expr)) stop("x and y should be in the same length")
			y <- expr
			expr <- as.matrix(expr)
		}else if (is.matrix(expr)){
			size <- dim(expr)
			y <- rowMeans(expr)
			if(n != size[1]) stop("x and dim(y)[1] should be in the same length")
		}
		# nonmonotonic or monotonic
		if(eq == 'Brain_Consens' || eq == 'BCV' || eq == 'Cedergreen' || eq == 'Biphasic' || eq == 'Hill_six') Hormesis <- TRUE else Hormesis <- FALSE

		dframe <- data.frame(x, y)
		fit <- tryCatch({
			if(requireNamespace("nls2", quietly = TRUE)){
				if(eq == "Weibull" || eq == "Logit" || eq == "Hill" || eq == "Hill_two"){
						fit <- nls2::nls2(fun, data = dframe, start = list(Alpha = param[1], Beta = param[2]), control = nls.control(maxiter = 1000), algorithm = algo)
				
				}else if(eq == "BCW" || eq == "BCL" || eq == "GL" || eq == 'Brain_Consens' || eq == "Hill_three" || eq == 'Weibull_three' || eq == 'Logit_three'){
						fit <- nls2::nls2(fun, data = dframe, start = list(Alpha = param[1], Beta = param[2], Gamma = param[3]), control = nls.control(maxiter = 1000), algorithm = algo)
				
				}else if(eq == 'BCV'|| eq == 'Cedergreen' || eq == "Hill_four" || eq == 'Weibull_four' || eq == 'Logit_four'){
						#fit <- nls2(fun, data = dframe, start = list(Alpha = param[1], Beta = param[2], Gamma = param[3], Delta = param[4]), control = nls.control(maxiter = 1000), algorithm = "grid-search")
						fit <- nls2::nls2(fun, data = dframe, start = list(Alpha = param[1], Beta = param[2], Gamma = param[3], Delta = param[4]), control = nls.control(maxiter = 1000), algorithm = algo)
				}else if(eq == 'Beckon' || eq == 'Biphasic'){
						fit <- nls2::nls2(fun, data = dframe, start = list(Alpha = param[1], Beta = param[2], Gamma = param[3], Delta = param[4], Epsilon = param[5]), control = nls.control(maxiter = 1000), algorithm = algo)
						#fit <- nls2::nls2(fun, data = dframe, start = list(Alpha = param[1], Beta = param[2], Gamma = param[3], Delta = param[4], Epsilon = param[5]))
				}else if(eq == 'Hill_six'){
						fit <- nls2::nls2(fun, data = dframe, start = list(Alpha = param[1], Beta = param[2], Gamma = param[3], Alpha_one = param[4], Beta_one = param[5], Gamma_one = param[6]), control = nls.control(maxiter = 1000), algorithm = algo)
				} else{
					stop('input right equaiton name')
				}
			} else{
				stop('input the right equation name')
			}		
					
			fitInfo <- summary(fit) # fitting information	
			yhat <- predict(fit, x) # y prediction
			sst <- sum((y - mean(y))^2) # total sum of squares
			sse <- sum((y - yhat)^2) # sum of squared errors
			r2 <- 1 - sse / sst # coefficient of determination
			adjr2 <- 1 - sse * (n - 1) / (sst * (n - m)) # adjusted coefficient of determination
			rmse <- sqrt(sse / (n - m)) # root-mean-square error
			mae <- sum(abs(y - yhat)) / n # mean absolute error
			#Spiess A-N, Neumeyer N. 2010. An evaluation of R2 as an inadequate measure for nonlinear models in pharmacological and biochemical research: A Monte Carlo approach. BMC Pharmacol. 10: 11.
			lnL <- 0.5 * (-n * (log(2 * pi) + 1 - log(n) + log(sse)))
			aic <- 2 * m - 2 * lnL # Akaike information criterion 
			aicc <- aic + 2 * m * (m + 1) / (n - m - 1)
			bic <- m * log(n) - 2 * lnL # Bayesian information criterion
			sta <- t(c(r2, adjr2, mae, rmse, aic, aicc, bic))
			colnames(sta) <- c('r2', 'adjr2', 'MAE', 'RMSE', 'AIC', 'AICc', 'BIC')
			paramHat <- t(as.matrix(summary(fit)$parameters[, 1]))	
			
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
			if(!missing(effv) && Hormesis == FALSE){
				## effect concentration and confidence intervals 
				ecx <- ECx(eq, paramHat, effv[1])
			}else{
				ecx = n.a.
			}
			
			if(!missing(effv) && Hormesis == TRUE){
				effv <- sort(effv)
				effv_pos <- effv[effv > 0]
				ecx <- nmECx(eq, paramHat, effv[1], minx)
			}
			
			message('done...')
			if(Hormesis == FALSE){
				list(p = paramHat, sta = sta, ecx = ecx)
			}else{
				list(p = paramHat, sta = sta, minx = minx, miny = miny,  ecx = ecx)
			}
			
		}, error = function(e){
			#message("try...")
			if(Hormesis == FALSE){
				list(p = rep(n.a., m), sta = rep(n.a., 7), ecx = n.a.)
			}else{
				list(p = rep(n.a., m), sta = rep(n.a., 7), minx = n.a., miny = n.a., ecx = n.a.)
			}
			
		}, finally = {
			#message('done...')
		})
	}	
	###################################################	
	staval <- staval
	if (is.vector(conc)) conc <- t(conc)
	if (is.vector(rspn)) rspn <- t(rspn)
	if(nrow(conc) != nrow(rspn) || ncol(conc) != ncol(rspn))
		stop('row(column) of conc and rspn mismatch')
	param_name <- c('Alpha', 'Beta', 'Gamma', 'Delta', 'Epsilon', 'Zeta')
	param <- switch(eq,
		Hill = staval$Hill,
		Hill_two = staval$Hill_two,
		Hill_three = staval$Hill_three,
		Hill_four = staval$Hill_four,
		Weibull = staval$Weibull,
		Weibull_three = staval$Weibull_three,
		Weibull_four = staval$Weibull_four,
		Logit = staval$Logit,
		Logit_three = staval$Logit_three,
		Logit_four = staval$Logit_four,
		BCW = staval$BCW,
		BCL = staval$BCL,
		GL = staval$GL,
		Brain_Consens = staval$Brain_Consens,
		BCV = staval$BCV,
		Biphasic = staval$Biphasic,
		Hill_six = staval$Hill_six
	)
		
	if (nrow(param) > highBar) param <- param[sample(nrow(param), bar), ]
	m <- ncol(param)
	n <- nrow(param)
	#n <- 1
	curve_num <- nrow(conc)
	if(eq == 'Brain_Consens' || eq == 'BCV' || eq == 'Cedergreen' || eq == 'Beckon' || eq == 'Biphasic' || eq == 'Hill_six') Hormesis <- TRUE else Hormesis <- FALSE

	for(i in seq(curve_num)){
		conc_i <- as.vector(conc[i, ])
		#if(mean(rspn) > 1) rspn_i <- as.vector(rspn[i, ]) / 100  else rspn_i <- as.vector(rspn[i, ])
		rspn_i <- as.vector(rspn[i, ])
		
		for(j in seq(n)){
			fit <- fitKenel(conc_i, rspn_i, eq = eq , param = param[j, ],  effv = effv, algo = "default")
			if(fit$p[1] != 1000001)	break
		}
		if(Hormesis == TRUE){
			fit_ith <- c(fit$p, fit$sta, fit$minx, fit$miny, fit$ecx)
		} else{
			fit_ith <- c(fit$p, fit$sta, fit$ecx)
		}
		if(i == 1) fit_sta <- t(fit_ith) else fit_sta <- rbind(fit_sta, fit_ith)
	}

	#if(isTRUE(rownames(conc))) rownames(fit_sta) <- rownames(conc) else rownames(fit_sta) <- paste('fit_', as.vector(seq(nrow(conc))), sep ='')
	
	if(Hormesis == TRUE){
		colnames(fit_sta) <- c(param_name[1 : m], c('r2', 'adjr2', 'MAE', 'RMSE', 'AIC', 'AICc', 'BIC','minx', 'miny', 'ecx'))
	} else{
		colnames(fit_sta) <- c(param_name[1 : m], c('r2', 'adjr2', 'MAE', 'RMSE', 'AIC', 'AICc', 'BIC', 'ecx'))
	}

	#rownames(fit_sta) <- rownames(conc)
	if(isTRUE(rownames(conc))) rownames(fit_sta) <- rownames(conc) else rownames(fit_sta) <- paste('fit_', as.vector(seq(nrow(conc))), sep ='')

	if(sav == TRUE) {
		svfile = paste(eq, "_", Sys.Date(), ".fit", sep = "")
		write.table(fit_sta, svfile, sep = "\t", quote = F)
	}
	list(sta = fit_sta)
}
