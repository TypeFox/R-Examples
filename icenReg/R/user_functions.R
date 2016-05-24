ic_sp <- function(formula, data, model = 'ph', weights = NULL, bs_samples = 0, useMCores = F, 
                  useGA = T, maxIter = 5000, baseUpdates = 5){
  if(missing(data)) data <- environment(formula)
	cl <- match.call()
	mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action", "offset"), names(mf), 0L)
	mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    x <- model.matrix(mt, mf, contrasts)
	if(is.matrix(x))	xNames <- colnames(x)
    else				xNames <- as.character(formula[[3]])
	if('(Intercept)' %in% colnames(x)){	
		ind = which(colnames(x) == '(Intercept)')
		x <- x[,-ind]
		xNames <- xNames[-ind]
	}
	
  useFullHess = FALSE  
    
    
	if(length(xNames) == 0 & bs_samples > 0){
		 cat('no covariates included, so bootstrapping is not useful. Setting bs_samples = 0')
		 bs_samples = 0
	}
	
  yMat <- as.matrix(y)[,1:2]
  if(is(y, 'Surv')){
    rightCens <- mf[,1][,3] == 0
  	yMat[rightCens,2] <- Inf
  	exact <- mf[,1][,3] == 1
  	yMat[exact, 2] = yMat[exact, 1]
  }
	storage.mode(yMat) <- 'double'
    
    if(sum(is.na(mf)) > 0)
    	stop("NA's not allowed. If this is supposed to be right censored (i.e. [4, NA] was supposed to be right censored at t = 4), replace NA with Inf")
        
    testMat <- cbind(x, 1)
    invertResult <- try(diag(solve(t(testMat) %*% testMat )), silent = TRUE)
    if(is(invertResult, 'try-error'))
	    stop('covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level')
	if(model == 'ph')	callText = 'ic_ph'
	else if(model == 'po')	callText = 'ic_po'
	else stop('invalid choice of model. Current optios are "ph" (cox ph) or "po" (proportional odds)')
	
	if(is.null(weights)) 				weights = rep(1, nrow(yMat))
	if(length(weights) != nrow(yMat))	stop('weights improper length')
	if(any(is.na(weights) > 0) )		stop('NAs not allowed in weights')
	if(any(weights < 0)	)				stop('negative weights not allowed')
	
	if(is.null(ncol(x)) ) recenterCovars = FALSE
	
  other_info <- list(useGA = useGA, maxIter = maxIter, 
                     baselineUpdates = baseUpdates, 
                     useFullHess = useFullHess, 
                     useEM = FALSE)  

  fitInfo <- fit_ICPH(yMat, x, callText, weights, other_info)
	dataEnv <- list()
	dataEnv[['x']] <- as.matrix(x, nrow = nrow(yMat))
	if(ncol(dataEnv$x) == 1) colnames(dataEnv[['x']]) <- as.character(formula[[3]])
	dataEnv[['y']] <- yMat
	seeds = as.integer( runif(bs_samples, 0, 2^31) )
	bsMat <- numeric()
	if(useMCores) `%mydo%` <- `%dopar%`
	else          `%mydo%` <- `%do%`
	i <- NULL  #this is only to trick R CMD check, 
	           #as it does not recognize the foreach syntax
	if(bs_samples > 0){
	    bsMat <- foreach(i = seeds, .combine = 'rbind') %mydo%{
	        set.seed(i)
	        sampDataEnv <- bs_sampleData(dataEnv, weights)
			    getBS_coef(sampDataEnv, callText = callText,
				           other_info = other_info)
	    	}
	 }
	    
	 if(bs_samples > 0){
	   	 names(fitInfo$coefficients) <- xNames
	   	 colnames(bsMat) <- xNames
	   	 incompleteIndicator <- is.na(bsMat[,1])
	   	 numNA <- sum(incompleteIndicator)
	   	 if(numNA > 0){
	    		if(numNA / length(incompleteIndicator) >= 0.1)
	    		cat('warning: ', numNA,
	    		    ' bootstrap samples (out of ', bs_samples, 
	    		    ') were dropped due to singular covariate matrix.',
              'Likely due to very sparse covariate. Be wary of these results.\n', sep = '')
	    		bsMat <- bsMat[!incompleteIndicator,]
	    	}
	      covar <- cov(bsMat)
        est_bias <- colMeans(bsMat) - fitInfo$coefficients 
        fitInfo$coef_bc <- fitInfo$coefficients - est_bias
    }else{ 
        bsMat <- NULL
        covar <- NULL
        coef_bc <- NULL
    }
    names(fitInfo$coefficients) <- xNames
    fitInfo$bsMat <- bsMat
    fitInfo$var <- covar
    fitInfo$call = cl
    fitInfo$formula = formula
    fitInfo$.dataEnv <- new.env()
    if(!missing(data)){ fitInfo$.dataEnv$data = data }
    fitInfo$par = 'semi-parametric'
    fitInfo$model = model
    fitInfo$reg_pars <- fitInfo$coefficients
    fitInfo$terms <- mt
    fitInfo$xlevels <- .getXlevels(mt, mf)
    if(fitInfo$iterations == maxIter){
      warning(paste0('Maximum iterations reached in ic_sp.', 
              '\n\nCommon cause of problem is when many observations are uncensored in Cox-PH model',
              '\nas this causes heavy numeric instability in gradient ascent step.', 
              '\n(note:proportional odds model is more numerically stable).',
              '\nICM step is still stable, so try increasing maxIter two fold and observe if', 
              '\ndifference in final llk is less than 0.1. If not, increase maxIter until this occurs.'))
    }
   return(fitInfo)
}


vcov.icenReg_fit <- function(object,...) object$var


getSCurves <- function(fit, newdata = NULL){
	if(inherits(fit, 'impute_par_icph'))	stop('getSCurves currently not supported for imputation model')
	if(inherits(fit, 'ic_par'))				stop('getSCurves does not support ic_par objects. Use getFitEsts() instead. See ?getFitEsts')
	etas <- get_etas(fit, newdata)
	grpNames <- names(etas)
	transFxn <- get_link_fun(fit)
	if(fit$par == 'semi-parametric' | fit$par == 'non-parametric'){
		x_l <- fit$T_bull_Intervals[1,]
		x_u <- fit$T_bull_Intervals[2,]
		x_l <- c(x_l[1], x_l)
		x_u <- c(x_l[1], x_u)
		Tbull_intervals <- cbind(x_l,  x_u)
		colnames(Tbull_intervals) <- c('lower', 'upper')
		s <- 1 - c(0, cumsum(fit$p_hat))
		ans <- list(Tbull_ints = Tbull_intervals, "S_curves" = list())
		
		for(i in 1:length(etas)){
			eta <- etas[i]
			ans[["S_curves"]][[grpNames[i] ]] <- transFxn(s, eta)
		}
		class(ans) <- 'sp_curves'
		return(ans)
	}
	else{
	  	stop('getSCurves only for semi-parametric model. Try getFitEsts')
	}
}


plot.icenReg_fit <- function(x, y, fun = 'surv', 
                             lgdLocation = 'topright', xlab = "time", ...){
	if(inherits(x, 'impute_par_icph'))	stop('plot currently not supported for imputation model')
  argList <- list(...)
  colors <- argList$col
	if(missing(y)) y <- argList$newdata	
	newdata <- y
  nRows <- 1
  if(!is.null(newdata)) nRows <- nrow(newdata)
  if(fun == 'surv'){ s_trans <- function(x){x}; yName = 'S(t)'}
  else if(fun == 'cdf'){ s_trans <- function(x){1-x}; yName = 'F(t)' }
  else stop('"fun" option not recognized. Choices are "surv" or "cdf"')

  addList <- list(xlab = xlab, ylab = yName)
  argList <- addListIfMissing(addList, argList)
  firstPlotList <- argList
  firstPlotList[['type']] <- 'n'
  firstPlotList[['x']] <- 1
  firstPlotList[['y']] <- 1
  
    
	if(x$par == 'semi-parametric' | x$par == 'non-parametric'){
		curveInfo <- getSCurves(x, y)
		allx <- c(curveInfo$Tbull_ints[,1], curveInfo$Tbull_ints[,2])
		dummyx <- range(allx, finite = TRUE)
		dummyy <- c(0,1)
	  firstPlotList[['xlim']] <- dummyx
	  firstPlotList[['ylim']] <- dummyy
    
		x_l <- curveInfo$Tbull_ints[,1]
		x_u <- curveInfo$Tbull_ints[,2]
		k <- length(x_l)
		ss <- curveInfo$S_curves

		do.call(plot, firstPlotList)
		
		if(is.null(colors))  colors <- 1:length(ss)
		if(length(colors) == 1) colors <- rep(colors, length(ss)) 
		for(i in 1:length(ss)){
			lines(x_l, s_trans(ss[[i]]), col = colors[i], type = 's')
			lines(x_u, s_trans(ss[[i]]), col = colors[i], type = 's')
			lines(c(x_l[k], x_u[k]), s_trans(c(ss[[i]][k], ss[[i]][k])), col = colors[i])
		}
		if(length(ss) > 1){
			grpNames <- names(ss)
			legend(lgdLocation, legend = grpNames, lwd = rep(1, length(grpNames) ), col = colors)
		}
	}
	else if(inherits(x, 'par_fit')){
    ranges <- matrix(nrow = nRows, ncol = 2)
		ranges[,1] <- getFitEsts(x, newdata = newdata, p = 0.05 )
    ranges[,2] <- getFitEsts(x, newdata = newdata, p = 0.95 )
    
    addList <- list(xlab = xlab, ylab = yName, 
                    xlim = range(as.numeric(ranges)), ylim = c(0,1))
    argList <- addListIfMissing(addList, argList)
    firstPlotList <- argList
    do.call(plot, firstPlotList)

    ranges[,1] <- getFitEsts(x, newdata = newdata, p = 0.005 )
		ranges[,2] <- getFitEsts(x, newdata = newdata, p = 0.995 )
		if(is.null(colors))  colors <- 1:nRows
		
		for(i in 1:nrow(ranges)){
			grid = ranges[i,1] + 0:100/100 * (ranges[i,2] - ranges[i,1])
			est.s <- 1 - getFitEsts(x, newdata = subsetData_ifNeeded(i, newdata), q = grid)
			lines(grid, s_trans(est.s), col = colors[i])
		}
		if(nrow(ranges) > 1){
			grpNames <- rownames(newdata)
			legend(lgdLocation, legend = grpNames, lwd = rep(1, length(grpNames) ), col = 1:ncol(ranges))
		}
	}
}


lines.icenReg_fit <- function(x, y, fun = 'surv', ...){
  argList <- list(...)
  colors <- argList$col
  if(missing(y)) y <- argList$newdata	
  newdata <- y
  nRows <- 1
  if(!is.null(newdata)) nRows <- nrow(newdata)
  if(fun == 'surv'){ s_trans <- function(x){x}; yName = 'S(t)'}
  else if(fun == 'cdf'){ s_trans <- function(x){1-x}; yName = 'F(t)' }
  else stop('"fun" option not recognized. Choices are "surv" or "cdf"')
  
#  addList <- list(xlab = xlab, ylab = yName)
#  argList <- addListIfMissing(addList, argList)

  if(x$par == 'semi-parametric' | x$par == 'non-parametric'){
    argList <- addIfMissing('s', 'type', argList)
    curveInfo <- getSCurves(x, y)
    allx <- c(curveInfo$Tbull_ints[,1], curveInfo$Tbull_ints[,2])
    dummyx <- range(allx, finite = TRUE)
    dummyy <- c(0,1)
    x_l <- curveInfo$Tbull_ints[,1]
    x_u <- curveInfo$Tbull_ints[,2]
    k <- length(x_l)
    ss <- curveInfo$S_curves
    if(is.null(colors))  colors <- 1:length(ss)
    if(length(colors) == 1) colors <- rep(colors, length(ss)) 
    for(i in 1:length(ss)){
      argList[['x']] <- x_l
      argList[['y']] <- s_trans(ss[[i]])
      argList[['col']] <- colors[i]
      do.call(lines, argList)     
      argList[['x']] <- x_u
      do.call(lines, argList)     
      argList[['x']] <- c(x_l[k], x_u[k])
      argList[['y']] <- s_trans(c(ss[[i]][k], ss[[i]][k]))
      do.call(lines, argList)
    }
  }
  else if(inherits(x, 'par_fit')){
    ranges <- matrix(nrow = nRows, ncol = 2)
    ranges[,1] <- getFitEsts(x, newdata = newdata, p = 0.05 )
    ranges[,2] <- getFitEsts(x, newdata = newdata, p = 0.95 )
    if(is.null(colors))  colors <- 1:nRows
    for(i in 1:nrow(ranges)){
      grid = ranges[i,1] + 0:100/100 * (ranges[i,2] - ranges[i,1])
      est.s <- 1 - getFitEsts(x, newdata = subsetData_ifNeeded(i, newdata), q = grid)
      argList[['x']] <- grid
      argList[['y']] <- s_trans(est.s)
      argList[['col']] <- colors[i]
      do.call(lines, argList)
    }
  }
}


OLD_lines.icenReg_fit <- function(x, y, fun = 'surv', lgdLocation = 'topright', xlab = "time",
                             colors = NULL, ...){
  if(missing(y)) y <- list(...)$newdata	
  newdata <- y
  nRows <- 1
  if(!is.null(newdata)) nRows <- nrow(newdata)
  if(fun == 'surv'){ s_trans <- function(x){x}; yName = 'S(t)'}
  else if(fun == 'cdf'){s_trans <- function(x){1-x}; yName = 'F(t)'}
  else stop('"fun" option not recognized. Choices are "surv" or "cdf"')
  
  if(x$par == 'semi-parametric'){
    curveInfo <- getSCurves(x, y)
    allx <- c(curveInfo$Tbull_ints[,1], curveInfo$Tbull_ints[,2])
    dummyx <- range(allx, finite = TRUE)
    dummyy <- c(0,1)
    
#    plot(dummyx, dummyy, xlab = xlab, ylab = yName, ..., type = 'n')
    x_l <- curveInfo$Tbull_ints[,1]
    x_u <- curveInfo$Tbull_ints[,2]
    k <- length(x_l)
    ss <- curveInfo$S_curves
    if(is.null(colors))  colors <- 1:length(ss)
    
    for(i in 1:length(ss)){
      lines(x_l, s_trans(ss[[i]]), col = colors[i], type = 's')
      lines(x_u, s_trans(ss[[i]]), col = colors[i], type = 's')
      lines(c(x_l[k], x_u[k]), s_trans(c(ss[[i]][k], ss[[i]][k])), col = colors[i])
    }
    if(length(ss) > 1){
      grpNames <- names(ss)
      legend(lgdLocation, legend = grpNames, lwd = rep(1, length(grpNames) ), col = colors)
    }
  }
  else if(inherits(x, 'par_fit')){
    ranges <- matrix(nrow = nRows, ncol = 2)
    ranges[,1] <- getFitEsts(x, newdata = newdata, p = 0.05 )
    ranges[,2] <- getFitEsts(x, newdata = newdata, p = 0.95 )
  #  plot(NA, xlim = range(as.numeric(ranges), finite = TRUE), ylim = c(0,1), xlab = xlab, ylab = yName)
    ranges[,1] <- getFitEsts(x, newdata = newdata, p = 0.005 )
    ranges[,2] <- getFitEsts(x, newdata = newdata, p = 0.995 )
    if(is.null(colors))  colors <- 1:nRows
    
    for(i in 1:nrow(ranges)){
      grid = ranges[i,1] + 0:100/100 * (ranges[i,2] - ranges[i,1])
      est.s <- 1 - getFitEsts(x, newdata = subsetData_ifNeeded(i, newdata), q = grid)
      lines(grid, s_trans(est.s), col = colors[i])
    }
    if(nrow(ranges) > 1){
      grpNames <- rownames(newdata)
      legend(lgdLocation, legend = grpNames, lwd = rep(1, length(grpNames) ), col = 1:ncol(ranges))
    }
  }
}



summary.icenReg_fit <- function(object,...)
	new('icenRegSummary', object)
summary.ic_npList <- function(object, ...)
  object

	
summaryOld.icenReg_fit <- function(object,...){
	sigfigs = 4
	fit <- object
	if(inherits(fit, 'impute_par_icph')) cat('\nMultiple Imputations Cox PH model for interval censored data\n')
	else{
		fullModelName <- if(fit$model == 'ph') "Cox PH" else "Proportional Odds"
		baseline = fit$par
		fitDisc <- paste0("\nModel: ", fullModelName, "\nBaseline: ",  baseline, "\n")
		cat(fitDisc)
	}
	if(!is.null(fit$var)){
		colNames <- c('Estimate', 'Exp(Est)', 'Std. Error', 'z-value', 'p')
		coefs <- fit$coefficients
		output <- matrix(nrow = length(coefs), ncol = length(colNames))
		se <- sqrt(diag(fit$var))
		for(i in seq_along(coefs)){
			output[i, 1] <- coefs[i]
			output[i, 2] <- exp(coefs[i])
			output[i, 3] <- se[i]
			output[i, 4] <- coefs[i]/se[i]
			output[i, 5] <- 2*(1 - pnorm(abs(output[i,4])))
		}
		colnames(output) <- colNames
		rownames(output) <- names(coefs)
		output <- signif(output, sigfigs)
		cat("Call = \n")
		print(fit$call)
		cat('\n')
		print(output)
		if(inherits(fit, 'ic_ph') | inherits(fit, 'ic_po')){
			cat('\nfinal llk = ', fit$llk, '\n')
			cat('Iterations = ', fit$iterations, '\n')
			cat('Bootstrap samples = ', nrow(fit$bsMat), '\n')
			if(nrow(fit$bsMat) < 100)
				cat('CAUTION: recommend more bootstrap samples for inference!\n')
		}
		
		if(inherits(fit, 'impute_par_icph')){
			cat('\nnumber of imputations = ', nrow(fit$imp_coef), '\n')
		}
		if(inherits(fit, 'par_fit')){
			cat('\nfinal llk = ', fit$llk, '\n')
			cat('Iterations = ', fit$iterations,'\n')
		}
	}
	else{
		colNames <- c('Estimate', 'Exp(Est)')
		coefs <- fit$coefficients
		output <- matrix(nrow = length(coefs), ncol = length(colNames))
		se <- sqrt(diag(fit$var))
		for(i in seq_along(coefs)){
			output[i, 1] <- coefs[i]
			output[i, 2] <- exp(coefs[i])
			}
		colnames(output) <- colNames
		rownames(output) <- names(coefs)
		output <- signif(output, sigfigs)
		cat("Call = \n")
		print(fit$call)
		print(output)
		cat('final llk = ', fit$llk, '\n')
		cat('Iterations = ', fit$iterations, '\n')	
		cat('Standard Errors not available. To get standard errors, rerun ic_ph with "bs_samples" > 0 (suggested at least 1000)')
	}
}






# Multiple Imputations Regression Models

impute_ic_ph <- function(formula, data, imps = 100, eta = 10^-10, rightCenVal = 10000, seed = NULL, useMCores = FALSE){
	cl <- match.call()
	mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    
    if (is.empty.model(mt)){
    	stop('no covariates included. Try using computeMLE in MLECens package')
    }
     x <- model.matrix(mt, mf, contrasts)
	if('(Intercept)' %in% colnames(x)){	
		ind = which(colnames(x) == '(Intercept)')
		x <- x[,-ind]
	}
	if(!is(y, 'Surv'))	stop('response must be a "Surv" item with type = "interval2"')
	
    yMat <- as.matrix(y)[,1:2]
    rightCens <- mf[,1][,3] == 0
	yMat[rightCens,2] <- Inf
	exact <- mf[,1][,3] == 1
	yMat[exact, 2] = yMat[exact, 1]
	storage.mode(yMat) <- 'double'
		
	param_y <- mf[[1]]
	l_e0 <- param_y[,1] < eta
	r_e0 <- param_y[,2] < eta
	param_y[l_e0,1] <- eta
	param_y[r_e0,2] <- eta
	param_data  <- data.frame(mf[,-1])
	param_formula <- formula
	param_formula[[2]] <- as.name('param_y')
	paramFit <- fullParamFit_exp(param_formula, mf, param_y, rightCenVal)
	
	
	fits_from_imputes <- list()
	obs_l <- yMat[,1]
	obs_r <- yMat[,2]
	use_n <- length(obs_l)
	one_vec <- rep(1, use_n)
	
	coxphFormula <- formula
	coxphFormula[[2]][[2]] <- as.name('sim_y')
	coxphFormula[[2]][[3]] <- as.name('one_vec')
	coxphFormula[[2]][[4]] <- 'right'
	
	mf_forImputes <- mf
	mf_forImputes <- data.frame(mf_forImputes)
	mf_forImputes[['one_vec']] <- one_vec

	if(!is.numeric(seed))
		seed <- round( runif(1, max = 10000000) )

	if(!useMCores){
		set.seed(seed)
		for(i in 1:imps){
			sampPars <- simPars_fromFit(paramFit)
			sim_y <- imputeCensoredData_exp(obs_l, obs_r, sampPars, maxVal = rightCenVal)
			mf_forImputes[['sim_y']] <- sim_y		
			fits_from_imputes[[i]] <- coxph(coxphFormula, data = mf_forImputes)
		}
	}
	
	if(useMCores){
		fits_from_imputes <- foreach(i = 1:imps) %dopar%{
			set.seed(i + seed)
			sampPars <- simPars_fromFit(paramFit)
			sim_y <- imputeCensoredData_exp(obs_l, obs_r, sampPars, maxVal = rightCenVal)
			mf_forImputes[['sim_y']] <- sim_y
			return(coxph(coxphFormula, data = mf_forImputes) ) 
		}
	}
	
	k <- length(paramFit$coefficients)
	coef <- numeric()
	ave_var <- matrix(0, nrow = k-1, ncol = k-1)
	for(i in 1:imps){
		thisFit <- fits_from_imputes[[i]]
		coef <- rbind(coef, thisFit$coefficients)
		ave_var <- ave_var + thisFit$var/imps	
	}
	var <- ave_var + cov(coef)
	mean_coef <- colMeans(coef)
	fit <- list(coefficients = mean_coef, var = var, imp_coef = coef, average_imp_covariace = ave_var, call = cl, formula = formula)
    class(fit) <- c('icenReg_fit', 'impute_par_icph')
	return(fit)
}

simIC_weib <- function(n = 100, b1 = 0.5, b2 = -0.5, model = 'ph', 
					   shape = 2, scale = 2, 
					   inspections = 2, inspectLength = 2.5,
					   rndDigits = NULL, prob_cen = 1){
	rawQ <- runif(n)
    x1 <- runif(n, -1, 1)
    x2 <- 1 - 2 * rbinom(n, 1, 0.5)
    nu <- exp(x1 * b1 + x2 * b2)
    
    if(model == 'ph')		adjFun <- function(x, nu) {1 - x^(1/nu)}
	else if(model == 'po') 	adjFun <- function(x, nu) {1 - x*(1/nu) / (x * 1/nu - x + 1)}
    adjQ <- adjFun(rawQ, nu)
    trueTimes <- qweibull(adjQ, shape = shape, scale = scale)
    
    obsTimes <- runif(n = n, max = inspectLength)
    if(!is.null(rndDigits))
    	obsTimes <- round(obsTimes, rndDigits)
    
    l <- rep(0, n)
    u <- rep(0, n)
    
    caught <- trueTimes < obsTimes
    u[caught] <- obsTimes[caught]
    l[!caught] <- obsTimes[!caught]
    
    if(inspections > 1){
    	for(i in 2:inspections){
		    oldObsTimes <- obsTimes
    		obsTimes <- oldObsTimes + runif(n, max = inspectLength)
		    if(!is.null(rndDigits))
    			obsTimes <- round(obsTimes, rndDigits)
    		caught <- trueTimes >= oldObsTimes  & trueTimes < obsTimes
    		needsCatch <- trueTimes > obsTimes
    		u[caught] <- obsTimes[caught]
    		l[needsCatch] <- obsTimes[needsCatch]
    	}
    }
    else{
    	needsCatch <- !caught	
    }
    u[needsCatch] <- Inf
    
    if(sum(l > u) > 0)	stop('warning: l > u! Bug in code')
    
    isCensored <- rbinom(n = n, size = 1, prob = prob_cen) == 1

    l[!isCensored] <- trueTimes[!isCensored]
    u[!isCensored] <- trueTimes[!isCensored]
    
    if(sum(l == Inf) > 0){
      allTimes <- c(l,u)
      allFiniteTimes <- allTimes[allTimes < Inf]
      maxFiniteTime <- max(allFiniteTimes)
      l[l == Inf] <- maxFiniteTime
    }
    return(data.frame(l = l, u = u, x1 = x1, x2 = x2))
}


simICPO_beta <- function(n = 100, b1 = 1, b2 = -1, inspections = 1, shape1 = 2, shape2 = 2, rndDigits = NULL){
	rawQ <- runif(n)
    x1 <- rnorm(n)
    x2 <- rbinom(n, 1, 0.5) - 0.5
    nu <- exp(x1 * b1 + x2 * b2)
    adjQ <- 1 - rawQ*(1/nu) / (rawQ * 1/nu - rawQ + 1)
    trueTimes <- qbeta(adjQ, shape1 = shape1, shape2 = shape2)
    
    inspectionError = 1 / (inspections + 1)
    obsTimes <- 1 / (inspections + 1) + runif(n, min = -inspectionError, max = inspectionError)
    if(!is.null(rndDigits))
    	obsTimes <- round(obsTimes, rndDigits)
    
    l <- rep(0, n)
    u <- rep(0, n)
    
    caught <- trueTimes < obsTimes
    u[caught] <- obsTimes[caught]
    l[!caught] <- obsTimes[!caught]
    
    if(inspections > 1){
    	for(i in 2:inspections){
		    oldObsTimes <- obsTimes
    		obsTimes <- i / (inspections+1) + runif(n, min = -inspectionError, max = inspectionError)
		    if(!is.null(rndDigits))
    			obsTimes <- round(obsTimes, rndDigits)
    		caught <- trueTimes >= oldObsTimes  & trueTimes < obsTimes
    		needsCatch <- trueTimes > obsTimes
    		u[caught] <- obsTimes[caught]
    		l[needsCatch] <- obsTimes[needsCatch]
    	}
    }
    else{
    	needsCatch <- !caught	
    }
    u[needsCatch] <- 1
    
    if(sum(l > u) > 0)	stop('warning: l > u! Bug in code')
    
    return(data.frame(l = l, u = u, x1 = x1, x2 = x2))
}





diag_covar <- function(object, varName, 
           data, model, weights = NULL,
           yType = 'meanRemovedTransform', 
           factorSplit = TRUE, 
           numericCuts, col, 
           xlab, ylab, main, 
           lgdLocation = NULL){
	if(!yType %in% c('survival', 'transform', 'meanRemovedTransform')) stop("yType not recognized. Options = 'survival', 'transform' or 'meanRemovedTransform'")
  if(missing(data)){
    if(!is(object, 'icenReg_fit')) stop("either object must be icenReg_fit, or formula with data supplied")
    data <- object$getRawData()
  }
  max_n_use <- nrow(data) #for backward compability. No longer need max_n_use
  subDataInfo <- subSampleData(data, max_n_use, weights)
	data <- subDataInfo$data
	weights <- subDataInfo$w
	if(is.null(weights)) weights <- rep(1, nrow(data))
	
	fullFormula <- getFormula(object)
	if(missing(model)) model <- object$model
	if(is.null(model)) stop('either object must be a fit, or model must be provided')
	if(missing(data))	data <- getData(object)
	if(missing(varName)){
		allVars <- getVarNames_fromFormula(fullFormula)
		nV <- length(allVars)
		k <- ceiling(sqrt(nV))
		if(k > 1) par(mfrow = c( ceiling(nV/k), k) )
		for(vn in allVars){
			useFactor <- length( unique((data[[vn]])) ) < 5
			diag_covar(object, vn, factorSplit = useFactor, model = model, data = data, yType = yType, weights = weights, lgdLocation = lgdLocation)
			}
		return(invisible(NULL))
	}

	if(model == 'ph')				s_trans <- function(x){ isOk <- x > 0 & x < 1
															ans <- numeric()
															ans[isOk] <- log(-log(x[isOk]) )
															ans
															}
	
	else if(model == 'po')		 	s_trans <- function(x){ isOk <- x > 0 & x < 1
															ans <- numeric()
															ans[isOk] <- log(x[isOk]/(1-x[isOk]))
															return(ans)
														   }
															

	allData <- data
	vals <- allData[[varName]]
	if(is.null(vals))	stop(paste('Cannot find variable', varName, 'in original dataset'))
	orgCall <- fullFormula
	reducCall <- removeVarFromCall(orgCall, varName)
	if(identical(orgCall, reducCall))	stop('varName not found in original call')
	
	spltInfo <- NULL
	if(factorSplit){
		factVals <- factor(vals)
		if(length(levels(factVals)) > 20) stop('Attempting to split on factor with over 20 levels. Try using numeric version of covariate and use numericCuts instead')
		spltInfo <- makeFactorSplitInfo(vals, levels(factVals))
	}
	else{
		if(missing(numericCuts))	numericCuts <- median(data[[varName]])
		spltInfo <- makeNumericSplitInfo(vals, numericCuts)
	}
	
	spltFits <- splitAndFit(newcall = reducCall, data = allData, varName = varName, 
							splitInfo = spltInfo, fitFunction = ic_sp, model = model, weights = weights)
		
	allX <- numeric()
	allY <- numeric()
	fitNames <- ls(spltFits)
	for(nm in fitNames){
		allX <- c(allX, as.numeric(spltFits[[nm]]$T_bull_Intervals) )
	}
	
	xlim <- range(allX, na.rm = TRUE, finite = TRUE)
	ylim <- sort( s_trans(c(0.025, 0.975)) )
	
	if(missing(col))
		col <- 1 + 1:length(fitNames)
	if(missing(xlab))	xlab = 't'
	if(missing(main)) 	main = varName
	if(missing(ylab)){
		if(model == 'ph'){
			ylab = 'log(-log(S(t)))'
			lgdLoc <- 'bottomright'
		}
		else if(model == 'po'){
			ylab = 'log(S(t)/(1 - S(t)))'
			lgdLoc <- 'bottomleft'	
		}
		else stop('model not recognized!')
	}
	t_vals <- xlim[1] + 1:999/1000 * (xlim[2] - xlim[1])
	estList <- list()
	meanVals <- 0
	if(yType == 'survival'){
		ylab = 'S(t)'
		lgdLoc <- 'bottomleft'
		s_trans <- function(x) x
		ylim = c(0,1)
	}
	if(yType == 'meanRemovedTransform'){
		ylab = paste('Mean Removed', ylab)
	}
		
	for(i in seq_along(fitNames)){
		if(yType == 'transform' | yType == 'survival'){
			nm <- fitNames[i]
			thisCurve <- getSCurves(spltFits[[nm]])
			theseTbls <- thisCurve$Tbull_ints
			thisS <- thisCurve$S_curves$baseline
			#thisCol <- col[i]
			#lines(theseTbls[,1], s_trans(thisS), col = thisCol, type = 's')
			#lines(theseTbls[,2], s_trans(thisS), col = thisCol, type = 's')
			estList[[nm]] <- list(x = theseTbls, y = s_trans(thisS) )
		}
		else if(yType == 'meanRemovedTransform'){	
			nm <- fitNames[i]
			estList[[nm]] <- s_trans(1 - getFitEsts(spltFits[[nm]], q = t_vals) )
			meanVals <- estList[[nm]] + meanVals
			}
	}

	if(yType == 'meanRemovedTransform'){
		meanVals <- meanVals/length(estList)
		ylim = c(Inf, -Inf)
		for(i in seq_along(estList)){
			theseLims <- range(estList[[i]] - meanVals, finite = TRUE, na.rm = TRUE)
			ylim[1] <- min(theseLims[1], ylim[1])
			ylim[2] <- max(theseLims[2], ylim[2])
		}
		
		yrange <- ylim[2] - ylim[1]
		ylim[2] = ylim[2] + yrange * 0.2
		ylim[1] = ylim[1] - yrange * 0.2
	}
		
	plot(NA, xlab = xlab, ylab = ylab, main = main, xlim = xlim, ylim = ylim)

	
	if(yType == 'meanRemovedTransform'){
		for(i in seq_along(estList)){
			lines(t_vals, estList[[i]] - meanVals, col = col[i])
		}
	}	
	else if(yType == 'transform' | yType == 'survival'){
		for(i in seq_along(estList)){
			xs <- estList[[i]]$x
			y  <- estList[[i]]$y
			lines(xs[,1], y, col = col[i], type = 's')
			lines(xs[,2], y, col = col[i], type = 's')
		}
	}
	if(!is.null(lgdLocation))	lgdLoc <- lgdLocation
	legend(lgdLoc, legend = fitNames, col = col, lwd = 1)
}




ic_par <- function(formula, data, model = 'ph', dist = 'weibull', weights = NULL){
  if(missing(data)) data <- environment(formula)
	cl <- match.call()
	mf <- match.call(expand.dots = FALSE)
#    m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(mf), 0L)
    m <- match(c("formula", "data", "subset", "na.action", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    x <- model.matrix(mt, mf, contrasts)
	if(is.matrix(x))	xNames <- colnames(x)
    else				xNames <- as.character(formula[[3]])
	if('(Intercept)' %in% colnames(x)){	
		ind = which(colnames(x) == '(Intercept)')
		x <- x[,-ind]
		xNames <- xNames[-ind]
	}
		
  yMat <- as.matrix(y)[,1:2]
  
  if(is(y, "Surv")){
    rightCens <- mf[,1][,3] == 0
	  yMat[rightCens,2] <- Inf
	
	  exact <- mf[,1][,3] == 1
	  yMat[exact, 2] = yMat[exact, 1]
  }
    storage.mode(yMat) <- 'double'
    
    if(sum(is.na(mf)) > 0)
    	stop("NA's not allowed. If this is supposed to be right censored (i.e. [4, NA] was supposed to be right censored at t = 4), replace NA with Inf")
        
    testMat <- cbind(x, 1)
    invertResult <- try(diag(solve(t(testMat) %*% testMat )), silent = TRUE)
    if(is(invertResult, 'try-error'))
	    stop('covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level')
	    
	callText <- paste(dist, model)

	if(is.null(weights))	weights = rep(1, nrow(yMat))
	if(length(weights) != nrow(yMat))	stop('weights improper length!')
	if(min(weights) < 0)				stop('negative weights not allowed!')
	if(sum(is.na(weights)) > 0)			stop('cannot have weights = NA')
	if(is.null(ncol(x))) recenterCovar = FALSE
   	fitInfo <- fit_par(yMat, x, parFam = dist, link = model, 
   	                   leftCen = 0, rightCen = Inf, uncenTol = 10^-6, 
   	                   regnames = xNames, weights = weights,
   	                   callText = callText)
	fitInfo$call = cl
	fitInfo$formula = formula
  fitInfo$.dataEnv = new.env()
  if(!missing(data)){ fitInfo$.dataEnv$data = data }
	fitInfo$par = dist
	fitInfo$model = model
  fitInfo$terms <- mt
  fitInfo$xlevels <- .getXlevels(mt, mf)

  return(fitInfo)
}

getFitEsts <- function(fit, newdata, p, q){
  if(missing(newdata)) newdata <- NULL
  etas <- get_etas(fit, newdata)
  
  if(missing(p))	p <- NULL
  if(missing(q))  q <- NULL
  if(!is.null(q)) {xs <- q; type = 'q'}
  else{ 
    type = 'p'
    if(is.null(p)) xs <- 0.5
    else		   xs <- p
  }
  
  if(length(etas) == 1){etas <- rep(etas, length(xs))}
  if(length(xs) == 1){xs <- rep(xs, length(etas))}
  if(length(etas) != length(xs) ) stop('length of p or q must match nrow(newdata) OR be 1')

  regMod <- fit$model
  
  if(inherits(fit, 'sp_fit'))	{
    scurves <- getSCurves(fit, newdata = NULL)
    baselineInfo <- list(tb_ints = scurves$Tbull_ints, s = scurves$S_curves$baseline)
    baseMod = 'sp'
  }
  if(inherits(fit, 'par_fit')){	
    baseMod <- fit$par
    baselineInfo <- fit$baseline
  }
  if(type == 'q')
    ans <- getSurvProbs(xs, etas, baselineInfo = baselineInfo, regMod = regMod, baseMod = baseMod)
  else if(type == 'p')
    ans <- getSurvTimes(xs, etas, baselineInfo = baselineInfo, regMod = regMod, baseMod = baseMod)
  return(ans)
}



getFitEsts_OLD <-function(fit, newdata, p, q){
	if(missing(newdata)) newdata <- NULL
	etas <- get_etas(fit, newdata)

	if(missing(p))	p <- NULL
	if(missing(q))  q <- NULL
	if(!is.null(q)) {xs <- q; type = 'q'}
	else{ 
		type = 'p'
		if(is.null(p)) xs <- 0.5
		else		   xs <- p
	}
	
	if(length(etas) == 1){etas <- rep(etas, length(xs))}
	if(length(xs) == 1){xs <- rep(xs, length(etas))}
	if(length(etas) != length(xs) ) stop('length of p or q must match nrow(newdata) OR be 1')
	
	if(inherits(fit, 'sp_fit'))	{
		scurves <- getSCurves(fit, newdata)
#		ans <- matrix(nrow = length(xs), ncol = length(etas))
		ans <- numeric()
#		colnames(ans) <- names(scurves$S_curves)
				
		for(i in 1:length(etas)){
# 			if(type == 'p') ans[,i] <- get_tbull_mid_q(xs, scurves[[2]][[i]], scurves[[1]])
# 			else ans[,i] 			<- get_tbull_mid_p(xs, scurves[[2]][[i]], scurves[[1]])
		  if(type == 'p') ans[i] <- get_tbull_mid_q(xs[i], scurves[[2]][[i]], scurves[[1]])
		  else ans[i] 			<- get_tbull_mid_p(xs[i], scurves[[2]][[i]], scurves[[1]])
		  
		}
		return(ans)
	}
	if(inherits(fit, 'par_fit')){	
		s_fun <- get_s_fun(fit)
		link_fun <- get_link_fun(fit)
#		ans <- matrix(nrow = length(xs), ncol = length(etas))
		ans <- numeric()
#		colnames(ans) <- names(etas)
#		rownames(ans) <- xs
		if(type == 'p'){
			optimReadyFun <- function(x, p, baselinePars, eta, s_fun, link_fun){
				s_o <- s_fun(x, baselinePars)
				s_c <- link_fun(s_o, eta)
				f <- 1 - s_c
				return( (f-p)^2 )
			}	
			
# 			for(i in 1:length(xs)){
# 				for(j in 1:length(etas)){
# 					upperBound <- findUpperBound(xs[i], s_fun, link_fun, fit, etas[j])
# 					ans[i,j] <- optimize(optimReadyFun, interval = c(0, upperBound), 
# 									p = xs[i], baselinePars = fit$baseline,
# 									etas[j], s_fun, link_fun, tol = 10^-6)$minimum
# 				
# 				}
# 			}
			upperBound = 1
			
			applyFUN <- function(i, xs, link_fun, s_fun, baseline, etas){
			  upperBound <- findUpperBound(upperBound,xs[i], s_fun, link_fun, fit, etas[i])
			   optimize(optimReadyFun, interval = c(0, upperBound), 
			                     p = xs[i], baselinePars = fit$baseline,
			                     etas[i], s_fun, link_fun, tol = 10^-6)$minimum
			}
			ans <- sapply(1:length(xs), FUN = applyFUN, 
			              xs = xs, etas = etas, 
			              link_fun = link_fun, s_fun = s_fun,
			              baseline = fit$baseline)
			
			
# 			for(i in 1:length(xs)){
# 			    upperBound <- findUpperBound(upperBound,xs[i], s_fun, link_fun, fit, etas[i])
# 			    ans[i] <- optimize(optimReadyFun, interval = c(0, upperBound), 
# 			                         p = xs[i], baselinePars = fit$baseline,
# 			                         etas[i], s_fun, link_fun, tol = 10^-6)$minimum
# 			    
# 			}
			return(ans)
		}
		
# 		for(i in 1:length(xs)){
# 			for(j in 1:length(etas))
# 				ans[i,j] <- 1-link_fun(s_fun(xs[i], fit$baseline), etas[j])
# 		}
  
		applyFUN <- function(i, xs, link_fun, s_fun, baseline, etas){ 1 - link_fun(s_fun(xs[i], baseline), etas[i])}
		ans <- sapply(1:length(xs), FUN = applyFUN, 
		              xs = xs, etas = etas, 
		              link_fun = link_fun, s_fun = s_fun,
		              baseline = fit$baseline)
# 		for(i in 1:length(xs)){
# 		    ans[i] <- 1-link_fun(s_fun(xs[i], fit$baseline), etas[i])
# 		}
		return(ans)
	}
	stop('getFitEsts not currently supported for this object')
}

diag_baseline <- function(object, data, model = 'ph', weights = NULL,
						  dists = c('exponential', 'weibull', 'gamma', 'lnorm', 'loglogistic', 'generalgamma'),
						  cols = NULL, lgdLocation = 'bottomleft',
						  useMidCovars = T){
	newdata = NULL
	if(useMidCovars) newdata <- 'midValues'
	formula <- getFormula(object)
	if(missing(data))	data <- getData(object)
	max_n_use = nrow(data)	#no longer necessary, for backward compatability				  	
	
	subDataInfo <- subSampleData(data, max_n_use, weights)
	sp_data <- subDataInfo$data
	weights <- subDataInfo$w

	sp_fit <- ic_sp(formula, data = sp_data, bs_samples = 0, model = model)
	plot(sp_fit, newdata)
	xrange <- range(getSCurves(sp_fit)$Tbull_ints, finite = TRUE)
	grid <- xrange[1] + 0:100/100 *(xrange[2] - xrange[1])
	if(is.null(cols)) cols <- 1 + 1:length(dists)
	for(i in seq_along(dists)){
		this_dist <- dists[i]
		par_fit <- ic_par(formula, data = data, model = model, dist = this_dist)
		y <- getFitEsts(par_fit, newdata = newdata, q = grid)
		lines(grid, 1 - y, col = cols[i])
	}
	legend(lgdLocation, legend = c('Semi-parametric', dists), col = c('black', cols), lwd = 1)
}

predict.icenReg_fit <- function(object, type = 'response',
                                newdata = NULL, ...)
      #imputeOptions = fullSample, fixedParSample, median
  {
  if(is.null(newdata)) newdata <- object$getRawData()
  if(type == 'lp')
    return( log(get_etas(object, newdata = newdata)))
  if(type == 'response')
    return(getFitEsts(fit = object, newdata = newdata))
  stop('"type" not recognized: options are "lp", "response" and "impute"')
}


imputeCens<- function(fit, newdata = NULL, imputeType = 'fullSample', numImputes = 5){
#  if(is(fit, 'sp_fit'))
#    stop("imputation not available for semi-parametric model. 
#         This is due to the lack of distributional theory for baseline distribution")
  if(is.null(newdata)) newdata <- fit$getRawData()
  yMat <- expandY(fit$formula, newdata, fit)
  p1 <- getFitEsts(fit, newdata, q = as.numeric(yMat[,1]) ) 
  p2 <- getFitEsts(fit, newdata, q = as.numeric(yMat[,2]) ) 
  ans <- matrix(nrow = length(p1), ncol = numImputes)
  storage.mode(ans) <- 'double'
  if(imputeType == 'median'){
    p_med <- (p1 + p2)/2
    ans <- getFitEsts(fit, newdata, p = p_med)
    isLow <- ans < yMat[,1]
    ans[isLow] <- yMat[isLow,1]
    isHi <- ans > yMat[,2]
    ans[isHi] <- yMat[isHi]
    return()
  }
  if(imputeType == 'fixedParSample'){
    for(i in 1:numImputes){
      p_samp <- runif(length(p1), p2, p1)
      theseImputes <- getFitEsts(fit, newdata, p = p_samp)
      isLow <- theseImputes < yMat[,1]
      theseImputes[isLow] <- yMat[isLow,1]
      isHi <- theseImputes > yMat[,2]
      theseImputes[isHi] <- yMat[isHi,2]
      ans <- fastMatrixInsert(theseImputes, ans, colNum = i)
    }
    return(ans)
  }
  if(imputeType == 'fullSample'){
    isSP <- is(fit, 'sp_fit')
    for(i in 1:numImputes){
      orgCoefs <- getSamplablePars(fit)
      if(!isSP){
        coefVar <- getSamplableVar(fit)
        sampledCoefs <- sampPars(orgCoefs, coefVar)
      }
      else{
        sampledCoefs <- getBSParSample(fit)
      }
      setSamplablePars(fit, sampledCoefs)
      p1 <- getFitEsts(fit, newdata, q = as.numeric(yMat[,1]) ) 
      p2 <- getFitEsts(fit, newdata, q = as.numeric(yMat[,2]) ) 
      p_samp <- runif(length(p1), p1, p2)
      theseImputes <- getFitEsts(fit, newdata, p = p_samp)
      isLow <- theseImputes < yMat[,1]
      theseImputes[isLow] <- yMat[isLow,1]
      isHi <- theseImputes > yMat[,2]
      theseImputes[isHi] <- yMat[isHi,2]
      fastMatrixInsert(theseImputes, ans, colNum = i)
      setSamplablePars(fit, orgCoefs)
    }
    return(ans)
  }
  stop('imputeType type not recognized.')
}

plot.sp_curves <- function(x, sname = 'baseline', xRange = NULL, ...){
  if(is.null(xRange))
    xRange <- range(c(x[[1]][,1], x[[1]][,2]), finite = TRUE)
  dotList <- list(...)
  addList <- list(xlim = xRange, ylim = c(0,1), x = NA)
  dotList <- addListIfMissing(addList, dotList)
  do.call(plot, dotList)
  lines(x, sname = sname, ...)
}

lines.ic_npList <- function(x, fitNames = NULL, ...){
  if(is.null(fitNames)){
    fitNames <- names(x$scurves)
    lines(x, fitNames, ...)
  }
  dotList <- list(...)
  cols <- dotList$col

  for(i in seq_along(fitNames)){
    thisName <- fitNames[i]
    dotList$col <- cols[i]
    dotList$x <- x$scurves[[thisName]]
    do.call(lines, dotList)
  }
}

plot.ic_npList <- function(x, fitNames = NULL, lgdLocation = 'bottomleft', ... ){
  addList <- list(xlim = x$xRange,
                  ylim = c(0,1),
                  xlab = 't', 
                  ylab = 'S(t)', 
                  x = NA)
  dotList <- list(...)
  dotList <- addListIfMissing(addList, dotList)
  do.call(plot, dotList)  
  grpNames <- names(x$fitList)
  cols <- dotList$col
  if(is.null(cols)) cols = 2:(length(grpNames) + 1)
  if(length(cols) != length(grpNames)) 
    stop('length of number of strata not equal to number of colors')
  dotList$col <- cols
  dotList$fitNames = fitNames
  dotList$x <- x
  do.call(lines, dotList)
  legend(lgdLocation, legend = grpNames, col = cols, lty = 1)
}

lines.sp_curves <- function(x, sname = 'baseline',...){
  firstTimeObs <- x[[1]][1,1]
  firstTimeAssume <- firstTimeObs
  if(firstTimeObs > 0)
    firstTimeAssume <- 0
  lines(c(firstTimeAssume, firstTimeObs), c(1,1), ...)
  lines(x[[1]][,1], x[[2]][[sname]], ..., type = 's')
  lines(x[[1]][,2], x[[2]][[sname]],   ..., type = 's')
  lastObs <- nrow(x[[1]])
  lastTimes <- x[[1]][lastObs,]
  if(lastTimes[2] == Inf) lastTimes[2] <- lastTimes[1]
  lastTimes[2] <- lastTimes[2] + (lastTimes[2] - firstTimeObs)
  lines(lastTimes, c(0,0), ... ) 
}

dGeneralGamma <- function(x, mu, s, Q){
  max_n <- getMaxLength(list(x, mu, s, Q) )
  x <- updateDistPars(x, max_n)
  mu <- updateDistPars(mu, max_n)
  s <- updateDistPars(s, max_n)
  Q <- updateDistPars(Q, max_n)
  
  ans <- .Call('dGeneralGamma', x, mu, s, Q)
  return(ans)
}

qGeneralGamma <- function(p, mu, s, Q){
  x <- p
  max_n <- getMaxLength(list(x, mu, s, Q) )
  x <- updateDistPars(x, max_n)
  mu <- updateDistPars(mu, max_n)
  s <- updateDistPars(s, max_n)
  Q <- updateDistPars(Q, max_n)
  
  ans <- .Call('qGeneralGamma', x, mu, s, Q)
  return(ans)
  
}

pGeneralGamma <- function(q, mu, s, Q){
  x <- q
  max_n <- getMaxLength(list(x, mu, s, Q) )
  x <- updateDistPars(x, max_n)
  mu <- updateDistPars(mu, max_n)
  s <- updateDistPars(s, max_n)
  Q <- updateDistPars(Q, max_n)
  
  ans <- .Call('qGeneralGamma', x, mu, s, Q)
  return(ans)
  
}


ic_np <- function(formula = NULL, data, maxIter = 1000, tol = 10^-10){
  if(is.null(formula)){ return(ic_npSINGLE(data, maxIter = maxIter, tol = tol)) }
  if(!inherits(formula, 'formula')) {
    #Covering when user ONLY provides data as first unlabeled argument
    data <- formula
    return(ic_npSINGLE(data, maxIter = maxIter, tol = tol))
  }
  
  if(missing(data)) data <- environment(formula)
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  yMat <- as.matrix(y)[,1:2]
  if(is(y, 'Surv')){
    rightCens <- mf[,1][,3] == 0
    yMat[rightCens,2] <- Inf
    exact <- mf[,1][,3] == 1
    yMat[exact, 2] = yMat[exact, 1]
  }
  storage.mode(yMat) <- 'double'
  
  
  
  formFactor <- formula[[3]]
  if( length(formFactor) != 1 ){ 
    stop('predictor must be either single factor OR 0 for ic_np')
  }
  if(formFactor == 0){ return(ic_npSINGLE(yMat, maxIter = maxIter, tol = tol)) }
  thisFactor <- data[[ as.character(formFactor) ]]
  if(!is.factor(thisFactor)){ stop('predictor must be factor') }
  
  theseLevels <- levels(thisFactor)
  fitList <- list()
  
  for(thisLevel in theseLevels){
    thisData <- yMat[thisFactor == thisLevel, ]
    if(nrow(thisData) > 0)
      fitList[[thisLevel]] <- ic_npSINGLE(thisData, maxIter = maxIter, tol = tol)
  }
  ans <- ic_npList(fitList)
  return(ans)
}

ic_npSINGLE <- function(data,  maxIter = 1000, tol = 10^-10){
  data <- as.matrix(data)
  if(ncol(data) != 2) stop("data should be an nx2 matrix or data.frame")
  if(any(data[,1] > data[,2]) ) stop(paste0("data[,1] > data[,2].",
                                          "This is impossible for interval censored data") )
  storage.mode(data) <- "double"
  mis <- findMaximalIntersections(data[,1], data[,2])
  fit <- .Call("EMICM", mis$l_inds, mis$r_inds, as.integer(maxIter))
  tbulls <- rbind(mis$mi_l, mis$mi_r)
  ans <- new('ic_np')
  #ans <- list(phat = fit[[1]], Tbull_ints = tbulls, llk = fit[[2]], iters = fit[[3]])
  ans$p_hat <- fit[[1]]
  ans$T_bull_Intervals <- tbulls
  ans$coefficients <- numeric()
  ans$llk <- fit[[2]]
  ans$iterations <- fit[[3]]
  ans$par <- 'non-parametric'
  ans$model = 'none'
  ans$var <- matrix(nrow = 0, ncol = 0) 
  dataEnv <- new.env()
  dataEnv[['data']] <- data
  ans[['.dataEnv']] <- dataEnv
  return(ans)
}