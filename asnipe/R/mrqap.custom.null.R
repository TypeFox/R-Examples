
mrqap.custom.null <- function (formula, random.y, intercept = TRUE, directed = "undirected", diagonal = FALSE, test.statistic = "t-value", tol = 1e-07) {
	
	form <- formula
	formula <- model.frame(formula, na.action=NULL)
	i <- attr(attr(formula, 'terms'), 'response')
	y <- as.matrix(formula[i])
	x_names <- attr(attr(formula, 'terms'),'term.labels')
	x <- array(NA,c(nrow(formula[[x_names[1]]]),ncol(formula[[x_names[1]]]),length(x_names)))
	for (i in 1:length(x_names)) {
		x[,,i] <- as.matrix(formula[[x_names[i]]])
	}
	
	if (is.null(x_names)) {
		if (intercept) {
			x_names <- c("intercept",paste("x",1:(dim(x)[3]),sep=""))
		} else {
			x_names <- paste("x",1:(dim(x)[3]),sep="")
		}
	} else {
		if (intercept) {
			x_names <- c("intercept",x_names)
		}
	}
	
	if (ncol(y) != ncol(random.y)) {
		stop("Random y matrices have wrong dimentsions!")
	} else {
		randomisations <- dim(random.y)[1]
	}
		
	getfit <- function(y, x, tol) {
        newy <- matrix(c(y),ncol=1)
        newx <- vector()
        for (i in 1:dim(x)[3]) newx <- cbind(newx, matrix(c(x[,,i]),ncol=1))
        keep <- is.na(newy) | apply(is.na(newx), 1, any)
    	list(qr(newx[!keep, ], tol = tol), newy[!keep])
    }
                
	n <- dim(y)[2]
	
    which_xs <- matrix(TRUE, n, n)
    if (!diagonal) diag(which_xs) <- FALSE
    if (directed == "undirected") which_xs[upper.tri(which_xs)] <- FALSE


	if (intercept) {
		newx <- array(NA,c(n,n,dim(x)[3]+1))
		newx[,,1] <- matrix(1,n,n)
		newx[,,2:(dim(x)[3]+1)] <- x
		x <- newx
	}
	nx <- dim(x)[3]	
	
	if (nx == 1) {
		stop("x must contain more than one predictor variable")
	}
	
	if (!diagonal) {
		for (i in 1:dim(x)[3]) diag(x[,,i]) <- NA
		diag(y) <- NA
	}
	
	if (directed=="undirected") {
		for (i in 1:dim(x)[3]) x[,,i][upper.tri(x[,,i])] <- NA
		y[upper.tri(y)] <- NA
	}
	
	get_fit_obs <- getfit(y, x, tol)
	fit_obs <- list()
	fit_obs$formula <- form
	fit_obs$coefficients <- qr.coef(get_fit_obs[[1]], get_fit_obs[[2]])
	fit_obs$fitted.values <- qr.fitted(get_fit_obs[[1]], get_fit_obs[[2]])
    fit_obs$residuals <- qr.resid(get_fit_obs[[1]], get_fit_obs[[2]])
    fit_obs$rank <- get_fit_obs[[1]]$rank
    fit_obs$n <- length(get_fit_obs[[2]])
    fit_obs$df.residual <- fit_obs$n - fit_obs$rank

	if (test.statistic == "beta") {
        fit_obs$test.statistic <- fit_obs$coefficients
    } else if (test.statistic == "t-value") {
        fit_obs$test.statistic <- fit_obs$coefficients/sqrt(diag(chol2inv(get_fit_obs[[1]]$qr)) * sum(fit_obs$residuals^2)/(fit_obs$n - fit_obs$rank))
	}

	test.statistic_rand <- matrix(0, randomisations, nx)
	for (i in 1:nx) {
		get_x_correlation <- getfit(x[,,i],x[,,c((1:nx)[-i]),drop=FALSE], tol = tol)
		x_residuals <- x[,,i]
		x_residuals[which_xs] <- qr.resid(get_x_correlation[[1]], get_x_correlation[[2]])[which_xs]
		if (directed == "undirected") x_residuals[upper.tri(x_residuals)] <- t(x_residuals)[upper.tri(x_residuals)]
		
		for (j in 1:randomisations) {
			rany <- random.y[j,,]
			if (!diagonal) {
				diag(rany) <- NA
			}
	
			if (directed=="undirected") {
				rany[upper.tri(rany)] <- NA
			}
			
			get_fit_rand <- getfit(rany, x, tol = tol)
			
			if (test.statistic == "beta") {
        			test.statistic_rand[j, i] <- qr.coef(get_fit_rand[[1]], get_fit_rand[[2]])[i]
    			} else if (test.statistic == "t-value") {
    				fit_rand <- list()
				fit_rand$coefficients <- qr.coef(get_fit_rand[[1]], get_fit_rand[[2]])
				fit_rand$fitted.values <- qr.fitted(get_fit_rand[[1]], get_fit_rand[[2]])
				fit_rand$residuals <- qr.resid(get_fit_rand[[1]], get_fit_rand[[2]])
				fit_rand$rank <- get_fit_rand[[1]]$rank
				fit_rand$n <- length(get_fit_rand[[2]])
				fit_rand$df.residual <- fit_obs$n - fit_obs$rank
        			test.statistic_rand[j, i] <- (fit_rand$coefficients/sqrt(diag(chol2inv(get_fit_rand[[1]]$qr)) * sum(fit_rand$residuals^2)/(fit_rand$n - fit_rand$rank)))[i]
			}
		}			
	}
	fit_obs$P.greater <- apply(sweep(test.statistic_rand, 2, fit_obs$test.statistic, ">="), 2, mean)
	fit_obs$P.lesser <- apply(sweep(test.statistic_rand, 2, fit_obs$test.statistic, "<="), 2, mean)
	fit_obs$P.values <- apply(sweep(abs(test.statistic_rand), 2, abs(fit_obs$test.statistic), ">="), 2, mean)
	fit_obs$AIC <- (length(x_names)-1) - log(abs(prod(diag(get_fit_obs[[1]]$qr))))*nrow(x)
	
	names(fit_obs$P.values) <- x_names
	
	class(fit_obs) <- "mrqap.dsp"
	
	return(fit_obs)
	
}