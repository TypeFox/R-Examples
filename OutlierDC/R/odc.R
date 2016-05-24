##########################################################################
#
#    Outlier Detection in Censored Data using Quantile Regression 
#
#      by
#      Soo-Heang Eo and HyungJun Cho
#      Deparment of Statistics 
#      Korea University
#
#      Feb 2013
#
##########################################################################

odc <- function(formula, data, 
	method = c("score", "boxplot","residual"), 
	rq.model = c("Wang", "PengHuang", "Portnoy"), 
	#bound = c("UB", "LB", "both" ),
	k_r = 1.5, k_b = 1.5, h = .05){

    ##########
    #Preparation
    call <- match.call()
    rq.model <- match.arg(rq.model)
	method <- match.arg(method)
    formula <- Formula(formula)
	#bound <- match.arg(bound)
	bound = "UB"

	if(is.data.frame(data)) data <- as.data.frame(data)
	mf1 <- model.frame(formula, data = data)
	X.mat <- model.matrix(formula, data = mf1)
	resp <- model.response(mf1)
	y <- resp[,1]
	status <- resp[,2]
    n <- nrow(data)
    rownames(data) <- 1:n

	# Set parallel computing in order to incease computing power
	# Start constructing data matrix
    # End constructing data
	cat("Please wait... \n")
	if(rq.model != "Wang") fit <- crq(formula, data = data, method = rq.model)

	result <- new("OutlierDC")
	result@call <- call
	result@formula <- formula
	
	betas <- matrix(NA, nrow = ncol(X.mat), ncol = 5)
	if(rq.model == "Wang"){
		betas[,1] <- LCRQ(y = y, x = X.mat[ ,-1], delta = status, tau = .10, h= h)
    	betas[,2] <- LCRQ(y = y, x = X.mat[ ,-1], delta = status, tau = .25, h= h)
    	betas[,3] <- LCRQ(y = y, x = X.mat[ ,-1], delta = status, tau = .50, h= h)
    	betas[,4] <- LCRQ(y = y, x = X.mat[ ,-1], delta = status, tau = .75, h= h)
    	betas[,5] <- LCRQ(y = y, x = X.mat[ ,-1], delta = status, tau = .90, h= h)
	}
	else {
    	betas[,1] <- coef(fit, taus = .10)
    	betas[,2] <- coef(fit, taus = .25)
    	betas[,3] <- coef(fit, taus = .50)
		betas[,4] <- coef(fit, taus = .75)
    	betas[,5] <- coef(fit, taus = .90)
	}
	
	fit.q10 <- (X.mat %*%  betas[,1])[,1]
	fit.q25 <- (X.mat %*%  betas[,2])[,1]
	fit.q50 <- (X.mat %*%  betas[,3])[,1]
	fit.q75 <- (X.mat %*%  betas[,4])[,1]
	fit.q90 <- (X.mat %*%  betas[,5])[,1]
	fitted.mat <- cbind(fit.q10, fit.q25, fit.q50, fit.q75, fit.q90)

	if(method == "score"){
	    # Calculate outlier score
		score <- ifelse(y >= fit.q50, 
			(y - fit.q50) / (fit.q75 - fit.q50), 
			(y - fit.q50) / (fit.q50 - fit.q25))
		outlier <- FALSE
		result@score <- score
		n.outliers <- as.integer(0)
	}
    if(method == "boxplot"){
    	# calculate fences
    	iqr <- fit.q75 -fit.q25
    	upper.fence <- fit.q75 + (k_b * iqr)
    	lower.fence <- fit.q25 - (k_b * iqr)

    	if(bound == "both"){
	    	outlier <- ifelse(y > fit.q50,
	    				y > upper.fence, y < lower.fence)
	    	outlier <- ifelse(is.na(outlier), FALSE, outlier)
    	}
    	else if(bound == "LB"){
    		# calculate only lower bound
	    	outlier <- y <= lower.fence
	    	outlier <- ifelse(is.na(outlier), FALSE, outlier)
    	}
    	else if(bound == "UB"){
    		# calculate only upper bound
	    	outlier <- y >= upper.fence
	    	outlier <- ifelse(is.na(outlier), FALSE, outlier)
    	}

		result@lower <- lower.fence
		result@upper <- upper.fence
		n.outliers <- sum(outlier)
    }
    else if(method == "residual"){
    	score <- y - fit.q50
    	score <- score
    	cutoff <- median(abs(score) / qnorm(0.75), na.rm = TRUE)

    	if (bound == "both")	outlier <- abs(score) > (k_r * cutoff)
    	else if (bound == "UB") outlier <- score > (k_r * cutoff)
    	else if (bound == "LB") outlier <- score < -1 * (k_r * cutoff)

    	result@score <- score # absolute value of residuals
    	result@cutoff  <- cutoff
    	n.outliers <- sum(outlier)
    }
	## Declare outliers using D statistics
	## manipulate results object
	rownames(betas) <- colnames(X.mat)
	colnames(betas) <- c("q10", "q25", "q50", "q75", "q90")
	betas <- as.data.frame(betas)   
	
	result@raw.data <- data
	result@coefficients <- betas
	result@fitted.mat <- fitted.mat
	result@method <- method
	result@rq.model <- rq.model
	result@outliers <- outlier	
	result@n.outliers <- n.outliers
	result@refined.data <- data[!outlier,, drop = FALSE]
	result@outlier.data <- data[outlier,, drop = FALSE]
	result@k_r <- k_r
	result@k_b <- k_b
	result@bound <- bound	
	cat("Done. \n")
	return(result)
}
#END################################################################s
