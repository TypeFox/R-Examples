build.dist.struct <-
function(z, X, exact = NULL, calip.option = 'propensity', calip.cov = NULL, caliper = 0.2, verbose = FALSE){
	
	cal.penalty <- 100
	if(is.null(exact)) exact = rep(1, length(z))	
	if(!(calip.option %in% c('propensity','user','none'))){
		stop('Invalid calip.option specified.')
	}
    if (is.vector(X)) X <- matrix(X, length(X), 1)
    if(!(length(z) == (dim(X)[1]))){
    	stop("Length of z does not match row count in X")
    }
    if(!(length(exact) == length(z))){
    	stop("Length of exact does not match length of z")
    }
    if(!(all((z == 1) | (z == 0)))){
    	stop("The z argument must contain only 1s and 0s")
    }
   	   
	if(is.data.frame(X) || is.character(X)){
		if(!is.data.frame(X)) X <- as.data.frame(X)
		X.chars <- which(laply(X, function(y) 'character' %in% class(y)))
		if(length(X.chars) > 0){
			if (verbose) print('Character variables found in X, converting to factors.')
			for(i in X.chars){
				X[,i] <- factor(X[,i])
				
			}
		}
	    #if some variables are factors convert to dummies
	     X.factors <-  which(laply(X, function(y) 'factor' %in% class(y)))
	     
   		#handle missing data
   		for(i in which(laply(X, function(x) any(is.na(x))))){
			if (verbose) print(paste('Missing values found in column', i ,'of X; imputing and adding missingness indicators'))
   			if(i %in% X.factors){
   				#for factors, make NA a new factor level
   				X[,i] <- addNA(X[,i])
   			}else{
   				#for numeric/logical, impute means and add a new indicator for missingness
   				X[[paste(colnames(X)[i],'NA', sep = '')]] <- is.na(X[,i])
   				X[which(is.na(X[,i])),i] <- mean(X[,i], na.rm = TRUE)
   			}
   		}
		for(i in rev(X.factors)){
	     	dummyXi <- model.matrix(as.formula(
	     		paste('~',colnames(X)[i], '-1')),data=X)
	     	X <- cbind(X[,-i], dummyXi)
	    }
	      
    }else{
    	#handle missing data
    	for(i in c(1:ncol(X))){
    		if(any(is.na(X[,i]))){
   				X <- cbind(X,is.na(X[,i]))
   				colnames(X)[ncol(X)] <- paste(colnames(X)[i],'NA', sep = '')
   				X[which(is.na(X[,i])),i] <- mean(X[,i], na.rm = TRUE)    
    		}
    	}
			
	}
    #get rid of columns that do not vary
	varying <- apply(X,2, function(x) length(unique(x)) > 1)
	if(!all(varying) && verbose) print('Constant-value columns found in X, they will not be used to calculate Mahalanobis distance.')
	X <- X[,which(varying),drop = FALSE]
	
        
    if (calip.option == 'propensity') {
        calip.cov <- glm.fit(cbind(rep(1, nrow(X)),X), z, family = binomial())$linear.predictors
        cal <- sd(calip.cov) * caliper
    }else if(calip.option == 'user'){
    		stopifnot(!is.null(calip.cov))
    		cal <- sd(calip.cov) * caliper
    }
    nobs <- length(z)
    rX <- as.matrix(X)
    for (j in 1:(dim(rX)[2])) rX[, j] <- rank(rX[, j])
    cv <- cov(rX)
    vuntied <- var(1:nobs)
    rat <- sqrt(vuntied/diag(cv))
    if(length(rat) == 1){
    		cv <- as.matrix(rat) %*% cv %*% as.matrix(rat)
	}else{
		cv <- diag(rat) %*% cv %*% diag(rat)
	}
    #library(MASS)
    icov <- ginv(cv)
    nums <- 1:nobs
    ctrl.nums <- 1:(sum(z == 0))
    treated <- nums[z == 1]
    
    #find distance between each treated and each control it will be connected to and store in a distance structure
    dist.struct <- list()
    for (i in c(1:length(treated))) {
        controls <- nums[(z == 0) & (exact == exact[treated[i]])]
        control.names <- ctrl.nums[exact[z == 0] == exact[treated[i]]]
        costi <- mahalanobis(rX[controls, ,drop=FALSE], rX[treated[i], ], icov, inverted = T)
        if (calip.option != 'none') {
        		calip.update <- rep(0, length(costi))
        		calip.update[abs(calip.cov[treated[i]] - calip.cov[controls]) - cal > 0] <- Inf
        		costi <- costi + calip.update
        	}
        names(costi) <- control.names
		dist.struct[[i]] <- costi[is.finite(costi)]	
    }

	if (sum(laply(dist.struct, length)) == 0) stop('All matches forbidden. Considering using a wider caliper?')
	return(dist.struct)
}
