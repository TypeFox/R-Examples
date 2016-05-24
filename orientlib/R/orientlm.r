
orientlm <- function(observed, leftformula, trueorient = rotmatrix(diag(3)), rightformula,
                     data = list(), subset, weights, na.action, iterations = 5) {
    if (missing(leftformula) && missing(rightformula)) stop('Must give left or right formula (or both)')
    
    mf <- match.call(expand.dots = FALSE)
    mf$observed <- mf$leftformula <- mf$trueorient <- mf$rightformula <- NULL
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    
    if (!missing(leftformula)) {
    	mf$formula <- leftformula
    	leftframe <- eval(mf, parent.frame())
    	leftmatrix <- model.matrix(leftformula, leftframe)
    	if (nrow(leftmatrix) < length(observed)) 
    	   leftmatrix <- rbind(leftmatrix, matrix(1,length(observed)-nrow(leftmatrix),1))
    } else leftmatrix <- NULL
    
    if (!missing(rightformula)) {
		mf$formula <- rightformula
		rightframe <- eval(mf, parent.frame())
		rightmatrix <- model.matrix(rightformula, rightframe)
		if (nrow(rightmatrix) < length(observed)) 
			rightmatrix <- rbind(rightmatrix, matrix(1,length(observed)-nrow(rightmatrix),1))
    } else rightmatrix <- NULL
    
    if (!missing(subset)) {
		observed <- observed[subset]
		if (length(trueorient) > 1) trueorient <- trueorient[subset]
		if (!missing(weights)) weights <- weights[subset]
    }
    if (missing(weights)) weights <- rep(1,length(observed))
    
    # Now fit the model Vi = A1^t Ui A2
    
    if (missing(rightformula)) {   # Left-sided model  Vi = A1^t Ui 
 		obs <- observed %*% trueorient^(-1)
 		fit <- lm( rotvector(obs)@x ~ leftmatrix - 1, weights = weights)
 		A1 <- t(nearest.SO3(array(t(predict(fit)), c(3,3,length(observed)))))
 		predict <- t(A1) %*% trueorient
 		return(list(leftfit = fit, rightfit= NULL, A1 = A1, A2 = NULL, predict=predict))
    } 
    if (missing(leftformula)) {   # right-sided model  Vi = Ui A2
	 	obs <- trueorient^(-1) %*% observed
	 	fit <- lm( rotvector(obs)@x ~ rightmatrix - 1, weights = weights)
	 	A2 <- nearest.SO3(array(t(predict(fit)), c(3,3,length(observed))))
	 	predict <- trueorient %*% A2
	 	return(list(leftfit = NULL, rightfit = fit, A1 = NULL, A2 = A2, predict=predict))
    }
    # two-sided model
    A2 <- rotmatrix(diag(3))
    for (i in 1:iterations) {
    	obs <- observed %*% t(A2) %*% t(trueorient)
    	leftfit <- lm( rotvector(obs)@x ~ leftmatrix - 1, weights = weights)
 		A1 <- t(nearest.SO3(array(t(predict(leftfit)), c(3,3,length(observed)))))
 		
 		obs <- t(trueorient) %*% A1 %*% observed
 		rightfit <- lm( rotvector(obs)@x ~ rightmatrix - 1, weights = weights)
	 	A2 <- nearest.SO3(array(t(predict(rightfit)), c(3,3,length(observed))))
 	}
 	predict <- t(A1) %*% trueorient %*% A2
 	return(list(leftfit=leftfit, rightfit=rightfit, A1 = A1, A2 = A2, predict=predict))
}
 
