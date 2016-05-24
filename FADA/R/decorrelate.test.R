decorrelate.test <- function(faobject, data.test) {
	nbclass <- length(unique(faobject$groups))
	mod <- faobject$mod.decorrelate.test
	proba.test <- predict(mod, data.test$x, type = "response")

	if (nbclass == 2) {
		proba.test <- matrix(c(1 - proba.test, proba.test), ncol = 2, byrow = FALSE)
	}
	proba.test <- matrix(proba.test, ncol = nbclass, byrow = FALSE)

	if (is.null(faobject$B)) {
		return(list(meanclass = faobject$meanclass, fa.training = faobject$fa.training, 
			fa.testing = data.test$x, Psi = faobject$Psi, B = faobject$B, factors.training = faobject$factors.training, 
			factors.testing = NULL, groups = faobject$groups, proba.training = faobject$proba.training, 
			proba.testing = proba.test, mod.decorrelate.test = faobject$mod.decorrelate))
	}

	sB <- scale(t(faobject$B), center = FALSE, scale = sqrt(faobject$Psi))
	G <- solve(diag(ncol(faobject$B)) + tcrossprod(sB))
	sB <- scale(t(faobject$B), center = FALSE, scale = faobject$Psi)

	testzclass <- lapply(1:nbclass, function(i) {
		tcrossprod(tcrossprod(scale(data.test$x, center = faobject$meanclass[, i], 
			scale = FALSE), sB), G)
	})
	testzclass <- lapply(1:nbclass, function(i) {
		sweep(testzclass[[i]], 1, proba.test[, i], "*")
	})
	testz <- Reduce("+", testzclass)

	fatest <- data.test$x - tcrossprod(testz, faobject$B)
	l <- list(meanclass = faobject$meanclass, fa.training = faobject$fa.training, 
		fa.testing = fatest, Psi = faobject$Psi, B = faobject$B, factors.training = faobject$factors.training, factors.testing = testz, groups = faobject$groups, proba.training = faobject$proba.training, proba.testing = proba.test, mod.decorrelate.test = faobject$mod.decorrelate, data.train= faobject$data.train, 
			maxnbfactors = faobject$maxnbfactors, min.err = faobject$min.err, EM = faobject$EM, maxiter = faobject$maxiter)	
	return(l)

}
