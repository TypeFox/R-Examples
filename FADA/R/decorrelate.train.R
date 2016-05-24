decorrelate.train <- function(data.train, nbf = NULL, maxnbfactors=12, diagnostic.plot = FALSE, 
min.err = 0.001, verbose = TRUE,EM = TRUE, maxiter = 15,...) {
	p <- ncol(data.train$x)
	n <- nrow(data.train$x)
	nbclass <- length(unique(data.train$y))
	cl <- sort(unique(data.train$y))
	if (!all(cl == c(1:nbclass))) {
		stop("Group variable must be 1,2, ...")
	}
	if (p <= 3) {
		stop("The number of variables must be at least 4.")
	}
	meanclass <- sapply(1:nbclass, function(i) {
		colMeans(data.train$x[data.train$y == i, ])
	})
	cdta <- data.train$x
	for (i in 1:nbclass) {
		cdta[data.train$y == i, ] <- sweep(data.train$x[data.train$y == i, ], 2, 
			meanclass[, i], "-")
	}
	if (is.null(nbf)) {
		nbf <- nbfactors(scale(cdta, center = FALSE, scale = TRUE), maxnbfactors = maxnbfactors, 
			EM = EM, diagnostic.plot = diagnostic.plot, jumps.nbfactor = 0.05)$optimalnbfactors
	}
	if (verbose) 
		print(paste("Number of factors:", nbf, "factors", sep = " "))
	coefpx <- LassoML(data.train,...)
	mod <- coefpx$mod
	if (nbf == 0) {
		return(list(meanclass = meanclass, fa.training = data.train$x, Psi = colVars(cdta), 
			B = NULL, factors.training = NULL, groups = data.train$y, proba.training = coefpx$proba.train, 
			mod.decorrelate.test = mod, data.train = data.train$x))
	} else {
		eps <- 1
		iter <- 0
		proba.train <- coefpx$proba.train
		if (verbose) 
			print("Objective criterion: ")
		while (eps > min.err & iter <= maxiter) {
			fa <- emfa(cdta, nbf = nbf, minerr = min.err, verbose = FALSE, EM = EM)
			sB <- scale(t(fa$B), center = FALSE, scale = sqrt(fa$Psi))
			G <- solve(diag(ncol(fa$B)) + tcrossprod(sB))
			sB <- scale(t(fa$B), center = FALSE, scale = fa$Psi)
			zclass <- lapply(1:nbclass, function(i) {
				sweep(tcrossprod(tcrossprod(scale(data.train$x, center = meanclass[, 
					i], scale = FALSE), sB), G), 1, proba.train[, i], "*")
			})
			z <- Reduce("+", zclass)
			fadta <- data.train
			fadta$x <- data.train$x - tcrossprod(z, fa$B)
			fameanclass <- sapply(1:nbclass, function(i) {
				colMeans(fadta$x[data.train$y == i, ])
			})
			coefpx <- LassoML(fadta,...)
			proba.train <- coefpx$proba.train
			v = colVars(cdta)
			eps <- colSums((fameanclass/sqrt(v) - meanclass/sqrt(v))^2)
			eps <- (mean(eps))
			iter <- iter + 1
			if (verbose) 
				print(eps)
			meanclass <- fameanclass
			cdta <- data.train$x
			for (i in 1:nbclass) {
				cdta[data.train$y == i, ] <- sweep(data.train$x[data.train$y == i, 
					], 2, meanclass[, i], "-")
			}
		}
		fa <- emfa(cdta, nbf = nbf, minerr = 1e-06, EM = EM)
		sB <- scale(t(fa$B), center = FALSE, scale = sqrt(fa$Psi))
		G <- solve(diag(ncol(fa$B)) + tcrossprod(sB))
		sB <- scale(t(fa$B), center = FALSE, scale = fa$Psi)
		zclass <- lapply(1:nbclass, function(i) {
			tcrossprod(tcrossprod(scale(data.train$x, center = meanclass[, i], 
				scale = FALSE), sB), G)
		})
		zclass <- lapply(1:nbclass, function(i) {
			sweep(zclass[[i]], 1, proba.train[, i], "*")
		})
		z <- Reduce("+", zclass)
		fadta <- data.train
		fadta$x <- data.train$x - tcrossprod(z, fa$B)
		l <- list(meanclass = meanclass, fa.training = fadta$x, Psi = fa$Psi, B = fa$B, 
			factors.training = z, groups = data.train$y, proba.training = proba.train, 
			mod.decorrelate.test = mod, data.train = data.train$x, 
			maxnbfactors = maxnbfactors, min.err = min.err, EM = EM, maxiter = maxiter)
		return(l)
	}
}
