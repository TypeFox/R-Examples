local.rbfdot <- function(data, n=5) {
	data <- as.matrix(data)
	nelts <- dim(data)[1]
	d <- dim(data)[2]
	## first, standardize data
	means <- apply(data[,1:d], 2, mean)
	sds <- apply(data[,1:d], 2, sd)
	sds[sds==0] <- 1 # manage NULL case
	data[,1:d] <- sweep(data[,1:d], 2, means, "-")
	data[,1:d] <- sweep(data[,1:d], 2, sds, "/")

	kern <- as.matrix(dist(data[,1:d]))

	# according to specc "local" implementation, use a median
	sigmas <- sapply(1:nelts, function(i) {
		dummy <- sort(kern[i,])
		return(median(dummy[1:n]))
	})
	# exterior product
	sigmamat <- as.matrix(sigmas) %*% t(as.matrix(sigmas))
	
	kern <- exp(- (kern^2) / sigmamat)
	return(kern)
}


speccalt <- function(kern, k=NA, maxk=20) {
	nelts <- dim(kern)[1]
	diag(kern) <- 0 # "zeroed" according to Zelnik-Manor and Perona 2004
	deg <- sapply(1:nelts, function(i) { # degree matrix
		return(sum(kern[i,]))
	})
  
	# hack : one "zero-based" (see "Lrw" in von Luxburg (2006)) for Bartlett procedure, 
	# the other for the actual clustering.
	L <- diag(1/sqrt(deg)) %*% kern %*% diag(1/sqrt(deg))
	L2 <- diag(nelts) - diag(1/deg) %*% kern
	eig <- eigen(L)
	eig2 <- eigen(L2)
	
	if(is.na(k)) {
		k <- bartlett(eig2$values, maxk=maxk)$k # automatic selection heuristic
	}
	
	Y <- eig$vectors[,1:k]
	sumsY <- rowSums(Y^2)
	Y <- t(sapply(1:nelts, function(i) Y[i,] / sqrt(sumsY[i])))
	
	return(kmeans(Y, k, 200)$cluster)	
	
}


bartlett <- function(eigvals, thres=0.95, maxk=20) {	 # adaptation of Bartlett statistical test for equal variances
	n <- length(eigvals) - 1
	m <- length(eigvals)
	# thres eigenvals to 10^{-12}
	eigvals[eigvals<10^{-12}] <- 10^{-12}
	stats <- numeric()
	tests <- numeric()
	# as we exclude the last eigenval, we start the tests at m=2 - which naively always has 0 tstat.
	for(p in 2:maxk) {
		q <- m-p
		nu <- (p-1) * (p+2) / 2
		# approx : with such laplacians, eigenvals up to rank k are close to 0,
		# so we approximate these by 1.
		fact <- n-q - (2 * p^2 + 2) / (6 * p) + q * (mean(eigvals[(q+1):n])^2) / ((1 - mean(eigvals[(q+1):n]))^2)
		Vk <- prod((p-1) * eigvals[(q+1):n] / sum(eigvals[(q+1):n]))
		tstat <- -fact * log(Vk)
		stats <- c(stats, tstat)
		tests <- c(tests, pchisq(tstat, nu))
	}
	
	# if thres is exceeded -> return the preceding index
	# if not exceeded, return the index preceding the max of the stat.
	if(any(tests > thres)) {
		thresind <- which(tests>thres)[1]
	} else {
		thresind <- which.max(tests)
	}
	
	return(list(k=thresind, stats=stats, tests=tests))
}




