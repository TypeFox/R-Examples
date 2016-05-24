	
splitEMInit <- function(y, g, conds, lib.size,
	lib.type, alg.type, fixed.lambda, equal.proportions, 
	prev.labels, prev.probaPost, init.runs, init.iter, verbose, s=NA) {

	## g is the new number of clusters IN ADDITION TO FIXED LAMBDA
	## NB: This function inspiried by init2.k() function from poisson.glm.mix package
	## (written by Panos Papastamoulis)

	n <- dim(y)[1]
	unique.labels <- unique(prev.labels)
	## Check whether any of the clusters has one or zero observations & remove from consideration
	tab <- table(prev.labels) 
	## Fix to make sure that we only consider clusters with at least 2 observations
	unique.labels <- names(tab)
	if(length(which(tab < 2)) > 0) {
		unique.labels <- names(tab)[-which(tab < 2)]
	}

	K <- g
	if(class(fixed.lambda) == "list") {
		K <- g + length(fixed.lambda);	
	}	
	prev.K <- K - 1
	
	## Calculate per-class entropy of previous model
	## and choose cluster with largest entropy to split
	perEntropy <- rep(NA, length(unique.labels))
	names(perEntropy) <- unique.labels
	for(k in 1:length(unique.labels)) {
		perEntropy[k] <- -sum(log(prev.probaPost[which(prev.labels == as.numeric(unique.labels[k])),as.numeric(unique.labels[k])]))
	}
	cluster.choose <- as.numeric(unique.labels[which(perEntropy == max(perEntropy))])
	index1 <- which(prev.labels == cluster.choose)
	
	## Random selection of observations within splitted component
	## Repeat this init.runs times
	LL.all <- rep(NA, init.runs)
	init.all <- vector("list", init.runs)
	for(iter in 1:init.runs) {
		u.numbers <- runif(length(index1))
		t <- matrix(0, nrow = n, ncol = K)
		t[,1:prev.K] <- prev.probaPost
		t[index1,g] <- t[index1,cluster.choose] * u.numbers
		t[index1,cluster.choose] <- t[index1,cluster.choose] * (1-u.numbers)

		## Smoothing t values
		epsilon <- 1e-10
		maxcut <- 1 - epsilon; mincut <- epsilon
		t <- apply(t, 2, pmax, mincut)
		t <- apply(t, 2, pmin, maxcut)
		t <- t/rowSums(t)

		## Initialize small EM-algorithm with new z values
		initialize <- probaPostInit(y = y, g = g, lib.size = lib.size,
			lib.type = lib.type, conds = conds, alg.type = alg.type,
			fixed.lambda = fixed.lambda,
			equal.proportions = equal.proportions, probaPost.init = t,
			init.iter = init.iter, verbose = verbose, s = s)
		LL.all[iter] <- initialize$log.like
		init.all[[iter]] <- initialize
	}
	
	init.select <- which(LL.all == max(LL.all, na.rm = TRUE))[1]
	final.init <- init.all[[init.select]]
	lambda.init <- final.init$lambda
	pi.init <- final.init$pi

	return(list(lambda.init = lambda.init, pi.init = pi.init))
}

