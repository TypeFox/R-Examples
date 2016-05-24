
PoisMixClus <- function(y, g, conds, lib.size = TRUE, lib.type = "TMM",  
	init.type = "small-em", init.runs = 1, init.iter = 10, alg.type = "EM", cutoff = 10e-6, 
	iter = 1000, fixed.lambda = NA, equal.proportions = FALSE, prev.labels = NA, 
	prev.probaPost = NA, verbose = FALSE, interpretation = "sum", EM.verbose = FALSE, wrapper=FALSE, s=NA) {

	## fixed.lambda should be a list of length (number of fixed clusters)
	## g gives the number of clusters IN ADDITION to the fixed clusters
	## 	specified by fixed.lambda
	## equal.proportions should be TRUE or FALSE

	## Next loop only run if PoisMixClus is called directly (otherwise called by PoisMixClusWrapper)
	if(wrapper==FALSE) {
		if(is.matrix(y) == FALSE & is.data.frame(y) == FALSE) 
			stop(paste(sQuote("y"), "must be a matrix"))
		if(min(y) < 0 | sum(round(y)) != sum(y)) 
			stop(paste(sQuote("y"), "must be a matrix made up of nonnegative counts"))
		if(min(rowSums(y)) == 0)
			stop(paste("at least one observation in", sQuote("y"), 
			"contains all 0's and must be removed from the data"))
		if(is.vector(conds) == FALSE | length(conds) != ncol(y))
			stop(paste(sQuote("conds"), "must be a vector the same length as the number of columns in", sQuote("y")))
		if(is.logical(lib.size) == FALSE)
			stop(paste(sQuote("libsize"), "must be", dQuote("TRUE"), "(PMM-II) or", 
				dQuote("FALSE"), "(PMM-I)"))
		if(length(init.type) > 1)
			stop(paste(sQuote("init.type"), "must be of length 1"))
		if(init.type != "small-em" & init.type != "kmeans" & init.type != "split.small-em") 
			stop(paste(sQuote("init.type"), "must be one of", dQuote("small-em"), "or", 
				dQuote("kmeans"), "or", dQuote("split.small-em")))
		if(alg.type != "EM" & alg.type != "CEM")
			stop(paste(sQuote("alg.type"), "must be one of", dQuote("EM"), "or", dQuote("CEM")))
		if(length(alg.type) > 1)
			stop(paste(sQuote("alg.type"), "must be one of", dQuote("EM"), "or", dQuote("CEM")))
		if(is.logical(verbose) == FALSE)
			stop(paste(sQuote("verbose"), "must be", dQuote("TRUE"), "or", dQuote("FALSE")))
		if(class(fixed.lambda) != "list" & is.na(fixed.lambda[1]) == FALSE)
			stop(paste(sQuote("fixed.lambda"), "must be", dQuote("NA") , "or a list."))
		if(is.vector(prev.labels) == FALSE & is.na(prev.labels[1]) == FALSE)
			stop(paste(sQuote("prev.labels"), "must be", dQuote("NA") , "or a vector of labels."))
		if(is.na(s[1]) == FALSE & length(s) != length(conds) & sum(s) != 1)
			stop(paste(sQuote("s"), "must be", dQuote("NA") , "or a vector of normalized library sizes summing to 1."))

		## Grouping columns of y in order of condition (all replicates put together)
		o.ycols <- order(conds)
		y <- y[,o.ycols]
		conds <- conds[o.ycols]
		conds.names <- unique(conds)
		d <- length(unique(conds))
		r <- as.vector(table(conds))
		if(length(rownames(y)) == 0) rn <- 1:nrow(y);
		if(length(rownames(y)) > 0) rn <- rownames(y);
		y <- as.matrix(y, nrow = nrow(y), ncol = ncol(y))
		rownames(y) <- rn;
		n <- dim(y)[1];cols <- dim(y)[2]
		w <- rowSums(y)
		## Only calculate s values if they are not provided
		if(is.na(s[1]) == TRUE) {
			s <- rep(NA, cols)
			if(lib.size == FALSE) {
				s <- rep(1, cols)
			}
			if(lib.size == TRUE) {
				if(lib.type == "TC") s <- colSums(y) / sum(y);
				if(lib.type == "UQ") s <- apply(y, 2, quantile, 0.75) / sum(apply(y, 2, quantile, 0.75));
				if(lib.type == "Med") s <- apply(y, 2, median) / sum(apply(y, 2, median));
				if(lib.type == "DESeq") {
					## Code from DESeq, v1.8.3
					loggeomeans <- rowMeans(log(y))
					s <- apply(y, 2, function(x) 
						exp(median((log(x)-loggeomeans)[is.finite(loggeomeans)])))
					s <- s / sum(s)
				}
				if(lib.type == "TMM") {
					f <- calcNormFactors(as.matrix(y), method = "TMM")
					s <- colSums(y)*f / sum(colSums(y)*f)
				} 
			}
		}
		s.dot <- rep(NA, d) 
		for(j in 1:d) {
			s.dot[j] <- sum(s[which(conds == (unique(conds))[j])])
		}
	}
	
	if(wrapper==TRUE) {
		conds.names <- unique(conds)
		d <- length(unique(conds))
		r <- as.vector(table(conds))
		n <- dim(y)[1];cols <- dim(y)[2]
		w <- rowSums(y)
		K <- g
		s.dot <- rep(NA, d) 
		for(j in 1:d) {
			s.dot[j] <- sum(s[which(conds == (unique(conds))[j])])
		}
	}

	K <- g
	if(class(fixed.lambda) == "list") {
		for(ll in 1:length(fixed.lambda)) {
			if(is.vector(fixed.lambda[[ll]]) == FALSE |
				length(fixed.lambda[[ll]]) != d)
				stop(paste(sQuote("fixed.lambda"), "must be", dQuote("NA") , 
					"or a list of length equal to the number of conditions."))
			if(length(which(fixed.lambda[[ll]] == 0)) > 0) {
				if(length(which(fixed.lambda[[ll]] == 1)) + 
					length(which(fixed.lambda[[ll]] == 0)) == length(fixed.lambda[[ll]])) {

					tmp <- 1/sum(s.dot[which(fixed.lambda[[ll]] == 1)])
					new <- fixed.lambda[[ll]]
					new[which(fixed.lambda[[ll]] == 1)] <- tmp
					message(cat("Fixed lambda\n", fixed.lambda[[ll]], "\n", "modified to\n",
						new, "\n", "to satisfy imposed parameter constraints.\n"))
					fixed.lambda[[ll]][which(fixed.lambda[[ll]] == 1)] <- tmp
				}
			}
			if(abs(as.numeric(fixed.lambda[[ll]] %*% s.dot)-1) > 10e-8)
				warning(paste("Check that constraint on lambda*s.dot is upheld for",sQuote("fixed.lambda")))
		}
		K <- g + length(fixed.lambda);	
	}
	diff <- 100 ## Convergence criterion
	index <- 0; go <- 1;

	## Inital values
	## init.type: "kmeans", "small-em", "split.small-em"

	if(init.type == "kmeans") {
		init.alg <- "kmeanInit";
		init.args <- list(y = y, g = K, conds = conds, lib.size = lib.size, 
			lib.type = lib.type, fixed.lambda = fixed.lambda,
			equal.proportions = equal.proportions, s = s)
	}
	if(init.type == "small-em") {
		init.alg <- "emInit"
		init.args <- list(y = y, g = g, conds = conds, lib.size = lib.size, 
			lib.type = lib.type, alg.type = alg.type, init.run = init.runs,
			init.iter = init.iter, fixed.lambda = fixed.lambda, 
			equal.proportions = equal.proportions, verbose = verbose, s = s)
	}
	if(init.type == "split.small-em") {
		init.alg <- "splitEMInit"
		init.args <- list(y = y, g = g, conds = conds, lib.size = lib.size,
			lib.type = lib.type, alg.type = alg.type,
			fixed.lambda = fixed.lambda,
			equal.proportions = equal.proportions, 
			prev.labels = prev.labels, prev.probaPost = prev.probaPost,
			init.iter = init.iter, init.runs = init.runs, verbose = verbose, s = s)
	}
	## Adding quote = TRUE to speed up do.call
	param.init <- do.call(init.alg, init.args, quote=TRUE)

	if(equal.proportions == FALSE) {
		pi <- pi.old <- param.init$pi.init
	}
	if(equal.proportions == TRUE) {
		pi <- pi.old <- rep(1/K, K)
	}
	lambda <- lambda.old <- param.init$lambda.init
	mean.calc <- mean.old <- PoisMixMean(y = y, g = K, conds = conds, 
		s = s, lambda = lambda)

	while(go == 1) {

		############
		## E-step ##
		############
		t <- probaPost(y, K, conds, pi, s, lambda)

		############
		## C-step ##
		############
		if(alg.type == "CEM") {
			## If two values of t_{ik} are map, 
			## arbitrarily choose the first
			partition <- unlist(apply(t, 1, 
				function(x) which(x == max(x, na.rm = TRUE))[1]))
			partition.mat <- matrix(0, nrow = n, ncol = K)
			for(i in 1:n) partition.mat[i,partition[i]] <- 1;
		}

		############
		## M-step ##
		############
		if(alg.type == "CEM") {
			if(equal.proportions == FALSE) {
				for(k in 1:K) {
					pi[k] <- length(which(partition == k))/n
				}
			}
			if(equal.proportions == TRUE) {
				pi <- rep(1/K, K)
			}
			denom <- colSums(partition.mat * w)
			if(class(fixed.lambda) != "list") {
				for(j in 1:d) {
					denom.bis <- denom * s.dot[j]
					num <- colSums(partition.mat * 
					rowSums(as.matrix(y[,which(conds == (unique(conds))[j])])))
					lambda[j,] <- num / denom.bis
				}
			}
			if(class(fixed.lambda) == "list") {
				for(ll in 1:length(fixed.lambda)) {
					lambda[,ll] <- fixed.lambda[[ll]]
				}
				for(j in 1:d) {
					denom.bis <- denom * s.dot[j]
					denom.bis <- denom.bis[-c(1:length(fixed.lambda))]
					num <- colSums(partition.mat * 
						rowSums(as.matrix(y[,which(conds == (unique(conds))[j])])))
					num <- num[-c(1:length(fixed.lambda))]
					lambda[j,-c(1:length(fixed.lambda))] <- num / denom.bis
				}
			}
		}

		if(alg.type == "EM") {
			if(equal.proportions == FALSE) {
				pi <- colSums(t)/n
			}
			if(equal.proportions == TRUE) {
				pi <- rep(1/K, K)
			}
			denom <- colSums(t * w)
			if(class(fixed.lambda) != "list") {
				denom.bis <- matrix(rep(denom, length(s.dot)) * rep(s.dot, each=K), byrow=T, ncol=K)
				num <- matrix(rowsum(matrix(y, ncol=n, nrow=length(conds), byrow=T), group=conds), 
					nrow=n, ncol=d, byrow=T)
				num <- matrix(unlist(lapply(lapply(1:d, function(x) t*num[,x]), colSums), recursive=FALSE,
					use.names=FALSE), nrow=d, ncol=K, byrow=T)
				lambda <- num/denom.bis
			}
			if(class(fixed.lambda) == "list") {
				for(ll in 1:length(fixed.lambda)) {
					lambda[,ll] <- fixed.lambda[[ll]]
				}
				## This loop could be improved for speed as above
				for(j in 1:d) {
					denom.bis <- denom * s.dot[j]
					denom.bis <- denom.bis[-c(1:length(fixed.lambda))]
					num <- colSums(t * 
						matrix(rep(rowSums(as.matrix(y[,which(conds == (unique(conds))[j])])),K), 
						ncol = K))
					num <- num[-c(1:length(fixed.lambda))]
					lambda[j,-c(1:length(fixed.lambda))] <- num / denom.bis
				}
			}	
		}

		#################
		## Convergence ##
		#################
		mean.calc <- PoisMixMean(y, g = K, conds, s, lambda)
		diff <- abs(logLikePoisMixDiff(y, mean.calc, pi, mean.old, pi.old))
		lambda.old <- lambda; pi.old <- pi; mean.old <- mean.calc;

		index <- index + 1
		if(verbose == TRUE) print(paste("Log-like diff:", diff))
		if(diff < cutoff) go <- 0;
		if(iter != FALSE & iter == index) go <- 0;
	}
	
	if(EM.verbose == TRUE) {
		cat("#####################################\n")
		cat("Number of EM iterations:", index, "\n")
		cat("Last log-likelihood difference:", diff, "\n")
		cat("#####################################\n")
	}

	#####################################
	## Final estimates of lambda and p ##
	#####################################
	names(pi) <- paste("Cluster", 1:K)
	colnames(lambda) <- paste("Cluster", 1:K)
	rownames(lambda) <- conds.names
	lambda.final <- lambda
	pi.final <- pi

	## Check to make sure one of the components is not degenerate
	if(min(pi) == 0 | is.nan(sum(lambda)) == TRUE) {
		probaPost <- NA
		labels <- NA
		BIC <- NA
		ICL <- NA
	}

	if(min(pi) > 0 | is.nan(sum(lambda)) == FALSE) {

		mean.calc <- PoisMixMean(y, g = K, conds, s, lambda)
		LL.tmp <- mylogLikePoisMix(y, mean.calc, pi)
		LL <- LL.tmp$ll
	
		######################
		## Determine labels ##
		######################
		t <- probaPost(y, K, conds, pi, s, lambda)
		## If two clusters have exactly identical map estimators,
		## arbitrarily choose the first one
		map <- unlist(apply(t, 1, function(x) which(x == max(x, 
			na.rm = TRUE))[1]))
		z <- matrix(0, nrow = n, ncol = K)
		for(i in 1:n) z[i,map[i]] <- 1;
		probaPost <- t
		labels <- map

		##############################
		## Calculate BIC, ICL       ##
		##############################
		if(equal.proportions == FALSE & class(fixed.lambda) != "list") {
#			np <- (K-1) + n + (d-1)*K 	# pi + w + lambda
			## CHANGE September 25, 2013: w is not considered as a parameter
			np <- (K-1) + (d-1)*K 	# pi + lambda
		}
		if(equal.proportions == TRUE & class(fixed.lambda) != "list") {
#			np <- n + (d-1)*K 	# w + lambda
			## CHANGE September 25, 2013: w is not considered as a parameter
			np <- (d-1)*K 		# lambda
		}
		if(equal.proportions == FALSE & class(fixed.lambda) == "list") {
#			np <- (K-1) + n + (d-1)*(K-length(fixed.lambda))	# pi + w + lambda not fixed
			## CHANGE September 25, 2013: w is not considered as a parameter
			np <- (K-1) + (d-1)*(K-length(fixed.lambda))		# pi + lambda not fixed
		}
		if(equal.proportions == TRUE & class(fixed.lambda) == "list") {
#			np <- n + (d-1)*(K-length(fixed.lambda))			# w + lambda not fixed
			## CHANGE September 25, 2013: w is not considered as a parameter
			np <- (d-1)*(K-length(fixed.lambda))				# lambda not fixed

		}
		BIC <- -LL + (np/2) * log(n)
#		entropy <- -2*sum(z*log(t), na.rm = TRUE)
		## CHANGE October 18, 2013: replace z with t in the entropy calculation for ICL
#		entropy <- -2*sum(t*log(t), na.rm = TRUE)
		## CHANGE July 25, 2014: typo in entropy (thanks, Melina Gallopin!)
		entropy <- -sum(t*log(t), na.rm = TRUE)
		ICL <- BIC + entropy
	}
	## Should cluster behavior (lambda) be interpretated wrt the gene means or sums?	
	if(interpretation == "mean") {
		s <- s * q
		w <- w / q
	}
	results <- list(lambda = lambda.final, pi = pi.final, labels = labels, 
		probaPost = probaPost, log.like = LL, BIC = -BIC, ICL = -ICL, 
		alg.type = alg.type, lib.size = lib.size, lib.type = lib.type, s = s,
		conds = conds, iterations = index, logLikeDiff = diff, model.selection = NA)

	class(results) <- "HTSCluster"
	return(results)
}
