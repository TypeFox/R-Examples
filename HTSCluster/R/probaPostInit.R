
probaPostInit <- function(y, g, conds, lib.size, lib.type, 
	alg.type = "EM", fixed.lambda, equal.proportions,
	probaPost.init, init.iter, verbose, s=NA) 
{

	## fixed.lambda should be a list of length (number of fixed clusters)
	## g gives the number of clusters IN ADDITION to the fixed clusters
	## 	specified by fixed.lambda
	## equal.proportions should be TRUE or FALSE

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
	if(alg.type != "EM" & alg.type != "CEM")
		stop(paste(sQuote("alg.type"), "must be one of", dQuote("EM"), "or", dQuote("CEM")))
	if(length(alg.type) > 1)
		stop(paste(sQuote("alg.type"), "must be one of", dQuote("EM"), "or", dQuote("CEM")))
	if(is.logical(verbose) == FALSE)
		stop(paste(sQuote("verbose"), "must be", dQuote("TRUE"), "or", dQuote("FALSE")))
	if(class(fixed.lambda) != "list" & is.na(fixed.lambda[1]) == FALSE)
		stop(paste(sQuote("fixed.lambda"), "must be", dQuote("NA") , "or a list."))
	if(is.matrix(probaPost.init) == FALSE)
		stop(paste(sQuote("z.init"), "must be a matrix of posterior probabilities."))

#	Already done in PoisMixClus and PoisMixClusWrapper
#	## Grouping columns of y in order of condition (all replicates put together)
#	o.ycols <- order(conds)
#	y <- y[,o.ycols]
#	conds <- conds[o.ycols]
	conds.names <- unique(conds)

	d <- length(unique(conds))
	r <- as.vector(table(conds))
	diff <- 100 ## Convergence criterion
	if(length(rownames(y)) == 0) rn <- 1:nrow(y);
	if(length(rownames(y)) > 0) rn <- rownames(y);
	y <- as.matrix(y, nrow = nrow(y), ncol = ncol(y))
	rownames(y) <- rn;
	n <- dim(y)[1];cols <- dim(y)[2]

	w <- rowSums(y)
	if(is.na(s[1]) == FALSE) {
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

	K <- g
	if(class(fixed.lambda) == "list") {
		K <- g + length(fixed.lambda);	
	}
	
	## Inital values using probaPost.init
	pi <- pi.old <- rep(NA, K)
	lambda <- lambda.old <- matrix(NA, nrow = d, ncol = K)
	t <- probaPost.init

	for(index in 0:init.iter) {

		if(index > 0) {
			############
			## E-step ##
			############
			t <- probaPost(y, K, conds, pi, s, lambda)
		}

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
				for(j in 1:d) {
					denom.bis <- denom * s.dot[j]
					num <- colSums(t * 
						matrix(rep(rowSums(as.matrix(y[,which(conds == (unique(conds))[j])])),K), 
					ncol = K))
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
					num <- colSums(t * 
						matrix(rep(rowSums(as.matrix(y[,which(conds == (unique(conds))[j])])),K), 
						ncol = K))
					num <- num[-c(1:length(fixed.lambda))]
					lambda[j,-c(1:length(fixed.lambda))] <- num / denom.bis
				}
			}	
		}
		lambda.old <- lambda; pi.old <- pi;
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
		probaPost <- NA; labels <- NA; BIC <- NA;	ICL <- NA
	}

	if(min(pi) > 0 | is.nan(sum(lambda)) == FALSE) {
		mean.calc <- PoisMixMean(y, g = K, conds, s, lambda)
		LL.tmp <- mylogLikePoisMix(y, mean.calc, pi)
		LL <- LL.tmp$ll
	}
	results <- list(lambda = lambda.final, pi = pi.final, log.like = LL)
	return(results)
}



