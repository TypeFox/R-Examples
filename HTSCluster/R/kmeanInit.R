kmeanInit <- function(y, g, conds, lib.size, lib.type, fixed.lambda,
	equal.proportions, s=NA) {

	if(is.matrix(y) == FALSE & is.data.frame(y) == FALSE) 
		stop(paste(sQuote("y"), "must be a matrix"))
	if(min(y) < 0 | sum(round(y)) != sum(y)) 
		stop(paste(sQuote("y"), "must be a matrix made up of nonnegative counts"))
	if(min(rowSums(y)) == 0)
		stop(paste("at least one observation in", sQuote("y"), "contains all 0's and must be removed from the data"))
	if(length(g) != 1)
		stop(paste(sQuote("g"), "(the number of clusters) must be a nonnegative integer"))
	if(g < 0 | round(g) != g) 
		stop(paste(sQuote("g"), "(the number of clusters) must be a nonnegative integer"))
	if(is.vector(conds) == FALSE | length(conds) != ncol(y))
		stop(paste(sQuote("conds"), "must be a vector the same length as the number of columns in", sQuote("y")))
	if(is.logical(lib.size) == FALSE)
		stop(paste(sQuote("libsize"), "must be", dQuote("TRUE"), "(PMM-II) or", 
			dQuote("FALSE"), "(PMM-I)"))
	n <- dim(y)[1];cols <- dim(y)[2];
	y <- as.matrix(y, nrow = n, ncol = cols)
	d <- length(unique(conds))
	r <- as.vector(table(conds))
	w <- rowSums(y)
	if(is.na(s[1]) == TRUE) {
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

	## Use K-means to create initial partition
	g.init <- kmeans(y / w, g)
	partition <- g.init$cluster
	partition.mat <- matrix(0, nrow = n, ncol = g)
	for(i in 1:n) partition.mat[i,partition[i]] <- 1;

	## Calculate lambda.init and p.init
	denom <- colSums(partition.mat * w)
	lambda.init <- matrix(NA, nrow = d, ncol = g)
	for(j in 1:d) {
		denom.bis <- denom * s.dot[j]
		num <- colSums(partition.mat *
			rowSums(as.matrix(y[,which(conds == (unique(conds))[j])])))
		lambda.init[j,] <- num / denom.bis
	}
	if(class(fixed.lambda) == "list") {
		index <- lapply(fixed.lambda, function(x) {
			colSums((lambda.init - x)^2)})
		index.order <- lapply(index, order)
		fixed <- c()
		## Pick which components correspond to fixed components
		## Assign fixed lambdas to "closest" component (Euclidean distance),
		## starting with the first fixed lambda
		for(ll in 1:length(fixed.lambda)) {
			while(index.order[[ll]][1] %in% fixed) {
				index.order[[ll]] <- index.order[[ll]][-1]
			}
			lambda.init[,index.order[[ll]][1]] <- fixed.lambda[[ll]]	
			fixed <- c(fixed, index.order[[ll]][1])
		}
		## Rearrange lambda so that fixed values are in the first columns
		lambda.init <- cbind(lambda.init[,fixed], lambda.init[,-fixed])
	}
	if(equal.proportions == FALSE) {
		pi.init <- as.vector(table(g.init$cluster)/n)
		## Rearrange pi so that it corresponds to lambda
		if(is.list(fixed.lambda) == TRUE) {
			pi.init <- c(pi.init[fixed], pi.init[-fixed])
		}
	}
	if(equal.proportions == TRUE) {
		pi.init <- rep(1/g, g)
	}
	return(list(pi.init = pi.init, lambda.init = lambda.init))
}

