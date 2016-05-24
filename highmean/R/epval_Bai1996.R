epval_Bai1996 <- function(sam1, sam2, n.iter = 1000, seeds){
	if(missing(seeds)) seeds <- NULL
	if(length(seeds) != n.iter){
		seeds <- NULL
		cat("The length of seeds does not match the specified n.iter.\n")
		cat("Seeds for each permutation/resampling iteration are assigned randomly.\n")
	}
	n1 <- dim(sam1)[1]
	n2 <- dim(sam2)[1]
	tau <- (n1 + n2)/(n1*n2)
	n <- n1 + n2 - 2
	p <- dim(sam1)[2]
	sam <- rbind(sam1, sam2)
	diff <- colMeans(sam1) - colMeans(sam2)
	XX <- t(diff) %*% diff
	bai1996.stat <- as.numeric(XX)

	bai1996.stat.perm <- numeric(n.iter)
	for(i in 1:n.iter){
		if(!is.null(seeds)) set.seed(seeds[i])
		perm <- sample(1:(n1 + n2))
		sam.perm <- sam[perm,]
		sam1.perm <- sam.perm[1:n1,]
		sam2.perm <- sam.perm[(n1 + 1):(n1 + n2),]
		diff.perm <- colMeans(sam1.perm) - colMeans(sam2.perm)
		XX.perm <- t(diff.perm)%*%diff.perm
		bai1996.stat.perm[i] <- as.numeric(XX.perm)
	}
	bai1996.pval <- (sum(bai1996.stat.perm >= bai1996.stat) + 1)/(n.iter + 1)
	names(bai1996.pval) <- "Bai1996"
	out <- NULL
	out$sam.info <- c("n1" = n1, "n2" = n2, "p" = p)
	out$cov.assumption <- "the two groups have same covariance"
	out$method <- "permutation"
	out$pval <- bai1996.pval
	return(out)
}