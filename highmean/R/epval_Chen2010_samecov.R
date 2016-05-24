epval_Chen2010_samecov <- function(sam1, sam2, perm.iter = 1000, seeds){
	n1 <- dim(sam1)[1]
	n2 <- dim(sam2)[1]
	tau <- (n1 + n2)/(n1*n2)
	n <- n1 + n2 - 2
	p <- dim(sam1)[2]
	sam <- rbind(sam1, sam2)
	diff <- colMeans(sam1) - colMeans(sam2)
	col.sum1 <- colSums(sam1)
	col.sum2 <- colSums(sam2)
	rm1 <- sum(col.sum1^2) - sum(diag(sam1 %*% t(sam1)))
	rm2 <- sum(col.sum2^2) - sum(diag(sam2 %*% t(sam2)))
	chen2010.stat <- rm1/(n1*(n1 - 1)) + rm2/(n2*(n2 - 1)) - 2*sum(col.sum1*col.sum2)/(n1*n2)
	chen2010.stat <- as.numeric(chen2010.stat)

	chen2010.stat.perm <- numeric(perm.iter)
	for(i in 1:perm.iter){
		if(!is.null(seeds)) set.seed(seeds[i])
		perm <- sample(1:(n1 + n2))
		sam.perm <- sam[perm,]
		sam1.perm <- sam.perm[1:n1,]
		sam2.perm <- sam.perm[(n1 + 1):(n1 + n2),]
		diff.perm <- colMeans(sam1.perm) - colMeans(sam2.perm)
		col.sum1.perm <- colSums(sam1.perm)
		col.sum2.perm <- colSums(sam2.perm)
		rm1.perm <- sum(col.sum1.perm^2) - sum(diag(sam1.perm%*%t(sam1.perm)))
		rm2.perm <- sum(col.sum2.perm^2) - sum(diag(sam2.perm%*%t(sam2.perm)))
		chen2010.stat.perm[i] <- rm1.perm/(n1*(n1 - 1)) + rm2.perm/(n2*(n2 - 1)) - 2*sum(col.sum1.perm*col.sum2.perm)/(n1*n2)
		chen2010.stat.perm[i] <- as.numeric(chen2010.stat.perm[i])
	}
	chen2010.pval <- (sum(chen2010.stat.perm >= chen2010.stat) + 1)/(perm.iter + 1)
	names(chen2010.pval) <- "Chen2010"
	out <- NULL
	out$sam.info <- c("n1" = n1, "n2" = n2, "p" = p)
	out$cov.assumption <- "the two groups have same covariance"
	out$method <- "permutation"
	out$pval <- chen2010.pval
	return(out)
}