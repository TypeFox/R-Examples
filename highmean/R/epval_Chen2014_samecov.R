epval_Chen2014_samecov <- function(sam1, sam2, perm.iter = 1000, seeds){
	n1 <- dim(sam1)[1]
	n2 <- dim(sam2)[1]
	p <- dim(sam1)[2]
	Sn <- ((n1 - 1)*cov(sam1) + (n2 - 1)*cov(sam2))/(n1 + n2 - 2)
	sam <- rbind(sam1, sam2)
	test.stat <- max_threshold(sam1, sam2, Sn)
	test.stat <- as.numeric(test.stat)

	test.stat.perm <- numeric(perm.iter)
	for(i in 1:perm.iter){
		if(!is.null(seeds)) set.seed(seeds[i])
		perm <- sample(1:(n1 + n2))
		sam.perm <- sam[perm,]
		sam1.perm <- sam.perm[1:n1,]
		sam2.perm <- sam.perm[(n1 + 1):(n1 + n2),]
		Sn.perm <- ((n1 - 1)*cov(sam1.perm) + (n2 - 1)*cov(sam2.perm))/(n1 + n2 - 2)
		test.stat.perm[i] <- max_threshold(sam1.perm, sam2.perm, Sn.perm)
		test.stat.perm[i] <- as.numeric(test.stat.perm[i])
	}
	pval <- (sum(test.stat.perm >= test.stat) + 1)/(perm.iter + 1)
	names(pval) <- "Chen2014"
	out <- NULL
	out$sam.info <- c("n1" = n1, "n2" = n2, "p" = p)
	out$cov.assumption <- "the two groups have same covariance"
	out$method <- "permutation"
	out$pval <- pval
	return(out)
}