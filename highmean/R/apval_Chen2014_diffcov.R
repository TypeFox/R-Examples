apval_Chen2014_diffcov <- function(sam1, sam2){
	n1 <- dim(sam1)[1]
	n2 <- dim(sam2)[1]
	p <- dim(sam1)[2]
	sam.cov1 <- cov(sam1)
	sam.cov2 <- cov(sam2)
	test.stat <- max_threshold_diffcov(sam1, sam2, sam.cov1, sam.cov2)
	test.stat <- as.numeric(test.stat)
	pval <- 1 - exp(-exp(-test.stat))
	names(pval) <- "Chen2014"
	out <- NULL
	out$sam.info <- c("n1" = n1, "n2" = n2, "p" = p)
	out$cov.assumption <- "the two groups have different covariances"
	out$method <- "asymptotic distribution"
	out$pval <- pval
	return(out)
}