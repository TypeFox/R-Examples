apval_Bai1996 <- function(sam1, sam2){
	n1 <- dim(sam1)[1]
	n2 <- dim(sam2)[1]
	tau <- (n1*n2)/(n1 + n2)
	n <- n1 + n2 - 2
	p <- dim(sam1)[2]
	diff <- colMeans(sam1) - colMeans(sam2)
	XX <- sum(diff^2)
	sam.cov <- ((n1 - 1)*cov(sam1) + (n2 - 1)*cov(sam2))/n
	trS <- sum(diag(sam.cov))
	tr.cov2 <- n^2/((n + 2)*(n - 1))*(sum(sam.cov^2) - trS^2/n)
	test.stat <- (tau*XX - trS)/sqrt(2*(n + 1)/n*tr.cov2)
	test.stat <- as.numeric(test.stat)
	pval <- 1 - pnorm(test.stat)
	names(pval) <- "Bai1996"
	out <- NULL
	out$sam.info <- c("n1" = n1, "n2" = n2, "p" = p)
	out$cov.assumption <- "the two groups have same covariance"
	out$method <- "asymptotic distribution"
	out$pval <- pval
	return(out)
}