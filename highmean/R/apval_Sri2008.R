apval_Sri2008 <- function(sam1, sam2){
	n1 <- dim(sam1)[1]
	n2 <- dim(sam2)[1]
	tau <- (n1*n2)/(n1 + n2)
	n <- n1 + n2 - 2
	p <- dim(sam1)[2]
	sam.cov <- ((n1 - 1)*cov(sam1) + (n2 - 1)*cov(sam2))/n
	sam.corr <- cov2cor(sam.cov)
	trR2 <- sum(sam.corr^2)
	c <- 1 + trR2/p^(3/2)
	diff <- colMeans(sam1) - colMeans(sam2)
	diag.sam <- diag(sam.cov)
	diag.sam[diag.sam <= 10^(-10)] <- 10^(-10)
	XX <- sum(diff^2/diag.sam)
	test.stat <- (tau*XX - n*p/(n - 2))/sqrt(2*(trR2 - p^2/n)*c)
	test.stat <- as.numeric(test.stat)
	pval <- 1 - pnorm(test.stat)
	names(pval) <- "Sri2008"
	out <- NULL
	out$sam.info <- c("n1" = n1, "n2" = n2, "p" = p)
	out$cov.assumption <- "the two groups have same covariance"
	out$method <- "asymptotic distribution"
	out$pval <- pval
	return(out)
}