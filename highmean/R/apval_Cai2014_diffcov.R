apval_Cai2014_diffcov <- function(sam1, sam2){
	n1 <- dim(sam1)[1]
	n2 <- dim(sam2)[1]
	p <- dim(sam1)[2]
	b <- 2*log(p) - log(log(p))
	diff <- colMeans(sam1) - colMeans(sam2)
	sam.cov1 <- cov(sam1)
	sam.cov2 <- cov(sam2)
	sam.diag1 <- diag(sam.cov1)
	sam.diag1[sam.diag1 <= 10^(-10)] <- 10^(-10)
	sam.diag2 <- diag(sam.cov2)
	sam.diag2[sam.diag2 <= 10^(-10)] <- 10^(-10)
	test.stat <- max(diff^2/(sam.diag1/n1 + sam.diag2/n2))
	test.stat <- as.numeric(test.stat)
	stan.stat <- test.stat - b
	pval <- 1 - exp(-exp(-stan.stat/2)/sqrt(3.14159))
	out <- NULL
	out$sam.info <- c("n1" = n1, "n2" = n2, "p" = p)
	out$cov.assumption <- "the two groups have different covariances"
	out$method <- "asymptotic distribution"
	out$pval <- pval
	return(out)
}