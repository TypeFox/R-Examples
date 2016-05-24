apval_Cai2014_samecov <- function(sam1, sam2){
	n1 <- dim(sam1)[1]
	n2 <- dim(sam2)[1]
	tau <- (n1*n2)/(n1 + n2)
	p <- dim(sam1)[2]
	b <- 2*log(p) - log(log(p))
	diff <- colMeans(sam1) - colMeans(sam2)
	sam.cov <- ((n1 - 1)*cov(sam1) + (n2 - 1)*cov(sam2))/(n1 + n2 - 2)
	sam.diag <- diag(sam.cov)
	sam.diag[sam.diag <= 10^(-10)] <- 10^(-10)
	test.stat <- tau*max(diff^2/sam.diag)
	test.stat <- as.numeric(test.stat)
	stan.stat <- test.stat - b
	pval <- 1 - exp(-exp(-stan.stat/2)/sqrt(3.14159))
	names(pval) <- "Cai2014"
	out <- NULL
	out$sam.info <- c("n1" = n1, "n2" = n2, "p" = p)
	out$cov.assumption <- "the two groups have same covariance"
	out$method <- "asymptotic distribution"
	out$pval <- pval
	return(out)
}