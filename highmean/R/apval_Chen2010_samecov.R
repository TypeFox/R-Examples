apval_Chen2010_samecov <- function(sam1, sam2){
	n1 <- dim(sam1)[1]
	n2 <- dim(sam2)[1]
	p <- dim(sam1)[2]
	n <- n1 + n2 - 2
	sam.cov <- ((n1 - 1)*cov(sam1) + (n2 - 1)*cov(sam2))/n
	trS <- sum(diag(sam.cov))
	tr.cov2 <- n^2/((n + 2)*(n - 1))*(sum(sam.cov^2) - trS^2/n)
	T1 <- sam1 %*% t(sam1)
	T2 <- sam2 %*% t(sam2)
	P1 <- (sum(T1) - sum(diag(T1)))/(n1*(n1 - 1))
	P2 <- (sum(T2) - sum(diag(T2)))/(n2*(n2 - 1))
	P3 <- -2*sum(sam1 %*% t(sam2))/(n1*n2)
	T <- P1 + P2 + P3
	test.stat <- T/sqrt((2/(n1*(n1 - 1)) + 2/(n2*(n2 - 1)) + 4/(n1*n2))*tr.cov2)
	test.stat <- as.numeric(test.stat)
	pval <- 1 - pnorm(test.stat)
	names(pval) <- "Chen2010"
	out <- NULL
	out$sam.info <- c("n1" = n1, "n2" = n2, "p" = p)
	out$cov.assumption <- "the two groups have same covariance"
	out$method <- "asymptotic distribution"
	out$pval <- pval
	return(out)
}