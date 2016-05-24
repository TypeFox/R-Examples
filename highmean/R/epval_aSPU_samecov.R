epval_aSPU_samecov <- function(sam1, sam2, pow = c(1:6, Inf), perm.iter = 1000, seeds){ 
	n1 <- dim(sam1)[1]
	n2 <- dim(sam2)[1]
	p <- dim(sam1)[2]
	pval <- aSPUperm(sam1, sam2, pow = pow, n.perm = perm.iter, seeds = seeds)
	out <- NULL
	out$sam.info <- c("n1" = n1, "n2" = n2, "p" = p)
	out$pow <- pow
	out$cov.assumption <- "the two groups have same covariance"
	out$method <- "permutation"
	out$pval <- pval
	return(out)
}