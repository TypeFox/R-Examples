epval_Cai2014_diffcov <- function(sam1, sam2, n.resam = 1000, sam.cov1, sam.cov2, cov1.est, cov2.est, cv.fold = 5, norm = "F", seeds, optim.bandwidth1, optim.bandwidth2, output.opt.bw1, output.opt.bw2){
	n1 <- dim(sam1)[1]
	n2 <- dim(sam2)[1]
	p <- dim(sam1)[2]

	T_CLX_0 <- stat_Cai2014(sam1, sam2)
	T_CLX_resam <- rep(NA, n.resam)
	for(b in 1:n.resam){
		if(!is.null(seeds)) set.seed(seeds[b])
		sam1.b <- rmnorm(n = n1, mean = rep(0, p), cov1.est)
		sam2.b <- rmnorm(n = n2, mean = rep(0, p), cov2.est)
		T_CLX_resam[b] <- stat_Cai2014(sam1.b, sam2.b)
	}

	p.CLX <- (sum(T_CLX_resam > T_CLX_0) + 1)/(n.resam + 1)

	pval <- p.CLX
	names(pval) <- "Cai2014"
	out <- NULL
	out$sam.info <- c("n1" = n1, "n2" = n2, "p" = p)
	if(output.opt.bw1) out$opt.bw1 <- optim.bandwidth1
	if(output.opt.bw2) out$opt.bw2 <- optim.bandwidth2
	out$cov.assumption <- "the two groups have different covariances"
	out$method <- "parametric bootstrap"
	out$pval <- pval
	return(out)
}