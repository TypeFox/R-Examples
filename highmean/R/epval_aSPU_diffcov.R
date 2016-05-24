epval_aSPU_diffcov <- function(sam1, sam2, pow = c(1:6, Inf), n.resam = 1000, sam.cov1, sam.cov2, cov1.est, cov2.est, cv.fold = 5, norm = "F", seeds, optim.bandwidth1, optim.bandwidth2, output.opt.bw1, output.opt.bw2){
	n1 <- dim(sam1)[1]
	n2 <- dim(sam2)[1]
	p <- dim(sam1)[2]

	if(!is.element(Inf, pow)) pow <- c(pow, Inf)
	pow <- sort(pow)

	T_SPU_0 <- stat_SPU(sam1, sam2, pow)

	T_SPU_resam <- matrix(NA, n.resam, length(pow))

	for(b in 1:n.resam){
		if(!is.null(seeds)) set.seed(seeds[b])
		sam1.b <- rmnorm(n = n1, mean = rep(0, p), cov1.est)
		sam2.b <- rmnorm(n = n2, mean = rep(0, p), cov2.est)
		T_SPU_resam[b,] <- stat_SPU(sam1.b, sam2.b, pow)
	}

	p.SPU <- rep(NA, length(pow))
	for(i in 1:length(pow)){
		p.SPU[i] <- (sum(abs(T_SPU_resam[,i]) > abs(T_SPU_0[i])) + 1)/(n.resam + 1)
		p.SPU.resam <- (n.resam + 1 - rank(abs(T_SPU_resam[,i])))/n.resam
		if(i == 1){
			T.aSPU.resam <- p.SPU.resam
		}else{
			T.aSPU.resam[which(T.aSPU.resam > p.SPU.resam)] <- p.SPU.resam[which(T.aSPU.resam > p.SPU.resam)]
		}
	}

	T_aSPU_0 <- min(p.SPU)
	p.aSPU <- (sum(T.aSPU.resam < T_aSPU_0) + 1)/(n.resam + 1)

	pvals <- c(p.SPU, p.aSPU)
	names(pvals) <- c(paste("SPU", as.character(pow), sep = "_"), "aSPU")
	out <- NULL
	out$sam.info <- c("n1" = n1, "n2" = n2, "p" = p)
	out$pow <- pow
	if(output.opt.bw1) out$opt.bw1 <- optim.bandwidth1
	if(output.opt.bw2) out$opt.bw2 <- optim.bandwidth2
	out$cov.assumption <- "the two groups have different covariances"
	out$method <- "parametric bootstrap"
	out$pval <- pvals
	return(out)
}