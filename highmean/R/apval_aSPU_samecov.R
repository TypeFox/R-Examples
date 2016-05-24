apval_aSPU_samecov <- function(sam1, sam2, pow = c(1:6, Inf), sam.cov, cov.est, cv.fold = 5, norm = "F", optim.bandwidth, output.opt.bw){
	n1 <- dim(sam1)[1]
	n2 <- dim(sam2)[1]
	tau <- (n1*n2)/(n1 + n2)
	p <- dim(sam1)[2]
	if(cv.fold != floor(cv.fold) | cv.fold <= 1) stop("cv.fold must be an integer greater than 1.")

	if(!is.element(Inf, pow)) pow <- c(pow, Inf)
	pow <- sort(pow)

	diff <- colMeans(sam1) - colMeans(sam2)
	pval <- numeric(length(pow) + 1)
	L <- numeric(length(pow))
	L.e <- numeric(length(pow) - 1)
	L.var <- numeric(length(pow) - 1)
	stan.L <- numeric(length(pow) - 1)
	for(i in 1:length(pow)){
		if(pow[i] != Inf){
			L[i] <- sum(diff^(pow[i]))
			L.e[i] <- SPU_E(pow[i], n1, n2, cov.est)
			L.var[i] <- SPU_Var(pow[i], n1, n2, cov.est)
			stan.L[i] <- (L[i] - L.e[i])/sqrt(L.var[i])
			if(pow[i] %% 2 == 1) pval[i] <- 2*(1 - pnorm(abs(stan.L[i])))
			if(pow[i] %% 2 == 0) pval[i] <- 1 - pnorm(stan.L[i])
		}else{
			diag.sam.cov <- diag(sam.cov)
			diag.sam.cov[diag.sam.cov <= 10^(-10)] <- 10^(-10)
			L[i] <- max(diff^2/diag.sam.cov)
			stan.L.inf <- tau*L[i] - (2*log(p) - log(log(p)))
			pval[i] <- pval.inf <- 1 - exp(-exp(-stan.L.inf/2)/sqrt(3.14159))
		}
	}
	f.pow <- pow[pow != Inf]
	odd.ga <- f.pow[f.pow %% 2 == 1]
	odd.ga.id <- which(f.pow %% 2 == 1)
	even.ga <- f.pow[f.pow %% 2 == 0]
	even.ga.id <- which(f.pow %% 2 == 0)
	n.odd.ga <- length(odd.ga)
	n.even.ga <- length(even.ga)
	R_O <- matrix(NA, n.odd.ga, n.odd.ga)
	R_E <- matrix(NA, n.even.ga, n.even.ga)
	diag(R_O) <- 1
	diag(R_E) <- 1
	for(s in odd.ga){
		for(t in odd.ga){
			if(s != t){
				L.cov <- SPU_Cov(s, t, n1, n2, cov.est)
				io <- which(odd.ga == s)
				jo <- which(odd.ga == t)
				i <- which(pow == s)
				j <- which(pow == t)
				R_O[io, jo] <- L.cov/sqrt(L.var[i]*L.var[j])
			}
		}
	}
	for(s in even.ga){
		for(t in even.ga){
			if(s != t){
				L.cov <- SPU_Cov(s, t, n1, n2, cov.est)
				ie <- which(even.ga == s)
				je <- which(even.ga == t)
				i <- which(pow == s)
				j <- which(pow == t)
				R_E[ie, je] <- L.cov/sqrt(L.var[i]*L.var[j])
			}
		}
	}
	TO <- max(abs(stan.L[odd.ga.id]))
	TE <- max(stan.L[even.ga.id])
	pval_O <- 1 - pmvnorm(lower = -rep(TO, n.odd.ga), upper = rep(TO, n.odd.ga),
		mean = rep(0, n.odd.ga), sigma = R_O)
	pval_E <- 1 - pmvnorm(lower = rep(-Inf, n.even.ga), upper = rep(TE, n.even.ga),
		mean = rep(0, n.even.ga), sigma = R_E)
	pval.min <- min(c(pval_O, pval_E, pval.inf))
	pval[length(pow) + 1] <- 1 - (1 - pval.min)^3
	names(pval) <- c(paste("SPU_", pow, sep = ""), "aSPU")

	out <- NULL
	out$sam.info <- c("n1" = n1, "n2" = n2, "p" = p)
	out$pow <- pow
	if(output.opt.bw) out$opt.bw <- optim.bandwidth
	names(L) <- pow
	out$spu.stat <- L
	names(L.e) <- names(L.var) <- pow[pow < Inf]
	out$spu.e <- L.e
	out$spu.var <- L.var
	colnames(R_O) <- rownames(R_O) <- odd.ga
	out$spu.corr.odd <- R_O
	colnames(R_E) <- rownames(R_E) <- even.ga
	out$spu.corr.even <- R_E
	out$cov.assumption <- "the two groups have same covariance"
	out$method <- "asymptotic distribution"
	out$pval <- pval
	return(out)
}