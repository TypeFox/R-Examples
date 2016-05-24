aSPUperm <- function(sam1, sam2, pow = c(1:6, Inf), n.perm = 1000, seeds){
	n1 <- dim(sam1)[1]
	n2 <- dim(sam2)[1]
	p <- dim(sam1)[2]
	Sn <- ((n1 - 1)*cov(sam1) + (n2 - 1)*cov(sam2))/(n1 + n2 - 2)
	sam <- rbind(sam1, sam2)
	diff <- colMeans(sam1) - colMeans(sam2)
	Ts <- rep(NA, length(pow))
	for(j in 1:length(pow)){
		if(pow[j] < Inf){
			Ts[j] <- sum(diff^pow[j])
		}else{
			diag.Sn <- diag(Sn)
			diag.Sn[diag.Sn <= 10^(-10)] <- 10^(-10)
			Ts[j] <- max(diff^2/diag.Sn)
		}
	}

	p.spu <- rep(NA, length(pow))
	Ts.perm <- matrix(NA, length(pow), n.perm)

	for(b in 1:n.perm){
		if(!is.null(seeds)) set.seed(seeds[b])
		perm <- sample(1:(n1 + n2))
		sam.perm <- sam[perm,]
		sam1.perm <- sam.perm[1:n1,]
		sam2.perm <- sam.perm[(n1 + 1):(n1 + n2),]
		Sn.perm <- ((n1 - 1)*cov(sam1.perm) + (n2 - 1)*cov(sam2.perm))/(n1 + n2 - 2)
		diff.perm <- colMeans(sam1.perm) - colMeans(sam2.perm)
		for(j in 1:length(pow)){
			if(pow[j] < Inf){
				Ts.perm[j, b] <- sum(diff.perm^pow[j])
			}
			if(pow[j] == Inf){
				diag.Sn.perm <- diag(Sn.perm)
				diag.Sn.perm[diag.Sn.perm <= 10^(-10)] <- 10^(-10)
				Ts.perm[j, b] <- max(diff.perm^2/diag.Sn.perm)
			}
		}
	}

	for(j in 1:length(pow)){
		p.spu[j] <- (sum(abs(Ts[j]) <= abs(Ts.perm[j,])) + 1)/(n.perm + 1)
		p.spu.perm <- (n.perm + 1 - rank(abs(Ts.perm[j,])))/n.perm
		if(j == 1){
			minp.perm <- p.spu.perm
		}else{
			minp.perm[which(minp.perm > p.spu.perm)] <- p.spu.perm[which(minp.perm > p.spu.perm)]
		}
	}

	p.aspu <- (sum(minp.perm <= min(p.spu)) + 1)/(n.perm + 1)
	pvs <- c(p.spu, p.aspu)
	names(pvs) <- c(paste("SPU", as.character(pow), sep = "_"), "aSPU")
	return(pvs) 
}  