stat_SPU <- function(sam1, sam2, pow){
	n1 <- dim(sam1)[1]
	n2 <- dim(sam2)[1]
	p <- dim(sam1)[2]
	sam.cov1 <- cov(sam1)
	sam.cov2 <- cov(sam2)
	diff <- colMeans(sam1) - colMeans(sam2)
	Ts <- rep(NA, length(pow))
	for(j in 1:length(pow)){
		if(pow[j] < Inf){
			Ts[j] <- sum(diff^pow[j])
		}else{
			diag1 <- diag(sam.cov1)
			diag1[diag1 <= 10^(-10)] <- 10^(-10)
			diag2 <- diag(sam.cov2)
			diag2[diag2 <= 10^(-10)] <- 10^(-10)
			Ts[j] <- max(diff^2/(diag1/n1 + diag2/n2))
		}
	}
	names(Ts) <- paste("SPU", as.character(pow), sep = "_")
	return(Ts)
}