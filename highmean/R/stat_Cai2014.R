stat_Cai2014 <- function(sam1, sam2){
	n1 <- dim(sam1)[1]
	n2 <- dim(sam2)[1]
	sam.cov1 <- cov(sam1)
	sam.cov2 <- cov(sam2)
	diff <- colMeans(sam1) - colMeans(sam2)
	diag1 <- diag(sam.cov1)
	diag1[diag1 <= 10^(-10)] <- 10^(-10)
	diag2 <- diag(sam.cov2)
	diag2[diag2 <= 10^(-10)] <- 10^(-10)
	test.stat <- max(diff^2/(diag1/n1 + diag2/n2))
	test.stat <- as.numeric(test.stat)
	return(test.stat)
}