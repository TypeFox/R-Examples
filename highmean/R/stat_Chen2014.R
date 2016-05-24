stat_Chen2014 <- function(sam1, sam2){
	sam.cov1 <- cov(sam1)
	sam.cov2 <- cov(sam2)
	test.stat <- max_threshold_diffcov(sam1, sam2, sam.cov1, sam.cov2)
	test.stat <- as.numeric(test.stat)
	return(test.stat)
}