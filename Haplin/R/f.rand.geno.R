## MENDELIAN SELECTION
f.rand.geno <- function(x1, x2){
	## PASTE TOGETHER x1 AND x2, IN RANDOM PERMUTATION
	.rand <- as.logical(rbinom(length(x1), size = 1, prob = 0.5))
	.x1 <- ifelse(.rand, x1, x2)
	.x2 <- ifelse(!.rand, x1, x2)
	.x <- paste(.x1, .x2, sep = ";")	
	.x
}
