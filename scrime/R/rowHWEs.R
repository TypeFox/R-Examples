rowHWEs <- function(x, levels=1:3, affy=FALSE, check=TRUE){
	mat.obs <- rowTables(x, levels=levels, affy=affy, check=check)
	n <- rowSums(mat.obs)
	p <- (mat.obs[,1] + 0.5 * mat.obs[,2])
	pb <- n-p 
	mat.exp <- mat.obs
	mat.exp[,1] <- p^2 
	mat.exp[,2] <- 2 * p * pb
	mat.exp[,3] <- pb^2 
	stats <-  n * (rowSums(mat.obs * mat.obs / mat.exp) - 1)
	rawp <- pchisq(stats, 1, lower.tail=FALSE)
	structure(list(stats=stats, rawp=rawp))
}

 
	