proj_grad <- 
function(matA, ll, x = rep(0, ncol(matA))){
## FISTA algorithm for least square on the positive orthant
	x_1 <- x
	y <- x
	t <- 1
	
	res <- rep(0,1000)
	L <- max(svd(matA)$d)^2
	for(k in 1:1000){
		x_1 <- x
		temp <- t(matA)%*%(matA %*% y - ll)
		x <- y - 1/2/L * temp
		x[x <0 ] <- 0
		t1 <- (1+sqrt(1+4*t^2))/2
		
		y <- x + (t - 1) / t1 * (x - x_1)
		t <- t1
		res[k] <- sum((matA %*% x - ll)^2)
	}
	x
}
