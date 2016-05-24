SPU_Var <- function(gamma, n1, n2, cov){
	p <- dim(cov)[1]
	if(gamma == 1){
		var <- (1/n1 + 1/n2)*(t(rep(1, p)) %*% cov %*% rep(1, p))
		var <- as.numeric(var)
	}else{
		P1 <- SPU_E(2*gamma, n1, n2, cov)
		P2 <- -(SPU_E(gamma, n1, n2, cov))^2
		c_d <- c_and_d(gamma, gamma)
		n.case <- dim(c_d)[1]
		P3 <- 0
		diags <- diag(cov)
		p <- length(diags)
		mat1 <- matrix(rep(diags, p), p, p, byrow = FALSE)
		mat2 <- matrix(rep(diags, p), p, p, byrow = TRUE)
		for(i in 1:n.case){
			c1 <- c_d$c1[i]
			c2 <- c_d$c2[i]
			c3 <- c_d$c3[i]
			d1 <- c_d$d1[i]
			d2 <- c_d$d2[i]
			d3 <- c_d$d3[i]
			mat <- mat1^(c1 + d1)*mat2^(c2 + d2)*cov^(c3 + d3)
			diag(mat) <- 0
			N <- (factorial(gamma))^2*sum(mat)
			D <- (n1^(c1 + c2 + c3)*n2^(d1 + d2 + d3)*factorial(c1)*factorial(c2)*
				factorial(d1)*factorial(d2)*factorial(c3)*factorial(d3)*
				2^(c1 + c2 + d1 + d2))
			P3 <- P3 + N/D
		}
		var <- P1 + P2 + P3
	}
	return(var)
}