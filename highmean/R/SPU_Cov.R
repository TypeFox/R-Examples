SPU_Cov <- function(gamma1, gamma2, n1, n2, cov){
	ga1 <- ifelse(gamma1 > gamma2, gamma1, gamma2)
	ga2 <- ifelse(gamma1 <= gamma2, gamma1, gamma2)
	P1 <- SPU_E(ga1 + ga2, n1, n2, cov)
	P2 <- -SPU_E(ga1, n1, n2, cov)*SPU_E(ga2, n1, n2, cov)
	c_d <- c_and_d(ga1, ga2)
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
		N <- factorial(ga1)*factorial(ga2)*sum(mat)
		D <- (n1^(c1 + c2 + c3)*n2^(d1 + d2 + d3)*factorial(c1)*factorial(c2)*
			factorial(d1)*factorial(d2)*factorial(c3)*factorial(d3)*
			2^(c1 + c2 + d1 + d2))
		P3 <- P3 + N/D
	}
	L.cov <- P1 + P2 + P3
	return(L.cov)
}