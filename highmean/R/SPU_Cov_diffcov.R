SPU_Cov_diffcov <- function(gamma1, gamma2, n1, n2, cov1, cov2){
	ga1 <- ifelse(gamma1 > gamma2, gamma1, gamma2)
	ga2 <- ifelse(gamma1 <= gamma2, gamma1, gamma2)
	p <- dim(cov1)[1]
	P1 <- SPU_E_diffcov(ga1 + ga2, n1, n2, cov1, cov2)
	P2 <- -SPU_E_diffcov(ga1, n1, n2, cov1, cov2)*SPU_E_diffcov(ga2, n1, n2, cov1, cov2)
	c_d <- c_and_d(ga1, ga2)
	n.case <- dim(c_d)[1]
	P3 <- 0
	diags1 <- diag(cov1)
	mat1.col <- matrix(rep(diags1, p), p, p, byrow = FALSE)
	mat1.row <- matrix(rep(diags1, p), p, p, byrow = TRUE)
	diags2 <- diag(cov2)
	mat2.col <- matrix(rep(diags2, p), p, p, byrow = FALSE)
	mat2.row <- matrix(rep(diags2, p), p, p, byrow = TRUE)
	for(i in 1:n.case){
		c1 <- c_d$c1[i]
		c2 <- c_d$c2[i]
		c3 <- c_d$c3[i]
		d1 <- c_d$d1[i]
		d2 <- c_d$d2[i]
		d3 <- c_d$d3[i]
		mat <- mat1.col^c1*mat1.row^c2*cov1^c3*mat2.col^d1*mat2.row^d2*cov2^d3
		diag(mat) <- 0
		N <- factorial(ga1)*factorial(ga2)*sum(mat)
		D <- (n1^(c1 + c2 + c3)*n2^(d1 + d2 + d3)*factorial(c1)*factorial(c2)*factorial(d1)*factorial(d2)*factorial(c3)*factorial(d3)*2^(c1 + c2 + d1 + d2))
		P3 <- P3 + N/D
	}
	L.cov <- P1 + P2 + P3
	return(L.cov)
}