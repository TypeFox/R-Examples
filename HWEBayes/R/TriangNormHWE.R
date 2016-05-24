TriangNormHWE <-
function(nvec){
	   if (length(nvec) != 3) stop("TriangNormHWE: Dimension of nvec not equal to 3\n")
	   n11 <- nvec[1]; n12 <- nvec[2]; n22 <- nvec[3]
	   a <- 2*n11+n12; b <- n12+2*n22
	   term1 <- pbeta(.5,a+2,b+1)*beta(a+2,b+1)
	   term2 <- beta(a+1,b+2)
	   term3 <- pbeta(.5,a+2,b+1)*beta(a+2,b+1)
	   logconst <- lfactorial(sum(nvec)) - sum(lfactorial(nvec)) + (n12+2)*log(2)
	   TriangNormHWE <- exp(logconst)*(term1+term2-term3)
	   TriangNormHWE
}

