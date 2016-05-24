NOEC <- function(x, expr, sigLev = 0.05){
	# NOEC and LOEC calculation using Dunnett's test
	## Dunnett, C.W., 1964. New tables for multiple comparisons with controls. Biometrics 30, 482-491
	n <- nrow(expr)
	m <- ncol(expr)
	C <-  sum(expr)^2 / ((n + 1) * m) 
	Tj <- rowSums(expr)
	SSB <- sum(Tj^2) / m - C
	SST <- sum(expr^2) - C
	SSW <- SST - SSB
	SW <- sqrt((SST - SSB) / ((n + 1) * (m - 1)))
	tj <- (Tj / m) / (SW * sqrt(1 / m + 1 / m))
	probF <- qf(1 - sigLev / 2, (n + 1) * (m - 1), n)
	noecSign <- sign(abs(tj) - probF)
	idx.one <- which(noecSign == -1)
	if(length(idx.one) == 0) noec = NULL else noec = x[idx.one[length(idx.one)]]
	idx.two <- which(noecSign == 1)
	if(length(idx.two) == 0) loec = NULL else loec = x[idx.two[1]]
	mat <- cbind(x, tj, probF, noecSign)
	list(mat = mat, no = noec, lo = loec)
}
