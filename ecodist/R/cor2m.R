# Correlate variables from 1 matrix with variables from another
# matrix. If TRIM, set rho<critical value(alpha) to 0. Computes
# this critical value as a t-test with n-2 df.

cor2m <- function(x,y, trim=TRUE, alpha=0.05) {
	xz <- scale(as.matrix(x))
	yz <- scale(as.matrix(y))
	n <- dim(x)[[1]]
	cc <- t(yz) %*% xz / (n -1.0)
	rownames(cc) <- colnames(y)
	colnames(cc) <- colnames(x)
	if (trim) {
		rt <- cc * sqrt((n-2)/(1-cc^2))
		cv <- qt(1-alpha, n-2)
		cc[abs(rt) < cv] <- 0
	}
	cc
}

	
