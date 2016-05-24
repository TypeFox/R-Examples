incRpca.block <- function (x, B, lambda, U, n0 = 0, f, 
	q = length(lambda), center, byrow=FALSE) 
{
	if (!is.matrix(x)) x <- as.matrix(x)
	if (byrow == TRUE) x <- t(x)
	if (!missing(center)) x <- x - center

	# Data dimension
	n <- ncol(x)
	d <- nrow(x)
		
	# Initialize PCA if no initial values provided
	if (missing(U)) {
		svdx <- svd(x[,1:B])
		U <- svdx$u[,1:min(B,q)]
		lambda <- svdx$d[1:min(B,q)]^2/B
		x <- x[,-(1:B)]
		n <- n - B
		n0 <- n0 + B } 

	# Check dimension compatibility of lambda, U, x
	stopifnot(NCOL(U) == length(lambda) && NROW(U) == d)
	
	# Initialize B and f if missing
	if (missing(B)) B <- n
	nblock <- n %/% B
	if (missing(f)) f <- B/(n0+(0:(nblock-1))*B)

	# Loop
	for (i in 1:nblock)
	{
		lambda <- lambda * (1-f[i])
		l <- length(lambda)
		qwork <- min(q,l+B)
		ind <- seq.int((i-1)*B+1,i*B)
		QR <- qr(cbind(U %*% diag(sqrt(lambda),l,l), 
				x[,ind] * sqrt(f[i]/B)))
		svdR <- svd(qr.R(QR), qwork, 0)	
		lambda <- svdR$d[1:qwork]^2
		U <- qr.Q(QR) %*% svdR$u
	}

	return(list(values=lambda, vectors=U))
}




