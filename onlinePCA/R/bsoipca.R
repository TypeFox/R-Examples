bsoipca <- function(x, q, U, B, center, byrow=FALSE)
{
	# Convert x to matrix as needed
	if (!is.matrix(x)) x <- as.matrix(x)

	# Transpose x as needed
	if (byrow == TRUE) x <- t(x)

	# Center x as needed
	if (!missing(center)) x <- x - center

	# Data dimension
	n <- ncol(x)
	d <- nrow(x)

	# Determine q if missing but U specified
	if (missing(q) && !missing(U)) q <- NCOL(U)

	# Initialize PCs if missing
	if (missing(U))
		U <- qr.Q(qr(matrix(rnorm(d*q),d,q)))

	# Check that dimensions of x and U are compatible
	stopifnot(NROW(U) == d)

	
	# Determine block size if missing
	if (missing(B)) {	
		nblock <- ceiling(log(d))
		B <- floor(n/nblock)
	} else {nblock <- n %/% B}
		
	# Loop 
	for (i in 1:nblock)
	{
		ind <- seq.int((i-1)*B+1,i*B)
		S <- x[,ind] %*% crossprod(x[,ind],U)	
		U <- qr.Q(qr(S/B))
	}
	
	return(U)
}
