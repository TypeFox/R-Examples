##
##  w i l k i n s o n . R  Wilkinson matrix
##


wilkinson <- function(n){
	if (length(n) != 1 || n != round(n))
		stop("Argument 'n' must be an integer.")
	if (n <= 0) return(c())

	m <- (n-1)/2.0
	r <- rep(1, n-1)
	Diag(abs(-m:m)) + Diag(r, 1) + Diag(r, -1)
}
