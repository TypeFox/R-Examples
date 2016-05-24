rref <- function(A, tol=sqrt(.Machine$double.eps),verbose=FALSE,
		fractions=FALSE){
	## A: coefficient matrix
	## tol: tolerance for checking for 0 pivot
	## verbose: if TRUE, print intermediate steps
	## fractions: try to express nonintegers as rational numbers
	## Written by John Fox
	if (fractions) {
            ## Added MASS to Depends in DESCRIPTION file
		##mass <- require(MASS)
		##if (!mass) stop("fractions=TRUE needs MASS package")
	}
	if ((!is.matrix(A)) || (!is.numeric(A)))
		stop("argument must be a numeric matrix")
	n <- nrow(A)
	m <- ncol(A)
	for (i in 1:min(c(m, n))){
		col <- A[,i]
		col[1:n < i] <- 0
		# find maximum pivot in current column at or below current row
		which <- which.max(abs(col))
		pivot <- A[which, i]
		if (abs(pivot) <= tol) next     # check for 0 pivot
		if (which > i) A[c(i, which),] <- A[c(which, i),]  # exchange rows
		A[i,] <- A[i,]/pivot            # pivot
		row <- A[i,]
		A <- A - outer(A[,i], row)      # sweep
		A[i,] <- row                    # restore current row
		if (verbose)
			if (fractions) print(MASS::fractions(A))
			else print(round(A,round(abs(log(tol,10)))))
	}
	for (i in 1:n)
		if (max(abs(A[i,1:m])) <= tol)
			A[c(i,n),] <- A[c(n,i),] # 0 rows to bottom
	if (fractions) fractions (A)
	else round(A, round(abs(log(tol,10))))
}
