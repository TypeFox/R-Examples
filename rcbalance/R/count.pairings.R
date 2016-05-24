count.pairings <- function(z, exact) {
	if (length(z) != length(exact)) {
		stop("z and exact must be the same length")
	}
	if (length(unique(z)) != 2 ) {
		stop("z must contain exactly 2 distinct values")
	}
	matchtab <- table(exact,z)
	sum(matchtab[,1]*matchtab[,2])
}