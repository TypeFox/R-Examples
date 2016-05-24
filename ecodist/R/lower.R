lower <- function(m)
{
# Takes the lower triangle of a matrix
# Does NOT check for symmetric matrix

m <- as.matrix(m)

if(ncol(m) != nrow(m))
	stop("Matrix not square.")
if(ncol(m) == 1) {
   warning("Matrix is 1x1.")
   m
}
else 
	m[col(m) < row(m)]
}
