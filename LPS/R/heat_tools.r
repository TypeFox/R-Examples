# Distances : 1 - Pearson (time efficient implementation)
# Author : Sylvain Mareschal <maressyl@gmail.com>
dist.COR <- function(input) {
	# Distance matrix
	mtx <- 1 - cor(t(input), method="pearson", use="na.or.complete")

	# Lower half matrix for dist object
	object <- as.double(mtx[ row(mtx) > col(mtx) ])

	# dist object
	attr(object, "Size") <- dim(mtx)[1]
	attr(object, "Labels") <- dimnames(mtx)[[1]]
	attr(object, "Diag") <- FALSE
	attr(object, "Upper") <- FALSE
	attr(object, "method") <- "1 - Pearson"
	class(object) <- "dist"

	return(object)
}

# Agglomeration : Ward
# Author : Sylvain Mareschal <maressyl@gmail.com>
hclust.ward <- function(input) {
	hclust(input, method="ward")
}

