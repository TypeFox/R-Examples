##
##  Matlab flipping matrices functions
##


flipdim <- function(a, dim=1) {
	if (!is.matrix(a))
		stop("Argument 'a' must ba a matrix.")
	if (!(dim %in% c(1,2)))
		stop("Argument 'dim' must be 1 or 2 (for rows or columns).")

	if (dim == 1) {
		a <- a[nrow(a):1, ]
	} else {
		a <- a[ ,ncol(a):1]
	}
	return(a)
}

flipud <- function(a) {
	flipdim(a, 1)
}

fliplr <- function(a) {
	flipdim(a, 2)
}

rot90 <- function(a, k=1) {
	if (!is.matrix(a))
		stop("Argument 'a' must ba a matrix.")
	if (floor(k) != ceiling(k))
		k <- 0

	switch(EXPR = 1 + (k %% 4),
		a,
		t(a[, seq(ncol(a), 1, by=-1)]),
		a[seq(nrow(a), 1, by=-1), seq(ncol(a), 1, by=-1)],
		{a <- t(a); a[, seq(ncol(a), 1, by=-1)]}
	)
}
