
## 
extend_simple_sparse_array <-
function(x, MARGIN = 0L)
{
    if (!is.simple_sparse_array(x))
	stop("'x' not of class 'simple_sparse_array'")
    k <- MARGIN < 0L
    MARGIN[k] <- -MARGIN[k] -1L
    k <- MARGIN[1L]
    ## extend
    D <- c(1L, x$dim)
    I <- cbind(1L, x$i)
    N <- x$dimnames
    if (!is.null(N))
	N <- c(list(NULL), N)
    if (k > 0L)
	if (k > length(D))
	    stop("'MARGIN' invalid")
	else {	## order
	    i <- order(c(k + 1L, seq.int(length(D) - 1L))) 
	    D <- D[i]
	    I <- I[,i]
	    if (!is.null(N))
		N <- N[i]
	}
    x <- simple_sparse_array(I, x$v, D, N)
    rm(I, D, N)
    while (length(MARGIN <- MARGIN[-1L])) {
	k <- MARGIN > k
	MARGIN[k] <- MARGIN[k] + 1L
	x <- extend_simple_sparse_array(x, MARGIN[[1L]])
    }
    x
}

##
abind_simple_sparse_array <-
function(..., MARGIN = 1L)
{
    if (length(MARGIN) != 1L || 
	MARGIN == 0L)
	stop("'MARGIN' invalid")
    args <- list(...)
    if (length(args))
	args <- args[!sapply(args, is.null)]
    if (!length(args))
	return(NULL)
    x <- as.simple_sparse_array(args[[1L]])
    if (MARGIN < 0L)
	x <- extend_simple_sparse_array(x, MARGIN)
    if (length(args) == 1L)
	return(x)
    for (y in args[-1L]) {
	y <- as.simple_sparse_array(y)
	if (MARGIN < 0L)
	    y <- extend_simple_sparse_array(y, MARGIN)
	m <- abs(MARGIN)
	if (length(y$dim) == length(x$dim) - 1L)
	    y <- extend_simple_sparse_array(y, -min(m, length(x$dim)))
	else
	    if (length(y$dim) - 1L == length(x$dim)) {
		x <- extend_simple_sparse_array(x, -min(m, length(y$dim)))
	    } else
		if (length(y$dim) != length(x$dim))
		    stop("lengths of 'dim' do not conform")
	D <- x$dim
	m <- min(length(D), m)
	if (!identical(y$dim[-m], D[-m]))
	    stop("common parts of 'dim' do not conform")
	if (vector(typeof(x$v), 1L) != vector(typeof(y$v), 1L))
	    stop("definitions of ZERO of 'v' do not conform")
	V <- c(x$v, y$v)
	I <- y$i
	I[, m] <- I[, m] + D[m]
	I <- rbind(x$i, I)
	N <- x$dimnames
	if (!is.null(N[[m]])) {
	    N[[m]] <- 
		c(  N[[m]],
		    if (!is.null(y$dimnames[[m]]))
			y$dimnames[[m]]
		    else
			rep("", y$dim[[m]])
		)
	    if (is.null(names(N)))
                names(N) <- names(y$dimnames)
	} else
	    if (!is.null(y$dimnames[[m]])) {
		if (is.null(N))
		    N <- y$dimnames
		else
		    if (is.null(names(N)))
			names(N) <- names(y$dimnames)
		N[[m]] <-
		    c(
			rep("", D[m]),
			y$dimnames[[m]]
		     )
	    }
	D[m] <- D[m] + y$dim[m]
	x <- simple_sparse_array(I, V, D, N)
    }
    x
}


###
