###

## For performance reasons the constructor does not
## check for multiple or 'zero' elements.
##
## Argument 'strict' provides a choice whether to 
## enforce these constraints, or to reduce 'multiples'
## to NA (unless they all are identical) and remove
## 'zeros'.
##
reduce_simple_sparse_array <- 
function(x, strict = FALSE, order = FALSE) 
{
    if (!.Call(R__valid_ssa, x))
	stop("'x' not of class 'simple_sparse_array'")
    I <- x$i
    if (length(i <- attributes(I)) > 1L)
	dim(I) <- i$dim
    rm(i)
    V <- .Call(R_unattr, x$v)
    if (length(V)) {
	## reduce multiple entries
	i <- .Call(R_match_matrix, I, NULL, NULL)
	if (length(i[[1L]]) > length(i[[2L]])) {
	    if (strict)
		stop("multiple entries")
	    I <- I[i[[2L]],, drop = FALSE]
	    i <- structure(
		i[[1L]],
		levels = seq_len(dim(I)[1L]),
		class  = "factor"
	    )
	    V <- split(V, i)
	    rm(i)
	    names(V) <- NULL
	    V <- sapply(V, function(x) 
		if (length(x) > 1L) {
		    x <- as.list(x)
		    if (all(sapply(x[-1L], identical, x[[1L]])))
			x[[1L]]
		    else
			NA
		} else 
		    x, 
		USE.NAMES = FALSE)
	    warning("NAs introduced by reduction")
	} else
	    rm(i)
	## remove 'zero' entries
	i <- which(V == vector(typeof(V), 1L))
	if (strict)
	    stop("zero entries")
	if (length(i)) {
	    i <- -i
	    V <- V[i]
	    I <- I[i,, drop = FALSE]
	}
	rm(i)
	## order entries
	if (order) {
	    i <- do.call("order", rev(.Call(R_split_col, I)))
	    if (!identical(i, seq_along(i))) {
		V <- V[i]
		I <- I[i,, drop = FALSE]
	    }
	    rm(i)
	}
    }
    D <- as.vector(x$dim)
    N <- x$dimnames
    N <-
    if (!length(N) ||
	(is.null(names(N)) &&
	  all(sapply(N, is.null))))
	NULL
    else
	lapply(N, as.vector)
    simple_sparse_array(I, V, D, N)
}

##
drop_simple_sparse_array <-
function(x) 
{
    if (!is.simple_sparse_array(x))
	stop("'x' not of class 'simple_sparse_array'")
    dx <- x$dim
    if (any(dx == 0L))
	return(vector(typeof(x$v), 0L))	    ## sanitize
    k <- which(dx == 1L)
    if (length(k) == length(dx))
	return(x$v)
    if (length(k)) {
	k <- -k
	x$i <- x$i[, k, drop = FALSE]
	x$dim <- dx[k]
	if (!is.null(x$dimnames))
	    x$dimnames <- x$dimnames[k]
    }
    x
}


## see simplify2array
simplify_simple_sparse_array <-
function(x, higher = TRUE)
{
    if (!is.simple_sparse_array(x))
	stop("'x' not of class 'simple_sparse_array'")
    V <- x$v
    if (is.atomic(V) || 
	  !length(V))
	return(x)
    i <- unique(unlist(lapply(V, length)))
    ## FIXME not implemented
    if (length(i) > 1L) 
	return(x)
    if (!i)
	return(x)
    if (i == 1L) {
	x$v <- unlist(V, recursive = FALSE)
	return(x)
    }
    I <- x$i
    D <- x$dim
    N <- x$dimnames
    if (higher &&
	length(d <- unique(lapply(V, dim))) == 1L &&
      !is.null(d <- unlist(d))) {
	i <- d
	n <- dimnames(V[[1L]])
    } else 
	if (!is.null(n <- names(V[[1L]])))
	    n <- list(n)
    V <- unlist(V, recursive = FALSE)
    ## FIXME not optimized
    for (k in rev(i)) {
	l <- dim(I)[1L]
	if (k > 1L)
	    I <- apply(I, 2L, rep, each = k)
	I <- cbind(rep(seq.int(k), l), I)
    }
    if (!is.null(N)) {
	if (!is.list(n))
	    n <- list(n)
	N <- c(n, N)
    } else
	if (!is.null(n)) 
	    N <- list(n, vector("list", length(D)))
    D <- c(i, D)
    simple_sparse_array(I, V, D, N)
}


###
