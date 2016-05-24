###

rollup <- 
function(x, MARGIN, INDEX, FUN, ...)
    UseMethod("rollup")

rollup.array <-
function(x, MARGIN, INDEX = NULL, FUN = sum, ..., DROP = FALSE) {
    if (is.character(MARGIN))
        MARGIN <- match(MARGIN, names(dimnames(x)))
    if (!all(match(MARGIN, seq_along(dim(x)), nomatch = 0L)))
        stop("'MARGIN' invalid")
    if (is.null(INDEX))
	INDEX <- vector("list", length(MARGIN)) 
    else {
	if (is.atomic(INDEX))
	    INDEX <- list(INDEX)
	if (length(INDEX) != length(MARGIN))
	    stop("'INDEX' invalid length")
    }
    names(INDEX) <- MARGIN
    FUN <- match.fun(FUN)
    d <- dim(x)
    n <- dimnames(x)
    if (is.null(n))
	n <- vector("list", length(d))
    i <- arrayInd(seq_along(x), .dim = d)
    for (k in MARGIN) {
	z <- INDEX[[as.character(k)]]
	z <-
	if (is.null(z))
	    structure(
		rep(1L, d[k]),
		levels = "1",
		class  = "factor"
	    )
	else {
	    if (length(z) != d[k])
		stop(gettextf("INDEX [%s] invalid length", k),
                     domain = NA)
	    factor(z)
	}
	i[, k] <- z[i[, k]]
	z <- levels(z)
	d[k]   <- length(z)
	n[[k]] <- z
	rm(z)
    }
    i <- .Call(R_vector_index, d, i)
    i <- structure(
	i,
	levels = seq_len(prod(d)),
	class  = "factor"
    )
    i <- split.default(x, i)
    names(i) <- NULL
    i <- lapply(i, FUN, ...)
    if (all(unlist(lapply(i, length)) == 1L))
	i <- unlist(i, recursive = FALSE, use.names = FALSE)
    ## NOTE see drop_simple_sparse_array
    if (DROP) {
	if (any(d == 0L))
	    return(i)
	k <- which(d == 1L)
	if (length(k) == length(d))
	    return(i)
	if (length(k)) {
	    k <- -k
	    d <- d[k]
	    n <- n[k]
	}
    }
    array(i, d, n)
}

rollup.matrix <- rollup.array

rollup.simple_sparse_array <-
function(x, MARGIN, INDEX = NULL, FUN = sum, ..., DROP = FALSE,
	    EXPAND = c("none", "sparse", "dense", "all")) {
    if (is.character(MARGIN)) 
	MARGIN <- match(MARGIN, names(dimnames(x)))
    if (!all(match(MARGIN, seq_along(dim(x)), nomatch = 0L)))
	stop("'MARGIN' invalid")
    if (is.null(INDEX))
	INDEX <- vector("list", length(MARGIN)) 
    else {
	if (is.atomic(INDEX))
	    INDEX <- list(INDEX)
	if (length(INDEX) != length(MARGIN))
	    stop("'INDEX' invalid length")
    }
    names(INDEX) <- MARGIN
    FUN <- match.fun(FUN)
    EXPAND <- match(
	match.arg(EXPAND), 
	eval(formals(rollup.simple_sparse_array)$EXPAND)
    )
    D <- dim(x)
    I <- x$i
    if (EXPAND > 1L) {
	if (EXPAND > 2L)
	    P <- array(1L, dim(I))
	T <- vector("list", length(D))
	for (k in seq_along(D)[-MARGIN])
	    T[[k]] <- rep(1L, D[k])
    }
    N <- dimnames(x)
    if (is.null(N))
	N <- vector("list", length(D))
    V <- x$v
    if (EXPAND < 4L &&
	!.Call(R__valid_v, V))
	stop("component 'v' contains 'ZERO' value(s)")
    for (k in MARGIN) {
	z <- INDEX[[as.character(k)]]
	z <-
	if (is.null(z))
	    structure(
		rep(1L, D[k]),
		levels = "1",
		class  = "factor"
	    )
	else {
	    if (length(z) != D[k])
		stop(gettextf("INDEX [%s] invalid length", k),
                     domain = NA)
	    factor(z)
	}
	l <- levels(z)
	D[k]   <- length(l)
	N[[k]] <- l
	i <- I[, k]
	if (EXPAND > 1L) {
	    if (EXPAND > 2L) {
		p <- .Call(R_part_index, z)
		T[[k]] <- attr(p, "table")
		P[, k] <- p[i]
		rm(p)
	    } else
		T[[k]] <- tabulate(z, length(l))
	}
	i <- z[i]
	rm(l, z)
	I[, k] <- i
	i <- is.na(i)
	i <- which(i)
	if (length(i)) {
	    i <- - i
	    I <- I[i,, drop = FALSE]
	    V <- V[i]
	    if (EXPAND > 2L)
		P <- P[i,, drop = FALSE]
	}
	rm(i)
    }
    if (EXPAND == 4L) {
	## NOTE see src/main/unique.c in the R
	##	source code.
	k <- prod(D)
	if (k > 1073741824L)
	    stop("number of cells %d too large for hashing", k)
	i <- .Call(R_vector_index, D, I)
	I <- arrayInd(seq_len(k), .dim = D)
	k <- .Call(R_vector_index, D, I)
	i <- match(i, k)
	rm(k)
    } else {
	i <- .Call(R_match_matrix, I, NULL, NULL)
	I <- I[i[[2L]],, drop = FALSE]
	i <-   i[[1L]]
    }
    i <- structure(
	i, 
	levels = seq_len(dim(I)[1L]), 
	class  = "factor"
    )
    if (EXPAND == 1L) {
	V <- split.default(V, i)
	rm(i)
	names(V) <- NULL
	V <- lapply(V, FUN, ...)
    } else {
	.pt <- proc.time()
	message(gettextf("processing %d cells ... ", dim(I)[1L]),
		appendLF = FALSE,
                domain = NA)
	i <- split.default(seq_along(i), i)
	names(i) <- NULL
	V <- mapply(function(i, z) {
		z <- I[z, ]
		z <- mapply("[", T, z)
		if (EXPAND > 2L) {
		    ## NOTE this consumes less computation time
		    ##	    and memory than
		    ## z <- array(vector(typeof(V),1L), z)
		    ## z[P[i,, drop = FALSE]] <- V[i]
		    z <- .Call(R_ini_array, z, P, V, i)
		    FUN(z, ...)
		} else
		    FUN(V[i], prod(z) - length(i), ...)
	    },
	    i,
	    seq_along(i),
	    SIMPLIFY = FALSE, USE.NAMES = FALSE
	)
	rm(i, T)
	if (EXPAND > 2L)
	    rm(P)
	message(sprintf("[%.2fs]\n", (proc.time() - .pt)[3L]),
		appendLF = FALSE,
                domain = NA)
    }
    if (all(unlist(lapply(V, length)) == 1L)) {
	V <- unlist(V, recursive = FALSE, use.names = FALSE)
	i <- V == vector(typeof(V), 1L)
	i <- which(i)
	if (length(i)) {
	    i <- - i
	    I <- I[i,, drop = FALSE]
	    V <- V[i]
	}
    }
    x <- simple_sparse_array(I, V, D, N)
    rm(I, V, D, N)
    if (DROP)
	x <- drop_simple_sparse_array(x)
    x
}


rollup.simple_triplet_matrix <- 
function(x, MARGIN, INDEX = NULL, FUN = sum, ...) {
    FUN <- match.fun(FUN)
    if (!identical(FUN, sum)) {
	if (!is.null(list(...)$DROP))
	    stop("'DROP' not supported")
	x <- rollup.simple_sparse_array(as.simple_sparse_array(x), 
	    MARGIN, INDEX, FUN, ...
	)
	return(as.simple_triplet_matrix(x))
    }
    if (is.character(MARGIN)) 
	MARGIN <- match(MARGIN, names(dimnames(x)))
    if (!all(match(MARGIN, seq_along(dim(x)), nomatch = 0L)))
	stop("'MARGIN' invalid")
    if (is.null(INDEX))
	INDEX <- vector("list", length(MARGIN))
    else {
	if (is.atomic(INDEX))
	    INDEX <- list(INDEX)
	if (length(INDEX) != length(MARGIN))
	    stop("'INDEX' invalid length")
    }
    names(INDEX) <- MARGIN
    for (k in MARGIN) {
	x <- switch(k,
	    t(rollup(t(x), 2L, INDEX[as.character(k)], FUN, ...)),
	    {
		z <- INDEX[[as.character(k)]]
		z <- 
		if (is.null(z))
		    structure(
			rep(1L, dim(x)[k]),
			levels = "1",
			class  = "factor"
		    )
		else {
		    if (length(z) != dim(x)[k])
			stop(gettextf("INDEX [%s] invalid length", k),
                             domain = NA)
		    factor(z)
		}
		.Call(R_row_tsums, 
		      x, z, 
		      if (is.null(list(...)$na.rm))
			  FALSE
		      else
			  as.logical(list(...)$na.rm), 
		      FALSE
		)
	    }
	)
    }
    x
}

##
rollup.default <-
function(x, MARGIN, INDEX = NULL, FUN = sum, ..., DROP = FALSE) {
    if (!length(dim(x)))
	stop("dim(x) must have a positive length")
    rollup(as.array(x), MARGIN, INDEX, FUN, ..., DROP = DROP)
}

###
