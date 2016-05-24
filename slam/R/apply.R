## CB 2013/12
colapply_simple_triplet_matrix <-
function(x, FUN, ...) {
    FUN <- match.fun(FUN)
    out <- .External(R_col_apply_stm, x, FUN, ...)
    if (length(out)) {
	if (all(unlist(lapply(out, length)) == 1L))
	    out <- unlist(out, recursive = FALSE, use.names = FALSE)
	names(out) <- colnames(x)
    } 
    else
	## NOTE we always supplie as matrix in case dimensions 
	##	must conform with further arguments.
	storage.mode(out) <-
	    typeof(FUN(as.matrix(x), ...))
    out
}

rowapply_simple_triplet_matrix <-
function(x, FUN, ...) {
    FUN <- match.fun(FUN)
    if (!is.simple_triplet_matrix(x))
	stop("'x' not of class simple_striplet_matrix")
    colapply_simple_triplet_matrix(t(x), FUN, ...)
}

## FIXME a workaround for a proper C implementation.
crossapply_simple_triplet_matrix <- 
function(x, y = NULL, FUN, ...) {
    FUN <- match.fun(FUN)
    if (is.null(y)) {
	if (!is.simple_triplet_matrix(x))
	    stop("'x' not of class simple_triplet_matrix")
	Y <- x
	out <- colapply_simple_triplet_matrix(x, function(x) {
	out <- colapply_simple_triplet_matrix(Y, FUN, x, ...)
	    Y <<- Y[, -1L]
	    out
	})
	out <- unlist(out, recursive = FALSE, use.names = FALSE)
	Y <- simple_triplet_zero_matrix(x$ncol)
	Y <- row(Y) >= col(Y)
	out[Y] <- out
	out <- matrix(out, nrow = x$ncol, ncol = x$ncol, byrow = TRUE,
	    dimnames = if (!is.null(colnames(x)))
		list(colnames(x), colnames(x))
	)
	out[Y] <- t(out)[Y]
	return(out)
    }
    if (is.simple_triplet_matrix(y)) {
	if (!is.simple_triplet_matrix(x))
	    return(
		t(crossapply_simple_triplet_matrix(y, as.matrix(x), 
		    function(y, x) FUN(x, y, ...)))
	    )
	if (x$nrow != y$nrow)
	    stop("the numer of rows of 'x' and 'y' do not conform")
	## Fix asymmetric performance.
	if (x$ncol > y$ncol)
	    return(
		t(crossapply_simple_triplet_matrix(y, x,
		    function(y, x) FUN(x, y, ...)))
	    )
	if (y$ncol > 0L &&
	    x$ncol > 0L) {
	    out <- colapply_simple_triplet_matrix(x, function(x)
		   colapply_simple_triplet_matrix(y, function(y)
			FUN(x, y, ...)))
	}
	else
	    out <- colapply_simple_triplet_matrix(x[, 0L], 
		FUN, as.matrix(y[, 0L]), ...) 
    }
    else {
	if (!is.simple_triplet_matrix(x))
	    stop("'x, y' not of class simple_triplet_matrix")
	y <- as.matrix(y)
	if (x$nrow != nrow(y))
	    stop("the numer of rows of 'x' and 'y' do not conform")
	if (ncol(y) > 0L &&
	     x$ncol > 0L) {
	    Y <- split(y, factor(col(y), levels = seq_len(ncol(y))))
	    out <- colapply_simple_triplet_matrix(x, function(x) {
		out <- lapply(Y, function(y)
		    FUN(x, y, ...))
		if (all(unlist(lapply(out, length)) == 1L))
		    out <- unlist(out, recursive = FALSE, use.names = FALSE)
		out
	    })
	    rm(Y)
	}
	else
	    out <- colapply_simple_triplet_matrix(x[, 0L], 
		FUN, y[, 0L, drop = FALSE], ...)
    }
    out <- unlist(out, recursive = FALSE, use.names = FALSE)
    out <- matrix(out, nrow = x$ncol, ncol = ncol(y), byrow = TRUE,
	dimnames = 
	    if (!is.null(colnames(x)) || !is.null(colnames(y)))
		list(colnames(x), colnames(y))
    )
    out
}

tcrossapply_simple_triplet_matrix <-
function(x, y = NULL, FUN, ...) {
    FUN <- match.fun(FUN)
    if (is.simple_triplet_matrix(x))
	crossapply_simple_triplet_matrix(t(x),
	    if (is.null(y))
		y
	    else
		if (is.simple_triplet_matrix(y))
		    t(y)
		else
		    t(as.matrix(y)),
	    FUN, ...
    )
    else
	if (is.simple_triplet_matrix(y))
	    crossapply_simple_triplet_matrix(t(as.matrix(x)), t(y), FUN, ...)
	else
	    stop("'x, y' not of class simple_triplet_matrix")
}

###
