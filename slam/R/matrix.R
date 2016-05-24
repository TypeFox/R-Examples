## A simple class for sparse (triplet) matrices.

## Mostly intended for being able to take advantage of LP solvers which
## allow for sparse specifictions of (possible rather large) constraint
## matrices.

simple_triplet_matrix <-
function(i, j, v, nrow = max(i), ncol = max(j), dimnames = NULL)
{
    stm <- list(i = as.integer(i), j = as.integer(j), v = v,
                nrow = as.integer(nrow), ncol = as.integer(ncol),
                dimnames = dimnames)
    class(stm) <- "simple_triplet_matrix"
    if(!.Call(R__valid_stm, stm))
	stop("failed to create a valid 'simple_sparse_matrix' object")
    stm
}

as.simple_triplet_matrix <-
function(x)
    UseMethod("as.simple_triplet_matrix")

as.simple_triplet_matrix.simple_triplet_matrix <- identity

as.simple_triplet_matrix.matrix <-
function(x)
{
    x <- unclass(x)
    if(!prod(dim(x)))
        return(simple_triplet_matrix(integer(), integer(), c(x),
                                     nrow = nrow(x), ncol = ncol(x),
                                     dimnames = dimnames(x)))
    ind <- which(is.na(x) | (x != vector(typeof(x), 1L)), arr.ind = TRUE)
    dimnames(ind) <- NULL
    simple_triplet_matrix(ind[, 1L], ind[, 2L], x[ind],
                          nrow = nrow(x), ncol = ncol(x),
                          dimnames = dimnames(x))
}

as.simple_triplet_matrix.default <-
function(x)
    as.simple_triplet_matrix(unclass(as.matrix(x)))

## Sparse matrix classes in package 'Matrix'.

as.simple_triplet_matrix.dgTMatrix <-
function(x)
{
    simple_triplet_matrix(x@i + 1L, x@j + 1L, x@x,
                          x@Dim[1L], x@Dim[2L], x@Dimnames)
}

as.simple_triplet_matrix.dgCMatrix <-
function(x)
{
    nc <- x@Dim[2L]
    simple_triplet_matrix(x@i + 1L, rep.int(seq_len(nc), diff(x@p)),
                          x@x,
                          x@Dim[1L], nc, x@Dimnames)
}

as.simple_triplet_matrix.dgRMatrix <-
function(x)
{
    nr <- x@Dim[1L]
    simple_triplet_matrix(rep.int(seq_len(nr), diff(x@p)), x@j + 1L,
                          x@x,
                          nr, x@Dim[2L], x@Dimnames)
}

## See Work/Matrix.R for S4 methods for coercing simple triplet matrices
## to Matrix objects.

## Sparse matrix classes in package 'SparseM'.

as.simple_triplet_matrix.matrix.coo <-
function(x)
    simple_triplet_matrix(x@ia, x@ja, x@ra,
                          x@dimension[1L], x@dimension[2L])

as.simple_triplet_matrix.matrix.csc <-
function(x)
{
    nc <- x@dimension[2L]
    simple_triplet_matrix(x@ja, rep.int(seq_len(nc), diff(x@ia)), x@ra,
                          x@dimension[1L], nc)
}

as.simple_triplet_matrix.matrix.csr <-
function(x)
{
    nr <- x@dimension[1L]
    simple_triplet_matrix(rep.int(seq_len(nr), diff(x@ia)), x@ja, x@ra,
                          nr, x@dimension[2L])
}

## Sparse matrix class in package 'spam'.

as.simple_triplet_matrix.spam <-
function(x)
{
    nr <- x@dimension[1L]
    simple_triplet_matrix(rep.int(seq_len(nr), diff(x@rowpointers)),
                          x@colindices, x@entries,
                          nr, x@dimension[2L])
}

as.matrix.simple_triplet_matrix <-
function(x, ...)
{
    nr <- x$nrow
    nc <- x$ncol
    y <- matrix(vector(typeof(x$v), nr * nc), nr, nc)
    y[cbind(x$i, x$j)] <- x$v
    dimnames(y) <- x$dimnames
    y
}

as.array.simple_triplet_matrix <-
function(x, ...)
    as.array(as.matrix.simple_triplet_matrix(x, ...))

as.simple_triplet_matrix.simple_sparse_array <-
function(x) {
    dx <- x$dim
    if(length(dx) == 1L) {
        simple_triplet_matrix(
	    i = x$i[, 1L],
            j = rep.int(1L, nrow(x$i)),
            v = x$v,
            nrow = dx,
            ncol = 1L,
            dimnames = 
		if (!is.null(x$dimnames))
		    c(x$dimnames, list(NULL))
		else
		    NULL
	)
    } 
    else 
    if(length(dx) == 2L) {
	simple_triplet_matrix(
	    i = x$i[,1L],
	    j = x$i[,2L],
	    v = x$v,
	    nrow = x$dim[1L],
	    ncol = x$dim[2L],
	    dimnames = x$dimnames
	)
    }
    else
	stop("Unsupported number of dimensions")
}


is.simple_triplet_matrix <-
function(x)
    inherits(x, "simple_triplet_matrix")

is.numeric.simple_sparse_array <-
is.numeric.simple_triplet_matrix <-
function(x)
    is.numeric(x$v)

Math.simple_triplet_matrix <-
function(x, ...)
{
    ## Functions in the Math group mapping 0 to 0:
    funs <- c("abs", "sign", "sqrt",
              "floor", "ceiling", "trunc", "round", "signif")
    if(is.na(match(as.character(.Generic), funs)))
        stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
                      .Generic, .Class),
             domain = NA)

    x$v <- get(.Generic)(x$v, ...)
    x
}

Ops.simple_triplet_matrix <-
function(e1, e2)
{
    ## Currently, we only implement the following (for numeric
    ## operands):

    ## * Unary plus and minus.
    ##
    ## * Addition, subtraction and multiplication of two compatible
    ##   simple triplet matrices (or operands coercible to these).
    ##   [Division by a simple triplet matrix typically involves
    ##   division by zero and hence is not provided.]
    ##
    ## * Multiplication and division of a simple triplet matrix x by a
    ##   number or a vector of length nrow(x) (allowing to conveniently
    ##   scale the rows of a numeric simple triplet matrix).
    ##
    ## * Non-equality comparison of a simple triplet matrix with 0.
    ##
    ## * Comparisons of the elements of a simple triplet matrix with a
    ##   number.
    ##
    ## More could be added (but note that the elements could have
    ## arbitrary modes).

    ## Drop zero-valued elements
    .reduce <- function(x) {
	ind <- which(!x$v)
	if(length(ind)) {
	    ind <- -ind
	    x$i <- x$i[ind]
	    x$j <- x$j[ind]
	    x$v <- x$v[ind]
	}
	x
    }

    op <- as.character(.Generic)

    if(nargs() == 1L) {
        if(op == "+") return(e1)
        if(op == "-") {
            e1$v <- - e1$v
            return(e1)
        }
        stop(gettextf("Unary '%s' not defined for \"%s\" objects.",
                      .Generic, .Class),
             domain = NA)
    }

    if(!(op %in% c("+", "-", "*", "/", "^",
                   "==", "!=", "<", "<=", ">", ">=")))
        stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
                      .Generic, .Class),
             domain = NA)

    ## Require numeric operands for the arithmetic operators.
    if(!is.numeric(e1) || !is.numeric(e2))
        stop("Not implemented.")

    if(op %in% c("==", "!=", "<", "<=", ">", ">=")) {
        if(length(e2) == 1L) {
            if(is.na(e2))
                stop("NA/NaN handling not implemented.")
	    names(e2) <- NULL
            ind <- if(do.call(.Generic, list(0, e2))) {
                ## This inverts the sparse storage advantage, and hence
                ## will typically be inefficient.  Need to find the row
                ## and column positions of the zero entries.
                m <- matrix(TRUE, e1$nrow, e1$ncol)
                m[cbind(e1$i, e1$j)] <- FALSE
                which(m, arr.ind = TRUE)
            } else 
		integer()
            e1$v <- do.call(.Generic, list(e1$v, e2))
	    e1 <- .reduce(e1)
            if(n <- NROW(ind)) {
                e1$i <- c(e1$i, ind[, 1L])
                e1$j <- c(e1$j, ind[, 2L])
                e1$v <- c(e1$v, rep.int(TRUE, n))
            }
            return(e1)
        }
        stop("Not implemented.")
    }

    if(op == "^") {
        ## Allow for taking (single) positive exponents.
        if(is.object(e2) || (length(e2) != 1L) ||
           !is.finite(e2) || (e2 <= 0))
            stop("Not implemented.")
	names(e2) <- NULL
        e1$v <- e1$v ^ e2
        return(e1)
    }

    .make_dimnames <- function(e1, e2) {
        if(is.null(rnms <- rownames(e1)))
            rnms <- rownames(e2)
        if(is.null(cnms <- colnames(e1)))
            cnms <- colnames(e2)
        if(is.null(rnms) && is.null(cnms))
            NULL
        else {
            out <- list(rnms, cnms)
            if(is.null(nms <- names(dimnames(e1))))
                nms <- names(dimnames(e2))
            names(out) <- nms
            out
        }
    }

    ## Obviously, the following could be generalized ...

    if(op == "*") {
	if(!is.object(e1)) {
	    e3 <- e2
	    e2 <- e1
	    e1 <- e3 
	}
        if(!is.object(e2)) {
            if(length(e2) == 1L) {
		if(!is.finite(e2))
		    return(as.simple_triplet_matrix(as.matrix(e1) * e2))
		names(e2) <- NULL
		e1$v <- e1$v * e2
		if(!e2)
		    e1 <- .reduce(e1)
                return(e1)
            }
            if(length(e2) == e1$nrow) {
		names(e2) <- NULL
		pos <- which(!is.finite(e2))
		if(length(pos)) {
		    ## replace with dense rows
		    ind <- match(e1$i, pos, nomatch = 0L) == 0L
		    e1$v <- c(e1$v[ind],
                              as.matrix(e1[pos, ]))
		    e1$i <- c(e1$i[ind],
                              rep.int(pos, e1$ncol))
		    e1$j <- c(e1$j[ind],
                              rep(seq_len(e1$ncol), each = length(pos)))
		}
                e1$v <- e1$v * e2[e1$i]
		if(any(!e2))
		    e1 <- .reduce(e1)
                ## Could add something like
                ##    if(is.null(e1$dimnames) &&
                ##       !is.null(nms <- names(e2))) {
                ##        e1$dimnames <- list(nms, NULL)
                ##    }
                ## but then multiplying a matrix and a vector does not
                ## seem to do this either ...
                return(e1)
            }
	    if(is.matrix(e2)) {
		if(!all(dim(e2) == c(e1$nrow, e1$ncol)))
		    stop("Incompatible dimensions.")
		pos <- which(!is.finite(e2))
		if(length(pos)) {
		    ## add zeros
		    pos <- pos[match(pos, e1$i + (e1$j - 1L) * e1$nrow,
				     nomatch = 0L) == 0L] - 1L
		    if(length(pos)) {
			e1$v <- c(e1$v, vector(typeof(e1$v), length(pos)))
			e1$i <- c(e1$i, pos %%  e1$nrow + 1L)
			e1$j <- c(e1$j, pos %/% e1$nrow + 1L)
		    }
		}
		e1$v <- e1$v * e2[cbind(e1$i, e1$j)]
		if (any(!e2))
		    e1 <- .reduce(e1)
                e1$dimnames <- .make_dimnames(e1, e2)
		return(e1)
	    }
	    stop("Not implemented.")
        }
        ## This leaves multiplying two simple triplet matrices.
        e1 <- as.simple_triplet_matrix(e1)
        e2 <- as.simple_triplet_matrix(e2)
        ## Check dimensions: currently, no recycling.
        if(((nr <- e1$nrow) != e2$nrow) || ((nc <- e1$ncol) != e2$ncol))
            stop("Incompatible dimensions.")
        if(length(e1$v) < length(e2$v)) {
            ## Swap e1 and e2 so that duplicated indices can be found
            ## more efficiently.
            e3 <- e1
            e1 <- e2
            e2 <- e3
        }
        ## Find duplicated indices.
        ## pos <- match(paste(e2$i, e2$j, sep = "\r"),
        ##              paste(e1$i, e1$j, sep = "\r"),
        ##              nomatch = 0L)
        pos <- .Call(R_match_matrix, cbind(e1$i, e1$j),
                                     cbind(e2$i, e2$j), 0L)[[2L]]
        ind <- which(pos > 0L)
        if(!all(is.finite(e1$v)) || !all(is.finite(e2$v))) {
	    ## Augment and reduce
	    e2$i <- c(e2$i[ind], e2$i[-ind], e1$i[-pos])
	    e2$j <- c(e2$j[ind], e2$j[-ind], e1$j[-pos])
	    e2$v <- c(e2$v[ind] * e1$v[pos],
                      vector(typeof(e2$v), 1L) * c(e2$v[-ind], e1$v[-pos]))
	    e2$dimnames <- .make_dimnames(e1, e2)
	    return(.reduce(e2))
	} else
	    return(simple_triplet_matrix(e2$i[ind], e2$j[ind],
                                         e2$v[ind] * e1$v[pos],
                                         nr, nc, .make_dimnames(e1, e2)))
    }

    ## This is slightly inefficent but special value handling is already
    ## in place.  Note v / 0 = v * 0^(-1) = v * Inf.
    if(op == "/") {
        if(!is.object(e2))
	    return(e1 * e2^(-1))
	e2 <- as.matrix(e2)
	if (!is.object(e1))
	    return(as.simple_triplet_matrix(e1 * e2^(-1)))
	return(e1 * e2^(-1))
    }

    ## This leaves adding and subtracting two simple triplet matrices.
    e1 <- as.simple_triplet_matrix(e1)
    e2 <- if(op == "+")
        as.simple_triplet_matrix(e2)
    else
        as.simple_triplet_matrix(-e2)
    ## Check dimensions: currently, no recycling.
    if((e1$nrow != e2$nrow) || (e1$ncol != e2$ncol))
        stop("Incompatible dimensions.")
    if(length(e1$v) < length(e2$v)) {
        ## Swap e1 and e2 so that duplicated indices can be found more
        ## efficiently.
        e3 <- e1
        e1 <- e2
        e2 <- e3
    }
    ## Find duplicated indices.
    ## pos <- match(paste(e2$i, e2$j, sep = "\r"),
    ##              paste(e1$i, e1$j, sep = "\r"),
    ##              nomatch = 0L)
    pos <- .Call(R_match_matrix, cbind(e1$i, e1$j),
                                 cbind(e2$i, e2$j), 0L)[[2L]]
    ind <- which(pos == 0L)
    ## Notice 0 + special value = special value.
    e1$v[pos] <- e1$v[pos] + e2$v[pos > 0L]
    e1$i <- c(e1$i, e2$i[ind])
    e1$j <- c(e1$j, e2$j[ind])
    e1$v <- c(e1$v, e2$v[ind])
    e1$dimnames <- .make_dimnames(e1, e2)
    .reduce(e1)
}

Summary.simple_triplet_matrix <-
function(..., na.rm = FALSE)
{
    v <- unlist(lapply(list(...),
                       function(e) {
                           v <- as.simple_triplet_matrix(e)$v
                           if(length(v) < prod(dim(e)))
                               v <- c(v, vector(typeof(v), 1L))
                           v
                       }),
                recursive = FALSE)
    do.call(.Generic, list(v, na.rm = na.rm))
}

dim.simple_triplet_matrix <-
function(x)
    c(x$nrow, x$ncol)

`dim<-.simple_triplet_matrix` <-
function(x, value)
{
    value <- as.integer(value)
    if((length(value) != 2L) || any(is.na(value)))
        stop("invalid dim replacement value")
    nr <- x$nrow
    nc <- x$ncol
    if(prod(value) != nr * nc)
        stop("invalid dim replacement value")

    pos <- nr * (x$j - 1L) + x$i - 1L

    nr <- value[1L]
    nc <- value[2L]

    x$i <- pos %%  nr + 1L
    x$j <- pos %/% nr + 1L
    x$nrow <- nr
    x$ncol <- nc
    x$dimnames <- NULL

    x
}

dimnames.simple_triplet_matrix <-
function(x)
    x$dimnames

`dimnames<-.simple_sparse_array` <-
`dimnames<-.simple_triplet_matrix` <- 
function(x, value)
{
    if(!is.null(value)) {
        ## NOTE that if length(value) < length(dim(x)) we
	##      have to assume that the dimensions with index
	##      seq_len(length(value)) are to be set. For 
	##	example, we are called with a list of length
	##      one if we call dimnames(x)[[1L]] <- value and
	##      dimnames(x) == NULL (because of [[<-)
	##
        if(!is.list(value) || length(value) > length(dim(x)))
            stop("Invalid dimnames.")
        if(!length(value))
            value <- NULL
        else {
            dnx <- vector("list", length(dim(x)))
	    len <- sapply(value, length)
	    ind <- which(len > 0L)
	    if (any(len[ind] != dim(x)[ind]))
		stop("Invalid component length.")
            dnx[ind] <- lapply(value[ind], as.character)
	    if (!is.null(names(value))) {
		ind <- seq_len(length(value))
		names(dnx)[ind] <- names(value)
	    }
        }
    }
    ## See the constructor (above).
    if(is.null(value))
        x$dimnames <- NULL
    else
        x$dimnames <- dnx
    x
}

`[.simple_triplet_matrix` <-
function(x, i, j, drop = FALSE)
{
    ## (Well, we certainly don't drop ...)

    ## (See e.g. `[.data.frame` for the trickeries of subscript methods:
    ## e.g., 
    ##   x[i = sample.int(nr, k), , drop = FALSE]
    ## counts 4 arguments (x, i, j and drop) where j is missing ...

    na <- nargs() - !missing(drop)
    if((na == 1L) ||
       (na == 2L) && missing(i) ||
       (na == 3L) && missing(i) && missing(j))
        return(x)
    
    nr <- x$nrow
    nc <- x$ncol

    if(na == 2L) {
        ## Single index subscripting.
        ## Mimic subscripting matrices: no named argument handling in
        ## this case.
        if(is.character(i))
            out <- vector(typeof(x$v))[rep.int(NA, length(i))]
        else if(!is.matrix(i)) {
            if(is.logical(i)) {
		if(nr * nc > 16777216L)
		  stop("Logical vector subscripting disabled for this object.")
                i <- which(rep(i, length.out = nr * nc))
	    }
            else if(!is.numeric(i))
                stop(gettextf("Invalid subscript type: %s.",
                              typeof(i)),
                     domain = NA)
	    else
		## 52-bit safe
		if(nr * nc > 4503599627370496)
		  stop("Numeric vector subscripting disabled for this object.")
	    ## Shortcut
	    if(!length(i))
		return(vector(mode = typeof(x$v), length = 0L))
	    if(is.double(i))
		i <- trunc(i)
            ## Let's hope we have a vector.
            ## What if we have both negatives and positives?
            if(all(i >= 0, na.rm = TRUE)) {
                i <- i[i > 0]
                out <- vector(mode = typeof(x$v), length = length(i))
		if(length(out)) {
		    is.na(i) <- i > nr * nc
		    is.na(out) <- is.na(i)
		    i <- match(i, (x$j - 1L) * nr + x$i, 0L)
		    out[i > 0L] <- x$v[i]
		}
            } else if(!any(is.na(i)) && all(i <= 0)) {
		if(nr * nc > 16777216L)
		  stop("Negative vector subscripting disabled for this object.")
                out <- vector(mode = typeof(x$v), nr * nc)
                out[(x$j - 1L) * nr + x$i] <- x$v
                out <- out[i]
            }
            else stop("Cannot mix positive and negative subscripts.")
        }
        else {
	    ## Shortcut
	    if(!nrow(i))
		return(vector(mode = typeof(x$v), length = 0L))
	    ## Ignore dimensions
	    if(ncol(i) != 2L)
		return(do.call("[.simple_triplet_matrix",
			       list(x = x, as.vector(i))))

            ## Note that negative values are not allowed in a matrix
            ## subscript.
	    if(is.double(i))
		i <- trunc(i)
            if(any(i < 0, na.rm = TRUE))
                stop("Invalid subscript.")
            ## Rows containing zero indices can be dropped.
            ## Rows with NA indices should give NA (at least for
            ## non-recursive x).
	    k <- .Call(R_all_row, i > 0, FALSE)
            i <- i[k, ,drop = FALSE]
            out <- vector(mode = typeof(x$v), length = nrow(i))
	    if(length(out)) {
		if (any(i > rep(c(nr, nc), each = nrow(i)), na.rm = TRUE))
		    stop("subscript out of bounds")
		k <- k[k]
		is.na(out) <- is.na(k)
		rm(k)
		## See duplicated.matrix
		## pos <- match(paste(i[, 1L], i[, 2L], sep = "\r"),
		##	     paste(x$i, x$j, sep = "\r"),
		##	     nomatch = 0L)
		storage.mode(i) <- "integer"
		i <- .Call(R_match_matrix, cbind(x$i, x$j), i, 0L)[[2L]]
		out[i > 0L] <- x$v[i]
	    }
        }
    }
    else {
        ## Two index subscripting is rather tricky, as it can also be
        ## used for rearranging and "recycling" rows and columns.  Let
        ## us not support the latter for now, so that selected rows and
        ## columns must be unique.
        pos <- NULL
        if(!missing(i)) {
	    if(any(is.na(i)))
		stop("NA indices not allowed.")
            pi <- seq_len(nr)
            if(is.logical(i)) {
                i <- rep(i, length.out = nr)
                nr <- sum(i)
                pos <- i[x$i]
            } else {
                if(is.character(i)) {
                    i <- match(i, rownames(x))
                    if(any(is.na(i)))
                        stop("Subscript out of bounds.")
                    if(any(duplicated(i)))
                        stop("Repeated indices currently not allowed.")
                } else if(is.numeric(i)) {
		    if(is.double(i))
			i <- trunc(i)
                    if(all(i >= 0)) {
                        i <- i[i > 0]
			if(any(i > nr))
			    stop("subscript out of bounds")
                        if(any(duplicated(i)))
                            stop("Repeated indices currently not allowed.")
                    } else if(all(i <= 0))
                        i <- pi[i]
                    else
                        stop("Cannot mix positive and negative subscripts.")
                } else {
                    stop(gettextf("Invalid subscript type: %s.",
                                  typeof(i)),
                         domain = NA)
                }
                nr <- length(i)
                pos <- match(x$i, i, 0L) > 0L
            }
            pi[i] <- seq_len(nr)
        }
        if(!missing(j)) {
	    if(any(is.na(j)))
		stop("NA indices not allowed.")
            pj <- seq_len(nc)
            if(is.logical(j)) {
                j <- rep(j, length.out = nc)
                nc <- sum(j)
                pos <- if(is.null(pos))
                    j[x$j]
                else
                    j[x$j] & pos
            } else {
                if(is.character(j)) {
                    j <- match(j, colnames(x))
                    if(any(is.na(j)))
                        stop("Subscript out of bounds.")
                    if(any(duplicated(j)))
                        stop("Repeated indices currently not allowed.")
                } else if(is.numeric(j)) {
		    if(is.double(j))
			j <- trunc(j)
                    if(all(j >= 0)) {
                        j <- j[j > 0]
			if(any(j > nc))
			    stop("subscript out of bounds")
                        if(any(duplicated(j)))
                            stop("Repeated indices currently not allowed.")
                    } else if(all(j <= 0))
                        j <- pj[j]
                    else
                        stop("Cannot mix positive and negative subscripts.")
                } else {
                    stop(gettextf("Invalid subscript type: %s.",
                                  typeof(j)),
                         domain = NA)
                }
                nc <- length(j)
                pos <- if(is.null(pos))
                    (match(x$j, j, 0L) > 0L)
                else
                    (match(x$j, j, 0L) > 0L) & pos
            }
            pj[j] <- seq_len(nc)
        }

        if(!is.null(dnx <- x$dimnames)) {
	    if (!missing(i)) {
		dnx[1L] <- list(dnx[[1L]][i])
		if (!length(dnx[[1L]]))
		    dnx[1L] <- list(NULL)
	    }
	    if (!missing(j)) {
		dnx[2L] <- list(dnx[[2L]][j])
		if (!length(dnx[[2L]]))
		    dnx[2L] <- list(NULL)
	    }
	    if (!length(dnx[[1L]]) && !length(dnx[[2L]]))
		dnx <- NULL
	}

        i <- if(missing(i)) x$i[pos] else pi[x$i[pos]]
        j <- if(missing(j)) x$j[pos] else pj[x$j[pos]]
        out <- simple_triplet_matrix(i, j, x$v[pos], nr, nc, dnx)
    }

    out
}

rbind.simple_triplet_matrix <-
function(..., deparse.level = 1L)
{
    args <- lapply(Filter(Negate(is.null), list(...)),
                   as.simple_triplet_matrix)
    ## Ignore 'deparse.level' ...
    out <- Reduce(function(x, y) {
        if((nc <- ncol(x)) != ncol(y))
            stop("Numbers of columns of matrices must match.")
        nr <- nrow(x)
        simple_triplet_matrix(c(x$i, y$i + nr),
                              c(x$j, y$j),
                              c(x$v, y$v),
                              nrow = nr + nrow(y), ncol = nc)
    }, args)
    ## Handle dimnames in one final step.
    rnms <- lapply(args, rownames)
    rnms <- if(!all(sapply(rnms, is.null))) {
	rnms <- mapply(function(rnm, n)
	    if(is.null(rnm))
		rep.int("", n)
	    else
		rnm,
	    rnms,
	    lapply(args, nrow),
	    SIMPLIFY = FALSE
	)
	
        do.call("c", rnms)
    }
    else
        NULL
    cnms <- Find(Negate(is.null), lapply(args, colnames))
    dimnames(out) <- list(rnms, cnms)
    out
}

cbind.simple_triplet_matrix <-
function(..., deparse.level = 1L)
{
    args <- lapply(Filter(Negate(is.null), list(...)),
                   as.simple_triplet_matrix)
    ## Ignore 'deparse.level' ...
    out <- Reduce(function(x, y) {
        if((nr <- nrow(x)) != nrow(y))
            stop("Numbers of rows of matrices must match.")
        nc <- ncol(x)
        simple_triplet_matrix(c(x$i, y$i),
                              c(x$j, y$j + nc),
                              c(x$v, y$v),
                              nrow = nr, ncol = nc + ncol(y))
    }, args)
    ## Handle dimnames in one final step.
    cnms <- lapply(args, colnames)
    cnms <- if(!all(sapply(cnms, is.null))) {
	cnms <- mapply(function(cnm, n)
	    if(is.null(cnm))
		rep.int("", n)
	    else
		cnm,
	    cnms,
	    lapply(args, ncol),
	    SIMPLIFY = FALSE
	)
        do.call("c", cnms)
    }
    else
        NULL
    rnms <- Find(Negate(is.null), lapply(args, rownames))
    dimnames(out) <- list(rnms, cnms)
    out
}

t.simple_triplet_matrix <-
function(x)
    simple_triplet_matrix(x$j, x$i, x$v, x$ncol, x$nrow, rev(x$dimnames))

duplicated.simple_triplet_matrix <-
function(x, incomparables = FALSE, MARGIN = 1L, fromLast = FALSE, ...)
{
    ## We could use the duplicated method for class matrix, but at the
    ## expense of going from sparse to dense ...
    if(!is.logical(incomparables) || incomparables)
        .NotYetUsed("incomparables != FALSE")
    if(MARGIN == 1L) {
        i <- x$i
        j <- x$j
        len <- x$nrow
    } else if(MARGIN == 2L) {
        i <- x$j
        j <- x$i
        len <- x$ncol
    } else
        stop("Invalid margin.")
    o <- order(i, j)
    y <- split(paste(j[o], x$v[o], sep = "\r"), i[o])
    tmp <- character(len)
    names(tmp) <- seq_along(tmp)
    tmp[names(y)] <- sapply(y, paste, collapse = "\r")
    duplicated(tmp, fromLast = fromLast)
}

unique.simple_triplet_matrix <-
function(x, incomparables = FALSE, MARGIN = 1L, fromLast = FALSE, ...)
{
    if(!is.logical(incomparables) || incomparables)
        .NotYetUsed("incomparables != FALSE")
    ind <- !duplicated(x, MARGIN = MARGIN, fromLast = fromLast)
    if(MARGIN == 1L)
        x[which(ind), ]
    else
        x[, which(ind)]
}

c.simple_triplet_matrix <-
function(..., recursive = FALSE)
{
    args <- list(...)
    ind <- which(sapply(args, inherits, "simple_triplet_matrix"))
    args[ind] <-
        lapply(args[ind],
               function(x) {
                   y <- vector(typeof(x$v), x$nrow * x$ncol)
                   y[x$i + (x$j - 1L) * x$nrow] <- x$v
                   y
               })
    do.call("c", args)
}

print.simple_triplet_matrix <-
function(x, ...)
{
    writeLines(gettextf("A %s simple triplet matrix.",
                        paste(dim(x), collapse = "x")))
    invisible(x)
}

mean.simple_triplet_matrix <-
function(x, ...)
{
    sum(x$v) / prod(dim(x))
}

aperm.simple_triplet_matrix <-
function(a, perm = NULL, ...)
{
    s <- c(1L, 2L)
    if(!is.null(perm)) {
        perm <- if(is.character(perm))
            match(perm, names(a$dimnames))
        else if(is.numeric(perm))
            match(perm, s)
        else NULL
        if(length(perm) != length(s) || any(is.na(perm)))
            stop("Invalid permutation.")
        if(all(perm == s))
            return(a)
    }
    ## Transpose.
    t.simple_triplet_matrix(a)
}

as.vector.simple_triplet_matrix <-
function(x, mode = "any")
    as.vector(as.matrix(x), mode)

## Utilities for creating special simple triplet matrices:

simple_triplet_zero_matrix <-
function(nrow, ncol = nrow, mode = "double")
    simple_triplet_matrix(integer(), integer(), vector(mode, 0L),
                          nrow, ncol)

simple_triplet_diag_matrix <-
function(v, nrow = length(v))
{
    v <- rep(v, length.out = nrow)
    i <- seq_len(nrow)
    simple_triplet_matrix(i, i, v, nrow, nrow)
}

