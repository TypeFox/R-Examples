## A simple class for sparse arrays.

## Not very useful yet: need at least a subscript method.
## (Unfortunately, additional methods such as for rowSums/colSums or
## apply, etc., are not straightforward to add in an S3 world ...)

simple_sparse_array <-
function(i, v, dim = NULL, dimnames = NULL)
{
    ## See examples
    storage.mode(i) <- "integer"
    if (!is.matrix(i))
	dim(i) <- c(length(i), 1L)
    ## <FIXME>
    ## Add some sanity checking eventually ...
    ## i should be a matrix of indices (non-"zero" entries).
    ## v should be a "vector" of non-zero values, with length equal to
    ## the number of rows of i.
    ## </FIXME>
    if(is.null(dim)) dim <- if(NROW(i)) apply(i, 2L, max) else c(0L, 0L)
    ## <FIXME>
    ## Add checks for dimnames: should be NULL or a list of entries
    ## which are either NULL or character vectors as long as the
    ## corresponding dim.
    ## </FIXME>
    ssa <- list(i = i, v = v, dim = as.integer(dim), dimnames = dimnames)
    class(ssa) <- "simple_sparse_array"
    ## Note that this should never be true as it implies that either
    ## the class is wrong or the container is malformed.
    if (!.Call(R__valid_ssa, ssa))
        stop("failed to create a valid 'simple_sparse_array' object")
    ssa
}

as.simple_sparse_array <-
function(x)
    UseMethod("as.simple_sparse_array")

as.simple_sparse_array.simple_sparse_array <- identity

as.simple_sparse_array.array <-
function(x)
{
    x <- unclass(x)
    dx <- dim(x)
    if(!prod(dx))
	return(simple_sparse_array(matrix(integer(), 0L, length(dx)), 
				   c(x), dx, dimnames(x)))
    ind <- which(is.na(x) | (x != vector(typeof(x), 1L)), arr.ind = TRUE)
    dimnames(ind) <- NULL
    simple_sparse_array(ind, x[ind], dx, dimnames(x))
}

as.simple_sparse_array.matrix <- as.simple_sparse_array.array

as.simple_sparse_array.simple_triplet_matrix <-
function(x) 
    simple_sparse_array(cbind(x$i, x$j), x$v, c(x$nrow, x$ncol), dimnames(x))

as.simple_sparse_array.default <-
function(x)
    as.simple_sparse_array(unclass(as.array(x)))

as.array.simple_sparse_array <-
function(x, ...)
{
    v <- x$v
    dim <- x$dim
    y <- array(vector(typeof(v), 1L), dim = dim,
               dimnames = x$dimnames)
    y[x$i] <- v
    y
}

is.simple_sparse_array <-
function(x)
    inherits(x, "simple_sparse_array")

Math.simple_sparse_array <-
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

Summary.simple_sparse_array <-
function(..., na.rm = FALSE)
{
    v <- unlist(lapply(list(...),
                       function(e) {
                           v <- as.simple_sparse_array(e)$v
                           if(length(v) < prod(dim(e)))
                               v <- c(v, vector(typeof(v), 1L))
                           v
                       }),
                recursive = FALSE)
    do.call(.Generic, list(v, na.rm = na.rm))
}

dim.simple_sparse_array <-
function(x)
    x$dim

`dim<-.simple_sparse_array` <-
function(x, value)
{
    value <- as.integer(value)
    if(!length(value) || any(is.na(value)))
        stop("invalid dim replacement value")
    dx <- dim(x)
    if(prod(value) != prod(dx))
        stop("invalid dim replacement value")

    x$i <- arrayInd(.Call(R_vector_index, x$dim, x$i), value)
    x$dim <- value
    x$dimnames <- NULL
    
    x
}

dimnames.simple_sparse_array <-
function(x)
    x$dimnames

## FIXME we now have drop_simple_sparse_array
`[.simple_sparse_array` <-
function(x, ...)
{
    ## Note that calling x[] with a simple sparse array x will call the
    ## subscript method with args x and missing ...
    na <- nargs()
    if((na == 1L) || (na == 2L) && missing(..1))
        return(x)

    nd <- length(x$dim)

    ## Note there is a limit to representing integer numbers as 
    ## doubles.
    spos <- function(i) {
	if(!nrow(i)) 
	    return(vector(mode = typeof(i), length = 0L))
        ## Scalar positions of array index matrices i in the usual row
        ## major ordering of arrays.
	if(ncol(i) > 1L) {
	    ## This may not work on systems with BLAS issues
	    ## as.vector(tcrossprod(c(1L, cumprod(x$dim[-nd])), i - 1L)) + 1L
	    1L + row_sums((i - 1L) * rep(c(1L, cumprod(x$dim)[-nd]), 
					 each = nrow(i)))
	} else
	    as.vector(i)
    }

    if(na == 2L) {
        i <- ..1
        ## Single index subscripting.
        if(is.logical(i))
            stop("Logical subscripting currently not implemented.")
        else if(is.character(i))
            stop("Character subscripting currently not implemented.")
        else if(!is.matrix(i)) {
	    ## 52-bit safe
	    if(prod(x$dim) > 4503599627370496)
	      stop("Numeric vector subscripting disabled for this object.")
	     ## Shortcut
	     if(!length(i)) 
		return(vector(mode = typeof(x$v), length = 0L))
            ## Let's hope we have a vector.
            ## What if we have both negatives and positives?
	    if(is.double(i))
		i <- trunc(i)
            if(all(i >= 0, na.rm = TRUE)) {
                i <- i[i > 0]
                out <- vector(mode = typeof(x$v), length = length(i))
		if(length(out)) {
		    ## Missing values.
		    is.na(i) <- i > prod(x$dim)
		    is.na(out) <- is.na(i)
		    i <- match(i, spos(x$i), 0L)
		    out[i > 0L] <- x$v[i]
		}
            } else if(!any(is.na(i)) && all(i <= 0)) {
		if(prod(x$dim) > 16777216L)
		  stop("Negative vector subsripting disabled for this object.")
                out <- vector(mode = typeof(x$v), prod(x$dim))
                out[spos(x$i)] <- x$v
		## NOTE this fails if NAs are introduced by 
		##	coercion to integer. 
                out <- out[i]
            }
            else stop("Cannot mix positive and negative subscripts.")
        }
        else {
	     ## Shortcut
	     if(!nrow(i)) 
		return(vector(mode = typeof(x$v), length = 0L))
	     ## Ignore dimensions. 
	     if(ncol(i) != nd) 
		return(do.call("[.simple_sparse_array", 
			       list(x = x, as.vector(i))))
            ## Note that negative values are not allowed in a matrix
            ## subscript.
	    if(is.double(i))
		i <- trunc(i)
            if(any(i < 0, na.rm = TRUE))
                stop("Invalid subscript.")
	    k <- .Call(R_all_row, i > 0, FALSE)
            i <- i[k, , drop = FALSE]
            out <- vector(mode = typeof(x$v), length = nrow(i))
	    if(length(out)) {
		if(any(i > rep(x$dim, each = nrow(i)), na.rm = TRUE))
		    stop("subscript out of bounds")
		## Missing values.
		k <- k[k]
		is.na(out) <- is.na(k)
		rm(k)
		## This is not really the fastest way to match rows, but is
		## there an obvious better one?
		## pos <- match(split(i, row(i)), split(x$i, row(x$i)), 0L)
		storage.mode(i) <- "integer"
		i <- .Call(R_match_matrix, x$i, i, 0L)[[2L]]
		out[i > 0L] <- x$v[i]
	    }
        }
    }
    else {
        if(na != (nd + 1L))
            stop("Incorrect number of dimensions.")
        ## Get indices.
	args <- vector("list", na - 1L)
	for(k in seq_along(args)) {
	    n <- as.name(sprintf("..%i", k))
	    if (!do.call(missing, list(n)))
		args[[k]] <- eval(n)
	}
        ## Ready to go.
        dx <- x$dim
        pos <- rep.int(TRUE, length(x$v))
        ind <- vector("list", length = nd)
        for(k in seq_len(nd)) {
            i <- args[[k]]              # Given indices.
            if(is.null(i)) {
		ind[[k]] <- seq_len(dx[k])
		next
	    }
            else if(!is.numeric(i))
                stop("Only numeric multi-index subscripting is implemented.")
            else {
		if (any(is.na(i)))
		    stop("NA indices currently not allowed")
		if(is.double(i))
		    i <- trunc(i)
                if(all(i >= 0)) {
                    i <- i[i > 0]
                    if(any(duplicated(i)))
                        stop("Repeated indices currently not allowed.")
		    if(any(i > dx[k]))
			stop("subscript out of bounds")
                } else if(all(i <= 0))
		    ## NOTE this fails if NAs are introduced by 
		    ##	    coercion to integer. 
                    i <- seq_len(dx[k])[i]
                else
                    stop("Cannot mix positive and negative subscripts.")
                ind[[k]] <- i
                dx[k] <- length(i)
                j <- match(x$i[, k], i, 0L)
                x$i[j > 0L, k] <- seq_along(i)[j]
                pos <- pos & (j > 0L)
            }
        }
        if(!is.null(dnx <- x$dimnames))
            dnx[] <- Map("[", dnx, ind)
        out <- simple_sparse_array(x$i[pos, , drop = FALSE], x$v[pos],
                                   dx, dnx)
    }

    out

}

## <TODO>
## Add duplicated and unique methods for simple sparse arrays along the
## lines of the corresponding methods for simple triplet matrices.
## </TODO>

print.simple_sparse_array <-
function(x, ...)
{
    writeLines(gettextf("A simple sparse array of dimension %s.",
                        paste(dim(x), collapse = "x")))
    invisible(x)
}

mean.simple_sparse_array <-
function(x, ...)
{
    sum(x$v) / prod(dim(x))
}

aperm.simple_sparse_array <-
function(a, perm = NULL, ...)
{
    s <- seq_along(a$dim)
    if(is.null(perm))
        perm <- rev(s)
    else {
        perm <- if(is.character(perm))
            match(perm, names(a$dimnames))
        else if(is.numeric(perm))
            match(perm, s)
        else NULL
        if(length(perm) != length(s) || any(is.na(perm)))
            stop("Invalid permutation.")
    }
    simple_sparse_array(a$i[, perm, drop = FALSE], a$v,
                        a$dim[perm], a$dimnames[perm])
    
}

as.vector.simple_sparse_array <-
function(x, mode = "any")
    as.vector(as.array(x), mode)

simple_sparse_zero_array <-
function(dim, mode = "double")
{
    ld <- length(dim)
    if (!ld)
	stop("'dim' must have positive length")
    simple_sparse_array(matrix(integer(), 0L, ld), vector(mode, 0L), dim)
}

