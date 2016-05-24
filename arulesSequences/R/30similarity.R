
## ceeboo 2008, 2016

setGeneric("similarity",
    function(x, y = NULL, ...) standardGeneric("similarity")) 

## FIXME includes timed sequences
setMethod("similarity", signature(x = "sequences"),
    function(x, y = NULL, method = c("jaccard", "dice", "cosine", "subset"),
        strict = FALSE)
    {
        method <- match(match.arg(method), eval(formals()$method)) - 1L

        if (!is.null(y)) {
            if (!is(y, "sequences"))
                stop("'y' not of class sequences")
            ## conform
            k <- match(y@elements, x@elements)
            n <- which(is.na(k))
            if (length(n)) {
                k[n] <- length(x@elements) + seq(length(n))
                x@data@Dim[1] <- x@data@Dim[1] + length(n)
                if (!strict)
                    x@elements <- c(x@elements, y@elements[n])
            }
            if (any(k != seq_len(length(k))))
                y@data <- .Call(R_recode_ngCMatrix, y@data, k)
            if (y@data@Dim[1] <  x@data@Dim[1])
                y@data@Dim[1] <- x@data@Dim[1]

            dimnames(y@data) <- list(NULL, y@sequenceInfo[["sequenceID"]])
            y <- y@data
        }
        dimnames(x@data) <- list(NULL, x@sequenceInfo[["sequenceID"]])

        .Call(R_similarity_sgCMatrix,
              x@data, y, if (strict) NULL else x@elements@items@data, method)
    }
)

setMethod("is.subset", signature(x = "sequences"),
    function(x, y = NULL, proper = FALSE) {
	## inefficient
	if (FALSE) {
	    s <- similarity(x, y, method = "subset") > 0;
	    if (proper) {
		if (!is.null(y))
		    r <- similarity(y, x, method = "subset") > 0
		else 
		    r <- s
		if (FALSE) {
		    ## FIXME R-2.7.0 bug
		    s <- s & !selectMethod("t", class(r))(r)
		} else {
		    j <- rep(seq(0, length.out = length(s@p) - 1L), 
			     diff(s@p))
		    k <- 
		    match(
			paste(s@i, j),
			if (is.null(y))
			    paste(j, s@i)
			else
			    paste(rep(seq(0, length.out = length(r@p) - 1L),
				  diff(r@p)), r@i),
			nomatch = 0L
		    ) == 0L
		    s@p <- cumsum(tabulate(j[k] + 2L, nbins = length(s@p)))
		    s@i <- s@i[k]
		    s@x <- s@x[k]
		}
	    }
	    s
	} else {
	    if (!is.null(y)) {
		if (!is(y, "sequences"))
		    stop("'y' not of class sequences")
		## conform
		k <- match(y@elements, x@elements)
		n <- which(is.na(k))
		if (length(n)) {
		    k[n] <- length(x@elements) + seq(length(n))
		    x@data@Dim[1] <- x@data@Dim[1] + length(n)
		    x@elements <- c(x@elements, y@elements[n])
		}
		if (any(k != seq_len(length(k))))
		    y@data <- .Call(R_recode_ngCMatrix, y@data, k)
		if (y@data@Dim[1] <  x@data@Dim[1])
		    y@data@Dim[1] <- x@data@Dim[1]
		s <- .Call(R_pnssuperset, y@data, x@data,
			   x@elements@items@data, proper, FALSE)
	    } else
		s <- .Call(R_pnssuperset, x@data, NULL,
			   x@elements@items@data, proper, FALSE)
	    ## FIXME
	    s <- new("ngCMatrix",
		p = c(0L, cumsum(sapply(s, length))),
		i = unlist(s) - 1L,
		Dim = c(length(x),
			if (is.null(y))
			    length(x) 
			else 
			    length(y)
			),
		Dimnames = list(
		    x@sequenceInfo[["sequenceID"]],
		    if (is.null(y))
			x@sequenceInfo[["sequenceID"]]
		    else
			y@sequenceInfo[["sequenceID"]]
		)
	    )
	    as(s, "lgCMatrix")
	}
    }
)

setMethod("is.superset", signature(x = "sequences"),
    function(x, y = NULL, proper = FALSE) {
	if (!is.null(y))
	    s <- is.subset(y, x, proper)
	else
	    s <- is.subset(x, y, proper)
	selectMethod("t", class(s))(s)
    }
)

##

setAs("dsCMatrix", "dist",
    function(from)
        .Call(R_as_dist_dsCMatrix, if (from@uplo == "L") from else t(from)))

###
