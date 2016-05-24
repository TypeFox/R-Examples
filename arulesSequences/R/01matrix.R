##
## sgCMatrix is designed to work with ptree but 
## otherwise does not deserve the qualifier Matrix.
##
## ceeboo 2007, 2008, 2016

# fixme: special cases

setClass("sgCMatrix",
    representation(
        p        = "integer",
        i        = "integer",
        Dim      = "integer",
        Dimnames = "list",
        factors  = "list"
    ),

    prototype(p        = integer(1), 
              i        = integer(0), 
              Dim      = integer(2), 
              Dimnames = list(NULL, NULL),
              factors  = list()),

    validity = function(object) .Call(R_valid_sgCMatrix, object)
)

setMethod("dim", signature(x = "sgCMatrix"),
    function(x) x@Dim)

setMethod("dimnames", signature(x = "sgCMatrix"),
    function(x) x@Dimnames)

## fixme: does not validate
setReplaceMethod("dimnames", signature(x="sgCMatrix"),
    function(x, value) {
        if (is.null(value))
            x@Dimnames <- list(NULL, NULL)
        else {
            if (!is.list(value) || length(value) != 2)
                stop("value must be list of length 2")

            value <- lapply(value, unlist)

            if ((l <- length(value[[1]])) && l != x@Dim[1])
                stop("length of value [1] invalid")
            if ((l <- length(value[[2]])) && l != x@Dim[2])
                stop("length of value [2] invalid")
    
            x@Dimnames <- lapply(value, 
		function(value) 
		    if (!is.null(value)) 
			as.character(value)
		    else
			value
	    )
        }
        x
    }
)

setMethod("[", signature(x = "sgCMatrix", i = "ANY", j = "ANY", drop = "ANY"),
    function(x, i, j, ..., drop) {
        if (!missing(i))
            x <- .Call(R_rowSubset_sgCMatrix, x, i) 
        if (missing(j))
            return(x)
        
        .Call(R_colSubset_ngCMatrix, x, j) 
    }
)

#

setAs("list", "sgCMatrix",
    function(from) {
        if (!length(from))
            return(new("sgCMatrix"))
        ## flatten non-atomic elements
        from <- lapply(from, lapply, paste, collapse=",")
        p <- sapply(from, length)
        names(p) <- NULL
        p <- cumsum(p)
        i <- factor(unlist(from, use.names = FALSE))
        new("sgCMatrix", p        = c(0L, p),
                         i        = c(i) - 1L,
                         Dim      = c(length(levels(i)), length(p)),
                         Dimnames = list(levels(i), names(from)))
    }
)

setAs("sgCMatrix", "list",
    function(from) {
        i <- .Call(R_asList_ngCMatrix, from, from@Dimnames[[1]])
        names(i) <- from@Dimnames[[2]]
        i
    }
)

setAs("ngCMatrix", "sgCMatrix",
    function(from) {
        class(from) <- "sgCMatrix"
        from
    }
)

setMethod("show", signature(object = "sgCMatrix"),
    function(object) {
        cat(object@Dim[1], "x", object@Dim[2],
            "sparse pseudo Matrix of class", class(object), "\n")
        ## fixme
        invisible(NULL)
    }
)

###
