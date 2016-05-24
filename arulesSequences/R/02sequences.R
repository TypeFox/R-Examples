
## class "sequences" represents a set of ordered sets of
## elements <X1,X2,...>, i.e. itemsets Xi.
##
## note that observed sequences should be represented by
## class "transactions". hence the class name serves as
## a blocker.
##
## ceeboo 2007, 2008, 2014, 2015

setClass("sequences",
    representation(
        elements     = "itemsets",
        data         = "sgCMatrix",
        sequenceInfo = "data.frame",
	tidLists     = "tidLists_or_NULL"
    ),
    contains = "associations",

    prototype(quality = data.frame(), sequenceInfo = data.frame()),

    validity = function(object) {
        if (dim(object@elements@items)[1] != dim(object@data)[1])
            return("slots 'elements' and 'data' do not conform")
        if (length(object@quality) &&
            dim(object@quality)[1] != dim(object@data)[2])
            return("slots 'data' and 'quality' do not conform")
        if (length(object@sequenceInfo) &&
            dim(object@sequenceInfo)[1] != dim(object@data)[2])
            return("slots 'data' and 'sequenceInfo' do not conform")
	if (length(object@tidLists) &&
	    dim(object@tidLists)[1] != dim(object@data)[2])
	    return("slots 'data' and 'tidLists' do not conform")

        TRUE
    }
)

# no method for itemsets

setMethod("nitems", signature(x = "sequences"), 
    function(x, itemsets = FALSE) {
        i <- .Call(R_rowSums_sgCMatrix, x@data) > 0
        if (!itemsets) {
            i <- .Call(R_colSubset_ngCMatrix, x@elements@items@data, i)
            i <- .Call(R_rowSums_ngCMatrix, i) > 0
        }
        sum(i)
    }
)

setMethod("dim", signature(x = "sequences"),
    function(x) c(x@data@Dim[2], x@elements@items@data@Dim))

setMethod("length", signature(x = "sequences"),
    function(x) dim(x@data)[2])

setMethod("size", signature(x = "sequences"),
    function(x, type = c("size","itemsets","length","items")) {
        type <- match.arg(type)
        switch(type, 
            size = .Call(R_colSums_ngCMatrix, x@data),
            itemsets = {
                s <- .Call(R_asList_ngCMatrix, x@data, NULL)
                names(s) <- NULL
                sapply(s, function(x) length(unique(x)))
            },
            length = {
                # fixme: inefficient
                s <- size(x@elements)
                s <- .Call(R_asList_ngCMatrix, x@data, s)
                names(s) <- NULL
                sapply(s, sum)
            },
            items = {
                # fixme: inefficient
                s <- .Call(R_asList_ngCMatrix, x@elements@items@data, NULL)
                s <- .Call(R_asList_ngCMatrix, x@data, s)
                names(s) <- NULL
                sapply(s, function(x) length(unique(unlist(x))))
            }
        )
    }
) 

setGeneric("ritems",
    function(x, ...) standardGeneric("ritems"))

setMethod("ritems", signature(x = "sequences"),
    function(x, type = c("min","max"), itemsets = FALSE) {
        type = match.arg(type)

        if (itemsets)
            s <- .Call(R_asList_ngCMatrix, x@data, NULL)
        else {
            s <- .Call(R_asList_ngCMatrix, x@elements@items@data, NULL)
            s <- .Call(R_asList_ngCMatrix, x@data, s)
            s <- lapply(s, unlist)
        }
        names(s) <- NULL
        f <- switch(type,
                min = function(x) min(x),
                max = function(x) max(x)
             )
        sapply(s, function(x) f(table(x)))
    }
)

#

setMethod("itemFrequency", signature(x = "sequences"),
    function(x, itemsets = FALSE, type = c("absolute", "relative")) {
        type <- match.arg(type)
        if (itemsets) 
            s <- .Call(R_rowSums_sgCMatrix, x@data)
        else {
            s <- .Call(R_asList_ngCMatrix, x@elements@items@data, NULL)
            s <- .Call(R_asList_ngCMatrix, x@data, s)
            s <- lapply(s, function(x) unique(unlist(x)))
            s <- factor(unlist(s), levels = 
                        seq_len(x@elements@items@data@Dim[1]))
            s <- as.vector(table(s))
        }
        switch(type, "absolute" = s,
                     "relative" = s / x@data@Dim[2])
    }
)

setGeneric("itemTable",
    function(x, ...) standardGeneric("itemTable"))

setMethod("itemTable", signature(x = "sequences"),
    function(x, itemsets = FALSE) {
        if (itemsets)
            t <- .Call(R_asList_ngCMatrix, x@data, NULL)
        else {
            t <- .Call(R_asList_ngCMatrix, x@elements@items@data, NULL)
            t <- .Call(R_asList_ngCMatrix, x@data, t)
            t <- lapply(t, unlist, recursive = FALSE)
        }
        t <- lapply(t, table)
        t <- unlist(t)
        t <- table(names(t), counts = t)
        names(dimnames(t))[1] <- 
            if (itemsets) "itemsets" else "items"
        ## convenience
        k <- sort(as.integer(rownames(t)))
        t <- t[k,, drop = FALSE]
        t
    }
)

## reduce is experimental

setMethod("[", signature(x = "sequences", i = "ANY", j = "ANY", drop = "ANY"),
    function(x, i, j, ..., reduce = FALSE, drop) {
        if (!missing(j))
            stop("incorrect number of dimensions (j not possible)")
        if (!missing(i)) {
            if (length(x@quality))
                x@quality <- x@quality[i,, drop = FALSE]
            if (length(x@sequenceInfo))
                x@sequenceInfo <- x@sequenceInfo[i,, drop = FALSE]
	    if (length(x@tidLists))
		x@tidLists <- x@tidLists[i,, drop = FALSE]
            x@data <- x@data[,i]
        }
        if (reduce) {
            i <- .Call(R_rowSums_sgCMatrix, x@data) > 0
            x@data <- x@data[i,]
            x@elements <- x@elements[i]
            i <- .Call(R_rowSums_ngCMatrix, x@elements@items@data) > 0
            x@elements@items <- x@elements@items[,i]
        }
        validObject(x, complete = TRUE)
        x
    }
)

setAs("list", "sequences",
    function(from) {
        if (!length(from))
            return(new("sequences"))
        ## FIXME we cannot provide a flatten option as
        ##       further arguments are not available
        if (!all(unlist(lapply(from, sapply, is.atomic))))
            stop("item(s) not atomic")

        i <- unlist(from, recursive = FALSE, use.names = FALSE)
        ## FIX Matrix mess
        i <- lapply(i, sort)
        p <- sapply(i, length)
        p <- cumsum(p)
        i <- unclass(factor(unlist(i)))
        l <- attr(i, "levels")
        attr(i, "levels") <- NULL
        i <- new("ngCMatrix", p   = c(0L, p),
                              i   = i - 1L, 
                              Dim = c(length(l), length(p)))

        s <- .Call(R_pnindex, i, NULL, FALSE)

        s <- .Call(R_colSubset_ngCMatrix, i, !duplicated(s))

        e <- new("itemMatrix", data     = s,
                               itemInfo = data.frame(labels = l,
						     stringsAsFactors = FALSE))

        e <- new("itemsets", items = e)

        i <- .Call(R_pnindex, s, i, FALSE)
        
        p <- sapply(from, length)
        names(p) <- NULL
        p <- cumsum(p)

        s <- new("sgCMatrix", p   = c(0L, p),
                              i   = i - 1L,
                              Dim = c(s@Dim[2], length(p)))

        new("sequences", elements     = e,
                         data         = s,
                         sequenceInfo = data.frame(sequenceID =
			    if (!is.null(names(from)))
				names(from),
			    stringsAsFactors = FALSE
			 ))
    }
)

setAs("sequences", "list",
    function(from) LIST(from, decode = TRUE))

setMethod("LIST", signature(from = "sequences"),
    function(from, decode = TRUE) {
        d <- if (decode) 
                 as.character(from@elements@items@itemInfo[['labels']])
             else 
                 NULL

        i <- .Call(R_asList_ngCMatrix, from@elements@items@data, d)
        i <- .Call(R_asList_ngCMatrix, from@data, i)
        if (decode)
            names(i) <- from@sequenceInfo[['sequenceID']]
        i
    }
)

#

.formatElements <- function(x, itemSep = ",", setStart = "{", setEnd = "}")
    paste(setStart, paste(x, collapse = itemSep), setEnd,
          collapse = "", sep = "")


setMethod("labels", signature(object = "sequences"),
    function(object, setSep = ",", seqStart = "<", seqEnd = ">", decode = TRUE, ...) {
        object <- LIST(object, decode)
        names(object) <- NULL
        sapply(lapply(object, lapply, .formatElements, ...), .formatElements,
               setSep, seqStart, seqEnd)
    }
)

setMethod("itemLabels", signature(object = "sequences"),
    function(object, itemsets = FALSE, ...) {
        l <- itemLabels(object@elements)
        if (itemsets) {
            l <- .Call(R_asList_ngCMatrix, object@elements@items@data, l)
            l <- sapply(l, .formatElements, ...)
        }
        l
    }
)

setReplaceMethod("itemLabels", signature(object = "sequences"),
    function(object, value) {
        itemLabels(object@elements@items) <- value
        object
    }
)

setAs("sequences", "data.frame",
    function(from) {
        if (!length(from))
            return(data.frame())
        d <- data.frame(sequence = labels(from))
        if (length(from@quality))
            d <- cbind(d, from@quality)
        if (length(from@sequenceInfo))
            d <- cbind(d, from@sequenceInfo)
        d
    }
)

##

setMethod("inspect", signature(x = "sequences"),
    function(x, setSep = ",", seqStart = "<", seqEnd = ">", decode = TRUE) {
        if (!length(x))
            return()

        q <- quality(x)
        x <- LIST(x, decode)

        i <- sapply(unlist(x, recursive = FALSE), length)
        i <- cumsum(i)
        
        s <- sapply(x, function(x) sum(sapply(x, length)))
        s <- cumsum(s)

        x <- unlist(x)
        
        tmp <- matrix(" ", 5, length(x))
        
        tmp[3,] <- x
    
        tmp[4,]   <- ","            # item separator
        tmp[4, i] <- "}"            # itemset closing delimiter
        tmp[5, i] <- setSep         # itemset seperator
        tmp[5, s] <- seqEnd         # sequence closing delimiter

        i  <- c(1, i[-length(i)]+1)
        o  <- c(1, s[-length(s)]+1)
    
        tmp[2, i] <- "{"            # itemset opening delimiter
        tmp[1, o] <- seqStart       # sequence opening delimiter

        tmp[1,] <- format(tmp[1,])

        tmp <- apply(tmp, 2, paste, collapse = "")
        
        out <- matrix("", 3+length(q), length(x)+1)
        
        out[2,] <- format(c("items", tmp), justify = "left")
        
        # fixme: original sequence numbers
        
        o <- o + 1                  # header

        out[1, o] <- seq(length(o))
        out[1,]   <- format(out[1,], justify = "right")

        s <- s + 1                  # header

        for (i in seq(q)) {
            tmp <- format(q[[i]], justify = "right")
            out[2+i, c(1, s)] <-
                format(c(names(q)[i], tmp), justify = "right")
        }
           
        out[2+length(q)+1,] <- "\n"
        
        cat("", out, "\n")
    }
)

##

setGeneric("sequenceInfo",
    function(object, ...) standardGeneric("sequenceInfo"))

setMethod("sequenceInfo", signature(object = "sequences"),
    function(object) object@sequenceInfo)

setGeneric("sequenceInfo<-",
    function(object, value) standardGeneric("sequenceInfo<-"))

setReplaceMethod("sequenceInfo", signature(object = "sequences"),
    function(object, value) {
        object@sequenceInfo <- value
        validObject(object)
        object
    }
)

setMethod("itemInfo", signature(object = "sequences"),
    function(object) itemInfo(object@elements))

setReplaceMethod("itemInfo", signature(object = "sequences"),
    function(object, value) {
        itemInfo(object@elements@items) <- value
        object
    }
)

setGeneric("itemsets",
    function(x, ...) standardGeneric("itemsets"))

setMethod("itemsets", signature(x = "sequences"),
    function(x) x@elements)

# fixme: prototype?

setClass("summary.sequences",
    representation(
        items    = "integer",
        elements = "integer",
        sizes    = "table",
        lengths  = "table",
	tidLists = "logical"
    ),
    contains = "summary.associations",

    prototype(
	length  = 0L, 
	sizes   = table(NULL), 
	lengths = table(NULL), 
	quality = table(NULL)
    )
)

setMethod("summary", signature(object = "sequences"),
    function(object, maxsum = 6) {
        if (!length(object))
            return(new("summary.sequences", 
		       ## see cspade
		       info     = object@info,
		       tidLists = FALSE))

        maxsum <- max(0, maxsum-1)

        i <- itemFrequency(object, type = "absolute")
        names(i) <- itemLabels(object)
        i <- sort(i, decreasing = TRUE)
        i <- c(head(i, maxsum), sum(tail(i, length(i)-maxsum)))
        names(i)[length(i)] <- "(Other)"

        s <- itemFrequency(object, itemsets = TRUE, type = "absolute")
        names(s) <- labels(object@elements, itemSep = ",")
        s <- sort(s, decreasing = TRUE)
        s <- c(head(s, maxsum), sum(tail(s, length(s)-maxsum)))
        names(s)[length(s)] <- "(Other)"

        q <- if (length(object@quality))
                 summary(object@quality) 
             else 
                 summary(NULL)

        new("summary.sequences", 
		length   = length(object),
                items    = i,
                elements = s,
                sizes    = table(sizes = size(object, "size")),
                lengths  = table(lengths = size(object, "length")),
                quality  = q,
                info     = object@info,
		tidLists = !is.null(object@tidLists))
    }
)

setMethod("show", signature(object = "summary.sequences"),
    function(object) {
        cat("set of", object@length, "sequences with\n")
        cat("\nmost frequent items:\n")
        print(object@items)

        cat("\nmost frequent elements:\n")
        print(object@elements)

        cat("\nelement (sequence) size distribution:\n")
        print(object@sizes)

        cat("\nsequence length distribution:\n")
        print(object@lengths)

        cat("\nsummary of quality measures:\n")
        print(object@quality)

	cat("\nincludes transaction ID lists:", object@tidLists, "\n")           
        if (length(object@info)) {
            info <- object@info
            if (is.language(info$data)) info$data <- deparse(info$data)
            cat("\nmining info:\n")
            print(data.frame(info, row.names = ""))
        }
        invisible(NULL)
    }
)

##
setMethod("is.closed", signature(x = "sequences"),
    function(x) {
        support <- quality(x)$support
	if (is.null(support))
	    stop("'x' does not contain support information")
	if (any(duplicated(x)))
	    stop("'x' not unique")
	if (FALSE) {
	    if (any(is.na(support)))
		stop("missing values not implemented")
	    m <- is.subset(x)
	    if (!all(m@x))
		stop("reduce not implemented")
	    m@x <- support[m@i + 1L] <=
	           support[rep(seq_len(length(m@p) - 1L), diff(m@p))]
	    m <- selectMethod("rowSums", class(m))(m) == 1L
	    names(m) <- NULL
	    m
	} else
	    .Call(R_pnsclosed, x@data, x@elements@items@data,
		  rank(support, na.last = "keep", ties.method = "min"),
		  FALSE)
    }
)

##

setMethod("is.maximal", signature(x = "sequences"),
    function(x) {
        u <- unique(x)
	if (FALSE) { 
	    m <- is.subset(u)
	    m <- selectMethod("rowSums", class(m))(m) == 1L
	    names(m) <- NULL
	} else 
	    m <- .Call(R_pnscount, u@data, u@data, 
		       u@elements@items@data, FALSE) == 1L
        i <- match(x, u)
        m[i]
    }
)

##

setMethod("duplicated", signature(x = "sequences"),
    function(x, incomparables = FALSE)  {
        i <- .Call(R_pnindex, x@data, NULL, FALSE)
        duplicated(i, incomparables)
    }
)

setMethod("unique", signature(x = "sequences"),
    function(x, incomparables = FALSE) x[!duplicated(x)])

## note that %in% does no longer dispatch to match unless
## the right hand operand is not of the same class

setMethod("%in%", signature(x = "sequences", table = "sequences"),
    function(x, table)
	match(x, table, nomatch = 0L) > 0L
)

setMethod("match", signature(x = "sequences", table = "sequences"),
    function(x, table, nomatch = NA_integer_, incomparables = NULL) {
        k <- match(x@elements, table@elements)
        n <- which(is.na(k))
        if (length(n)) {
            k[n] <- table@data@Dim[1] + seq(length(n))
            table@data@Dim[1] <- table@data@Dim[1] + length(n)
        }
        if (any(k != seq_len(length(k))))
            x@data <- .Call(R_recode_ngCMatrix, x@data, k)
        if (x@data@Dim[1] <  table@data@Dim[1])
            x@data@Dim[1] <- table@data@Dim[1]
        i <- .Call(R_pnindex, table@data, x@data, FALSE)
        match(i, seq(length(table)), nomatch =nomatch, 
                                     incomparables = incomparables)
    }
)

## the semantics are wrong with respect to
## the usual behavior.

setMethod("%in%", signature(x = "sequences", table = "character"),
    function(x, table) {
        k <- x@elements@items %in% table
        .Call(R_colSums_ngCMatrix, x@data[k,]) > 0
    }
)

setMethod("%pin%", signature(x = "sequences", table = "character"),
    function(x, table) {
        if (length(table) > 1)
            stop("'table' contains more than one item label pattern")
        k <- x@elements@items %pin% table
        .Call(R_colSums_ngCMatrix, x@data[k,]) > 0
    }
)

setMethod("%ain%", signature(x = "sequences", table = "character"),
    function(x, table) {
        p <- match(table, itemLabels(x@elements@items))
        if (any(is.na(p)))
            stop("table contains an unknown item label")
        p <- size(x@elements@items[,p])
        k <- p > 0
        p <- .Call(R_asList_ngCMatrix, x@data[k,], p[k])
        sapply(p, sum) >= length(table)
    }
)

setGeneric("%ein%",
    function(x, table) standardGeneric("%ein%"))

setMethod("%ein%", signature(x = "sequences", table = "character"),
    function(x, table) {
        p <- x@elements@items %ain% table
        p <- .Call(R_asList_ngCMatrix, x@data, p)
        sapply(p, any)
    }
)

##

.combineMeta <- function(x, y, name) {
    if (length(slot(x, name))) {
        if (length(slot(y, name)))
            slot(x, name) <- rbind(slot(x, name), slot(y, name))
        else {
            k <- rbind(NA, slot(x, name))
            slot(x, name) <-
                rbind(slot(x, name), k[rep(1, length(y)),, drop = FALSE])
        }
    } else
    if (length(slot(y, name))) {
        k <- rbind(NA, slot(y, name))
        slot(x, name) <-
            rbind(k[rep(1, length(x)),, drop = FALSE],  slot(y, name))
    }
    x
}

setMethod("c", signature(x = "sequences"),
    function(x, ..., recursive = FALSE) {
        args <- list(...)
        if (recursive)
            args <- unlist(args)
        for (y in args) {
            if (!is(y, "sequences"))
                stop("can only combine sequences")
            info <- y@info
            if (length(info)) {
                k <- match(names(info), names(x@info))
                k <- mapply(identical, info, x@info[k])
                info <- info[k]
            }
            x@info <- info
            x <- .combineMeta(x, y, "quality")
            x <- .combineMeta(x, y, "sequenceInfo")
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
            x@data <- .Call(R_cbind_ngCMatrix, x@data, y@data)
	    ## see arules
	    x@tidLists <- NULL
        }
        validObject(x, complete = TRUE)
        x
    }
)

## avoid uneccessary variables.

setMethod("subset", signature(x = "sequences"),
    function(x, subset) {
        if (missing(subset))
            return(x)
        i <- eval(substitute(subset),
                  envir = c(x@quality, x@sequenceInfo))
        x[i]
    }
)

setMethod("tidLists", signature(x = "sequences"),
    function(x) x@tidLists
)


###
