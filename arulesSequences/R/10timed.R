
## extends class squences for timing information
## and provides methods for analysis of it. note
## that the class is more storage-efficient than
## class transactions but the latter can contain
## additional event information.
##
## ceeboo 2007, 2008, 2014, 2015

setClass("timedsequences",
    representation(
        time     = "ngCMatrix",
        timeInfo = "data.frame"
    ),
    contains = "sequences",

    prototype(timeInfo = data.frame()),

    validity = function(object) {
        if (object@data@Dim[2] != object@time@Dim[2])
            stop("slots 'time' and 'data' do not conform")
        if (length(object@timeInfo) && 
            length(object@timeInfo[[1]]) != object@time@Dim[1])
            return("slots 'time' and 'timeInfo' do not conform")

        TRUE
    }
)

## fixme: due to a bug in ngCMatrix-class duplicate
##        row indexes pass the validity test.
setAs("transactions", "timedsequences",
    function(from) {
	if (!length(from))
	    return(new("timedsequences"))
        if (!all(c("sequenceID", "eventID") %in% names(transactionInfo(from))))
            stop("transactionInfo: missing sequenceID or eventID")
        if (any(is.na(transactionInfo(from)[['sequenceID']])))
            stop("sequenceID: missing values")
        if (any(is.na(transactionInfo(from)[['eventID']])))
            stop("eventID: missing values")

        t <- .as_integer(transactionInfo(from)[['eventID']])
        if (is.factor(t)) {
            warning("'eventID' is a factor")

            k <- order(transactionInfo(from)[['sequenceID']])
            if (any(k != seq_len(length(k)))) {
                warning("transactions not ordered")

                from <- from[k]
                t <- t[k]
            }

            i <- data.frame(labels  = levels(t), 
                            eventID = seq_len(length(levels(t))),
			    stringsAsFactors = FALSE)
        } else {
            k <- order(transactionInfo(from)[['sequenceID']], t)
            if (any(k != seq_len(length(k)))) {
                warning("transactions not ordered")

                from <- from[k]
                t <- t[k]
            }

            t <- factor(t)
            i <- data.frame(labels  = levels(t),
                            eventID = as.integer(levels(t)),
			    stringsAsFactors = FALSE)
        }

        s <- as(from , "itemMatrix")

        e <- new("itemsets", items = unique(s))

        s <- .Call(R_pnindex, e@items@data, s@data,  FALSE)
        s <- tapply(s, transactionInfo(from)[['sequenceID']], list)

        n <- data.frame(sequenceID = names(s),
			stringsAsFactors = FALSE)
        dim(s) <- NULL
        
        p <- cumsum(sapply(s, length))
        s <- unlist(s)

        s <- new("sgCMatrix", p   = c(0L, p),
                              i   = s - 1L,
                              Dim = c(length(e), length(p)))

        s <- new("sequences", data         = s,
                              elements     = e,
                              sequenceInfo = n)

        t <- new("ngCMatrix", p   = s@data@p,
                              i   = c(t) - 1L,
                              Dim = c(length(levels(t)), length(s)))

        new("timedsequences", s,
                              time     = t,
                              timeInfo = i)
    }
)

setMethod("dim", signature(x = "timedsequences"),
    function(x) c(x@data@Dim[2], x@elements@items@data@Dim, x@time@Dim[1])) 

setGeneric("timeInfo",
    function(object, ...) standardGeneric("timeInfo"))

setMethod("timeInfo", signature(object = "timedsequences"),
    function(object) object@timeInfo)

setGeneric("timeInfo<-", 
    function(object, value) standardGeneric("timeInfo<-"))

setReplaceMethod("timeInfo", signature(object = "timedsequences"),
    function(object, value) {
        object@timeInfo <- value
        validObject(object)
        object
    }
)

setGeneric("times",
    function(x, ...) standardGeneric("times"))

.timefun <- function(type)
    switch(type,
        times  = function(x) x,
        gaps   = function(x) diff(x),
        mingap = function(x) if (length(x <- diff(x))) min(x) else x,
        maxgap = function(x) if (length(x <- diff(x))) max(x) else x,
        span   = function(x) diff(range(x)),

        stop("type not implemented")
    )

## note (1) returns NA for sequences of length one
##      (2) simplifies if possible

setMethod("times", signature(x = "timedsequences"),
    function(x, type = c("times", "gaps", "mingap", "maxgap", "span")) {
        type <- match.arg(type)

        t <- .Call(R_asList_ngCMatrix, x@time, x@timeInfo[['eventID']])
        t <- sapply(t, .timefun(type))
        if (is.list(t)) {
            is.na(t) <- !sapply(t, length)
            t <- sapply(t, eval)
        }
        t
    }
)

setGeneric("timeFrequency", 
    function(x, ...) standardGeneric("timeFrequency"))

# note: (1) for span, mingap, and maxgap the total number is 
#       the number of sequences minus the number of sequences
#       with one element only. (2) for gaps the total number
#       is the total number of elements minus the number of
#       sequences. fixme: provide relative?

setMethod("timeFrequency", signature(x = "timedsequences"),
    function(x, type = c("times", "gaps", "mingap", "maxgap", "span")) {
        type = match.arg(type)

        if (type == "times")
            t <- NULL
        else
            t <- x@timeInfo[['eventID']]
        t <- .Call(R_asList_ngCMatrix, x@time, t)
        t <- lapply(t, .timefun(type))
        t <- table(unlist(t))
        if (type == "times")
            names(t) <- x@timeInfo[['labels']][as.integer(names(t))]
        t
    }
)

setGeneric("timeTable",
    function(x, ...) standardGeneric("timeTable"))

setMethod("timeTable", signature(x = "timedsequences"),
    function(x, type = c("times","gaps", "mingap", "maxgap", "span"), itemsets = FALSE) {
        type <- match.arg(type)

        if (itemsets) 
            i <- .Call(R_asList_ngCMatrix, x@data, NULL)
        else {
            i <- .Call(R_asList_ngCMatrix, x@elements@items@data, NULL)
            i <- .Call(R_asList_ngCMatrix, x@data, i)
        }
        if (type == "times")
            t <- NULL
        else
            t <- x@timeInfo[['eventID']]
        t <- .Call(R_asList_ngCMatrix, x@time, t)

        f <- .timefun(type)

        t <- mapply(function(t, i) {
            t <- rep(t, sapply(i, length))
            tapply(t, unlist(i), f, simplify = FALSE)
        }, t, i, USE.NAMES = FALSE, SIMPLIFY = FALSE)
        t <- unlist(t, recursive = FALSE)
        i <- rep(names(t), sapply(t, length))

        t <- table(i, unlist(t))
        if (type == "times") 
            colnames(t) <- x@timeInfo[['labels']][as.integer(colnames(t))]
        names(dimnames(t))[1] <- if (itemsets) "itemsets" else "items"
        names(dimnames(t))[2] <- type
        t
    }
)

setMethod("LIST", signature(from = "timedsequences"),
    function(from, decode = TRUE) {
        d <- if (decode)
                 as.character(from@elements@items@itemInfo[['labels']])
             else
                 NULL

        i <- .Call(R_asList_ngCMatrix, from@elements@items@data, d)
        i <- .Call(R_asList_ngCMatrix, from@data, i)
        t <- .Call(R_asList_ngCMatrix, from@time, from@timeInfo[['eventID']])
        i <- mapply(function(t, i) {
            names(i) <- t
            i
        }, t, i, SIMPLIFY = FALSE)
        if (decode)
            names(i) <- from@sequenceInfo[['sequenceID']]
        i
    }
)

setMethod("labels", signature(object = "timedsequences"),
    function(object, timeStart = "[", timeEnd = "]", setSep = ",", seqStart = "<", seqEnd = ">", decode = TRUE, ...) {
        object <- LIST(object, decode)
        names(object) <- NULL
        sapply(lapply(object, lapply, .formatElements, ...), function(x)
               paste(seqStart, paste(x, timeStart, names(x), timeEnd, sep = "",
                     collapse = setSep), seqEnd, sep = "", collapse = ""))
    }
)

## fixme
setMethod("inspect", signature(x = "timedsequences"),
    function(x, setSep = ",", seqStart = "<", seqEnd = ">", decode = TRUE) {
        if (!length(x))
            return()

        if (!all(s <- size(x))) {
            if (!any(s))
                return()
            x <- x[s > 0]
        }

        t <- .Call(R_asList_ngCMatrix, x@time, x@timeInfo[['labels']])
        t <- unlist(t)

        p <- sequenceInfo(x)

        q <- quality(x)
        x <- LIST(x, decode)

        i <- sapply(unlist(x, recursive = FALSE, use.names = FALSE), length)
        i <- cumsum(i)
        
        s <- sapply(x, function(x) sum(sapply(x, length)))
        s <- cumsum(s)

        x <- unlist(x, use.names = FALSE)
        
        tmp <- matrix(" ", 5, length(x))
        
        tmp[3,] <- x
    
        tmp[4,]   <- ","            # item separator
        tmp[4, i] <- "}"            # itemset closing delimiter
        tmp[5, i] <- setSep         # itemset seperator
        tmp[5, s] <- seqEnd         # sequence closing delimiter

        j <- c(1, i[-length(i)]+1)
        o <- c(1, s[-length(s)]+1)
    
        tmp[2, j] <- "{"            # itemset opening delimiter
        tmp[1, o] <- seqStart       # sequence opening delimiter


        tmp[1,] <- format(tmp[1,])

        tmp <- apply(tmp, 2, paste, collapse = "")
        
        out <- matrix("", 4+length(p)+length(q), length(x)+1)
        
        out[2,] <- format(c("items", tmp), justify = "left")
      
        i <- i + 1                  # header

        out[3, c(1,i)] <- c("eventID", t)
        out[3,]        <- format(out[3,], justify = "right")

        # fixme: original sequence numbers
        
        o <- o + 1                  # header

        out[1, o] <- seq(length(s))
        out[1,]   <- format(out[1,], justify = "right")

        s <- s + 1                  # header

        for (i in seq(p)) {
            tmp <- format(p[[i]], justify = "right")
            out[3+i, c(1, s)] <-
                format(c(names(p)[i], tmp), justify = "right")
        }

        for (i in seq(q)) {
            tmp <- format(q[[i]], justify = "right")
            out[3+length(p)+i, c(1, s)] <-
                format(c(names(q)[i], tmp), justify = "right")
        }
           
        out[3+length(p)+length(q)+1,] <- "\n"
        
        cat("", out, "\n")
    }
)

setClass("summary.timedsequences", 
    representation(
        times        = "integer",
        sequenceInfo = "data.frame",
        itemInfo     = "data.frame",
        timeInfo     = "data.frame"
    ),
    contains = "summary.sequences"
)

setMethod("summary", signature(object = "timedsequences"),
    function(object, maxsum = 6) {
        if (!length(object))
            return(new("summary.timedsequences", length = 0L))

        m <- max(0, maxsum-1)

        t <- timeFrequency(object, type = "times")
        t <- sort(t, decreasing = TRUE)
        t <- c(head(t, m), "(Other)" = sum(tail(t, -m)))
 
        new("summary.timedsequences", times        = t,
                                      sequenceInfo = head(sequenceInfo(object), 3),
                                      itemInfo     = head(itemInfo(object), 3),
                                      timeInfo     = head(timeInfo(object), 3),
        ## FIXME R-2.7.0 bug
           selectMethod("summary","sequences")(as(object, "sequences"), maxsum))
    }
)

setMethod("show", signature(object = "summary.timedsequences"),
function(object) {
        cat("set of", object@length, "timedsequences with\n")
        if (object@length) {
            cat("\nmost frequent times:\n")
            print(object@times)

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

            if (length(timeInfo)) {
                cat("\nincludes extended time information - examples:\n")
                print(object@timeInfo)
            }

            if (length(object@itemInfo)) {
                cat("\nincludes extended item information - examples:\n")
                print(object@itemInfo)
            }

            if (length(sequenceInfo)) {
                cat("\nincludes extended sequence information - examples:\n")
                print(object@sequenceInfo)
            }
        }
        invisible(NULL)
    }
)

## fixme
setMethod("[", signature(x = "timedsequences", i = "ANY", j = "ANY", drop = "ANY"),
    function(x, i, j, k, ..., reduce = FALSE, drop) {
        if (!missing(k)) {
            if (is.character(k))
                k <- x@timeInfo[['labels']] %in% k
            else
                k <- x@timeInfo[['labels']] %in% x@timeInfo[['labels']][k]

            y <- x
            x@time <- .Call(R_rowSubset_ngCMatrix, x@time, k)
            x@timeInfo <- x@timeInfo[k,, drop = FALSE]

            k <- y@time@i %in% (which(k) - 1L)

            x@data@i <- x@data@i[k]
            x@data@p <- x@time@p
            validObject(x, complete = TRUE)
        }
        if (!missing(j))
            stop("incorrect number of dimensions (j not possible)")
        if (!missing(i)) 
            x <- new("timedsequences", 
                     as(x, "sequences")[i, reduce = reduce],
                     time     = .Call(R_colSubset_ngCMatrix, x@time, i),
                     timeInfo = x@timeInfo)
        if (reduce) {
            if (missing(i))
                x <- new("timedsequences", 
                         as(x, "sequences")[reduce = reduce],
                         time     = x@time,
                         timeInfo = x@timeInfo)
            i <- .Call(R_rowSums_ngCMatrix, x@time) > 0
            x@time <- .Call(R_rowSubset_ngCMatrix, x@time, i)
            x@timeInfo <- x@timeInfo[i,, drop = FALSE]
            validObject(x)
        }
        x
    }
)

## fixme: sorting in natural order
setMethod("c", signature(x = "timedsequences"),
    function(x, ..., recursive = FALSE) {
        args <- list(...)
        if (recursive)
            args <- unlist(args)
        for (y in args) {
            if (!inherits(y, "timedsequences"))
                stop("can only combine timedsequences")

            u <- unique(c(x@timeInfo[['labels']], y@timeInfo[['labels']]))
            u <- .as_integer(u)

            if (is.factor(u))
                x@timeInfo <- data.frame(labels  = as.character(u),
                                         eventID = seq_len(length(u)),
					 stringsAsFactors = FALSE)
            else {
                u <- sort(u)
                x@timeInfo <- data.frame(labels  = as.character(u),
                                         eventID = u,
					 stringsAsFactors = FALSE)

                k <- match(x@timeInfo[['labels']], u)
                if (any(diff(k) < 0))
                    stop("'x' invalid time order")
                if (any(k != seq(length(k))))
                    x@time <- .Call(R_recode_ngCMatrix, x@time, k)
                if (x@time@Dim[1] <  length(u))
                    x@time@Dim[1] <- length(u)
            }

            k <- match(y@timeInfo[['labels']], u)
            if (any(diff(k) < 0))
                stop("'y' invalid time order")
            if (any(k != seq(length(k)))) 
                y@time <- .Call(R_recode_ngCMatrix, y@time, k)
            if (y@time@Dim[1] <  length(u))
                y@time@Dim[1] <- length(u)

            x@time <- .Call(R_cbind_ngCMatrix, x@time, y@time)
            validObject(x@time)

            x <- new("timedsequences", c(as(x, "sequences"),
                                         as(y, "sequences")),
                     time     = x@time,
                     timeInfo = x@timeInfo)
        }
	x
    }
)

##

setGeneric("timesets",
    function(object, ...) standardGeneric("timesets"))

setMethod("timesets", signature(object = "timedsequences"),
    function(object)
        new("itemMatrix", data        = object@time,
                          itemInfo    = object@timeInfo,
                          itemsetInfo = object@sequenceInfo)
)

setAs("timedsequences", "transactions",
    function(from) {
        i <- from@data@i + 1L
        i <- from@elements@items[i]

        t <- from@time@i + 1L
        t <- data.frame(
	    sequenceID = .as_integer(rep(from@sequenceInfo[['sequenceID']],
					 size(from))),
            eventID    = .as_integer(from@timeInfo[['labels']][t]))

        i <- new("transactions", as(i, "itemMatrix"))
	transactionInfo(i) <- t
	i
    }
)

##

setGeneric("firstOrder",
    function(x, ...) standardGeneric("firstOrder"))

## fixme
setMethod("firstOrder", signature(x = "timedsequences"),
    function(x, times = FALSE) {
        if (times) {
            t <- .Call(R_firstOrder_sgCMatrix, x@time)
            rownames(t) <- colnames(t) <- x@timeInfo[['labels']]
            t
        } else 
            .Call(R_firstOrder_sgCMatrix, x@data)
    }
)

###
