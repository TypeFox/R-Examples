
## ceeboo 2008, 2014, 2015

## helper
as_sequences_transactions <-
function(from) {
    if ( inherits(from, "timedsequences"))
	return(as(from, "transactions"))
    if (!inherits(from, "sequences"))
	stop("'from' not of class sequences")

    ## use order indexes of events
    n <- diff(from@data@p)
    t <- lapply(n, seq_len)
    t <- unlist(t, use.names = FALSE)

    t <- data.frame(
	    sequenceID = rep(seq_along(n), n),
            eventID    = t)

    i <- from@data@i + 1L
    i <- from@elements@items[i]

    i <- new("transactions", as(i, "itemMatrix"))
    transactionInfo(i) <- t
    i
}

support.idlists <-
function(x, transactions, type = "absolute", parameter = NULL, verbose = FALSE) {
    if (!inherits(x, "sequences"))
	stop("'x' not of class sequences")
    if (!inherits(transactions, c("transactions", "sequences")))
        stop("'transactions' not of class transactions or sequences")

    if (is.null(parameter))
	parameter <- list()

    parameter <- as(parameter, "SPparameter")
    ## not used
    parameter@support <- NA_real_
    parameter@maxsize <- NA_integer_
    parameter@maxlen  <- NA_integer_

    ## if (length(parameter@maxwin))
    ##     stop("'maxwin' is not supported")

    if (verbose) {
	cat("using method:", "idlists", "\n")
	t0 <- proc.time()
	cat("\nparameter specification:\n")
	cat(.formatSP(parameter), sep = "\n")
	## FIXME control
	cat("\npreprocessing ... ")
    }

    y <- transactions
    if (inherits(y, "sequences"))
	y <- as_sequences_transactions(y)
    else
	if (!all(c("sequenceID", "eventID") %in% names(transactionInfo(y))))
	    stop("transactionInfo: missing 'sequenceID' and/or 'eventID'")

    ## L1 order optimization
    i <- itemFrequency(y, type = "absolute")
    i <- order(i)
    if (any(i != seq_along(i))) {
	y@data <-
	    .Call(R_recode_ngCMatrix, y@data, order(i))
	y@itemInfo <- 
	    data.frame(labels = y@itemInfo[['labels']][i],
		       stringsAsFactors = FALSE)
	if (verbose)
	    cat("L1 ")
    }
    ## conform
    k <- match(itemLabels(x@elements@items), itemLabels(y))
    n <- which(is.na(k))
    ## Prepend unmatched optimization
    if (length(n)) {
	i <- seq_along(i) + length(n)
	y@data <-
	    .Call(R_recode_ngCMatrix, y@data, i)
	y@itemInfo <-
	    data.frame(labels = c(x@elements@items@itemInfo[['labels']][n],
				  y@itemInfo[['labels']]),
		       stringsAsFactors = FALSE)
	k <- k + length(n)
	k[n] <- seq_along(n)
	if (verbose)
	    cat("P ")
    }
    if (any(k != seq_along(k))) {
	x@elements@items@data <-
	    .Call(R_recode_ngCMatrix, x@elements@items@data, k)
	x@elements@items@itemInfo <- y@itemInfo
    }
    rm(n)

    if (type == "tidLists")
	k <- labels(x)
    else
	rm(k)
    x <- LIST(x, decode = FALSE)

    sid <- .as_integer(transactionInfo(y)[['sequenceID']])
    eid <- .as_integer(transactionInfo(y)[['eventID']])
    if (is.factor(eid))
	warning("'eventID' is a factor")
	
    y <- selectMethod("t", class(y@data))(y@data)

    if (verbose) {
	t1 <- proc.time()
	cat("[", t1[1]-t0[1], "s]\n", sep = "")
    }

    x <- .Call(R_ilscount, x, y, sid, eid,
	if (length(parameter@mingap)) parameter@mingap else NULL,
	if (length(parameter@maxgap)) parameter@mingap else NULL,
	if (length(parameter@maxwin)) parameter@maxwin else NULL,
	verbose,
	(type == 'tidLists')
    )

    switch(type,
	relative = {
	    n <- length(
		if (is.factor(sid))
		    levels(sid)
		else
		    unique(sid)	## FIXME
	    )
	    x / n
	},
	absolute = x,
	tidLists = {
	    x <- as(x, "tidLists")
	    x@itemInfo <-
		data.frame(labels = k, stringsAsFactors = FALSE)
	    k <- 
		if (is.factor(sid))
		    levels(sid)
		else
		    as.character(unique(sid))
	    x@transactionInfo <- 
		data.frame(sequenceID = k, stringsAsFactors = FALSE)
	    x
	}
    )
}

support.ptree <- 
function(x, transactions, type = "absolute", verbose = FALSE) {
    if (!inherits(x, "sequences"))
	stop("'x' not of class sequences")
    if (!inherits(transactions, c("transactions", "sequences")))
        stop("'transactions' not of class transactions or sequences")

    if (verbose) {
        t0 <- proc.time()
        cat("preprocessing ... ")
    }

    ## FIXME this is inefficient as we use the
    ##       ordering information only.
    y <- transactions
    if (inherits(y, "transactions"))
	y <-  as(y, "timedsequences")

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
    rm(k, n)

    if (verbose) {
        t1 <- proc.time()
        cat("[", t1[1]-t0[1], "s]\n", sep = "")
    }

    x <- .Call("R_pnscount", x@data, y@data, 
	x@elements@items@data, 
	verbose)

    switch(type,
        relative = x / y@data@Dim[2],
        absolute = x
    )
}

setMethod("support", signature(x = "sequences"),
    function(x, transactions, type = c("relative", "absolute"), control = NULL)
    {
        type <- match.arg(type)
        verbose <- control[['verbose']]
	if (is.null(verbose))
	    verbose <- FALSE
	else
	if (!is.logical(verbose))
	    stop("'verbose' invalid range")

	parameter <- control[['parameter']]
	if (is.null(parameter))
	    support.ptree(x, transactions, type, verbose)
	else {
	    if (!is.list(parameter))
		stop("'parameter' invalid range")
	    support.idlists(x, transactions, type, parameter, verbose)
	}
    }
)


##
setMethod("supportingTransactions", signature(x = "sequences"), 
    function(x, transactions, ...) {
	if (!inherits(transactions, c("transactions", "sequences")))
            stop("'transactions' not of class transactions or sequences")
	y <- transactions
	if (inherits(y, "transactions"))
	     y <- as(y, "timedsequences")
        s <- is.subset(x, y)
        s <- new("ngCMatrix", 
            i = s@i, 
            p = s@p, 
            Dim = s@Dim 
        )
	## The implementation in arules is no longer
	## exported, thus R-2.7.0 bug.
	s <- selectMethod("t", class(s))(s)
        new("tidLists", 
            data = s,
            itemInfo = data.frame(labels = labels(x), stringsAsFactors = FALSE),
            transactionInfo = data.frame(sequenceID = 
		sequenceInfo(y)[['sequenceID']], stringsAsFactors = FALSE)
        )
    }
)

###
