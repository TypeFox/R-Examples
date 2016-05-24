
##
## data interfaces to cSPADE
##
## ceeboo 2007, 2008, 2012, 2014, 2015, 2016

.as_integer <- function(x) {
    ## preserve factor
    if (typeof(x) != "integer") {
        ## must be atomic
        x <- factor(x)
        l <- suppressWarnings(as.integer(levels(x)))
        ## implicit coercion
        if (!any(is.na(l)) && all(l == levels(x)))
            x <- l[c(x)]
    }
    x
}

read_baskets <- function(con, sep = "[ \t]+", info = NULL, iteminfo = NULL,
			      encoding = "unknown") {
    x <- readLines(con, encoding = encoding)
    x <- sub("^[ \t]+", "", x)
    x <- strsplit(x, split = sep)
    if (!is.null(info)) {
        i <- info
        info <- lapply(seq(length(info)), function(k) sapply(x, "[", k))
        names(info) <- i
        x <- lapply(x, "[", -seq(length(info)))
        # fixme: warning
        x <- lapply(x, unique)
    }
    x <- as(x, "transactions")
    if (!is.null(info)) {
	if (!is.null(info[['sequenceID']])) {
	    info[['sequenceID']] <- .as_integer(info[['sequenceID']])
	    if (is.integer(info[['sequenceID']]))
		if (any(info[['sequenceID']] < 1L))
		    warning("sequenceID not positive")
	}
	if (!is.null(info[['eventID']])) {
            info[['eventID']] <- .as_integer(info[['eventID']])
	    if (is.factor(info[['eventID']]))
		warning("'eventID' is a factor")
	    else
		if (any(info[['eventID']] < 1L))
		    warning("eventID not positive")
	    if (!is.null(info[['sequenceID']]))
		if (any(order(info[['sequenceID']], info[['eventID']]) != 
			seq_along(info[['sequenceID']])))
		    warning("'sequenceID' and/or 'eventID not ordered")
	}
	if (TRUE) {
	    i <- sapply(info, is.character)
	    info[i] <- lapply(info[i], type.convert)
	}
        transactionInfo(x) <- data.frame(info, stringsAsFactors = FALSE)
    }
    if (!is.null(iteminfo)) {
        if (!is.data.frame(iteminfo))
            stop("'iteminfo' not a data frame")
        labels <- itemLabels(x)
        if (!all(labels %in% rownames(iteminfo)))
            stop("the row names of 'iteminfo' do not match the item labels")
        iteminfo <- iteminfo[labels,, drop = FALSE]
        if ("labels" %in% names(iteminfo))
            iteminfo[['labels']] <- as.character(iteminfo[['labels']])
        else
            iteminfo <- cbind(x@itemInfo, iteminfo)
        itemInfo(x) <- iteminfo
    }
    x
}

## currently internal only

read_spade <- 
function(con = "", decode = FALSE, labels = NULL, transactions = NULL,
		    class = NULL) {
    if (con == "")
        con <- stdin()
    else 
    if (is.character(con)) {
        con <- file(con, "r")
        on.exit(close(con))
    }
    if (!inherits(con, "connection")) 
        stop("'con' must be a character string or connection.")

    n <- readLines(con, 1)
    if (!length(n))
        stop("the number of lines is zero")
    n <- as.integer(strsplit(n, " ")[[1]][5])
        
    x <- readLines(con)
    # control not implemented (see the -t option)
    if (FALSE) {
	k <- grep("^PRUNE", x)
	if (length(k))
	    x <- x[-k]
    }
    if (!length(x))
        return(new("sequences", info = list(nsequences = n)))

    x <- strsplit(x, split = " -- ")
   
    # NOTE 1) position 1 contains the support count.
    #      2) the following K positions contain the 
    #         support counts of a partition (see the 
    #         -c option).
    #      3) the following positions represent pairs
    #         of SID FRQ identifying the containing data 
    #         sequences and their support counts (see
    #         the -y option).

    c <- strsplit(sapply(x, "[", 2), split = " ")
    if (!is.null(transactions)) {
	k <- lapply(c, function(x, i)
		## see NOTE 3)
		x <- matrix(x[i], nrow = 2L)[1L, ],
		## see NOTE 1) + 2)
		-seq_len(max(1L, length(levels(class))) + 1L)
	)
	k <- as(k, "tidLists")
	s <- transactionInfo(k)[['labels']]
	t <- transactionInfo(transactions)[['sequenceID']]
	k@transactionInfo <- data.frame(sequenceID =
	    if (is.factor(t))
		   levels(t)[as.integer(s)]
	    else
		s,
	    stringsAsFactors = FALSE
	)
	transactions <- k
	rm(k, s, t)
    }
    c <- lapply(seq_len(length(levels(class)) + 1L), function(k)
	    as.integer(sapply(c, "[", k)))
   
    # split into a list of lists (sequences) each 
    # containing a vector of character (itemsets)
   
    x <- lapply(strsplit(sapply(x, "[", 1), split = " -> "), strsplit, " ")
    if (decode)
        x <- lapply(x, lapply, as.integer)

    if (!length(x))
        stop("the number of sequences parsed is zero")

    x <- as(x, "sequences")
    names(c) <- c("support", levels(class))
    c <- mapply("/", c, c(n, if (!is.null(class)) table(class)), 
		SIMPLIFY = FALSE)
    x@quality <- data.frame(c, check.names = FALSE)
    x@info <- list(nsequences = n)

    k <- which(size(x) == 1)
    if (length(k) == length(x@elements)) {
        i <- x@data[,k]@i + 1L
        k[i] <- k
        quality(x@elements) <- x@quality[k,, drop = FALSE]
    } else
        stop("the data is incomplete")

    if (!is.null(labels)) {
        k <- as.integer(as.character(x@elements@items@itemInfo[['labels']]))
        itemLabels(x@elements@items) <- as.character(labels[k])
    }
    if (!is.null(transactions)) {
	transactions@itemInfo <- data.frame(labels = 
	    labels(x), stringsAsFactors = FALSE)
	x@tidLists <- transactions
    }
    validObject(x)
    x
}

## write data in text format for later
## processing by bin/makebin

write_cspade <- function(x, con) {
    if (!inherits(x, "transactions"))
        stop("'x' not of class transactions")

    r <- .Call(R_asList_ngCMatrix, x@data, NULL)
    r <- sapply(r, paste, collapse = " ")
    
    sid <- .as_integer(transactionInfo(x)[['sequenceID']])
    if (is.integer(sid))
	if (any(sid < 1L))
	    stop("'sequenceID' not positive")
    eid <- .as_integer(transactionInfo(x)[['eventID']])
    if (is.factor(eid))
        warning("'eventID' is a factor")
    else
	if (any(eid < 1L))
	    stop("'eventID' not positive")
    if (any(order(sid, eid) != seq_along(sid)))
        stop("'sequenceID' and/or 'eventID' not ordered")

    r <- rbind(as.character(as.integer(sid)),
               as.character(as.integer(eid)),
               as.character(size(x)), r)
    r <- apply(r, 2, paste, collapse = " ")
    
    writeLines(r, con)
}

## Write class(ification) file (see option -c)
##
## NOTE the offsets must be the same as in the
##      asc file
write_class <- function(x, con) {
    if (!is.factor(x))
	stop("'x' not a factor")
    s <- .as_integer(names(x))
    if (!length(s))
	s <- seq_along(x)
    else
	names(x) <- NULL
    if (any(is.na(x) | is.na(s)))
	stop("'x' invalid")
    writeBin(con = con,
        object = c(
            length(levels(x)),  # number of classes
            c(rbind(
                s,              # SID
                c(x) - 1L       # CLASS
            ))
        )
    )
}

## write data directly in binary format for
## later processing by bin/exttpose

makebin <- function(x, file) {
    if (!inherits(x, "transactions"))
        stop("'x' not of class transactions")

    sid <- .as_integer(transactionInfo(x)[['sequenceID']])
    eid <- .as_integer(transactionInfo(x)[['eventID']])
    if (is.factor(eid))
        warning("'eventID' is a factor")

    x <- as(x, "ngCMatrix")
    attr(x, "sid") <- sid
    attr(x, "eid") <- eid

    .Call(R_makebin, x, file)
}

## cSPADE wrapper
##
## note that we assume 1MB = 2^10 x 2^10 = 4^10 for the 
## computation of the number of database partitions.

cspade <- 
function(data, parameter = NULL, control = NULL, tmpdir = tempdir()) {

    ## workaround
    if (.Platform$OS  == "windows" && 
        .Platform$GUI == "Rgui")
	system2 <- function(command, args = character(), stdout = "", ...) {
	    if (is.character(stdout) && nzchar(stdout)) {
		args   <- c(args, ">", stdout)
		stdout <- NULL
	    }
	    args <- c("/c", shQuote(command), args)
	    command <- Sys.getenv("COMSPEC")
	    ## bail out
	    if (!nzchar(command))
		stop("environment variable 'COMSPEC' not set")
	    base::system2(command, args = args, stdout = stdout, ...)
	}

    if (!inherits(data, "transactions"))
        stop("'data' not of class transactions")
    if (!all(c("sequenceID", "eventID") %in% names(transactionInfo(data))))
        stop("transactionInfo: missing 'sequenceID' and/or 'eventID'")
    ## optional
    class <- transactionInfo(data)[['classID']]
    if (!is.null(class)) {
	names(class) <- transactionInfo(data)[['sequenceID']]
	class <- class[!duplicated(names(class))]
	class <- factor(class)
    }
    if (!all(dim(data))) 
        return(new("sequences"))
    parameter <- as(parameter, "SPparameter")
    control   <- as(control ,  "SPcontrol")

    if (control@verbose) {
        t1 <- proc.time()
        cat("\nparameter specification:\n")
        cat(.formatSP(parameter), sep = "\n")
        cat("\nalgorithmic control:\n")
        cat(.formatSP(control), sep = "\n")
        cat("\npreprocessing ...")
    }

    exe <- "bin"
    if (.Platform$r_arch != "")
	exe <- file.path(exe, .Platform$r_arch)
    exe <- system.file(exe, package = "arulesSequences")

    file <- tempfile(pattern = "cspade", tmpdir)
    on.exit(unlink(paste(file, "*", sep = ".")))

    ## preprocess
    opt <- ""
    nop <- ceiling((dim(data)[1] + 2 * length(data@data@i))
	 * .Machine$sizeof.long / 4^10 / 5)
    if (length(control@memsize)) {
        opt <- paste("-m", control@memsize)
        nop <- ceiling(nop * 32 / control@memsize)
    }
    if (length(control@numpart)) {
        if (control@numpart < nop)
            warning("'numpart' less than recommended")
        nop <- control@numpart
    }
    ## workaround
    out <- paste(file, "stdout", sep = ".")
    ## deprecated
    if (FALSE) {
	asc <- paste(file, "asc", sep = ".")
	write_cspade(data, con = asc)
	if (system2(file.path(exe, "makebin"), args = c(
		asc, paste(file, "data", sep = "."))) ||
	    system2(file.path(exe, "getconf"), args = c(
		"-i", file, "-o", file), stdout = out)
	) stop("system invocation failed")
	file.append("summary.out", out)
    } else
	makebin(data, file)
    if (system2(file.path(exe, "exttpose"), args = c(
            "-i", file, "-o", file, "-p", nop, opt, "-l -x -s",
            parameter@support), stdout = out)
       ) stop("system invocation failed")
    file.append("summary.out", out)

    if (!is.null(class))
	write_class(class, paste(file, "class", sep = "."))

    ## options
    if (length(parameter@maxsize))
        opt <- paste(opt, "-Z", parameter@maxsize, collapse = "")
    if (length(parameter@maxlen))
        opt <- paste(opt, "-z", parameter@maxlen,  collapse = "")
    if (length(parameter@mingap))
        opt <- paste(opt, "-l", parameter@mingap,  collapse = "")
    if (length(parameter@maxgap))
        opt <- paste(opt, "-u", parameter@maxgap,  collapse = "")
    if (length(parameter@maxwin))
        opt <- paste(opt, "-w", parameter@maxwin,  collapse = "")

    if (!length(control@bfstype) || !control@bfstype)
        opt <- paste(opt, "-r", collapse = "")
    if (control@tidLists)
	opt <- paste(opt, "-y", collapse = "")
    if (!is.null(class))
	opt <- paste(opt, "-c", collapse = "")

    if (control@verbose) {
        t2 <- proc.time()
        du <- sum(file.info(list.files(path = dirname(file),
            pattern = basename(file), full.names = TRUE))$size)
        cat(paste("", nop, "partition(s),",
            round(du / 4^10, digits = 2), "MB"))
        cat(paste(" [",format((t2-t1)[3], digits =2, format = "f"),
                  "s]", sep = ""))
        cat("\nmining transactions ...")
    }

    out <- paste(file, "out", sep = ".")
    ## workaround
    tmp <- paste(file, "tpose", sep = ".")
    if (nop > 1L)
	tmp <- paste(tmp, "P0", sep = ".")
    if (!file.exists(tmp))
	local({
	    n <- readBin(paste(file, "conf", sep = "."), "integer") 
	    cat(file = out, 
		sprintf("MINSUPPORT %s out of %s sequences\n", 
			ceiling(parameter@support * n), n)
	    )
	})
    else
	if (system2(file.path(exe, "spade"), args = c(
	    "-i", file, "-s", parameter@support, opt, "-e", nop, "-o"), 
	    stdout = out)
	) stop("system invocation failed")

    if (control@verbose) {
        t3 <- proc.time()
        du <- file.info(out)$size
        cat(paste("", round(du / 4^10, digits = 2), "MB"))
        cat(paste(" [",format((t3-t2)[3], digits =2, format = "f"),
                  "s]", sep = ""))
        cat("\nreading sequences ...")
    }

    out <- read_spade(con = out, labels = itemLabels(data),
	transactions = 
	    if (control@tidLists)
		data,
	class = class
    ) 

    out@info <- c(
        data = match.call()$data,
        ntransactions = length(data),
        out@info,
        support = parameter@support
    )

    if (control@verbose) {
        t4 <- proc.time()
        cat(paste(" [",format((t4-t3)[3], digits =2, format = "f"),
                  "s]", sep = ""))
        cat("\n\ntotal elapsed time: ", (t4-t1)[3], "s\n", sep ="")
    }
    if (!control@summary)
        unlink("summary.out")

    out
}

###
