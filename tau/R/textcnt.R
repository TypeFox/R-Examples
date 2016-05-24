## ceeboo 2008

## <NOTE>
## Currently the approach to counting of word sequences is a bad
## workaround. As R has a character cache this could be implemented
## using pointer comparisons and prefix trees, but the latter would be
## inefficient compared to using hash tables. However, the R API does
## not provide the necessary interfaces to its hash table code.
## </NOTE>

textcnt <- 
function(x, n = 3L, split = "[[:space:][:punct:][:digit:]]+",
         tolower = TRUE, marker = "_", words = NULL, lower = 0L,
         method = c("ngram", "string", "prefix", "suffix"),
         recursive = FALSE, persistent = FALSE, useBytes = FALSE,
	 perl = TRUE, verbose = FALSE, decreasing = FALSE)
{
    if (is.list(x)) {
        if (recursive) {
            if (persistent) {
                for (z in x[-length(x)])
                    textcnt(z, n, split, tolower, marker, words, lower,
                            method, recursive, persistent, perl, useBytes, 
			    verbose)
                x <- x[length(x)]
                persistent <- FALSE
            } else
                return(lapply(x, textcnt, n, split, tolower, marker, words,
                              lower, method, recursive, persistent, useBytes,
			      perl, verbose))
        }
    } else {
        if (is.null(x)) return (x)
        x <- list(x)
    }
    n <- as.integer(n)
    if (n < 0L)
	stop("'n' invalid value")
    ## shortcut
    if (n < 1L)
	return(NULL)
    method <- match.arg(method)
    if (!is.null(split))
	x <- 
	if (!is.null(formals(strsplit)$useBytes))
	    lapply(lapply(x, strsplit, split, perl = perl, useBytes = useBytes), unlist)
	else
	    lapply(lapply(x, strsplit, split, perl = perl), unlist)
    if (!useBytes && 
        tolower)
        x <- lapply(x, tolower)
    if (!is.null(words))
        x <- .Call(R_copyTruncate, x, words)
    if (method == "ngram") {
	if (length(grep(marker, unlist(x, use.names = FALSE), 
			useBytes = TRUE)))
	    stop("'marker' contained in 'x'")
	## <NOTE>
	## As the C code generates all n-grams up to the specified
	## length we have to remove padding counts later. [2010/2]
	## </NOTE>
	if (marker == "\2") {
	    ## pad
	    m <- paste(rep(marker, n - 1L), collapse = "")
	    x <- lapply(x, function(x) gsub("$(?<!^)", m, x, perl = TRUE,
					    useBytes = useBytes))  
	}
        ## add marker at both ends
        x <- lapply(x, function(x) gsub("^(?!$)|$(?<!^)", marker, x,
                                        perl = TRUE, useBytes = useBytes))
        x <- .Call(R_utf8CountNgram, x, n, lower, verbose, persistent, 
				       useBytes)
	if (length(x) && marker == "\2") {
	    ## Determine suffixes
	    i <- grep("^\2{2,}$", names(x), useBytes = useBytes)
	    ## Remove
	    x <- x[-i]
	    ## Reduce to prefix count
	    x["\2"] <- x["\2"] / (n + 1L)
	    ## Replace with default marker
	    names(x) <- gsub("\2",formals(textcnt)$marker, names(x),
			     useBytes = useBytes)
	}
	## <NOTE>
	## Because of its design the C code can only adjust the prefix
	## counts. Here we handle the suffix counts. Obviously, a suffix
	## without the marker is a prefix to that suffix and thus the
	## former contains the counts of the latter. Note that this
	## includes the special case of a string that is both, a prefix
	## and a suffix. [2009/9]
	## </NOTE>
	if (length(x) && marker == "\1") {
	    ## Determine suffixes.
	    i <- grep("\1$", names(x), useBytes = useBytes)
	    ## Match suffixes, including the marker.
	    m <- gsub("\1$(?<!^\1)", "", names(x)[i], perl = TRUE,
		      useBytes = useBytes)
	    m <- match(m, names(x), nomatch = 0L)
	    i <- i[m > 0]
	    m <- m[m > 0]
	    x[m] <- x[m] - x[i]	    ## reduce counts
	    x <- x[x > 0]	    ## remove
	    ## Replace with the default marker.
	    names(x) <- gsub("^\1|\1$", formals(textcnt)$marker, names(x),
                             useBytes = useBytes)
	}
    } else {
        ## preprocess 
	if (method == "string" && n > 1L) 
	    x <- lapply(x, function(x) {
		x <- .Call(R_removeBlank, x)
		x <- unlist(lapply(.Call(R_copyToNgram, x, n), paste,
		    	           collapse = " "))
		if (is.null(x))
		    x <- ""
		x
	    })
        x <- .Call(R_utf8CountString, x, n, lower,
                   match(method, c("string", "prefix", "suffix")) - 1L, 
                   verbose, persistent, useBytes)
    }
    if (!is.null(x)) {
	## Note ties are not reordered (see help sort) which
        ## means they are in prefix-tree order.
        if (decreasing)
            sort(x, decreasing = TRUE) -> x
	attr(x, "useBytes") <- useBytes
        class(x) <- "textcnt"
    }
    x
}

## util
format.textcnt <- 
function(x, ...)
    data.frame(
	frq = c(x),
	rank = rank(c(x)),
        bytes = nchar(names(x), type = "bytes", allowNA = TRUE),
	Encoding = Encoding(names(x)),
        row.names = names(x),
        stringsAsFactors = FALSE
    )


###
