factors <- function(e) {
    if (is.name(e)) list(e)
    else switch(as.character(e[[1]]),
    	"*" = c(factors(e[[2]]),factors(e[[3]])),
    	"(" = factors(e[[2]]),
    	list(e))
}

term2table <- function(rowterm, colterm, env, n) {
    rowargs <- factors(rowterm)
    colargs <- factors(colterm)
    allargs <- c(rowargs, colargs)
    rowsubset <- TRUE
    colsubset <- TRUE
    pctdenom <- NULL
    pctsubset <- TRUE
    values <- NULL
    summary <- NULL
    arguments <- NULL
    format <- NA
    justification <- NA
    for (i in seq_along(allargs)) {
        e <- allargs[[i]]
        fn <- ""
        if (is.call(e) && (fn <- as.character(e[[1]])) == ".Format")
            format <- e[[2]]
        else if (fn == "Justify") 
            justification <- as.character(e[[if (length(e) > 2) 3 else 2]])
        else if (fn == "Percent") {
            env1 <- new.env(parent=env)
            percent <- function(x, y) 100*length(x)/length(y)
            env1$Equal <- env1$Unequal <- function(...) sys.call()
            env1$Percent <- function(denom="all", fn=percent) {
              if (is.null(summary)) {
                if (identical(denom, "all")) summary <<- function(x) fn(x, values)
                else if (identical(denom, "row")) summary <<- function(x) fn(x, values[rowsubset])
                else if (identical(denom, "col")) summary <<- function(x) fn(x, values[colsubset])
                else if (is.call(denom) && deparse(denom[[1]]) %in% c("Equal", "Unequal")) { 
                    summary <<- local({
                        pctdenom <<- sapply(as.list(denom), deparse, width.cutoff = 500L)
             		pctsubset <<- pctdenom[1] == "Equal"
                    	function(x) {
			    fn(x, values[pctsubset])
			}})
		} else if (is.logical(denom)) summary <<- function(x) fn(x, values[denom])
		else summary <<- function(x) fn(x, denom)
                summaryname <<- "Percent"
              } else
    	        stop("Summary fn not allowed with Percent")
            }
            eval(e, env1)
	} else if (fn == "Arguments") {
	    if (is.null(arguments)) 
	      arguments <- e
	    else
	      stop("Duplicate Arguments list: ", deparse(arguments), " and ", deparse(e))
        } else if (fn != "Heading" && !identical(e, 1)) {
    	    arg <- eval(e, env)
    	    asis <- inherits(arg, "AsIs")
    	    if (asis || is.vector(arg) || inherits(arg, "labelledSubset")) {
    	    	if (missing(n)) n <- length(arg)
    	    	else if (n != length(arg))
    	    	    stop("Argument ", deparse(e), " is not length ", n)
    	    }
    	    
    	    if (!asis && is.logical(arg)) {
    	    	arg <- arg & !is.na(arg)
    	    	if (i <= length(rowargs))
    	    	    rowsubset <- rowsubset & arg
    	    	else
    	    	    colsubset <- colsubset & arg
    	    } else if (asis || is.atomic(arg)) {
    	    	if (is.null(values)) {
    	    	    values <- arg
    	    	    valuename <- e
    	    	} else 
    	    	    stop("Duplicate values: ", deparse(valuename),
    	    	         " and ", deparse(e))
    	    } else if (is.function(arg)) {
    	    	if (is.null(summary)) {
    	    	    summary <- arg
    	    	    summaryname <- e
    	    	} else
    	    	    stop("Duplicate summary fn: ", deparse(summaryname),
    	    	         " and ", deparse(e))
    	    } else 
    	    	stop("Unrecognized entry ", deparse(e))
    	}
    }
    if (!is.null(pctdenom)) { # We need a second pass to find the subsets
	for (i in seq_along(allargs)) {
	    e <- allargs[[i]]
	    fn <- ""
	    if (is.call(e)) 
	    	fn <- as.character(e[[1]])
	    if (!(fn %in% c(".Format", "Justify", "Percent", "Heading"))
	        && !identical(e, 1)) {
	        arg <- eval(e, env)
	        asis <- inherits(arg, "AsIs")
	        if (!asis && is.logical(arg)) {
	            if (inherits(arg, "labelledSubset")) {
	    		argexpr <- attr(arg, "label")
	    	    	arg <- arg & !is.na(arg)
	    	    	if (pctdenom[1] == "Equal" && argexpr %in% pctdenom[-1])
	    	    	    pctsubset <- pctsubset & arg
	    	    	else if (pctdenom[1] == "Unequal" && argexpr %in% pctdenom[-1])
	    	    	    pctsubset <- pctsubset | !arg
	    	    } else
	    	    	pctsubset <- pctsubset & !is.na(arg) & arg
	        }
	    }
	}
    }
    	
    if (missing(n))
    	stop("Length of ", deparse(rowterm), "~", deparse(colterm),
    	     " indeterminate")    
    if (is.null(summary)) {
	if (!is.null(arguments))
	    stop(deparse(arguments), " specified without summary function.")
        summary <- length
    }
    if (is.null(values) && is.null(arguments)) values <- rep(NA, n)
    subset <- rowsubset & colsubset 
    if (is.null(arguments)) 
	value <- summary(values[subset])
    else {
	arguments[[1]] <- summary
	for (i in seq_along(arguments)[-1]) {
	    arg <- eval(arguments[[i]], env)
	    if (length(arg) == n) 
		arg <- arg[subset]
	    arguments[[i]] <- arg
	}
	if (!is.null(values)) {
	    named <- !is.null(names(arguments))
	    for (i in rev(seq_along(arguments)[-1])) {
		arguments[[i+1]] <- arguments[[i]]
		if (named) names(arguments)[i+1] <- names(arguments)[i]
	    }
	    arguments[[2]] <- values[subset]
	    if (named) names(arguments)[2] <- ""
	}
	value <- eval(arguments, env)
    }
    if (length(value) != 1)
	warning("Summary statistic is length ", length(value), call. = FALSE)
    structure(list(value), n=n, format=format, 
                   justification=justification)
}

# This moves column names into their own column
moveColnames <- function(labels, do_move = (names != "")) {
    attrs <- attributes(labels)
    names <- colnames(labels)
    colnamejust <- attrs$colnamejust
    justification <- attrs$justification
    for (i in rev(seq_along(do_move))) {
    	if (do_move[i]) {
	    before <- seq_len(i-1)
	    after <- seq_len(ncol(labels) - i + 1) + i - 1
    	    labels <- cbind(labels[,before,drop=FALSE], "", labels[,after,drop=FALSE])
	    labels[1,i] <- names[i]
	    names <- c(names[before], "", "", names[after][-1])
	    if (length(colnamejust)) {
	    	colnamejust <- c(colnamejust[before], NA_character_, colnamejust[after])
	    }
	    if (length(justification)) 
	    	justification <- cbind(justification[,before,drop=FALSE], 
	    	    NA_character_, justification[,after,drop=FALSE])
	}
    }
    attrs$colnamejust <- colnamejust
    attrs$justification <- justification
    attrs$dim <- dim(labels)
    attrs$dimnames[[2]] <- names
    attributes(labels) <- attrs
    labels
}
    	
getLabels <- function(e, rows=TRUE, justify=NA, head=NULL, suppress=0) {
    op <- ""
    justification <- NULL
    colnamejust <- character(0)
    Heading <- head
    result <- if (rows) matrix(NA, ncol=0, nrow=1)
    	      else      matrix(NA, ncol=1, nrow=0)
    nrl <- ncl <- leftjustify <- leftheading <- leftsuppress <- 
           leftjustification <- leftcolnamejust <- 
           nrr <- ncr <- rightjustify <- rightheading <- rightsuppress <-
           rightjustification <- rightcolnamejust <- NULL
    getLeft <- function() {
	nrl <<- nrow(leftLabels)
	ncl <<- ncol(leftLabels)
	leftjustify <<- attr(leftLabels, "justify")
        leftheading <<- attr(leftLabels, "Heading")
	leftsuppress <<- attr(leftLabels, "suppress")
	leftjustification <<- attr(leftLabels, "justification")
	leftcolnamejust <<- attr(leftLabels, "colnamejust")
    }
    getRight <- function() {
	nrr <<- nrow(rightLabels)
	ncr <<- ncol(rightLabels)
	rightjustify <<- attr(rightLabels, "justify")
        rightheading <<- attr(rightLabels, "Heading")
	rightsuppress <<- attr(rightLabels, "suppress")
	rightjustification <- attr(rightLabels, "justification")
	rightcolnamejust <- attr(rightLabels, "colnamejust")
    }
    if (is.call(e) && (op <- as.character(e[[1]])) == "*")  {
        leftLabels <- getLabels(e[[2]], rows, justify, head, suppress)
        getLeft()
	# Heading and justify are the settings to carry on to later terms
	# justification is the matrix of justifications for
	# each label
	righthead <- Heading <- leftheading
	suppress <- leftsuppress
	if (!is.null(leftjustify))
	    justify <- leftjustify

	rightLabels <- getLabels(e[[3]], rows, justify, righthead, suppress)
	getRight()
	Heading <- rightheading
	suppress <- rightsuppress
	
	if (!is.null(rightjustify))
	    justify <- rightjustify

	if (rows) {
	    result <- justification <- matrix(NA_character_, nrl*nrr, ncl + ncr)
	    colnames(result) <- c(colnames(leftLabels), colnames(rightLabels))
	    colnamejust <- c(leftcolnamejust, rightcolnamejust)
	    for (i in seq_len(nrl)) {
		j <- 1 + (i-1)*nrr
		k <- seq_len(ncl)
		result[j, k] <- leftLabels[i,]
		if (!is.null(leftjustification))
		    justification[j, k] <- leftjustification[i,]
		j <- (i-1)*nrr + seq_len(nrr)
		k <- ncl+seq_len(ncr)
		result[j, k] <- rightLabels
		if (!is.null(rightjustification))
		    justification[j, k] <- rightjustification
	    }
	} else {
	    result <- justification <- matrix(NA_character_, nrl + nrr, ncl*ncr)
	    for (i in seq_len(ncl)) {
		j <- seq_len(nrl)
		k <- 1 + (i-1)*ncr
		result[j, k] <- leftLabels[,i]
		if (!is.null(leftjustification))
		    justification[j,k] <- leftjustification[,i]
		j <- nrl+seq_len(nrr)
		k <- (i-1)*ncr + seq_len(ncr)
		result[j, k] <- rightLabels
		if (!is.null(rightjustification))
		    justification[j,k] <- rightjustification
	    }
	}
    } else if (op == "+") {
        leftLabels <- getLabels(e[[2]], rows, justify, NULL, suppress)
        getLeft()
	Heading <- leftheading

	rightLabels <- getLabels(e[[3]], rows, justify, NULL, suppress)
	getRight()
	Heading <- rightheading
	suppress <- rightsuppress

	if (rows) {
	    # Here we have a problem:  we need to stack two things, each of which
	    # may have column names.  We use the following rule:
	    #  - if the column names for both things match, just use them.
	    #  - if the left one has a name, and the right doesn't, use the left name
	    #  - if both have names that don't match, add them as extra column(s)
	    #  - if the right has a name, and the left doesn't, treat as unmatched names
	    
	    leftnames <- colnames(leftLabels)
	    if (is.null(leftnames)) leftnames <- rep("", ncl)
	    rightnames <- colnames(rightLabels)
	    if (is.null(rightnames)) rightnames <- rep("", ncr)
	    if (!identical(rightnames, rep("", ncr)) &&
	        !identical(leftnames, rightnames)) {
	        rightLabels <- moveColnames(rightLabels)
		# some properties may have changed; just get them again
	        getRight()
	        leftLabels <- moveColnames(leftLabels)
	        getLeft()
		Heading <- rightheading
	        rightnames <- rep("", ncr)
	        leftnames <- rep("", ncl)
	    }    
	    cols <- max(ncl, ncr)
	    # Pad all to same width
	    if (ncl < ncr) {
	        pad <- rep("", ncr-ncl)
	    	leftnames <- c(pad, leftnames)
	    	if (length(leftcolnamejust)) {
	    	  pad[TRUE] <- NA_character_
	    	  leftcolnamejust <- c(pad, leftcolnamejust)
	    	}
	    	pad <- matrix("", nrl, ncr-ncl)
	    	leftLabels <- cbind(pad, leftLabels)
	    	if (!is.null(leftjustification)) {
	    	  pad[TRUE] <- NA_character_
	    	  leftjustification <- cbind(pad, leftjustification)
	    	}
	    	ncl <- ncr
	    } else if (ncl > ncr) {
	        pad <- rep("", ncl-ncr)
	    	rightnames <- c(pad, rightnames)
	    	if (length(rightcolnamejust)) {
	    	  pad[TRUE] <- NA_character_
	    	  rightcolnamejust <- c(pad, rightcolnamejust)
	    	}
	    	pad <- matrix("", nrr, ncl-ncr)
	    	rightLabels <- cbind(pad, rightLabels)
	    	if (!is.null(rightjustification)) {
	    	  pad[TRUE] <- NA_character_
	    	  rightjustification <- cbind(pad, rightjustification)
	    	}
	    	ncr <- ncl
	    }
	    result <- matrix("", nrl + nrr, cols)
	    justification <- matrix(NA_character_, nrl + nrr, cols)
	    colnames <- rep("", cols)
	    colnamejust <- rep(NA_character_, cols)
	    j <- seq_len(nrl)
	    result[j, ] <- leftLabels
	    colnames <- leftnames
	    if (length(leftcolnamejust))
	        colnamejust <- leftcolnamejust	    
	    if (length(leftjustification))
		justification[j, ] <- leftjustification

	    j <- nrl+seq_len(nrr)
	    result[j, ] <- rightLabels
	    
	    if (!is.null(rightjustification))
		justification[j, ] <- rightjustification 
	    if (!is.null(head)) {
		colnames[1] <- head
		head <- NULL
	    }
	    colnames(result) <- colnames
	} else {
	    nrows <- max(nrl, nrr)
	    result <- matrix("", nrows, ncl + ncr)
	    justification <- matrix(NA_character_, nrows, ncl + ncr)
	    j <- (nrows-nrl) + seq_len(nrl)
	    k <- seq_len(ncl)
	    result[j, k] <- leftLabels
	    if (!is.null(leftjustification))
		justification[j, k] <- leftjustification
	    j <- (nrows-nrr) + seq_len(nrr)
	    k <- ncl+seq_len(ncr)
	    result[j,k] <- rightLabels
	    if (!is.null(rightjustification))
		justification[j, k] <- rightjustification
	    if (!is.null(head)) {
		result <- rbind(rep(NA_character_, ncol(result)),
				result)
		result[1,1] <- head
		justification <- rbind(justification[1,], justification)
	    }
	}
    } else if (op == "(") {
    	return(getLabels(e[[2]], rows, justify, head, suppress))
    } else if (op == ".Format") {
    	result <- if (rows) matrix(NA, ncol=0, nrow=1)
    		  else      matrix(NA, ncol=1, nrow=0)
    } else if (op == "Heading") {
    	if (length(e) > 1) {
    	    override <- TRUE
    	    if (length(e) > 2) {
    	    	override <- as.logical(e[[3]])
    	    	if (is.na(override))
    	    	    stop("Second argument in ", deparse(e), " must be logical",
    	    	         call. = FALSE)
    	    }
    	    if (suppress <= 0 && (is.null(Heading) || override)) {
    	    	if (is.character(e[[2]]))
    	    	    Heading <- e[[2]]
    	    	else
    	    	    Heading <- deparse(e[[2]])
    	    	suppress <- 0
    	    } else 
    	    	suppress <- suppress - 1
    	} else
    	    suppress <- suppress + 1
    } else if (op == "Justify") {
    	justify <- as.character(e[[2]])
    } else if (op == "Arguments") {
	#suppress <- suppress + 1
    } else if (suppress > 0) {  # The rest just add a single label; suppress it
    	suppress <- suppress - 1
    } else if (!is.null(head)) {
    	result <- matrix(head, 1,1, dimnames=list(NULL, ""))
    	Heading <- NULL
    } else if (op == "Percent") {
        result <- matrix("Percent", 1,1, dimnames=list(NULL, ""))
    } else if (identical(e, 1)) 
    	result <- matrix("All", 1,1, dimnames=list(NULL, ""))
    else
    	result <- matrix(deparse(e), 1,1, dimnames=list(NULL, ""))
    if (is.null(justification))
    	justification <- matrix(justify, nrow(result), ncol(result))
    stopifnot(identical(dim(result), dim(justification)))
    structure(result, justification = justification, 
    	colnamejust = colnamejust, justify = justify,
   	Heading = Heading, suppress = suppress)
}

expandExpressions <- function(e, env) {
    if (is.call(e)) {
	if ((op <- as.character(e[[1]])) %in% c("*", "+", "~", "(", "=") ) {
	    e[[2]] <- expandExpressions(e[[2]], env)
	    if (length(e) > 2)
		e[[3]] <- expandExpressions(e[[3]], env)
	} else if (op == "Format" || op == ".Format" 
	        || op == "Heading" || op == "Justify"
	        || op == "Percent" || op == "Arguments")
	    e
	else {
	    v <- eval(e, envir=env)
	    if (is.language(v)) 
		e <- expandExpressions(v, env)
	}
    }
    e
}

collectFormats <- function(table) {
    formats <- list()
    recurse <- function(e) {
    	if (is.call(e)) {
    	    if ((op <- as.character(e[[1]])) %in% c("*", "+", "~", "(") ) {
    	    	e[[2]] <- recurse(e[[2]])
    	    	if (length(e) > 2)
	    	    e[[3]] <- recurse(e[[3]])
    	    } else if (op == c("Format")) {
    	    	if (length(e) == 2 && is.null(names(e)))
    	    	    formats <<- c(formats, list(e[[2]]))
    	    	else {
    	    	    e[[1]] <- as.name("format")
    	    	    formats <<- c(formats, list(e))
    	    	}    
    	    	e <- call(".Format", length(formats))
    	    }
    	}
    	e
    }
    result <- recurse(table)
    structure(result, fmtlist=formats)
}

checkDenomExprs <- function(e, subsetLabels) {
    if (is.call(e)) 
	if ((op <- as.character(e[[1]])) %in% c("*", "+", "~", "(", "=") ) {
	    checkDenomExprs(e[[2]], subsetLabels)
	    if (length(e) > 2)
		checkDenomExprs(e[[3]], subsetLabels)
	} else if (op == "Percent") {
	    e <- match.call(Percent, e)[["denom"]]
	    if (is.call(e) && as.character(e[[1]]) %in% c("Equal", "Unequal"))
		for (i in seq_along(e)[-1])
		    if (!(deparse(e[[i]]) %in% subsetLabels))
			stop("In ", deparse(e), ",\n", 
			    dQuote(deparse(e[[i]])), " is not a subset label.  Legal labels are\n", 
			    paste(subsetLabels, collapse=", "), call. = FALSE)
	}
}

collectSubsets <- function(e) {
    result <- c()
    if (is.call(e)) {
	if ((op <- as.character(e[[1]])) %in%  c("*", "+", "~", "(", "=") ) {
	    result <- c(result, collectSubsets(e[[2]]))
	    if (length(e) > 2)
		result <- c(result, collectSubsets(e[[3]]))
	} else if (op == "labelSubset") 
	    result <- c(result, match.call(labelSubset, e)[["label"]])
    }
    result
}

# This both expands factors and rewrites "bindings"

expandFactors <- function(e, env) {
    op <- ""
    if (is.call(e) && (op <- as.character(e[[1]])) %in% c("*", "+") ) 
    	call(op, expandFactors(e[[2]],env),expandFactors(e[[3]],env))
    else if (op == "(")
    	expandFactors(e[[2]],env)
    else if (op == "=") {
    	rhs <- expandFactors(e[[3]], env)
    	if (is.call(rhs) && as.character(rhs[[1]]) == "*"
    	 && is.call(rhs[[2]]) && as.character(rhs[[c(2,1)]]) == "Heading") {
    	    rhs[[c(2,2)]] <- as.name(deparse(e[[2]]))
    	    rhs
    	} else
    	    call("*", call("Heading", as.name(deparse(e[[2]]))), rhs)
    } else if (op == ".Format" || op == "Heading" || 
             op == "Justify" || op == "Percent" || op == "Arguments")
    	e
    else {
    	v <- eval(e, envir=env)
    	if (is.factor(v) & !inherits(v, "AsIs")) 
    	    e <- Factor(v, expr=e, override=FALSE)
    	e
    }
}


# A sum of products is a list whose elements are atoms or products of atoms.

sumofprods <- function(e) {
    if (!is.language(e)) return(list(e))
    if (is.expression(e)) e <- e[[1]]
    if (is.name(e)) result <- list(e)
    else {
	chr <- as.character(e[[1]])
	if (chr == "+") result <- c(sumofprods(e[[2]]),sumofprods(e[[3]]))
    	else if (chr == "*") {
	    left <- sumofprods(e[[2]])
	    right <- sumofprods(e[[3]])
	    result <- list()
	    for (i in 1:length(left)) 
	       for (j in 1:length(right)) 
		    result <- c(result, list(call("*", left[[i]], right[[j]])))
    	} else if (chr == "(") result <- sumofprods(e[[2]])
    	else if (chr == "~") stop("nested formulas not supported")
    	else result <- list(e)

    }    	
    result
}

tabledims <- function(e) {
    if (identical(e,1)) return(list(1))
    if (!is.language(e)) stop('Need an expression')
    if (is.expression(e)) e <- e[[1]]
    if (is.name(e)) result <- list(e)
    else {
        result <- list()
	chr <- as.character(e[[1]])
	if (chr == "~") {
	    if (length(e) == 2)
		result <- c(result, tabledims(e[[2]]))
	    else
	    	result <- c(result, tabledims(e[[2]]),tabledims(e[[3]]))
    	} else result <- list(e)
    }    	
    if (length(result) > 2) stop("only 2 dim tables supported")
    result
}

tabular <- function(table, ...) 
    UseMethod("tabular")
    
tabular.default <- function(table, ...) {
    tabular.formula(as.formula(table, env=parent.frame()), ...)
}

tabular.formula <- function(table, data=NULL, n, suppressLabels=0, ...) {
    formula <- table
    if (length(list(...)))
	warning(gettextf("extra argument(s) %s will be disregarded",
			 paste(sQuote(names(list(...))), collapse = ", ")),
		domain = NA)
    if (missing(n) && inherits(data, "data.frame"))
    	n <- nrow(data)
    if (is.null(data))
    	data <- environment(table)
    else if (is.list(data))
    	data <- list2env(data, parent=environment(table))
    else if (!is.environment(data))
    	stop("data must be a dataframe, list or environment")
    	
    table <- expandExpressions(table, data)
    table <- collectFormats(table)
    dims <- tabledims(table)
    if (length(dims) == 1) dims <- c(list(quote((` `=1))), dims)
    dims[[1]] <- expandFactors(dims[[1]], data)
    rlabels <- getLabels(dims[[1]], rows=TRUE, suppress=suppressLabels)
    suppressLabels <- attr(rlabels, "suppress")
    justify <- attr(rlabels, "justify")
    dims[[2]] <- expandFactors(dims[[2]], data)
    clabels <- getLabels(dims[[2]], rows=FALSE, justify=justify,
			 suppress=suppressLabels)
    
    # Check if the Percent calls name nonexistent terms
    subsetLabels <- unique(c(collectSubsets(dims[[1]]), collectSubsets(dims[[2]])))
    checkDenomExprs(dims[[1]], subsetLabels)
    checkDenomExprs(dims[[2]], subsetLabels)    
    
    rows <- sumofprods(dims[[1]])
    cols <- sumofprods(dims[[2]])

    result <- NULL
    formats <- NULL
    justifications <- NULL
    for (i in seq_along(rows)) {
    	row <- NULL
    	rowformats <- NULL
    	rowjustification <- NULL
    	for (j in seq_along(cols)) {
    	    # term2table checks that n matches across calls
    	    term <- term2table(rows[[i]], cols[[j]], data, n)
    	    n <- attr(term, "n")
    	    format <- attr(term, "format")
    	    justification <- attr(term, "justification")
    	    row <- cbind(row, term)
    	    rowformats <- cbind(rowformats, format)
    	    rowjustification <- cbind(rowjustification, justification)
    	}
    	result <- rbind(result, row)
    	formats <- rbind(formats, rowformats)
    	justifications <- rbind(justifications, rowjustification)
    }
    structure(result, formula=formula, rowLabels=rlabels, colLabels=clabels, table=table,
    	      formats = formats, justification = justifications, 
    	      class = "tabular")
}

justify <- function(x, justification="c", width=max(nchar(x))) {
    justification <- rep(justification, len=length(x))
    change <- justification %in% c("c", "l", "r")
    if (!any(change)) return(x)
    y <- x[change]
    justification <- justification[change]
    y <- sub("^ *", "", y)
    y <- sub(" *$", "", y)
    width <- rep(width, len=length(x))
    width <- width[change]
    lens <- nchar(y)
    ind <- justification == "c"
    if (any(ind)) {
    	left <- (width[ind] - lens[ind]) %/% 2
    	right <- width[ind] - lens[ind] - left
    	y[ind] <- sprintf("%*s%s%*s", left, "", y[ind], right, "")
    }
    ind <- justification == "l"
    if (any(ind)) {
    	right <- width[ind] - lens[ind]
    	y[ind] <- sprintf("%s%*s", y[ind], right, "")
    }
    ind <- justification == "r"
    if (any(ind)) {
    	left <- width[ind] - lens[ind]
    	y[ind] <- sprintf("%*s%s", left, "", y[ind])
    }
    x[change] <- y
    x
}

latexNumeric <- function(chars, minus=TRUE, leftpad=TRUE, rightpad=TRUE,
			 mathmode=TRUE) {
    regexp <- "^( *)([-]?)([^ -][^ ]*)( *)$"
    leadin <- sub(regexp, "\\1", chars)
    sign <- sub(regexp, "\\2", chars)
    rest <- sub(regexp, "\\3", chars)
    tail <- sub(regexp, "\\4", chars)
    
    if (minus && any(neg <- sign == "-")) {
    	if (any(leadin[!neg] == ""))
    	    leadin <- sub("^", " ", leadin)
    	leadin[!neg] <- sub(" ", "", leadin[!neg])
    	sign[!neg] <- "\\phantom{-}"
    }
    if (leftpad && any(ind <- leadin != "")) 
    	leadin[ind] <- paste("\\phantom{", 
    	                     gsub(" ", "0", leadin[ind]),
    	                     "}", sep="")
    if (rightpad && any(ind <- tail != ""))
    	tail[ind] <- paste("\\phantom{",
    			   gsub(" ", "0", tail[ind]),
    			   "}", sep="")
    if (mathmode)
    	paste("$", leadin, sign, rest, tail, "$", sep="")
    else
    	paste(leadin, sign, rest, tail, sep="")
}

format.tabular <- function(x, digits=4, justification="n", 
                           latex=FALSE, html=FALSE, 
                           leftpad = TRUE, rightpad = TRUE, minus = TRUE, mathmode = TRUE, ...) {
    if (latex && html) stop("Only one of 'latex' and 'html' may be requested")
    result <- unclass(x) 
    formats <- attr(x, "formats")
    fmtlist <- attr(attr(x, "table"), "fmtlist")
    justify <- attr(x, "justification")
    justify[is.na(justify)] <- justification
    ischar <- sapply(result, is.character)
    chars <- matrix(NA_character_, nrow(result), ncol(result))
    chars[ischar] <- unlist(result[ischar])
    lengths <- lapply(result, length)
    
    for (i in seq_len(ncol(result))) {
        ind <- col(result) == i & is.na(formats) & !ischar & lengths == 1
	if (any(ind)) {
            x <- do.call(c, result[ind])
    	    chars[ind] <- format(x, digits=digits, ...)
    	    if (is.numeric(x)) {
 		if (latex)
    	    	    chars[ind] <- latexNumeric(chars[ind], leftpad = leftpad, rightpad = rightpad, minus = minus, mathmode = mathmode)
    	    	else if (html)
    	    	    chars[ind] <- htmlNumeric(chars[ind], leftpad = leftpad, rightpad = rightpad, minus = minus)
    	    }
    	}
    }
    for (i in seq_along(fmtlist)) {
    	ind <- !is.na(formats) & formats == i
    	if (any(ind)) {        
            call <- fmtlist[[i]]
            isformat <- identical(call[[1]], as.name("format"))
            if (isformat) skip <- ischar | (lengths != 1)
            else skip <- ischar & FALSE
	    last <- length(call)
       	    x <- do.call(c, result[ind & !skip])
       	    call[[last+1]] <- x
       	    names(call)[last+1] <- "x"
       	    chars[ind & !skip] <- eval(call, parent.frame())
       	    if (isformat) {
       	    	if (latex) {
       	    	    if (is.numeric(x))
       	    	    	chars[ind] <- latexNumeric(chars[ind], leftpad = leftpad, rightpad = rightpad, minus = minus, mathmode = mathmode)
       	    	    else
       	    	    	chars[ind] <- texify(chars[ind])
       	    	} else if (html) {
       	    	    if (is.numeric(x))
       	    	    	chars[ind] <- htmlNumeric(chars[ind], leftpad = leftpad, rightpad = rightpad, minus = minus)
       	    	    else
       	    	    	chars[ind] <- htmlify(chars[ind])
       	    	}
       	    }
       	}
    }
    if (!latex && !html)
    	for (i in seq_len(ncol(result))) 
    	    chars[,i] <- justify(chars[,i], justify[,i])
    chars
}
