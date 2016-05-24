`[.tabular` <- function(x, i, j, ..., drop=FALSE) {
	
    if (drop) return(unclass(x)[i,j, drop=TRUE])

    attrs <- attributes(x)
    x <- unclass(x)[i,j, drop=FALSE]
    attrs$justification <- attrs$justification[i,j, drop=FALSE]
    attrs$formats <- attrs$formats[i,j, drop=FALSE]
    attrs$dim <- dim(x)
    attrs$dimnames <- dimnames(x)

    rowLabels <- attrs$rowLabels
    rattrs <- attributes(rowLabels)
    rattrs$justification <- rattrs$justification[i,, drop=FALSE]
    rattrs$dim <- dim(rattrs$justification)

    # This is the tricky bit
    orig <- row(rowLabels)
    for (k in seq_len(nrow(orig))[-1]) 
	orig[k,] <- ifelse(is.na(rowLabels[k,]), orig[k-1,], orig[k,])
    orig <- orig[i,, drop=FALSE]
    rowLabels <- rowLabels[cbind(c(orig), c(col(orig)))]
    dim(rowLabels) <- dim(orig)
    for (k in seq_len(nrow(orig))[-1])
	rowLabels[k,] <- ifelse(orig[k,] == orig[k-1,], NA, rowLabels[k,])
    rattrs$dimnames[1] <- list(NULL)

    # Do the same for column labels
    colLabels <- attrs$colLabels
    cattrs <- attributes(colLabels)
    cattrs$justification <- cattrs$justification[,j, drop=FALSE]
    cattrs$dim <- dim(cattrs$justification)
    orig <- col(colLabels)
    for (k in seq_len(ncol(orig))[-1]) 
	orig[,k] <- ifelse(is.na(colLabels[,k]), orig[,k-1], orig[,k])
    orig <- orig[,j, drop=FALSE]
    colLabels <- colLabels[cbind(c(row(orig)), c(orig))]
    dim(colLabels) <- dim(orig)
    for (k in seq_len(ncol(orig))[-1])
	colLabels[,k] <- ifelse(orig[,k] == orig[,k-1], NA, colLabels[,k])
    cattrs$dimnames[2] <- list(NULL)

    # Put it all back together
    attributes(rowLabels) <- rattrs
    attrs$rowLabels <- rowLabels
    attributes(colLabels) <- cattrs
    attrs$colLabels <- colLabels
    attributes(x) <- attrs
    x
}

cbind.tabular <- function(..., deparse.level = 1) {
    args <- list(...)
    if (!length(args)) return(NULL)
    result <- NULL
    while (length(args)) {
	x <- args[[1]]
	args <- args[-1]
	if (is.null(x)) next
	if (is.null(result)) {
	    attrs <- attributes(x)
	    result <- unclass(x)
	    fmtlist <- attr(attrs$table, "fmtlist")
	} else {	
	    xattrs <- attributes(x)	
	    if (nrow(result) != nrow(x) || !identical(attrs$rowLabels, xattrs$rowLabels) )
		stop("Cannot cbind if tables have different rows")
	    result <- cbind(result, unclass(x))

	    attrs$justification <- cbind(attrs$justification, xattrs$justification)
	    attrs$formats <- cbind(attrs$formats, xattrs$formats + length(fmtlist))
	    attrs$dim <- dim(result)
	    attrs$dimnames <- dimnames(result)
	
	    fmtlist <- c(fmtlist, attr(xattrs$table, "fmtlist"))
	    attrs$table <- c(attrs$table, xattrs$table)
	
	    # rowLabels are fine
	
	    colLabels <- attrs$colLabels
	    cattrs <- attributes(colLabels)
	    xcolLabels <- xattrs$colLabels
	    xcattrs <- attributes(xcolLabels)
	    colLabels <- cbind(colLabels, xcolLabels)
	    cattrs$dim <- dim(colLabels)
	    cattrs$dimnames <- dimnames(colLabels)
	    cattrs$justification <- cbind(cattrs$justification, xcattrs$justification)

	    attributes(colLabels) <- cattrs	
	    attrs$colLabels <- colLabels
	    attr(attrs$table, "fmtlist") <- fmtlist
	}
    }
    if (!is.null(result)) 
	attributes(result) <- attrs
    result
}

rbind.tabular <- function(..., deparse.level = 1) {
    args <- list(...)
    if (!length(args)) return(NULL)
    result <- NULL
    while (length(args)) {
	x <- args[[1]]
	args <- args[-1]	
	if (is.null(x)) next
	if (is.null(result)) {
	    attrs <- attributes(x)
	    result <- unclass(x)
	    fmtlist <- attr(attrs$table, "fmtlist")    
	} else {
	    xattrs <- attributes(x)	
	    if (ncol(result) != ncol(x) || !identical(attrs$colLabels, xattrs$colLabels) )
		stop("Cannot rbind if tables have different columns")
	    result <- rbind(result, unclass(x))

	    attrs$justification <- rbind(attrs$justification, xattrs$justification)
	    attrs$formats <- rbind(attrs$formats, xattrs$formats + length(fmtlist))
	    attrs$dim <- dim(result)
	    attrs$dimnames <- dimnames(result)
	
	    fmtlist <- c(fmtlist, attr(xattrs$table, "fmtlist"))
	    attrs$table <- c(attrs$table, xattrs$table)
	
	    rowLabels <- attrs$rowLabels
	    rattrs <- attributes(rowLabels)
	    xrowLabels <- xattrs$rowLabels
	    xrattrs <- attributes(xrowLabels)
	    rowLabels <- rbind(rowLabels, xrowLabels)
	    rattrs$dim <- dim(rowLabels)
	    rattrs$dimnames <- dimnames(rowLabels)
	    rattrs$justification <- rbind(rattrs$justification, xrattrs$justification)

	    # colLabels are fine
	
	    attributes(rowLabels) <- rattrs	
	    attrs$rowLabels <- rowLabels
	    attr(attrs$table, "fmtlist") <- fmtlist
	}
    }	
    if (!is.null(result)) 
	attributes(result) <- attrs
    result
}
