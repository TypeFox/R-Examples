as.matrix.tabular <- function(x, format=TRUE, 
			rowLabels=TRUE, colLabels=TRUE,
			justification="n", ...) {
    if (isTRUE(format)) {
    	chars <- format(x, ...)
	justify <- attr(x, "justification")
	rlabels <- attr(x, "rowLabels")
	rlabels[is.na(rlabels)] <- ""
	rjustify <- attr(rlabels, "justification")
	clabels <- attr(x, "colLabels")
	clabels[is.na(clabels)] <- ""
	cjustify <- attr(clabels, "justification")
	colnamejust <- attr(rlabels, "colnamejust")
	colnamejust[is.na(colnamejust)] <- justification
	corner <- matrix("", nrow(clabels), ncol(rlabels))
	corjustify <- matrix(NA, nrow(clabels), ncol(rlabels))
	for (i in seq_len(ncol(rlabels))) {
	    corner[nrow(clabels),i] <- colnames(rlabels)[i]
	    corjustify[nrow(clabels),i] <- rjustify[1,i]
	}
	if (rowLabels) {
	    if (colLabels) {
		result <- rbind(cbind(corner, clabels),
			cbind(rlabels, chars))
	        justify <- rbind(cbind(corjustify, cjustify),
	                cbind(rjustify, justify))
	    } else {
	    	result <- cbind(rlabels, chars)
	    	justify <- cbind(rjustify, justify)
	    }
	} else {
	    if (colLabels) {
	    	result <- rbind(clabels, chars)
	    	justify <- rbind(cjustify, justify)
	    } else 
	    	result <- chars
	}
	justify[is.na(justify)] <- justification
	for (i in seq_len(ncol(result)))
	    result[,i] <- justify(result[,i], justify[,i])
	dimnames(result) <- NULL
    } else {
    	result <- x
    	dim <- dim(x)
    	if (is.character(format)) format <- get(format, parent.frame)
    	if (is.function(format))
    	    result <- format(result)
        dim(result) <- dim
    }
    result
}


write.csv.tabular <- function(x, file="", 
    justification = "n", row.names=FALSE, ...) {
    options <- list(...)
    wtoptions <- names(options) %in% names(formals(write.table))

    result <- do.call(as.matrix, c(list(x, justification = justification),
    			        options[!wtoptions]))
    colnames(result) <- rep("", ncol(result))
    
    do.call(write.csv, c(list(result, file=file,
    	row.names=row.names), options[wtoptions]))
}

write.table.tabular <- function(x, file="", 
    justification = "n", row.names=FALSE, col.names=FALSE, ...) {
    options <- list(...)
    wtoptions <- names(options) %in% names(formals(write.table))
    result <- do.call(as.matrix, c(list(x, justification = justification),
    			        options[!wtoptions]))
    colnames(result) <- rep("", ncol(result))
    rownames(result) <- rep("", nrow(result))
    do.call(write.table, c(list(result, file=file,
    	row.names=row.names, col.names=col.names), options[wtoptions]))
}

print.tabular <- function(x, justification = "n", ...) {
   # chars <- format(x, justification = justification, ...)
    
    result <- as.matrix(x, justification = justification, ...)

    rownames(result) <- rep("", nrow(result))
    colnames(result) <- rep("", ncol(result))
    print(noquote(result))
    invisible(x)
}

