texify <- function(x) {
    x <- gsub("\\", "\\textbackslash{}", x, fixed = TRUE)
    Hmisc::latexTranslate(x)
}

latex.tabular <- function(object, file="", options=NULL, ...) {
    if (file == "")
    	out <- ""
    else {
    	out <- file(file, open="wt")
    	on.exit(close(out))
    }
    if (!is.null(options)) {
    	saveopts <- do.call(table_options, options)
    	on.exit(table_options(saveopts), add=TRUE)
    }
    opts <- table_options()
    
    mycat <- function(...) cat(..., file=out)
    
    chars <- format(object, latex = TRUE, minus = opts$latexminus, 
                               leftpad = opts$latexleftpad, 
                               rightpad = opts$latexrightpad,...) # format without justification
    
    vjust <- attr(object, "justification")
    vjustdefs <- rep(opts$justification, length.out=ncol(object))
    ind <- !is.na(vjust) & vjust != rep(vjustdefs, each=nrow(vjust))
    chars[ind] <- sprintf("\\multicolumn{1}{%s}{%s}",
    			  vjust[ind], chars[ind])
    
    rowLabels <- attr(object, "rowLabels")
    nleading <- ncol(rowLabels)
    rowLabels[is.na(rowLabels)] <- ""
    rjust <- attr(rowLabels, "justification")
    rjustdefs <- rep(opts$rowlabeljustification, length.out=nleading)
    ind <- !is.na(rjust) & rjust != rep(rjustdefs, each=nrow(rowLabels))
    rowLabels[ind] <- sprintf("\\multicolumn{1}{%s}{%s}",
    			      rjust[ind], rowLabels[ind])
    rlabels <- apply(rowLabels, 1, paste, collapse = " & ")
    colnamejust <- attr(rowLabels, "colnamejust")
    colnamejust <- rep(colnamejust, length.out=nleading)
    ind <- is.na(colnamejust)
    colnamejust[ind] <- rjustdefs[ind]
    ind <- colnamejust != rjustdefs
    colnames(rowLabels)[ind] <- sprintf("\\multicolumn{1}{%s}{%s}", 
    		colnamejust[ind], colnames(rowLabels)[ind])
    clabels <- attr(object, "colLabels")
    leadin <- paste(rep("&", max(nleading - 1, 0)), collapse=" ")
    cjust <- attr(clabels, "justification")
    ind <- is.na(cjust)
    cjust[ind] <- rep(vjustdefs, each=nrow(cjust))[ind]

    clines <- character(nrow(clabels))
    vdefs <- c(vjustdefs, "STOP") # match length to row
    for (i in seq_len(nrow(clabels))) {
    	row <- c(clabels[i,], "STOP") # add sentinel
    	rowjust <- c(cjust[i,], "STOP")
    	result <- leadin
    	titlerules <- ""
    	label <- ""
    	just <- vdefs[1]
    	ncols <- 0
    	for (j in seq_along(row)) {
    	    if (!is.na(row[j])) {
    	    	if (ncols > 0) {
    	    	    if (ncols > 1 || just != vdefs[j]) {
    	    	        result <- c(result, 
    	    	            sprintf("& \\multicolumn{%d}{%s}{%s}", ncols, just, label))
    	    	        if (ncols > 1 && i < nrow(clabels) && !is.null(opts$titlerule))
    	    	            titlerules <- sprintf("%s%s{%d-%d}", titlerules, opts$titlerule,
    	    	            			  firstcol, firstcol+ncols-1)
    	    	    } else
    	    	    	result <- c(result, paste("&", label))
    	    	}
    	    	ncols <- 1
    	    	firstcol <- j+nleading
    	    	label <- row[j]
    	    	just <- rowjust[j]
    	    } else
    	    	ncols <- ncols + 1
    	}
    	clines[i] <- paste( paste(result, collapse=" "), "\\\\", titlerules)
    	if (nleading == 0) # Remove extra leading &
    	    clines[i] <- sub("&", "", clines[i], fixed=TRUE)
    }
    clabels <- clines
    # Replace the leadin of the last line with the row label headings
    if (nchar(leadin))
    	clabels[length(clabels)] <- 
    	    sub(leadin, paste(colnames(rowLabels), collapse=" & "),
    	        clabels[length(clabels)], fixed = TRUE)
    else
    	clabels[length(clabels)] <-
    	    paste(paste(colnames(rowLabels), collapse=" & "), 
    	    	  clabels[length(clabels)])
    if (opts$doBegin)
        mycat("\\begin{", opts$tabular, "}{", 
    	  paste(rjustdefs, collapse=""),
    	  paste(vjustdefs, collapse=""), 
    	  "}\n", sep="")
    if (opts$doHeader) {
	mycat(opts$toprule, "\n", sep="")
	mycat(clabels, sep="\n")
	mycat(opts$midrule, "\n", sep="")
    }
    if (opts$doBody) {
	chars <- apply(chars, 1, paste, collapse=" & ")
	if (nleading > 0)
	    chars <- paste(" & ", chars, sep="")
	chars <- paste(chars, "\\\\")
	mycat(paste(rlabels, chars), sep="\n")
    }
    if (opts$doFooter)
    	mycat(opts$bottomrule, "\n")
    if (opts$doEnd)
    	mycat("\\end{", opts$tabular, "}\n", sep="")
    structure(list(file=file, style=character(0)), class="latex")
}
