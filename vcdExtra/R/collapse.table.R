# collapse a contingency table or ftable by re-assigning levels of table variables
# revised to accept an array also

collapse.table <- function(table, ...) {
	nargs <- length(args <- list(...))
	if (!nargs) 
		return(table)
	if (inherits(table, "ftable"))
		table <- as.table(table)
	if (inherits(table, "array"))
		table <- as.table(table)
	if (inherits(table, "table")) {
		tvars <- names(dimnames(table))
		table <- as.data.frame.table(table)
		freq <- table[,"Freq"]
	}
	else stop("Argument must be a table, array or ftable object")
	
	names <- names(args)
	for (i in 1:nargs) {
		vals <- args[[i]]
		nm <- names[[i]]
		if(any(nm==tvars)) levels(table[[nm]]) <- vals
		else warning(nm, " is not among the table variables.")
	}
#    term <- paste(tvars, collapse = '+')
#    form <- as.formula(paste("freq ~", term))
#    cat("term: ", term, "\n")
	xtabs(as.formula(paste("freq ~", paste(tvars, collapse = '+'))), data=table)
}
