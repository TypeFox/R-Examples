# find sequential indices for observations in a data frame
# corresponding to the unique combinations of the levels
# of a given model term from a model object or a data frame
#

dataIndex <- function(x, term){
	data <- if (is.data.frame(x)) x else
			if (inherits(x, "lm")) model.frame(x) else stop("Not a data frame or model object")

    factors <- data[, sapply(data, is.factor), drop=FALSE]
    term.factors <- unlist(strsplit(term, ":"))
    if (any(which <- !term.factors %in% colnames(factors))) 
        stop(paste(term.factors[which], collapse=", "), " not in the model")
    n.factors <- length(term.factors)
    factor.values <- factors[,term.factors, drop=FALSE]
    rows <- nrow(levs <- unique(factor.values))
    levs <- apply(levs, 1, paste, collapse=":")

		m <- match(term.factors,colnames(factors))
		data.levs <- apply(factors[m], 1, paste, collapse=":")
    match(data.levs, levs) 
}

