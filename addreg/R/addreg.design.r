addreg.design <- function(terms, data, allref, design.ref) {
    data.new <- data
	termlabels <- attr(terms,"term.labels")
    design.type <- sapply(allref, attr, "type")
	ref.vector <- as.vector(design.ref, mode = "integer")
	for(i in seq_len(length(design.ref))) {
		varname <- gsub("`","",termlabels[i])
        varref <- allref[[termlabels[i]]][[ref.vector[i]]]
		if(design.type[i] == 1) {
			cont.min <- min(data[[varname]], na.rm = TRUE)
			cont.max <- max(data[[varname]], na.rm = TRUE)
			if(varref == 1) data.new[[varname]] <- data[[varname]] - cont.min
			else data.new[[varname]] <- cont.max - data.new[[varname]]
		} else if(design.type[i] == 2) {
			data.new[[varname]] <- relevel(factor(data[[varname]]), ref = varref)
			contrasts(data.new[[varname]]) <- contr.treatment(levels(data.new[[varname]]), base = 1)
		} else if(design.type[i] == 3) {
			data.new[[varname]] <- factor(data[[varname]])
			contrasts(data.new[[varname]]) <- contr.isotonic(levels(data.new[[varname]]), perm = varref)
		}
	}
	attr(data.new, "terms") <- terms
	X <- model.matrix(terms, data.new)
	X
}