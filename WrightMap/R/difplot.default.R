difplot.default <- function(data, grouptype = NULL, group = NULL, item.names = NULL, ylim = c(-1, 
	1), ylab = NULL, ...) {
		
	table <- data

	if (is.null(group)) {
		if (is.null(grouptype)) {
			columns <- table[grep("^n_", names(table)) + 1]
			columns <- sapply(columns, unique,simplify = FALSE)
			columns <- sapply(columns, length)
			columns <- sort(columns)
			grouptype <- names(columns[1])
		}

		group <- table[grouptype][1, ]

	} else if (is.null(grouptype)) {
		grouptype <- names(table)[grep(paste("\\b", group, "\\b", sep = ""), table)[1]]
	}
	
	m <- table[grouptype] == group

	if (is.null(item.names)) {
		columns <- table[grep("\\bn_", names(table)) + 1]
		columns <- columns[-grep(paste("\\b", grouptype, "\\b", sep = ""), names(columns))]
		item.names <- columns[m, 1]
	}
	
	if(is.null(ylab))
		ylab <- paste0("Estimates of ",grouptype," effects (",group,")")
	
	plotCI(table$est[m], table$error[m], item.names, ylim = ylim, ylab = ylab, ...)
}