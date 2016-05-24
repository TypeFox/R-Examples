difplot.CQmodel <- function(data, table.name = NULL, grouptype = NULL, group = NULL, ...) {
	RMP <- data$RMP
	if (is.null(table.name)) {
		if (is.null(grouptype)) {
			if (is.null(group)) {
				table.at <- grep("(?<!^step)\\*(?!step$)", names(RMP), perl = TRUE)[1]
				table <- RMP[[table.at]]
				if (is.null(table)) 
					stop("No table of the form [factor]*[factor] found")
			} else {
				tables <- RMP[grep(paste("\\b", group, "\\b", sep = ""), RMP)]
				table.at <- grep("(?<!\\bstep)\\*(?!step\\b)", names(tables), perl = TRUE)[1]
				table <- tables[[table.at]]
				if (is.null(table)) 
					stop(paste("No table of the form [factor]*[factor] containing the group \"", group, "\" found", 
						sep = ""))
			}

		} else {
			table.at <- grep(paste("(?<!\\bstep)\\*", grouptype, "\\b|\\b", grouptype, "\\*(?!step\\b)", sep = ""), 
				names(RMP), perl = TRUE)[1]
			table <- RMP[[table.at]]
			if (is.null(table)) 
				stop(paste("No table of the form [factor]*", grouptype, " found", sep = ""))

		}
	} else {
		table.name2 <- gsub("\\*", "\\\\*", table.name)
		table.at <- grep(paste("^", table.name2, "$", sep = ""), names(RMP))[1]
		table <- RMP[[table.at]]
		if(is.null(table))
			stop(paste("No table named \'",table.name,"\'found",sep=" "))
	}
	difplot(table, grouptype, group, ...)
}