## ==============
## Displaying PST
## ==============

setMethod("print", "PSTr", function (x, max.level = NULL, digits = 1, give.attr = FALSE, 
    nest.lev = 0, indent.str = "", stem = "--"
	## , ...
	) {

    pasteLis <- function(lis, dropNam, sep = " = ") {
        lis <- lis[!(names(lis) %in% dropNam)]
        fl <- sapply(lis, format, digits = digits)
        paste(paste(names(fl), fl, sep = sep), collapse = ", ")
    }

	at <- attributes(x)

	path <- at[["path"]]
	prob <- at[["prob"]]
	hgt <- at[["order"]]
	n <- at[["n"]]

	le <- length(which.child(x))
	leaf <- x@leaf

	istr <- sub(" $", "`", indent.str)

	for (g in 1:length(leaf)) {
		cat(istr, stem, sep = "")

		## left <- if (leaf[g]) {"("} else {"["}
		left <- "("
		## right <- if (leaf[g]) {")"} else {"]"}
		right <- ")"
        
        if (give.attr) {
           		if (nzchar(at <- pasteLis(at, c("prob", "order", "path")))) { at <- paste(",", at) }
        }
        cat(left, path, right, "-[ ", sep="")
		## " (k=", format(hgt, digits = digits), "), ", 
		cat("p=(", sep="")
		cat(format(prob[g,], digits = digits, scientific=FALSE, nsmall=digits), sep=",")
		cat(") - n=",n[g], if (give.attr) at, " ]", 
				## if (!is.null(max.level) && nest.lev == max.level) " ..", "\n", 
				sep = "")

		if (is.null(which.child(x))) {
			if (leaf[g]) cat("--| \n") else cat("--() \n", sep="")
		} else {
			cat("\n")
		}

		istr <- sub(" $", " ", indent.str)
		stem <- "  "
	}

	if (!all(leaf)) {
		if (is.null(max.level) || nest.lev < max.level) {
			for (i in which.child(x)) {
				print(x[[i]], nest.lev = nest.lev + 1, 
					indent.str = paste(indent.str, if (i < le) " |" else "  "), 
					max.level = max.level, digits = digits, 
					give.attr = give.attr)
			}
		}
	}
}
)

