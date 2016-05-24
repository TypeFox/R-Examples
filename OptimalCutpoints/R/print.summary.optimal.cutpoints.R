print.summary.optimal.cutpoints <-
function(x, ...) {
	cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
	for(i in 1: length(x$p.table)) {
		if(!is.null(x$levels.cat)) {
			cat("\n*************************************************\n")
			cat(names(x$p.table)[i])
			cat("\n*************************************************\n")
		}
		cat("\nArea under the ROC curve (AUC): ", x$p.table[[i]][["AUC_CI"]], "\n\n")
		for (j in 1:(length(x$p.table[[i]]) - 1)) {
			cat(paste("CRITERION: " , names(x$p.table[[i]])[j], "\n", sep = ""))
			cat(paste("Number of optimal cutoffs: ", length(x$p.table[[i]][[j]]), "\n\n", sep = ""))
			if(length(x$p.table[[i]][[j]]) != 0) {
				for (k in 1:length(x$p.table[[i]][[j]])) {
					print(x$p.table[[i]][[j]][[k]], quote = FALSE, right = TRUE, na.print = "-")
					cat("\n")
				}
			}
		}
	}
	invisible(x)
}
