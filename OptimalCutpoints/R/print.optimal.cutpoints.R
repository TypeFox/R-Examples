print.optimal.cutpoints <-
function(x, digits = max(3L, getOption("digits") - 3L), ...) {
	cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
	methods <- x[x$methods]
	levels.cat <- if(is.null(x$levels.cat)) {"Global"} else {x$levels.cat}
	for (i in 1:length(levels.cat)) {
		if(length(levels.cat) > 1) {
			cat(paste("\nOptimal cutoffs - ", levels.cat[i], ":",sep = ""),"\n")
		} else {
			cat("\nOptimal cutoffs:\n")
		}
		res <- vector("list", length(methods))
		for(j in 1:length(methods)) {
			n.cutpoints <- length(methods[[j]][[i]][["optimal.cutoff"]][[1]])
			if(n.cutpoints != 0) {
				res[[j]] <- methods[[j]][[i]][["optimal.cutoff"]][[1]]		
			}
		}
		names(res) <- names(methods)
		n.max <- max(unlist(lapply(res, length)))
		m <- matrix(ncol = length(methods), nrow = n.max, dimnames = list(1:n.max, names(methods)))
		for(j in 1:length(methods)) {
			for(k in 1:length(res[[j]])) {
				m[k,j] <- sprintf(paste("%.",digits,"f", sep = ""), res[[j]][k])
			}
		}
		print(m, quote = FALSE, right = TRUE, na.print = "-", row.names = FALSE)
		cat("\nArea under the ROC curve (AUC): ", paste(round(methods[[1]][[i]][["measures.acc"]][["AUC"]][1], 3), " (", round(methods[[1]][[i]][["measures.acc"]][["AUC"]][2], 3),"",", ", round(methods[[1]][[i]][["measures.acc"]][["AUC"]][3], 3),")", sep = ""), "\n")
	}
}
