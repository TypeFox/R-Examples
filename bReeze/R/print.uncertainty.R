print.uncertainty <- 
function(x, ...) {
### summarising uncertainty object information
	
	cat("\n\tUncertainty\n\n")
	ucm <- x$uncertainty.meth
	cat("Uncertainties of applied methods:\n")
	ucm.units <- data.frame(t(names(ucm)))
	ucm.units[,] <- "[%]"
	obj <- as.data.frame(lapply(ucm, as.character))
	names(ucm) <- "uncertainty"
	names(ucm.units) <- names(obj) <- names(ucm)
	row.names(ucm.units) <- " "
	print(rbind(ucm.units, ucm), quote=FALSE)
	pe <- x$prob.exceedance
	cat("\nProbability of exceedance:\n")
	pe.units <- data.frame(t(names(pe)))
	pe.units[,1] <- paste0("[", attr(pe$probability, "unit"), "]")
	pe.units[,2] <- paste0("[", attr(pe$aep, "unit"), "]")
	obj <- as.data.frame(lapply(pe, as.character))
	names(pe)[2] <- "AEP"
	names(pe.units) <- names(obj) <- names(pe)
	row.names(pe.units) <- " "
	row.names(pe) <- row.names(pe)
	print(rbind(pe.units, pe), quote=FALSE)
	cat("\ncall: uncertainty(aep=", attr(x, "call")$aep, ", uc.values=c(", paste(attr(x, "call")$uc.values, collapse=", "), "), uc.names=c(\"", paste(attr(x, "call")$uc.names, collapse="\", \""), "\"), prob=c(", paste(attr(x, "call")$prob, collapse=", "), "), digits=c(", paste(attr(uc, "call")$digits, collapse=", "), "), print=", attr(x, "call")$print, ")\n\n", sep="")
}
