messageSGP <- function(tmp.message, domain=NULL, appendLF=TRUE) {
	if (!is.call(tmp.message)) {
		base::message(tmp.message)
	}
	PrintLogMessage(tmp.message)
	invisible()
}
PrintLogMessage <- function(tmp.message, domain=NULL) {
	# print log message to file
	dir.create("Logs", showWarnings = FALSE)
	logfile <- paste("Logs/SGP-", packageVersion("SGP"), "_", Sys.Date(), ".txt", sep="")

	if (is.call(tmp.message)) {
		tmp.message2 <- c(paste("\n\n\t", as.character(tmp.message)[1], "(\n\t\t", sep=""), paste(names(tmp.message)[-1], as.character(tmp.message)[-1], sep=" = ", collapse="\n\t\t"), ")\n\n")		
		cat(tmp.message2, file = logfile, append=TRUE)
	} else cat(tmp.message, "\n", file=logfile, sep="", append=TRUE)
}
