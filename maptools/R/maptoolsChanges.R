maptools <- function(changes=FALSE) {
	.DESC <- packageDescription("maptools")
	cat(.DESC[["Package"]], ", version ", .DESC[["Version"]],
	 	", ", .DESC[["Date"]], "\n", sep="")
	.CH <- NULL
	if (changes) {
		cat("\n")
		file <- system.file("changes", package = "maptools")
		.CH <- scan(file, list(version="character", 
			changes="character"), sep="\t", quiet=TRUE)
		for (i in length(.CH$changes):1) {
			cat(unlist(strsplit(.CH$changes[i], " ")), fill=TRUE, 
				labels=.CH$version[i])
			cat("\n")
		}
	}
	invisible(.CH)
}

