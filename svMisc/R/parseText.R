parseText <- function (text, firstline = 1, srcfilename = NULL,
encoding = "unknown")
{
	## Parse R instructions provided as a string and return the expression if it
	## is correct, or a 'try-error' object if it is an incorrect code, or NA if
	## the (last) instruction is incomplete
  	text <- paste(text, collapse = "\n")
	## if firstline is higher than 1, "align" code by prepending empty codes
	firstline <- as.integer(firstline)[1]
	if (firstline > 1)
		text <- paste(c(rep("", firstline - 1), text), collapse = "\n")
	if (is.null(srcfilename)) srcfilename <- "<text>"
	res <- tryCatch(parse(text = text, srcfile = srcfilecopy(srcfilename, text),
		encoding = encoding), error = identity)

	if (inherits(res, "error")) {
		## Check if this is incomplete code
		msg <- conditionMessage(res)
		
		## Incomplete string
		if (regexpr(gettext("INCOMPLETE_STRING", domain = "R"), msg) > 0)
			return(NA)
		## Incomplete instruction
		if (regexpr(gettext("end of input", domain = "R"), msg) > 0)
			return(NA)	
		
		## This should be incorrect R code
		## Rework the message a little bit... keep line:col position in front
		
		## TODO: from SciViews-K-dev:
		## This reformats the message as it would appear in the CLI:
		#errinfo <-
		#	strsplit(sub("(?:<text>:)?(\\d+):(\\d+): +([^\n]+)\n([\\s\\S]*)$", "\\1\n\\2\n\\3\n\\4", msg, perl=T), "\n", fixed=TRUE)[[1]]

		#errpos <- as.numeric(errinfo[1:2])
		#err <- errinfo[-(1:3)]
		#rx <- sprintf("^%d:", errpos[1])
		#errcode <- sub(rx, "", err[grep(rx, err)])
		#err <- simpleError(sprintf("%s in \"%s\"", errinfo[3], errcode))
		
		## -or-
		
		err <- res
		err$message <- res <- sub("^<.*>:", "", msg)
		## Call is from instructions in "text"... but from the corresponding line
		err$call <- strsplit(text, "\n")[[1]][as.integer(
			sub("^[^0-9]*([0-9]+):.*$", "\\1", res))]
		
		## ... until here...
		
		## Return a try-error object to remain compatible with previous versions
		## TODO: from SciViews-K-dev:
		#res <- .makeMessage(res)
		class(res) <- "try-error"
		attr(res, 'error') <- err
	}

    res
}
