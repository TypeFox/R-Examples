callTip <- function (code, only.args = FALSE, location = FALSE,
description = FALSE, methods = FALSE, width = getOption("width"))
{
	code <- attr(completion(code, types = NA, description = FALSE), "fguess")
	if (is.null(code) || !length(code) || code == "")
		return("")

	## Get the corresponding calltip
	ctip <- argsTip(code, only.args = only.args, width =NULL)  # Reflow later!
	if (is.null(ctip)) return("")
	
	## Do we need to append an indication of where this function is located?
	if (isTRUE(location)) {
 		where <- res <- eval(parse(text = paste("getAnywhere(", code, ")",
			sep = "")))$where[1]
		if (!is.na(where) && where != ".GlobalEnv")
			ctip <- paste(ctip, " [", sub("^package:", "", where), "]", sep = "")
	}
	## Reflow the tip now
	if (!is.null(width))
		ctip <- paste(strwrap(ctip, width = width, exdent = 4), collapse = "\n")
	
	## Do we add the description of this function?
	if (isTRUE(description)) {
		desc <- descFun(code)
		if (!is.null(desc) && length(desc) && desc != "") {
			if (!is.null(width))
				desc <- paste(strwrap(desc, width = width), collapse = "\n")
			ctip <- paste(ctip, "\n\n", desc, sep = "")
		}
	}
	
	## Do we add a short mention of available methods if the function is generic?
	if (isTRUE(methods)) {
		mets <- listMethods(code)
		if (length(mets)) {
			## How many 25 char strings can we put on width and 5 lines max?
			## Note: we use two space each time as separator, except for last
			## line => take this into account in the calculation
			if (is.null(width)) nitems <- 3 else nitems <- (width + 2) %/% 27
			if (nitems < 1) nitems <- 1
			
			## Make sure the list is not too long: restrict to nitems * 5 entries
			if (length(mets) > nitems * 5) mets <- c(mets[1:(nitems * 5)], "...")
			
			## Make sure each method description is not longer than 25 characters
			n <- nchar(mets)
			## Cut entries that are too long
			tooLong <- n > 25
			mets[tooLong] <- paste(substr(mets[tooLong], 1, 22), "...", sep = "")
				
			## Paste strings together
			mets <- paste(format(mets, width = 25), c(rep("  ", nitems - 1), "\n"),
				collapse = "", sep = "")
			## Add this info to the calltip
			ctip <- paste(ctip,
				"\n\nGeneric function with methods for the following classes:\n", mets,
				sep = "")
		}
	}
	ctip
}
