evalServer <- function (con, expr, send = NULL)
{
	## Evaluate expr on the R server, and return its value
	## con as returned by socketConnection(port = 8888)
	## send is optional. If supplied, expr must be a single unquoted object name.
	## Then send is evaluated on the client and the result is assigned
	## to that object on the server.
	## Robust flushing and dumping is just for windows. Linux is probably fine
	## without but no harm to leave in for now since binary mode will moot this.
	x <- substitute(expr)
	if (!missing(send) && (length(x) != 1 || mode(x) != "name"))
		stop("When send is supplied, expr must be a target variable name (unquoted) on the server to assign the result of the send expr to.")
	if (!is.character(x)) x <- deparse(x)

	readLines(con)  # Flush input stream just in case previous call failed to clean up
	if (missing(send)) {
		cat('..Last.value <- try(eval(parse(text = "', x,
			'"))); .f <- file(); dump("..Last.value", file = .f); flush(.f); seek(.f, 0); cat("\\n<<<startflag>>>", gsub("<pointer: [0-9a-fx]+>", "NULL", readLines(.f)), "<<<endflag>>>\\n", sep = "\\n"); close(.f); rm(.f, ..Last.value); flush.console()\n',
			file = con, sep = "")
		## It is important that one line only is written, so that other clients
		## don't mix in with these lines.
	} else {
		.f <- file()
		on.exit(close(.f))
		..Last.value <- send
		## dump() can stop prematurely if file=con, but also good to remove the /n
		## from dump()'s output before sending (to avoid possible conflicts with
		## other clients)
		dump("..Last.value", file <- .f)
		flush(.f)
		seek(.f, 0)
		cat(readLines(.f), ';', x,
			' <- ..Last.value; rm(..Last.value); cat("\\n<<<endflag>>>\\n"); flush.console()\n',
			file = con, sep = "")
	}
	objdump <- ""
	endloc <- NULL
	while (!length(endloc)) {
		obj <- readLines(con, n = 1000, warn = FALSE)
		## Wait for data to come back. Without this sleep, you get 20-30 calls
		## to readLines before data arrives.
		if (!length(obj)) {
			Sys.sleep(0.01)
			next
		}
		endloc <- grep("<<<endflag>>>", obj)
		if (length(endloc)) obj <- obj[0:(endloc[length(endloc)] - 1)]
		## This is more robust than paste'ing together a potentially very
		## large single string
		objdump <- c(objdump, obj)
	}
	if (!missing(send)) {
		if (!all(objdump == "")) stop(objdump)
		return(TRUE)
	}
	startloc <- grep("<<<startflag>>>", objdump)
	if (!length(startloc))
		stop("Unable to find <<<startflag>>>")
	## The startflag is because sometimes (strangely rarely) seek, flush and dump
	## can write return value to stdout which do not source.
	objdump <- objdump[-(1:startloc[length(startloc)])]
	## Fix any output buffer wrap issues. There are line breaks mid number
	## sometimes which don't source.
	## This is why warn = FALSE appears above in the call to readLines since it
	## warns about these noncomplete lines otherwise.
	nospace <- grep("[^ ]$", objdump)
	nospace <- nospace[nospace < length(objdump)]
	for (i in rev(nospace)) {  # Robust to consecutive lines to be joined
		objdump[i] <- paste(objdump[i], objdump[i + 1], sep = "")
		objdump[i + 1] <- ""
	}
	objcon <- textConnection(objdump)
	on.exit(close(objcon))
	source(objcon, local = TRUE, echo = FALSE, verbose = FALSE)
	return(..Last.value)
}
