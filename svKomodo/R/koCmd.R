koCmd <- function (cmd, data = NULL, async = FALSE, host = getOption("ko.host"),
port = getOption("ko.port"), kotype = getOption("ko.kotype"), timeout = 2,
type = c("js", "rjsonp", "output"), pad = NULL, ...)
{

    type <- match.arg(type)
	if (is.null(host)) host <- "localhost"  # Default value
	if (is.null(port)) port <- 7052         # Idem
	if (is.null(kotype)) kotype <- "file"	# Idem
	cmd <- gsub("\n", "\\\\n", cmd)
	cmd <- paste(cmd, collapse = " ")
    if (is.na(cmd) || is.null(cmd) || length(cmd) == 0) {
		warning("No command supplied in cmd argument")
		return("")
    }
    ## Do we need to paste data in the command?
	if (!is.null(data)) {
		"rework" <- function(data) {
			data <- as.character(data)
			data <- gsub("\n", "\\\\\\\\n", data)
			data <- paste(data, collapse = "\\\\n")
			return(data)
		}

		n <- names(data)
		if (is.null(n)) {
			## We assume that we replace '<<<data>>>'
			cmd <- gsub("<<<data>>>", rework(data), cmd)
		} else {	# Named data
			## We replace each <<<name>>> in turn
			for (i in 1:length(n))
				cmd <- gsub(paste("<<<", n[i], ">>>", sep = ""),
					rework(data[[n[i]]]), cmd)
		}
	}
	## What type of data do we send?
	cmd <- switch(type,
		js = paste("<<<js>>>", cmd, sep = ""),
		rjsonp = paste("<<<rjsonp>>>", pad, "(",
			paste(toRjson(cmd, ...), collapse = " "), ")", sep = ""),
		cmd)
		
	otimeout <- getOption("timeout")
	options(timeout = timeout)  # Default timeout is 120 seconds
	## Do we use file or socket server?
	if (kotype == "file") {
		## File is .sv<port> in the temporary dir (/tmp for Linux and Mac OS X)
		tempdir <- "/tmp"
		if (.Platform$OS.type == "windows") tempdir <- dirname(tempdir())
		svfile <- file.path(tempdir, paste(".sv", port, sep = ""))
		svoutfile <- paste(svfile, "out", sep = ".")
		## Make sure the output file is deleted
		unlink(svoutfile)
		## Send data to Komodo
		## Need a time-of-day fingerprint in ms
		fp <- as.integer((as.numeric(Sys.time()) %% 24*60*60) * 1000)
		fp <- gsub(" ", "0", format(fp, width = "8", scientific = FALSE))
		cmd <- paste("<<<", fp, ">>>", cmd, sep = "")
		svcon <- file(svfile, open = "w", blocking = TRUE, encoding = "UTF-8")
		try(cat(paste(cmd, collapse = "\n"), file = svcon), silent = TRUE)
		try(close(svcon), silent = TRUE)
		## Wait for output file... or max timeout (add 500ms to this timeout)
		timeout <- Sys.time() + timeout + 0.5
		while (!file.exists(svoutfile) && Sys.time() < timeout)
			Sys.sleep(0.2) # To avoid using too much resources
		if (file.exists(svoutfile)) {
			res <- try(readLines(svoutfile), silent = TRUE)
			Encoding(res) <- "UTF-8"
		} else res <- character(0)	
	} else { # This must be a socket server
		tryCatch({
				con <- socketConnection(host = host, port = port,
					blocking = TRUE)
				writeLines(cmd, con)
				res <- readLines(con)
				close(con)
			}, warning = function (e) {
				stop(simpleError("Komodo socket server is not available!",
				quote(koCmd)))
		})
#		ret <- try(writeLines(cmd, con), silent = TRUE)
#		if (!inherits(ret, "try-error"))
#			res <- try(readLines(con), silent = TRUE)
#		try(close(con), silent = TRUE)
	}
	options(timeout = otimeout)
    return(res)
}
