## A SciViews R server using HTTP R help server and JSONP for communcation
## Copyright (c) 2010, Ph. Grosjean (phgrosjean@sciviews.org)
## Use a HTTP request like this in the client:
## http://127.0.0.1:8888/custom/SciViews?msg&callback
## We must return something like (in a correct RJSONp object):
## <callback>({"result":{"String 1", "String 2", "..."}, "options":"<options>",
## "name":"<server.name>","port":"<port>"})
## Another (simpler) way to call it is by using
## http://127.0.0.1:8888/custom/SciViews?msg
## and in this case, the client must manage the simple string returned

## Get list of all names of clients that already connected to the http server
HttpClientsNames <- function ()
	sub("^HttpClient_", "", ls(envir = TempEnv(), pattern = "^HttpClient_"))

## Get or change the port of the http server
HttpServerPort <- function (port)
{
	if (!missing(port)) {
		port <- as.integer(round(port[1]))
		## This port is stored in 'ko.serve' option
		options(ko.serve = port)
		## If the server is running on another port, restart it now
		curport <- getNamespace("tools")$httpdPort
		if (curport > 0 && curport != port) startHttpServer(port = port)
		return(port)
	} else {  # Get the server port
		port <- getOption("ko.serve")
		if (is.null(port)) port <- 8888 else port <- as.integer(round(port[1]))
		return(port)
	}
}

## Get or change the name of the HTTP server
HttpServerName <- function (name)
{
	if (!missing(name)) {
		if (!is.character(name)) stop("'name' must be a string!")
		name <- as.character(name)[1]
		## This name is stored in the option R.id
		options(R.id = name)
		return(name)
	} else {  # Get the server name
		name <- getOption("R.id")
		if (is.null(name)) name <- "R"
		return(name)
	}
}

## Get or change http server options
parHttp <- function (client, ...)
{
	if (missing(client)) client <- "default" else
		client <- as.character(client)[1]
	
	## Set or get parameters for a given HTTP client
	serverport <- HttpServerPort()
	
	## No attempt is made to make sure this client exists
	sc <- paste("HttpClient", client, sep = "_")
	if (!exists(sc, envir = TempEnv(), inherits = FALSE,
		mode = "environment")) {
		## Create a new environment with default values
		e <- new.env(parent = TempEnv())
		e$client <- client
		e$serverport <- serverport
		e$prompt <- ":> "    # Default prompt
		e$continue <- ":+ "  # Default continuation prompt
		e$code <- ""         # Current partial code for multiline mode
		e$last <- ""         # String to add at the end of evaluations
		e$echo <- FALSE      # Don't echo commands to the console
		e$flag <- FALSE      # Do not flag pieces of code (not used yet!)
		e$multiline <- TRUE  # Allow for multiline code
		e$bare <- TRUE       # Always start in "bare" mode
		## Note: in bare mode, all other parameters are inactive!
		## and assign it to SciViews:TempEnv
		assign(sc, e, envir = TempEnv())
	} else e <- get(sc, envir = TempEnv(), mode = "environment")
	
	## Change or add parameters if they are provided
	args <- list(...)
	if (l <- length(args)) {
		change.par <- function (x, val, env) {
			if (is.null(x)) return(FALSE)  # Do nothing without a valid name
			if (is.null(val)) {
				suppressWarnings(rm(list = x, envir = env))  # Remove it
				return(TRUE)
			}
			env[[x]] <- val  # Add or change this variable in the environment
			return(TRUE)
		}
		n <- names(args)
		res <- rep(TRUE, l)
		for (i in seq_len(l)) res[i] <- change.par(n[i], args[[i]], e)
		if (any(!res)) warning("Non named arguments are ignored")
	}
	
	## If serverport has changed, update it now
	if (e$serverport != serverport) e$serverport <- serverport
	
	## Return e invisibly
	return(invisible(e))
}

## Stop the SciViews and R HTTP server and eliminate all tracks
stopHttpServer <- function (remove.clients = FALSE)
{
	## Eliminate the SciViews custom process function for HTTP server
	e <- getNamespace("tools")$.httpd.handlers.env
	if ("SciViews" %in% ls(envir = e)) rm(list = "SciViews", envir = e)
	
	## Do we also remove persistent data for clients?
	if (isTRUE(remove.clients))
		rm(list = ls(envir = TempEnv(), pattern = "^HttpClient_"),
			envir = TempEnv())
	
	## Stop the HTTP deamon
	try(tools::startDynamicHelp(FALSE), silent = TRUE)
}

## (Re)start HTTP help server on the choosen port
## TODO: allowing asking and returning results in the RJSON object
## TODO: conversion to UTF-8 encoding of the returned string
startHttpServer <- function (port = HttpServerPort(),
name = HttpServerName())
{
	if (!is.character(name)) stop("'name' must be a string!")
	name <- as.character(name)[1]
	
	## The port on which we want to run it
	if (!is.numeric(port[1]) || port[1] < 1)
		stop("'port' must be a positive integer!")
	port <- as.integer(round(port[1]))
	## The port on which the server currently runs
	curport <- getNamespace("tools")$httpdPort
	
	## Can we run the server?
	if (curport == -1L || nzchar(Sys.getenv("R_DISABLE_HTTPD")))
		stop("R http server is disabled or cannot start")
	
	## If it is currently running, stop it now
	if (curport != 0L) {
		if (curport != port)
			warning("R http server currently running on port ", curport,
				" and is restarted on port ", port, immediate. = TRUE)
		curport <- stopHttpServer()
	}
	
	## Start the http server on the right port
	if (curport == 0L) {
		oports <- getOption("help.ports")
		(on.exit(options(help.ports = oports)))
		options(help.ports = port)
		curport <- tools::startDynamicHelp()
	} else stop("Unable to start the http server")
	
	## Is the HTTP server running on the right port now?
	if (curport == port) {
		## Set the name of the HTTP server (for easier identification)
		HttpServerName(name)
		
		## Install the SciViews function that will process our requests
		e <- getNamespace("tools")$.httpd.handlers.env
		e[["SciViews"]] <- function (path, query, body, ...) {
			## Analyze the query: command + callback
			#cat(query, "\n", sep = " -- ")
			msg <- query[1]

#			## Strings are supposed to be send in UTF-8 format
#			Encoding(msg) <- "UTF-8"
#			msg <- enc2native(msg)

			l <- length(query)
			if (l == 1) callback <- NULL else {
				callback <- query[l]
#				Encoding(callback) <- "UTF-8"
			}
			
			## The HTTP request message cannot be too long.
			## So, for submission of very long R code, this mechanism
			## is not appropriate. Here we use a specially formatted msg
			## indicating that we should read code from a file instead.
			if (regexpr("^SOURCE=", msg) > 0) {
				srcfile <- sub("^SOURCE=", "", msg)
				on.exit(try(unlink(srcfile), silent = TRUE))
				if (!file.exists(srcfile) || inherits(msg <-
					try(readLines(srcfile, warn = FALSE, encoding = "UTF-8"),
						silent = TRUE), "try-error")) {
					res <- paste(
						gettext("Error: missing or unreadable source file"),
						" '", srcfile, "'\n", sep = "")
					cat(res)
					if (is.null(callback)) {
						return(NULL)
					} else return(Rjsonp(NULL, callback))
				} else msg <- paste(msg, collapse = "\n")
			}
			
			## Get the server name and port, and R encoding
			servername <- HttpServerName()
			serverport <- HttpServerPort()
			
			## Process the command in a similar way as processSocket() does
			## in the svSocket package... but return a RJSONP object if callback
			## is not NULL.
			## We use a custom function here to create this object faster than
			## by converting an R object to RJSON.
			Rjsonp <- function (res, callback) {
				## If no echo, return only a basic RJSONP object
				if (!returnResults || is.null(res)) {
					obj <- paste(callback,
						'(list("result" := NA, ',
						'"options" := list("echo" := FALSE), "name" := "',
						servername, '", "port" := ', serverport, '))', sep = "")
				} else {
					## Return a more consistent RJSONP object
					## Format main client options as a RJSON object
					options <- paste('list("echo" := ', pars$echo,
						', "bare" := ', pars$bare,
						', "partial" := ', (pars$code != ""), ')', sep = "")
					
					## Replace \n by \\n, etc. in res
					#res <- gsub("\n", "\\n", res, fixed = TRUE)
					res <- encodeString(res, quote = '"')
					
					## Check encoding and provide it if it is not UTF-8
					## No, provide it all the time!
					cs <- localeToCharset()[1]
					if (cs != "UTF-8") {
						encode <- paste (', "encoding" := "', cs, '"', sep = "")
					} else encode <- ""
					
					## Format the answer as a RJSONP object and return it
					obj <- paste(callback, '(list("result" := c(',
						paste(shQuote(res, type = "cmd"), collapse = ", "),
						'), "options" := ', options,
						', "name" := "', servername,
						'", "port" := ', serverport, encode, '))', sep = "")
					## Encode this string as UTF-8
					obj <- enc2utf8(obj)
				}
				#cat(obj, "\n")
				return(list(obj))
			}
			
			## Do we receive an <<<id=myID>>> sequence (name of the client)?
			if (regexpr("^<<<id=[a-zA-Z0-9]+>>>", msg) > 0) {
				## Get the identifier
				client <- sub("^<<<id=([a-zA-Z0-9]+)>>>.*$", "\\1", msg)
				## ... and eliminate that sequence
				msg <- sub("^<<<id=[a-zA-Z0-9]+>>>", "", msg)
			} else {
				## The client name is simply 'default'
				client <- "default"
			}
			
			## Do we receive <<<esc>>>? => break (currently, only breaks
			## multiline mode)
			if (substr(msg, 1, 9) == "<<<esc>>>") {
				pars <- parHttp(client, code = "")  # Reset multiline code
				msg <- substr(msg, 10, 1000000)
			}
			
			## Replace <<<n>>> by \n (for multiline code)
			msg <- gsub("<<<n>>>", "\n", msg)
			
			## Replace <<<s>>> by the corresponding client id and server port
			msg <- gsub("<<<s>>>", paste('"', client, '", ', serverport,
				sep = ""), msg)
			
			hiddenMode <- FALSE
			returnResults <- TRUE
			## If msg starts with <<<Q>>> or <<<q>>>, then disconnect server
			## before or after evaluation of the command, respectively
			## Since we always disconnect AFTER with http server, these options
			## have no effect here. They are used with the socket server only
			## If msg starts with <<<e>>>, evaluate command in the console and
			## disconnect
			## If msg starts with <<<h>>> or <<<H>>>, evaluate in hidden mode
			## and disconnect
			startmsg <- substr(msg, 1, 7)
			if (startmsg == "<<<Q>>>") {
				msg <- substr(msg, 8, 1000000)
				returnResults <- FALSE
			} else if (startmsg == "<<<q>>>") {
				msg <- substr(msg, 8, 1000000)
				parHttp(client, last = "")
			} else if (startmsg == "<<<e>>>") {
				msg <- substr(msg, 8, 1000000)
				## We just configure the server correctly
				parHttp(client, bare = FALSE, echo = TRUE, prompt = ":> ",
					continue = ":+ ", multiline = TRUE, last = "")
				## Add a command to the command history
				#timestamp("my R command", "", "", quiet = TRUE)
			} else if (startmsg == "<<<h>>>") {
				msg <- substr(msg, 8, 1000000)
				## Do not echo command on the server (silent execution)
				hiddenMode <- TRUE
				parHttp(client, bare = TRUE, last = "")
			} else if (startmsg == "<<<H>>>") {
				msg <- substr(msg, 8, 1000000)
				## Do not echo command on the server
				hiddenMode <- TRUE
				returnResults <- FALSE
				parHttp(client, bare = TRUE)
			} else if (startmsg == "<<<u>>>") {
				msg <- substr(msg, 8, 1000000)
				## Silent execution, nothing is returned to the client
				## (but still echoed to the server)
				hiddenMode <- FALSE
				returnResults <- FALSE
				parHttp(client, bare = TRUE)
			}
			
			## Get parameters for the client
			pars <- parHttp(client)
			if (Bare <- pars$bare) {
				Prompt <- ""
				Continue <- ""
				Echo <- FALSE
			} else {
				Prompt <- pars$prompt
				Continue <- pars$continue
				Echo <- pars$echo
			}
			## TODO: do we still need this?
			## Eliminate last carriage return
			msg <- sub("(.*)[\n][^\n]*$", "\\1", msg)
			if (!hiddenMode) {
				if (Echo) {
					## Note: command lines are now echoed directly in captureAll()
					## => no need of this any more!
					if (pars$code == "") Pre <- Prompt else Pre <- Continue
					#cat(Pre, msg, "\n", sep = "")
				}
				## Add previous content if we were in multiline mode
				if (pars$code != "") msg <- paste(pars$code, msg, sep = "\n")
				pars$code <- ""  # This changes the original data too!
			}
			
			## Parse the R code
			expr <- parseText(msg)
			## Is it a wrong code?
			if (inherits(expr, "try-error")) {
			    res <- paste(ngettext(1, "Error: ", "", domain = "R"),
			    sub("^[^:]+: ([^\n]+)\n[0-9]+:(.*)$", "\\1\\2", expr), sep = "")
			    if (Echo) cat(res)
			    if (is.null(callback)) {
					ret <- paste(res, pars$last, Prompt, sep = "")
					## Encode this as UTF-8
					ret <- enc2utf8(ret)
					return(ret)
				} else {
					return(Rjsonp(paste(res, pars$last, Prompt, sep = ""),
						callback))
				}
			}
			## Is it incomplete code?
			if (!is.expression(expr)) {
				## Is multiline mode allowed?
				if (!Bare && pars$multiline) {
					pars$code <- msg
					if (is.null(callback)) {
						if (returnResults) {
							ret <- paste(pars$last, Continue, sep = "")
							## Encode this as UTF-8
							ret <- enc2utf8(ret)
							return(ret)
						} else return(NULL)	
					} else {
						if (returnResults) {
							return(Rjsonp(paste(pars$last, Continue, sep = ""),
								callback))
						} else return(Rjsonp(NULL, callback))
					}
				} else {  # Multimode not allowed
				    res <- paste(
						gettext("Error: incomplete command in single line mode"),
						"\n", sep = "")
					if (Echo) cat(res)
					if (is.null(callback)) {
						if (returnResults) {
							ret <- paste(res, pars$last, Prompt, sep = "")
							## Encode this as UTF-8
							ret <- enc2utf8(ret)
							return(ret)
						} else return(NULL)
					} else {
						if (returnResults) {
							return(Rjsonp(paste(res, pars$last, Prompt,
								sep = ""), callback))
						} else return(Rjsonp(NULL, callback))
					}
				}
			}
			## Freeze parameters (unlinks from the environment)
			pars <- as.list(pars)
			## Is it something to evaluate?
			if (length(expr) < 1) {
				if (is.null(callback)) {
					ret <- paste(pars$last, Prompt, sep = "")
					## Encode this as UTF-8
					ret <- enc2utf8(ret)
					return(ret)
				} else {
					return(Rjsonp(paste(pars$last, Prompt, sep = ""), callback))
				}
			}
			## Correct code,... we evaluate it
			results <- captureAll(expr, echo = Echo, split = Echo)
			## Should we run taskCallbacks?
			# Note: these are installed in svKomodo package
			if (!hiddenMode) {
				h <- getTemp(".svTaskCallbackManager", default = NULL,
					mode = "list")
				if (!is.null(h)) h$evaluate()
			}
			## Collapse and add last and the prompt at the end
			results <- paste(results, collapse = "\n")
			#if (Echo) cat(results)
			if (!returnResults) {
				if (is.null(callback)) {
					return(NULL)
				} else {
					return(Rjsonp(NULL, callback))
				}
			}
			Prompt <- if (pars$bare) "" else pars$prompt
			results <- paste(results, pars$last, Prompt, sep = "")
			## Return the results in plain text, or RJSONP object
			if (is.null(callback)) {
				results <- enc2utf8(results)
				return(results)
			} else {
				return(Rjsonp(results, callback))
			}
		}
	}	
	return(invisible(curport))
}

.onLoad <- function (lib, pkg)
{
	## Try starting the Http server
	#try(startHttpServer(), silent = TRUE)
}

.onUnload <- function (libpath)
{
	## Make sure that the SciViews Http server is closed
	stopHttpServer(TRUE)
}

.packageName <- "svHttp"
