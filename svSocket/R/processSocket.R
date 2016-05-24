processSocket <- function (msg, socket, serverport, ...)
{
	## This is the default R function that processes a command send by a socket
    ## client. 'msg' is assumed to be R code contained in a string

	## Strings are supposed to be send in UTF-8 format
#	Encoding(msg) <- "UTF-8"
#	msg <- enc2native(msg)

	## Special case of a HEAD HTTPrequest
	## TODO: rework this, return a header and force socket disconnection on the server side
	## HEAD should return only the header with
	## HTTP/1.1 200 OK
	## Server: R socket server
	## Connection: close
	## Content-type: text/palin;charset=UTF-8
	## (followed by an empty line)
	## Get does the same, but process the R command and returns processed results
	## in the body, i.e, after the empty line
	## Respond immediately
	if (regexpr("^(HEAD|GET) /custom/SciViews\\?.* HTTP/.*$", msg) > 0 ) {
		## Get the code to execute
		code <- sub("^(HEAD|GET) /custom/SciViews\\?(.*) HTTP/.*$", "\\2", msg)
		code <- URLdecode(code)
		## Do we receive a <<<id=myID>>> sequence?
		if (regexpr("^<<<id=[a-zA-Z0-9]+>>>", code) > 0) {
			## Get the identifier
			client <- sub("^<<<id=([a-zA-Z0-9]+)>>>.*$", "\\1", code)
			## ... and eliminate that sequence
			code <- sub("^<<<id=[a-zA-Z0-9]+>>>", "", code)
		} else {
			## The client name is simply the socket name
			client <- socket
		}
		## Do some replacements
		## Replace <<<n>>> by \n (for multiline code)
		code <- gsub("<<<n>>>", "\n", code)
		## Replace <<<s>>> by the corresponding client id and server port
		code <- gsub("<<<s>>>", paste('"', client, '", ', serverport, sep = ""),
			code)
		## We ignore and eliminate the other <<<xxx>>> sequences
		code <- gsub("<<<[^>]+>>>", "", code)
		
		## Parse and execute this code
		expr <- parseText(code)
		results <- try(captureAll(expr, echo = FALSE, split = FALSE),
			silent = TRUE)

		## Make sure to return something different than "" (this is used to
		## sense if the R server is socket or HTTP)
		if (is.null(results) || is.na(results) || !length(results) ||
			results == "") results <- " "
		
		## Return captured results
#		return(enc2utf8(results))
		return(results)
	}

	## These are other key: value lines send by HTTP clients... just ignore
	if (regexpr("^[-a-zA-Z]+: ", msg) > 0) return("")

    ## Do we receive a <<<id=myID>>> sequence?
	if (regexpr("^<<<id=[a-zA-Z0-9]+>>>", msg) > 0) {
		## Get the identifier
		client <- sub("^<<<id=([a-zA-Z0-9]+)>>>.*$", "\\1", msg)
		## ... and eliminate that sequence
		msg <- sub("^<<<id=[a-zA-Z0-9]+>>>", "", msg)
	} else {
		## The client name is simply the socket name
		client <- socket
	}

	## Do we receive <<<esc>>>? => break (currently, only break multiline mode)
    if (substr(msg, 1, 9) == "<<<esc>>>") {
		## Reset multiline code and update clientsocket
        pars <- parSocket(client, serverport, clientsocket = socket, code = "")
        msg <- substr(msg, 10, 1000000)
    }

    ## Replace <<<n>>> by \n (for multiline code)
    msg <- gsub("<<<n>>>", "\n", msg)

    ## Replace <<<s>>> by the corresponding client identifier and server port
    msg <- gsub("<<<s>>>", paste('"', client, '", ', serverport, sep = ""), msg)

    hiddenMode <- FALSE
	returnResults <- TRUE
	## If msg starts with <<<Q>>> or <<<q>>>, then disconnect server before or
    ## after evaluation of the command, respectively
    ## If msg starts with <<<e>>>, evaluate command in the console and disconnect
    ## If msg starts with <<<h>>> or <<<H>>>, evaluate in hidden mode and disconnect
    startmsg <- substr(msg, 1, 7)
    if (startmsg == "<<<Q>>>") {
		msg <- substr(msg, 8, 1000000)
        ## Indicate to the client that he can disconnect now
        closeSocketClients(sockets = socket, serverport = serverport)
		returnResults <- FALSE
    } else if (startmsg == "<<<q>>>") {
		msg <- substr(msg, 8, 1000000)
        ## Remember to indicate disconnection at the end
        parSocket(client, serverport, clientsocket = socket, last = "\n\f")
    } else if (startmsg == "<<<e>>>") {
        msg <- substr(msg, 8, 1000000)
        ## We just configure the server correctly
        parSocket(client, serverport, clientsocket = socket, bare = FALSE,
			echo = TRUE, prompt = ":> ", continue = ":+ ", multiline = TRUE,
			last = "\n\f")
        ## Add a command to the command history
        #timestamp("my R command", "", "", quiet = TRUE)
    } else if (startmsg == "<<<h>>>") {
		msg <- substr(msg, 8, 1000000)
		## Do not echo command on the server (silent execution)
		hiddenMode <- TRUE
		parSocket(client, serverport, clientsocket = socket, bare = TRUE,
			last = "\n\f")
    } else if (startmsg == "<<<H>>>") {
		msg <- substr(msg, 8, 1000000)
		## Do not echo command on the server (silent execution with no return)
        closeSocketClients(sockets = socket, serverport = serverport)
		hiddenMode <- TRUE
		returnResults <- FALSE
		parSocket(client, serverport, clientsocket = socket, bare = TRUE)
	} else if (startmsg == "<<<u>>>") {
		msg <- substr(msg, 8, 1000000)
		## Silent execution, nothing is returned to the client
		## (but still echoed to the server)
		hiddenMode <- FALSE
		returnResults <- FALSE
		parSocket(client, serverport, clientsocket = socket, bare = TRUE)
	}

    ## Get parameters for the client
    pars <- parSocket(client, serverport)
    if (Bare <- pars$bare) {
        Prompt <- ""
        Continue <- ""
        Echo <- FALSE
    } else {
        Prompt <- pars$prompt
        Continue <- pars$continue
        Echo <- pars$echo
    }
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
#        return(enc2utf8(paste(res, pars$last, Prompt, sep = "")))
		return(paste(res, pars$last, Prompt, sep = ""))
    }
    ## Is it incomplete code?
    if (!is.expression(expr)) {
        ## Is multiline mode allowed?
        if (!Bare && pars$multiline) {
            pars$code <- msg
            if (returnResults) {
#				return(enc2utf8(paste(pars$last, Continue, sep = "")))
				return(paste(pars$last, Continue, sep = ""))
			} else return("")
        } else {  # Multimode not allowed
            res <- paste(gettext("Error: incomplete command in single line mode"),
                "\n", sep = "")
            if (Echo) cat(res)
			if (returnResults) {
#				return(enc2utf8(paste(res, pars$last, Prompt, sep = "")))
				return(paste(res, pars$last, Prompt, sep = ""))
			} else return("")
        }
    }
	## Freeze parameters (unlinks from the environment)
	pars <- as.list(pars)
    ## Is it something to evaluate?
    if (length(expr) < 1) {
#		return(enc2utf8(paste(pars$last, Prompt, sep = "")))
		return(paste(pars$last, Prompt, sep = ""))
	}
    ## Correct code,... we evaluate it
	## Something like this should allow for real-time echo in client, but it is too slow
	## and it outputs all results at the end...
	#results <- captureAll(expr, split = Echo, file = socketClientConnection(socket))
	results <- captureAll(expr, echo = Echo, split = Echo)
	## Should we run taskCallbacks?
	## Note: these are installed in svKomodo package
	if (!hiddenMode) {
		h <- getTemp(".svTaskCallbackManager", default = NULL, mode = "list")
		if (!is.null(h)) h$evaluate()
	}
    ## Collapse and add last and the prompt at the end
    results <- paste(results, collapse = "\n")
	#if (Echo) cat(results)
    if (!returnResults) return("")
	Prompt <- if (pars$bare) "" else pars$prompt
    results <- paste(results, pars$last, Prompt, sep = "")
#    return(enc2utf8(results))
	return(results)
}
