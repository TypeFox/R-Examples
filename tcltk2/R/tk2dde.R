### tk2dde.R - R functions to interface the dde tcl package provided in the
### standard Tcl installation under R (and in the ActiveState version) for Windows
### Copyright (c), Philippe Grosjean (phgrosjean@sciviews.org)
### Licensed under LGPL 3 or above
###
### Changes:
### - 2007-01-01: fisrt version (for tcltk2_1.0-0)
###
### To do:
### - Revise this code and secure it better where it could be (arguments, Tcl)
### - Add the dde eval command to evaluate a script in a different Tcl interpreter

.tk2dde.require <- function ()
{
	if (.Platform$OS.type != "windows")
		stop("This is a Windows-specific function!")
    ## Make sure tcl/tk dde is operational
	if (!capabilities("tcltk"))
		stop("This version of R cannot use Tcl/Tk!")
	res <- tclRequire("dde", warn = TRUE)
	if (inherits(res, "tclObj")) res <- tclvalue(res)
	if (res[1] == FALSE)
		stop("Unable to find the 'dde' Tcl/tk package!")
	return(res)  # The package version number

}

tk2dde <- function (topic = NULL)
{
    ## Initialize a tcltk dde server with name 'TclEval|topic'
    .tk2dde.require()

    ## If topic is NULL, just get my server name
    if (is.null(topic)) return(tclvalue(.Tcl("dde servername {}")))

    ## Otherwise topic must be character
    topic <- topic[1]
    if (!is.character(topic) || topic == "")
		stop("'topic' must be a non null character string!")

    ## Verify if I am not already registered under this topic
    if (tclvalue(.Tcl("dde servername {}")) == topic) return(0)	# OK

    ## Check that this server name does not exist yet
    if (length(grep(paste("[{]TclEval ", topic, "[}]", sep = ""),
		as.character(.Tcl("dde services TclEval {}")))) > 0)
    	return(1)	# This server name already exists => return 1 and don't set!

    ## Register me as a dde server with this topic name
    .Tcl(paste("dde servername", topic))
    ## Check that the server is set correctly
	## (if not, return 2 to warn that a problem occurred)
    if (tclvalue(.Tcl("dde servername {}")) == topic) return(0) else return(2)
}

tk2dde.exec <- function (service, topic, command, async = FALSE)
{
    ## Execute a command in the 'service|topic' dde server
    .tk2dde.require()

    if (!is.character(service) || !is.character(topic) || !is.character(command))
		stop("'service', 'topic' and 'command' must be character strings!")
    if (async[1] == TRUE) async <- "-async" else async <- ""

    ## Execute the command in a try(), to nicely catch the error
    ## class is "try-error" if an error occurs, otherwise, returns ""
    res <- (try(tclvalue(.Tcl(paste("dde execute ", async, " ",
		as.character(service[1]), " ", as.character(topic[1]), " ",
		as.character(command[1]), sep = "")))))
    return(res)
}

tk2dde.poke <- function (service, topic, item, data)
{
    ## Set a value (data) to 'item' in the 'service|topic' dde server's app
    .tk2dde.require()

    if (!is.character(service) || !is.character(topic))
		stop("'service' and 'topic' must be character strings!")
    if (!is.character(item))
		stop("'item' must be character strings!")
	## In Tcl, if 'data' is a character string, enclose it in curly braces
    data <- paste("{", paste(as.character(data), collapse = "\n"), "}", sep = "")

    ## For some reasons, dde poke does not seem to work with a TclEval serve...
	## use dde execute instead
    if (service == "TclEval") {
        Cmd <- paste("{set ", as.character(item[1]), " ", data, "}", sep = "")
		## This would not work with all kind of data!!!
        ## Also, if it is a vector, matrix, or array, it does not work properly!
        return(tk2dde.exec(service, topic, Cmd, async = TRUE))
    }

    ## Poke the data within a try(), to nicely catch the error
    ## class is "try-error" if an error occurs, otherwise, returns ""
    res <- (try(as.character(.Tcl(paste("dde poke", as.character(service[1]),
		as.character(topic[1]), as.character(item[1]), data)))))
    return(res)
}

tk2dde.request <- function (service, topic, item, binary = FALSE)
{
    ## Get the value for 'item' in 'service|topic' dde server
    .tk2dde.require()

    if (!is.character(service) || !is.character(topic))
		stop("'service' and 'topic' must be character strings!")
    if (!is.character(item))
		stop("'item' must be character strings!")
    if (binary[1] == TRUE) binary <- "-binary" else binary <- ""

    ## Request the value in a try(), to nicely catch the error
    ## class is "try-error" if an error occurs, otherwise, returns ""
    res <- (try(as.character(.Tcl(paste("dde request ", binary, " ",
		as.character(service[1]), " ", as.character(topic[1]), " ",
		as.character(item[1]), sep = "")))))
    return(res)
}

tk2dde.services <- function (service = "", topic = "")
{
    ## List the 'service|topic' dde currently available
    .tk2dde.require()

    ## Check arguments
    if (!is.character(service) || !is.character(topic))
		stop("'service' and 'topic' must be character strings!")
    service <- as.character(service[1])
    if (service == "") service <- "{}"	# This is an empty string in Tcl
    topic <- as.character(topic[1])
    if (topic == "") topic <- "{}"		# This is an empty string in Tcl

    ## Get the list of all 'service|topic' dde servers currently running
    res <- as.character(.Tcl(paste("dde services", service, topic)))
    return(res)
}
