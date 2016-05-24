guiCallTip <- function (code, file = NULL, onlyargs = FALSE, width = 60,
location = FALSE)
{
    ## This is an interface to CallTip for external programs
    ## Clear ::SciViewsR_CallTip
    .Tcl("set ::SciViewsR_CallTip {}")

    ## Using a callback, all args are strings => convert
    if (length(file) == 0 || file == "" || file == "NULL") file <- NULL
    onlyargs <- as.logical(onlyargs[1])
    width <- as.integer(width[1])

    ## Get the call tip
	ctip <- callTip(code, only.args = onlyargs, location = location)

    ## Possibly break long lines at reasonables widths
    if (onlyargs) Exdent <- 0 else Exdent <- 4
    if (!is.null(width) && !width < 1)
	   ctip <- paste(strwrap(ctip, width = width, exdent = Exdent), collapse = "\n")

	## Copy the result to a Tcl variable
    .Tcl(paste("set ::SciViewsR_CallTip {", ctip, "}", sep = ""))

    if (!is.null(file)) {  # Copy it also to the clipboard or a file
        ## If file = clipboard and this is Windows, copy it to the clipboard
        if (file == "clipboard") {
            if (.Platform$OS.type == "windows") {
                writeClipboard(ctip)
            } else {
                stop("'clipboard' not supported yet on platforms different than Windows!")
            }
        } else { # copy the call tip to the file
            cat(ctip, file = file)
        }
    }
	invisible(ctip)
}

"guiComplete" <-
function (code, file = NULL, sep = "|")
{
    # This is an interfacte to Complete for external programs
    # Clear ::SciViewsR_Complete
    .Tcl("set ::SciViewsR_Complete {}")

    # Using a callback, all args are strings => convert
    if (length(file) == 0 || file == "" || file == "NULL") file <- NULL
    sep = sep[1]

    # Get the completion list
	clist <- completion(code, sep = sep)

	# Copy the result to a Tcl variable
    .Tcl(paste("set ::SciViewsR_Complete {", clist, "}", sep = ""))

    if (!is.null(file)) { # Copy it also to the clipboard or a file
        # if file = clipboard and this is Windows, copy it to the clipboard
        if (file == "clipboard") {
            if (.Platform$OS.type == "windows") {
                writeClipboard(clist)
            } else {
                stop("'clipboard' not supported yet on platforms different than Windows!")
            }
        } else {  # Copy the completion list to the file
            cat(clist, file = file)
        }
    }
	invisible(clist)
}

guiDDEInstall <- function ()
{
    ## Register a dde server for R and install callbacks for serveur functions

    ## Make sure tcl/tk dde is operational
    if (.Platform$OS.type != "windows")
		return("DDE not installed: this is not Windows!")
	if (!capabilities("tcltk"))
		return("DDE not installed: this version of R cannot use Tcl/Tk!")
	tclRequire("dde", warn = TRUE)
	## Should be installed by default with the tcltk package under Windows

    ## Register a "SciViewsR" server
    topic <- "SciViewsR"
    ## Verify if I am not already registered under this topic
    if (!tclvalue(.Tcl("dde servername {}")) == topic) {
        ## Check that this server name does not exist yet
        if (length(grep(paste("[{]TclEval ", topic, "[}]", sep = ""),
			tclvalue(.Tcl("dde services TclEval {}")))) > 0) {
            mes <- "DDE not installed: server name already in use!"
			return(invisible(mes))
		}
        ## Register me as a dde server with this topic name
        .Tcl(paste("dde servername", topic))
        ## Check that the server is set correctly (if not, return an error)
        if (!tclvalue(.Tcl("dde servername {}")) == topic) {
            mes <- "DDE not installed: unknown error while starting the server!"
			return(invisible(mes))
		}
    }

    ## Install callbacks for guiXXXX functions, for DDE clients to access them
    ## guiCallTip()... Take care: must be adapted if you change guiCallTip()!
    res <- .Tcl.args(guiCallTip)
    .Tcl(paste("proc guiCallTip {code {file \"\"} {onlyargs FALSE}",
		" {width 60} {location FALSE} }", gsub("%", "$", res), sep = ""))

    ## guiComplete()... Take care: must be adapted if you change guiComplete()!
    res <- .Tcl.args(guiComplete)
    .Tcl(paste("proc guiComplete {code {file \"\"}",
		" {sep |} }", gsub("%", "$", res), sep = ""))

    ## Done
    return(invisible(""))  # OK!
}
