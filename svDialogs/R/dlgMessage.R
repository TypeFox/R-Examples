## Simplified versions of dlgMessage()
msgBox <- function (message)
{
	require(svDialogs)
	dlgMessage(message = message)$res
}

okCancelBox <- function (message)
{
	require(svDialogs)
	return(dlgMessage(message = message, type = "okcancel")$res == "ok")
}

## Define the S3 method
dlgMessage <- function (message, type = c("ok", "okcancel", "yesno",
"yesnocancel"), ..., gui = .GUI)
{
	if (!gui$startUI("dlgMessage", call = match.call(), default = "ok",
		msg = "Displaying a modal message dialog box",
		msg.no.ask = "A modal message dialog box was by-passed"))
		return(invisible(gui))
	
	## Check and rework main arguments and place them in gui$args
    if (missing(message)) message <- "[Your message here...]"
    message <- paste(as.character(message), collapse = "\n")
    type <- match.arg(type)
	gui$setUI(args = list(message = message, type = type))
	
	## ... and dispatch to the method
	UseMethod("dlgMessage", gui)
}

## Used to break the chain of NextMethod(), searching for a usable method
## in the current context
dlgMessage.gui <- function (message, type = c("ok", "okcancel", "yesno",
"yesnocancel"), ..., gui = .GUI)
{
	msg <- paste("No workable method available to display a message dialog box using:",
		paste(guiWidgets(gui), collapse = ", "))
	gui$setUI(status = "error", msg = msg, widgets = "none")
	stop(msg)
}

## The pure textual version used a fallback in case no GUI could be used
dlgMessage.textCLI <- function (message, type = c("ok", "okcancel", "yesno",
"yesnocancel"), ..., gui = .GUI)
{
	gui$setUI(widgets = "textCLI")	
	## Display the message and wait for user action
	if (gui$args$type == "ok") {
        readline(paste(gui$args$message, "\n(hit ENTER to continue) ", sep = ""))
		res <- "ok"
    } else {
        ## Use a non-graphical select.list() for the others
        choices <- switch(gui$args$type,
            okcancel = c("ok", "cancel"),
            yesno = c("yes", "no"),
            yesnocancel = c("yes", "no", "cancel")
        )
        res <- select.list(choices, title = gui$args$message, graphics = FALSE)
        if (res == "" && type != "yesno") res <- "cancel"
        if (res == "") res <- "no" # Selection of 0 with yes/no => no
	}
	gui$setUI(res = res, status = NULL)
	return(invisible(gui))
}

## The native version of the message box
dlgMessage.nativeGUI <- function (message, type = c("ok", "okcancel", "yesno",
"yesnocancel"), ..., gui = .GUI)
{
	gui$setUI(widgets = "nativeGUI")
	## A simple message box
    ## type can be 'ok' (info), 'okcancel', 'yesno', 'yesnocancel' (question)
    ## This dialog box is always modal
    ## Returns invisibly a character with the button that was pressed
	res <- switch(Sys.info()["sysname"],
		Windows = .winDlgMessage(gui$args$message, gui$args$type),
		Darwin = .macDlgMessage(gui$args$message, gui$args$type),
		.unixDlgMessage(gui$args$message, gui$args$type)
	)
	
	## Do we need to further dispatch?
	if (is.null(res)) NextMethod("dlgMessage", gui) else {
		gui$setUI(res = res, status = NULL)
		return(invisible(gui))
	}
}

## Windows version
.winDlgMessage <- function (message, type = c("ok", "okcancel", "yesno",
"yesnocancel"))
{
	res <- winDialog(type = type, message = message)
	## Rework result to match the result from the other functions
	if (type == "ok") return(invisible("ok")) else return(tolower(res))
}

## Mac OS X version	
.macDlgMessage <- function (message, type= c("ok", "okcancel", "yesno",
"yesnocancel"))
{
    ## Display a modal message with native Mac dialog box
	#if (.Platform$GUI == "AQUA") app <- "(name of application \"R\")" else
	## This works from Mac OS X 10.5 Leopard:
	if (.Platform$GUI == "AQUA") app <- "(name of application id \"Rgui\")" else
		if (.isJGR()) app <- "\"JGR\"" else app <- "\"Terminal\""
	type <- match.arg(type)
    buttons <- switch(type,
		ok = "\"OK\"",
		okcancel = "\"Cancel\",\"OK\"",
		yesno = "\"No\",\"Yes\"",
		yesnocancel = ",\"Cancel\",\"No\",\"Yes\"",
		stop("type can only be 'ok'n 'okcancel', 'yesno', 'yesnocancel'"))
	if (type == "ok") {
		beep <- " -e 'beep'"
		icon <- "caution"
		title <- "\"Information\""
		more <- " default button 1"
	} else {
		beep <- ""
		icon <- "note"
		title <- "\"Question\""
		if (type == "yesnocancel")
            more <- " default button 3 cancel button 1" else
			if (type == "yesno") more <- " default button 2" else
			more <- " default button 2 cancel button 1"
	}
	cmd <- paste("exit `osascript", beep, " -e 'tell application ", app,
		" to set dlg to display dialog \"", message, "\" with title ", title,
		more, " with icon ", icon, " buttons {", buttons,
		"}' -e 'if button returned of dlg is \"No\" then 2' 2> /dev/null`",
		sep = "")
	res <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE, wait = TRUE)
	## Decrypt result
	if (type == "ok")
		if (res > 0) return(NULL) else return(invisible("ok"))
	if (res == 2) return("no")
	if (res == 1) return("cancel")
	if (type == "okcancel") return("ok") else return("yes")
}

## Linux/Unix version
.unixDlgMessage <- function (message, type = c("ok", "okcancel", "yesno",
"yesnocancel"))
{
    ## zenity must be installed on this machine!
    if (Sys.which("zenity") == "") return(NULL)
    type <- match.arg(type)
    if (type == "ok") {
        alarm()
        msg <- paste("zenity --info --text=\"", message,
            "\" --title=\"Information\"", sep = "")
        res <- system(msg)
        if (res > 0) return(NULL) else return(invisible("ok"))
    } else if (type == "yesnocancel") {
        type <- "yesno"
        confirm <- TRUE
    } else confirm <- FALSE
    ## Now, we have only "okcancel" or "yesno"
    if (type == "okcancel") {
        msg <- paste("zenity --question --text=\"", message,
            "\" --ok-label=\"OK\" --cancel-label=\"Cancel\" --title=\"Question\"",
            sep = "")
        results <- c("ok", "cancel")
    } else {
        msg <- paste("zenity --question --text=\"", message,
            "\" --ok-label=\"Yes\" --cancel-label=\"No\" --title=\"Question\"",
            sep = "")
        results <- c("yes", "no")
    }
    res <- system(msg)
    if (res > 1) return(NULL) else res <- results[res + 1]
    ## Do we ask to continue (if was yesnocancel)?
    if (confirm) {
        conf <- system(paste("zenity --question --text=\"Continue?\"",
			"--ok-label=\"OK\" --cancel-label=\"Cancel\" --title=\"Confirm\""))
        if (conf == 1) return("cancel")
    }
    return(res)
}
