## A Tcl/Tk version of the message box
dlgMessage.tcltkGUI <- function (message, type = c("ok", "okcancel", "yesno",
"yesnocancel"), ..., gui = .GUI)
{
    gui$setUI(widgets = "tcltkGUI")
	## TODO: use custom dialog boxes, not tkmessageBox() cf. ugly in Linux
	## and buggy in Mac!!!
    ## Use tkmessageBox for now
    tkDefault <- switch(gui$args$type,
        "ok" = "ok",
        "okcancel" = "ok",
        "yesno" = "yes",
        "yesnocancel" = "yes",
        stop("'type' must be \"ok\", \"okcancel\", \"yesno\" or \"yesnocancel\"")
    )
    if (type == "ok") {
        alarm()
        tkIcon <- "info"
        tkTitle <- "Information"
    } else {
        tkIcon <- "question"
        tkTitle <- "Question"
    }
    res <- tkmessageBox(message = gui$args$message, title = tkTitle,
		type = gui$args$type, default = tkDefault, icon = tkIcon)
    res <- tolower(tclvalue(res))
	## Under Windows, we have to bring the R Console to top manually
	## TODO: this does not work in RTerm... but it is a minor ennoyance
	if (.Platform$OS.type == "windows") bringToTop(-1)
    gui$setUI(res = res, status = NULL)
	return(invisible(gui))
}
