## You have to run yad to get this dialog box!
## The form is defined as a list
## TXT, H, RO, NUM, CHK, CB, CBE, FL, SFL, DIR, CDIR, FN, MFL, DT, CLR, BTN or LBL
## TXT = simple text entry (default type)
## H = hidden text (password)
## RO = read-only text
## NUM = null of positive integers with up/down arrows
## CHK = checkbox, return TRUE or FALSE!
## CB = read-only combobox
## CBE = editable combobox
## FL = open file item
## SFL = save file item
## DIR = select directory
## CDIR = select or create directory
## FN = select font and size
## MFL = select mutliple files
## DT = date (use --date-format=%F for YYYY-MM-DD)
## CLR = select a color (but no palette!)
## BTN = display a button... content is shell command to execute!
## LBL = display a label
## ex:
#Form <- list(
#	"Name:TXT" = "some text",
#	"Age:NUM" = 25,
#	"Sex:CB" = c("male", "female"),
#	"Married:CHK"=TRUE
#)

## TODO: implement --image= argument of yad

## Define the S3 method
dlgForm <- function (form, title = "Fill the form", message = NULL,
columns = 1, strip.type = TRUE, ..., gui = .GUI)
{
	if (!gui$startUI("dlgForm", call = match.call(),
		default = lapply(form, "[", 1),
		msg = "Displaying a modal form dialog box",
		msg.no.ask = "A modal form dialog box was by-passed"))
		return(invisible(gui))
	
	## Check and rework main arguments and place them in gui$args
	if (!is.list(form) || is.null(names(form)))
		stop("'form' must be a named list")
	if (!length(title)) title <- "Fill the form" else
		title <- as.character(title)[1]
	if (length(message)) message <- paste(message, sep = "", collapse = "\n")
	columns <- as.integer(columns)
	if (!length(columns)) columns <- 1 else columns <- columns[1]
	if (columns < 1) columns <- 1
	strip.type <- isTRUE(as.logical(strip.type))
	gui$setUI(args = list(form = form, title = title, message = message,
		columns = columns, strip.type = strip.type))
	
	## ... and dispatch to the method
	UseMethod("dlgForm", gui)
}

## Used to break the chain of NextMethod(), searching for a usable method
## in the current context
dlgForm.gui <- function (form, title = "Fill the form", message = NULL,
columns = 1, strip.type = TRUE, ..., gui = .GUI)
{
	msg <- paste("No workable method available to display a form dialog box using:",
		paste(guiWidgets(gui), collapse = ", "))
	gui$setUI(status = "error", msg = msg, widgets = "none")
	stop(msg)
}

## The pure textual version used a fallback in case no GUI could be used
dlgForm.textCLI <- function (form, title = "Fill the form", message = NULL,
columns = 1, strip.type = TRUE, ..., gui = .GUI)
{
	gui$setUI(widgets = "textCLI")
	## Ask for the input for each form field using readline()
	cat("==", gui$args$title, "==\n")
	if (length(gui$args$message)) cat(gui$args$message, "\n", sep = "")
	form <- gui$args$form
	if (is.null(form)) return(list())
	res <- form
	## Special treatment for :CHK type
	isCHK <- grepl(":CHK$", names(form))
	names(form) <- sub(":[A-Z]+$", "", names(form))
	## Do we strip type?
	if (strip.type) res <- form
	for (i in 1:length(form)) {
		if (isCHK[i]) {
			if (isTRUE(as.logical(form[[i]])[1]))
				def <- "[Y]/n" else def <- "[N]/y"
			ans <- readline(paste(names(form)[i], " ", def, ": ", sep = ""))
			if (ans == "") ans <- def
			## Put TRUE if answer starts with y or Y, FALSE otherwise
			res[[i]] <- (grepl("^[yY]", ans))
		} else {
			def <- as.character(form[[i]])
			def[1] <- paste("[", def[1], "]", sep = "")
			res[[i]] <- readline(paste(names(form)[i],
				" ", paste(def, collapse = "/"), ": ", sep = ""))
			## TODO: how to enter multi-items here?
			## TODO: eliminate surrounding single or double quotes
			if (res[[i]] == "") res[[i]] <- form[[i]][1] else
				res[[i]] <- as(res[[i]], class(form[[i]])[1])
		}
	}
	## OK, Cancel or Redo?
	## TODO: the redo feature...
	ans <- readline("==\n[OK]/Cancel ")
	if (ans != "" && !grepl("^[oO]", ans)) res <- list()
	gui$setUI(res = res, status = NULL)
	return(invisible(gui))
}

## The native version of the input box
dlgForm.nativeGUI <- function (form, title = "Fill the form", message = NULL,
columns = 1, strip.type = TRUE, ..., gui = .GUI)
{
	gui$setUI(widgets = "nativeGUI")
	## A simple text input box using native window
	## Return either a string, or character(0) if 'Cancel' clicked
	res <- switch(Sys.info()["sysname"],
		Windows = .winDlgForm(gui$args$form, gui$args$title, gui$args$message,
			gui$args$columns, gui$args$strip.type),
		Darwin = .macDlgForm(gui$args$form, gui$args$title, gui$args$message,
			gui$args$columns, gui$args$strip.type),
		.unixDlgForm(gui$args$form, gui$args$title, gui$args$message,
			gui$args$columns, gui$args$strip.type)
	)
	
	## Do we need to further dispatch?
	if (is.null(res)) NextMethod("dlgForm", gui) else {
		gui$setUI(res = res, status = NULL)
		return(invisible(gui))
	}
}

## Windows version
.winDlgForm <- function (form, title, message, columns, strip.type)
{
	## Not yet => return NULL
	return(NULL)
}

## Mac OS X version
.macDlgForm <- function (form, title, message, columns, strip.type)
{
	## Not yet => return NULL
	return(NULL)
}

## Linux/Unix version
.unixDlgForm <- function (form, title, message, columns, strip.type) 
{
    ## yad must be installed on this machine!
    if (Sys.which("yad") == "") return(NULL)
    ## Avoid displaying warning message in case user clicks on Cancel
    owarn <- getOption("warn")
    on.exit(options(warn = owarn))
    options(warn = -1)
    ## Use yad to display the form dialog box
    if (!length(message)) {
		msg <- paste("yad --form --title=\"", title,
			"\" --on-top --center --skip-taskbar --separator=\"@@@\"",
			" --item-separator=\"&&&\" --date-format=%F --columns=", columns,
			sep = "")
	} else { # Message is provided
		msg <- paste("yad --form --title=\"", title, "\" --text=\"", message,
			"\" --on-top --center --skip-taskbar --separator=\"@@@\"",
			" --item-separator=\"&&&\" --date-format=%F --columns=", columns,
			sep = "")
	}
	## Add the definition of the various fields
	fields <- paste("--field=", shQuote(names(form)), sep = "", collapse = " ")
	## Add the default values
    default <- paste(shQuote(sapply(form, paste, collapse = "&&&")), collapse = " ")
	## Display the dialog box
	res <- system(paste(msg, fields, default), intern = TRUE)
    ## Did the user cancelled the dialog box
	if (!length(res)) return(list())
	## Reformat the result
	res <- strsplit(res, "@@@", fixed = TRUE)[[1]]
	## Replace results in initial form
	for (i in 1:length(form))
		form[[i]] <- as(strsplit(res[i], "&&&", fixed = TRUE)[[1]],
			class(form[[i]])[1])
	## Do we strip type?
	if (strip.type) names(form) <- sub(":[A-Z]+$", "", names(form))
	return(form)
}
