### tk2ico.R - Functions to interface the winico Tcl package under Windows
### Copyright (c), Philippe Grosjean (phgrosjean@sciviews.org)
### Licensed under LGPL 3 or above
###
### Changes:
### - 2012-04-28: winico.dll support is dropped; now using 'ico' package
### - 2007-01-01: first version (for tcltk2_1.0-0)

### TODO: implement ::ico::getIconByName, ::ico::getFileIcon & ::ico::writeIcon
### TODO: gif files are acceptable too for tk2ico.set(), example:
### Image <- tclVar()
### tcl("image", "create", "photo", Image, file = "myfile.gif")
### tcl("wm", "iconphoto", tt, Image) instead of tk2ico.set

.tk2ico.require <- function ()
{
	if (.Platform$OS.type != "windows")
		stop("This is a Windows-specific function!")
    ## Make sure tcl/tk dde is operational
	if (!capabilities("tcltk"))
		stop("This version of R cannot use Tcl/Tk!")
	if (!is.tk()) stop("Tk is required")
	res <- tclRequire("ico", warn = TRUE)
	if (inherits(res, "tclObj")) res <- tclvalue(res)
	if (res[1] == FALSE)
		stop("Unable to find the 'ico' Tcl/tk package!")
	return(res)  # The package version number
}

tk2ico.create <- function (iconfile, res = 0, size = 16)
{
	if (length(iconfile) != 1) stop("'iconfile' must be of length one!")
	if (!file.exists(iconfile <- as.character(iconfile)))
		stop(gettextf("File '%s' not found!", iconfile))

	## Just use tk2ico.lod() with different default args)
	tk2ico.load(file = iconfile, res = res, size =size)

	## This is old code for winico.dll
	#cmd <- paste("winico createfrom {", iconfile, "}", sep = "")
	#res <- try(icon <- .Tcl(cmd), silent = TRUE)
	#if (inherits(res, "try-error"))  # Tcl error is unreadable, put another one!
	#	stop("Error creating the icon resource; probably wrong 'iconfile'")

	#if (inherits(icon, "tclObj")) class(icon) <- c(class(icon), "tclIcon")
	#return(icon)
}

tk2ico.destroy <- function (icon)
{
	if (!is.tk()) return(NULL)
	if (.Platform$OS.type != "windows") return(NULL)

	if (!inherits(icon, "tclIcon"))
		stop("'icon' is not a \"tclIcon\" object!")

	res <- tclvalue(.Tcl(paste("catch {image delete ", icon, "}", sep = "")))
	return(res == "0")  # Return "0" if OK, "1" otherwise
}

tk2ico.hicon <- function (icon)
{
	## This is not supported any more, was a winico.dll-specific function!
	stop("This function is not supported any more in tcltk2 >= 1.2-0")
	
	#if (!is.tk()) return(NULL)
	#if (.Platform$OS.type != "windows") return(NULL)

	#if (!inherits(icon, "tclIcon"))
	#	stop("'icon' is not a \"tclIcon\" object!")

	#res <- try(hicon <- tcl("winico", "hicon", icon), silent = TRUE)
	#if (inherits(res, "try-error"))  # Tcl error is unreadable, put another one!
	#	stop("Error getting the icon handle for a \"tclIcon\" object!")

	#return(hicon)
}

tk2ico.info <- function (icon, convert = TRUE)
{
    ## This is not supported any more, was a winico.dll-specific function!
	stop("This function is not supported any more in tcltk2 >= 1.2-0")
	
	#if (!is.tk()) return(NULL)
	#if (.Platform$OS.type != "windows") return(NULL)

	#if (!inherits(icon, "tclIcon"))
	#	stop("'icon' is not a \"tclIcon\" object!")

	#res <- try(info <- as.character(tcl("winico", "info", icon)), silent = TRUE)
    #if (inherits(res, "try-error"))  # Tcl error message is unreadable!
    #	stop("Impossible to retrieve icon resource information!")

	#if (isTRUE(convert)) {  # Rework and transform into a data frame
    #    info <- strsplit(info, "-")
    #    info <- matrix(unlist(info), ncol = 8, byrow = TRUE)[, -1]
    #    info <- unlist(strsplit(info, " "))
    #    info <- info[info != ""]
    #    info <- info[(1:(length(info) / 2)) * 2]
    #    info <- matrix(info, ncol = 7)
    #    info <- data.frame(pos = as.integer(info[, 1]),
    #    	width = as.integer(info[, 2]), height = as.integer(info[, 3]),
    #    	geometry = I(info[, 4]), bpp = as.integer(info[, 5]),
    #    	hicon = I(info[, 6]), ptr = I(info[, 7]))
    #}
    #return(info)
}

tk2ico.list <- function (file = "shell32.dll")
{
    ## Make sure that the 'ico' package is loaded 
    .tk2ico.require()
	
	if (length(file) != 1)
		stop("'file' must be of length one!")

	## If the file is not found directly, try using Sys.which()
	if (!file.exists(file)) {
		File <- Sys.which(file)
		if (!file.exists(File))
			stop("file '", file, "' not found")
	} else File <- file

	cmd <- paste("::ico::icons {", File, "}", sep = "")
	res <- try(iconlist <- .Tcl(cmd), silent = TRUE)
	if (inherits(res, "try-error")) # Tcl error message is unreadable!
		stop("Unable to list the icon resources in 'file'!")
	return(as.character(iconlist))
}

tk2ico.sizes <- function (file = "shell32.dll", res = "application")
{
    ## Make sure that the 'ico' package is loaded 
    .tk2ico.require()
	
	if (length(file) != 1)
		stop("'file' must be of length one!")
	if (length(res) != 1)
		stop("'res' must be of length one!")
	## For compatibility reasons, res can be "application", "asterisk", "error",
	## "exclamation", "hand", "question", "information", "warning" or "winlogo"
	## but need to be changed into corresponding ID
	res <- as.character(res)[1]
	res <- switch(res,
		application = "154",
		asterisk = "173",
		error = "28",
		exclamation = "154",
		hand = "29",
		question = "263",
		information = "1001",
		warning = "200",
		winlogo = "47",
		res)
	
	## If the file is not found directly, try using Sys.which()
	if (!file.exists(file)) {
		File <- Sys.which(file)
		if (!file.exists(File))
			stop("file '", file, "' not found")
	} else File <- file

	cmd <- paste("::ico::iconMembers {", File, "} ", res,  sep = "")
	res <- try(iconsizes <- .Tcl(cmd), silent = TRUE)
	if (inherits(res, "try-error")) # Tcl error message is unreadable!
		stop("Unable to list sizes for the icon resource 'res' in 'file'!")
	iconsizes <- as.character(iconsizes)
	iconsizes <- unique(sub("^[^ ]+ ([0-9]+) .+$", "\\1", iconsizes))
	return(as.integer(iconsizes))
}

## For future use, when bpp will be added as argument
#tk2ico.bpp <- function (file = "shell32.dll", res = "application")
#{
#    ## Make sure that the 'ico' package is loaded 
#    .tk2ico.require()
#	
#	if (length(file) != 1)
#		stop("'file' must be of length one!")
#	if (length(res) != 1)
#		stop("'res' must be of length one!")
#	## For compatibility reasons, res can be "application", "asterisk", "error",
#	## "exclamation", "hand", "question", "information", "warning" or "winlogo"
#	## but need to be changed into corresponding ID
#	res <- as.character(res)[1]
#	res <- switch(res,
#		application = "154",
#		asterisk = "173",
#		error = "28",
#		exclamation = "154",
#		hand = "29",
#		question = "263",
#		information = "1001",
#		warning = "200",
#		winlogo = "47",
#		res)
#	
#	## If the file is not found directly, try using Sys.which()
#	if (!file.exists(file)) {
#		File <- Sys.which(file)
#		if (!file.exists(File))
#			stop("file '", file, "' not found")
#	} else File <- file
#
#	cmd <- paste("::ico::iconMembers {", File, "} ", res,  sep = "")
#	res <- try(iconbpp <- .Tcl(cmd), silent = TRUE)
#	if (inherits(res, "try-error")) # Tcl error message is unreadable!
#		stop("Unable to list sizes for the icon resource 'res' in 'file'!")
#	iconbpp <- as.character(iconbpp)
#	iconbpp <- unique(sub("^.+ ([0-9]+)$", "\\1", iconbpp))
#	return(as.integer(iconbpp))
#}

tk2ico.load <- function (file = "shell32.dll", res = "application", size = 16)
{
    ## Make sure that the 'ico' package is loaded 
    .tk2ico.require()
	
	if (length(file) != 1)
		stop("'file' must be of length one!")
	if (length(res) != 1)
		stop("'res' must be of length one!")
	## For compatibility reasons, res can be "application", "asterisk", "error",
	## "exclamation", "hand", "question", "information", "warning" or "winlogo"
	## but need to be changed into corresponding ID
	res <- as.character(res)[1]
	res <- switch(res,
		application = "154",
		asterisk = "173",
		error = "28",
		exclamation = "154",
		hand = "29",
		question = "263",
		information = "1001",
		warning = "200",
		winlogo = "47",
		res)

	## If the file is not found directly, try using Sys.which()
	if (!file.exists(file)) {
		File <- Sys.which(file)
		if (!file.exists(File))
			stop("file '", file, "' not found")
	} else File <- file

	## The old winico code!
	#cmd <- paste("winico load ", res, " {", file, "}", sep = "")
	cmd <- paste("::ico::getIcon {", File, "} ", res, " -res ", size, sep = "")
	res <- try(icon <- .Tcl(cmd), silent = TRUE)
	if (inherits(res, "try-error")) # Tcl error message is unreadable!
		stop("Unable to load the icon resource, 'file' or 'res' is wrong!")

	if (inherits(icon, "tclObj")) class(icon) <- c(class(icon), "tclIcon")
	return(icon)
}

`tk2ico.pos<-` <- function (icon, value)
{
    ## This is not supported any more, was a winico.dll-specific function!
	stop("This function is not supported any more in tcltk2 >= 1.2-0")
	
	#if (!is.tk()) return(NULL)
	#if (.Platform$OS.type != "windows") return(NULL)

	#if (!inherits(icon, "tclIcon"))
	#	stop("'icon' is not a \"tclIcon\" object!")
	#if (length(value) != 1 || !is.numeric(value))
	#	stop("'value' must be numeric and of length one!")
	#value <- round(value)

	#res <- tclvalue(.Tcl(paste("catch {winico pos ", icon, " ", value, "}",
	#	sep = "")))
	#if (res != "0") stop("Error while changing default position of the icon!")
	#return(icon)
}

## iconfile can be either an .ico file, or an .exe
tk2ico.setFromFile <- function (win, iconfile)
{
	## This is the simplest way to set a tk window icon
	tkwm.iconbitmap(win, iconfile)
}

tk2ico.set <- function (win, icon, pos = NULL, type = c("all", "small", "big"))
{
    ## Note that pos and typer are not supported any more (were winico.dll features)!
	if (!missing(pos)) warning("'pos' argument deprecated from tcltk2 1.2-2")
	if (!missing(type)) warning("'type' argument deprecated from tcltk2 1.2-2")
	## Integer for win is not supported any more
	if (inherits(win, "integer"))
		stop("Integers for argument win are not supported any more in tcltk > 1.2-0")
			
	## Make sure that the 'ico' package is loaded 
    .tk2ico.require()
	
	if (!inherits(win, "tkwin") || length(win) < 1)
		stop("'win' is not a \"tkwin\" object")
	if (!inherits(icon, "tclIcon"))
		stop("'icon' is not a \"tclIcon\" object!")
	
	## Change the icon of a Tk window
	tcl("wm", "iconphoto", win, icon)
	
	## This is old code that used winico.dll
	#if (!is.tk()) return(NULL)
	#if (.Platform$OS.type != "windows") return(NULL)

	## win is either a tkwin object, or the HWND of a foreign window
	#if (!inherits(win, c("tkwin", "integer")) || length(win) < 1)
	#	stop("'win' is not a \"tkwin\" object, or an integer (Window handle)!")
	#if (is.integer(win)) win <- win[1]
	#if (!inherits(icon, "tclIcon"))
	#	stop("'icon' is not a \"tclIcon\" object!")
	#if (!is.null(pos) && (length(pos) != 1 || !is.numeric(pos)))
	#	stop("'pos' must be numeric and of length one, or NULL!")
	#type <- match.arg(type)

	#if (type == "all") {  # We search for highest quality icons
    #	## Determine which icons are better in the resource
	#	info <- tk2ico.info(icon, convert = TRUE)
    #	if (nrow(info) == 1) {  # Only one icon in the ressource
    #    	resSmall <- tclvalue(.Tcl(paste("catch {winico setwindow ", win, " ",
	#			icon, " small}", sep = "")))
	#		resBig <- tclvalue(.Tcl(paste("catch {winico setwindow ", win, " ",
	#			icon, " big}", sep = "")))
    #	} else {  # There are several icons in the resource
    #    	## Are there 16x16 icons?
    #    	pos16 <- info[info$geometry == "16x16", ]$pos
    #    	if (length(pos16) > 0) pos16 <- max(pos16) else pos16 <- 0
    #    	resSmall <- tclvalue(.Tcl(paste("catch {winico setwindow ", win,
	#			" ", icon, " small ", pos16, "}", sep = "")))
    #    	## Are there 32x32 icons?
    #    	pos32 <- info[info$geometry == "32x32", ]$pos
    #    	if (length(pos32) > 0) pos32 <- max(pos32) else
    #        	pos32 <- max(info$pos)
    #        resBig <- tclvalue(.Tcl(paste("catch {winico setwindow ", win, " ",
	#			icon, " big ", pos32, "}", sep = "")))
    #    }
    #    ## Compute res, a vector of two logical values, with "0" == success
    #    res <- c((resSmall == "0"), (resBig == "0"))
    #    names(res) <- c("small", "big")

    #} else {  # Other type than 'all'
    #    if (is.null(pos)) pos <- ""
    #    res <- tclvalue(.Tcl(paste("catch {winico setwindow ", win, " ", icon,
	#		" ", type, " ", pos, "}", sep = "")))
	#	res <- (res == "0")  # Because "0" if OK, and "1" otherwise
    #}
    #return(res)
}

tk2ico.taskbar.add <- function(icon, pos = 0, text = tk2ico.text(icon),
leftmenu = NULL, rightmenu = NULL)
{
	## This is not supported any more, was a winico.dll-specific function!
	stop("This function is not supported any more in tcltk2 >= 1.2-0")
	
	#if (!is.tk()) return(NULL)
	#if (.Platform$OS.type != "windows") return(NULL)

	#if (!inherits(icon, "tclIcon"))
	#	stop("'icon' is not a \"tclIcon\" object!")

	#if (length(pos) != 1 || !is.numeric(pos))
	#	stop("'pos' must be numeric and of length one, or NULL!")

	#if (is.null(text)) text <- "" else
	#	text <- paste (" -text {", paste(text, collapse = "\n"), "}", sep = "")

	#if (!inherits(leftmenu, c("tkwin", "NULL")))
	#	stop("'leftmenu' must be a \"tkwin\" object or NULL!")
	#if (!inherits(rightmenu, c("tkwin", "NULL")))
	#	stop("'rightmenu' must be a \"tkwin\" object or NULL!")

	#if (!is.null(leftmenu) || !is.null(rightmenu)) {
	#	hicon <- tclvalue(tk2ico.hicon(icon))

	#	if (is.null(leftmenu)) {
	#		leftcmd <- ""; leftset <- ""
	#	} else {
	#		## The command to trigger the left-click menu
	#		leftcmd <- paste('if { $msg == "WM_LBUTTONUP" } { $::leftmenu',
	#			hicon, ' post $x $y }\n', sep = "")
	#		## ... and the appropriate name for our left-click menu
	#		leftset <- paste('set leftmenu', hicon, ' ', leftmenu$ID, sep = "")
	#	}

	#	if (is.null(rightmenu)) {
	#		rightcmd <- ""; rightset <- ""
	#	} else {
	#		## The command to trigger the right-click menu
	#		rightcmd <- paste('if { $msg == "WM_RBUTTONUP" } { $::rightmenu',
	#			hicon, ' post $x $y }\n', sep = "")
	#		## ... and the appropriate name for our right-click menu
	#		rightset <- paste('set rightmenu', hicon, ' ', rightmenu$ID, sep = "")
	#	}

	#	## Create the proc that will handle mouse clicks on our taskbar icon
	#	cmd <- paste('catch { ', leftset, '\n', rightset, '\n',
	#	'proc taskbarcallback', hicon, ' { hicon msg ico x y } {\n',
	#		leftcmd, rightcmd, '}}', sep = "")
	#	if (tclvalue(.Tcl(cmd)) != "0")
	#		stop("Error while creating the callback for this icon!")

	#	## Finally define the callback call
	#	callback <- paste(' -callback "taskbarcallback', hicon, ' ', hicon,
	#		' %m %i %X %Y"', sep = "")
	#} else callback <- ""

	## Install the taskbar icon
	#res <- tclvalue(.Tcl(paste("catch {winico taskbar add ", icon, " -pos ",
	#	round(pos), text, callback, "}", sep = "")))

	#return(res == "0")  # "0" if success, "1" if error
}

tk2ico.taskbar.delete <- function (icon)
{
	## This is not supported any more, was a winico.dll-specific function!
	stop("This function is not supported any more in tcltk2 >= 1.2-0")
	
	#if (!is.tk()) return(NULL)
	#if (.Platform$OS.type != "windows") return(NULL)

	#if (!inherits(icon, "tclIcon"))
	#	stop("'icon' is not a \"tclIcon\" object!")

	## First delete the Tcl procs that handle our taskbar icon events
	#hicon <- tclvalue(tk2ico.hicon(icon))
	#.Tcl(paste('catch {',
	#	'unset -nocomplain taskbarcallback', hicon, '\n',
	#	'unset -nocomplain leftmenu', hicon, '\n',
	#	'unset -nocomplain rightmenu', hicon, '\n',
	# 	'}', sep = ""))

	## Then delete the taskbar icon
	#res <- tclvalue(.Tcl(paste("catch {winico taskbar delete ", icon, "}",
	#	sep = "")))
	#return(res == "0")  # return "0" if OK, "1" otherwise
}

tk2ico.taskbar.modify <- function (icon, pos = NULL, text = NULL)
{
	## This is not supported any more, was a winico.dll-specific function!
	stop("This function is not supported any more in tcltk2 >= 1.2-0")
	
	#if (!is.tk()) return(NULL)
	#if (.Platform$OS.type != "windows") return(NULL)

	#if (!inherits(icon, "tclIcon"))
	#	stop("'icon' is not a \"tclIcon\" object!")

	#if (!is.null(pos)) {
	#	if (length(pos) != 1 || !is.numeric(pos))
	#		stop("'pos' must be numeric and of length one, or NULL!")
	#	pos <- paste(" -pos", round(pos))
	#}
	#if (is.null(text)) text <- "" else
	#	text <- paste (" -text {", paste(text, collapse = "\n"), "}", sep = "")

	#cmd <- paste("catch {winico taskbar modify ", icon, pos, text, "}", sep = "")
    #return(tclvalue(.Tcl(cmd)) == "0")   # "0" if OK, "1" otherwise
}

tk2ico.text <- function (icon)
{
	## This is not supported any more, was a winico.dll-specific function!
	stop("This function is not supported any more in tcltk2 >= 1.2-0")
	
	#if (!is.tk()) return(NULL)
	#if (.Platform$OS.type != "windows") return(NULL)

	#if (!inherits(icon, "tclIcon"))
	#	stop("'icon' is not a \"tclIcon\" object!")

	#res <- try(text <- tclvalue(tcl("winico", "text", icon)), silent = TRUE)
	#if (inherits(res, "try-error"))  # Tcl error is unreadable, put another one!
	#	stop("Error getting the text associated with an icon!")

    #return(res)
}

`tk2ico.text<-` <- function (icon, value)
{
	## This is not supported any more, was a winico.dll-specific function!
	stop("This function is not supported any more in tcltk2 >= 1.2-0")
	
	#if (!is.tk()) return(NULL)
	#if (.Platform$OS.type != "windows") return(icon)

	#if (!inherits(icon, "tclIcon"))
	#	stop("'icon' is not a \"tclIcon\" object!")
	#if (length(value) < 1)
	#	stop("'value' must not be empty or NULL!")

	#value <- paste(value, collapse = "\n")
    #res <- tclvalue(.Tcl(paste("catch {winico text ", icon, " {", value, "}}",
	#	sep = "")))
	#if (res != "0") stop("Error while changing text of the icon!")
    #return(icon)
}
