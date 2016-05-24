### tcltk2-Internal.R - Hidden functions for tcltk2
### Copyright (c), Philippe Grosjean (phgrosjean@sciviews.org)

.onLoad <- function(libname, pkgname) {
	libdir <- file.path(libname, pkgname, "tklibs")

	## A slightly modified version of addTclPath() that works also within SciViews
	addTclPath <- function (path = ".") {
		if (.Platform$OS.type == "windows") 
		    path <- gsub("\\\\", "/", path)
		a <- tclvalue(tcl("set", "::auto_path"))
		paths <- strsplit(a, " ", fixed = TRUE)[[1L]]
		if (!path %in% paths) 
		    tcl("lappend", "::auto_path", path)
	}
    res <- addTclPath(libdir)	# extend the Tcl/Tk path

    ## Make sure that Tcl/Tk locale is the same one as current R locale
	lang <- getLanguage()
	if (lang != "") {  # Set the same language for Tcl/Tk
		res <- tclRequire("msgcat")
	    if (inherits(res, "tclObj")) tcl("::msgcat::mclocale", lang)
	}

    if (is.tk()) {
		## Here is how we could install the supplementary material in Tcl/Tk
	
		## This is for a better management of scrollbars in listbox, text, canvas
		tclRequire("autoscroll") # Version 1.1
		tcl("source", file.path(libdir, "scrolledWidget.tcl"))
		
		##tclRequire("choosefont")      # Version 0.2
		##tclRequire("ctext")			# Version 3.1
		##tclRequire("cursor")       	# Version 0.1
		##tclRequire("mclistbox")    	# Version 1.2
		##tclRequire("swaplist")    	# Version 0.2
		##tclRequire("tablelist")    	# Version 5.5
		##Not provided any more -> tclRequire("Tktable")   		# Version 2.9

		## The following code is not implemented as Tcl package... just source it
		tcl("source", file.path(libdir, "notebook1.3", "notebook.tcl"))
	    tcl("source", file.path(libdir, "tree1.7", "tree.tcl"))

		## Do we try to load the tile widgets? (only if Tcl./Tk < 8.5)
		if (as.numeric(.Tcl("set ::tcl_version")) < 8.5) {
###				tcl("source", file.path(libdir, "fonts.tcl"))
				## Define fonts used in Tk (note: must be done AFTER loading tile!)
				## Default values for system fonts are calculated by tile...
				## but they should be computed from the system, actually
				## We collect back those values calculated by tile and possibly override
				## them with better values
###				tk2font.setstyle(system = TRUE, default.styles = TRUE, text = TRUE)
				### TODO: reflect possible changes to other graphical toolkits (how?)
		} else {	# There is a bug in mclistbox with Tcl/Tk 8.5
			## Patch by Christiane Raemsch, slightly modified by Ph. Grosjean
			## This is essentially the listbox procedure, but with an additional
			## focus argument required by mclistbox
			.Tcl('proc ::tk::ListboxBeginSelect {w el {focus 0}} {
				variable ::tk::Priv
				if {[$w cget -selectmode] eq "multiple"} {
					if {[$w selection includes $el]} {
						$w selection clear $el
					} else {
						$w selection set $el
					}
				} else {
					$w selection clear 0 end
					$w selection set $el
					$w selection anchor $el
					set Priv(listboxSelection) {}
					set Priv(listboxPrev) $el
				}
				event generate $w <<ListboxSelect>>
				if {$focus && [winfo exists $w]} {
					focus $w
				}
			}')
		}
		
		## Load additional ttk themes - No: load only on demand!
		## Not done any more on startup, done on demand in tk2theme() now
		#try(tclRequire("ttk::theme::plastik"), silent = TRUE)
		#try(tclRequire("ttk::theme::keramik"), silent = TRUE)
		#try(tclRequire("ttk::theme::keramik_alt"), silent = TRUE)
		#try(tclRequire("ttk::theme::clearlooks"), silent = TRUE)
		#try(tclRequire("ttk::theme::radiance"), silent = TRUE)
		
		## Which ttk theme should we use?
		## If the user specified a default theme, use it
		if (!.loadTheme()) {
			## ...otherwise, try to guess the best default value
			themes <- try(tk2theme.list(), silent = TRUE)
			if (!inherits(themes, "try-error")) {
				if ("aqua" %in% themes) { # This must be aquaTk on a Mac
					try(tk2theme("aqua"), silent = TRUE)
				} else if ("vista" %in% themes) { # This must be Vista or Win 7
					try(tk2theme("vista"), silent = TRUE)
				} else if ("xpnative" %in% themes) { # This must be XP
					try(tk2theme("xpnative"), silent = TRUE)
				} else if ("winnative" %in% themes) { # This must be a pre-XP windows
					try(tk2theme("winnative"), silent = TRUE)
				} else if (.isUbuntu()) {
					try(tk2theme("radiance"), silent = TRUE)
					## Special treatment for Ubuntu: change fonts to Ubuntu and Ubuntu mono
					## and use white text on black for tooltips
					tkfont.configure("TkDefaultFont", family = "Ubuntu", size = 11)
					tkfont.configure("TkMenuFont", family = "Ubuntu", size = 11)
					tkfont.configure("TkCaptionFont", family = "Ubuntu", size = 10)
					tkfont.configure("TkSmallCaptionFont", family = "Ubuntu", size = 9)
					tkfont.configure("TkTooltipFont", family = "Ubuntu", size = 9)
					tkfont.configure("TkMenuFont", family = "Ubuntu", size = 11)
					tkfont.configure("TkHeadingFont", family = "Ubuntu", size = 12)
					tkfont.configure("TkIconFont", family = "Ubuntu", size = 11)
					tkfont.configure("TkTextFont", family = "Ubuntu", size = 11)
					tkfont.configure("TkFixedFont", family = "Ubuntu Mono", size = 11)
					res <- tclRequire("tooltip")
					if (inherits(res, "tclObj")) {
						.Tcl(paste("set ::tooltip::labelOpts [list -highlightthickness 0",
							"-relief solid -bd 1 -background black -fg white]"))
					}
				} else { # A modern "default" theme that fit not too bad in many situations
					try(tk2theme("clearlooks"), silent = TRUE)
				}
			}
		}
		## Save default font as TkSysDefaultFont
		tk2font.set("TkSysDefaultFont", tk2font.get("TkDefaultFont"))
	}
	
	## Windows only
    if (.Platform$OS.type == "windows") {
		try(tclRequire("dde"), silent = TRUE)       # Version 1.2.2
        ## Not loaded automatically!
        #tclRequire("registry")  # Version 1.1.3
		## Support for winico.dll is drop from version 1.2-1!
    #    if (nzchar(r_arch <- .Platform$r_arch))
	#		tcl("load", file.path(libname, pkgname, "libs", r_arch, "Winico06.dll"))
	#	else
	#		tcl("load", file.path(libname, pkgname, "libs", "Winico06.dll"))
		## Also register the DDE server as TclEval|R
        try(tk2dde("R"), silent = TRUE)
    }
}

.onUnload <- function (libpath)
{
    # PhG: was .Last.lib()
	## Remove all currently scheduled tasks
	tclTaskDelete(id = NULL)
}

.saveTheme <- function ()
	cat(tk2theme(), "\n", sep = "", file = "~/.Rtk2theme")
	
.loadTheme <- function () {
	if (file.exists("~/.Rtk2theme")) {
		theme <- try(readLines("~/.Rtk2theme")[1], silent = TRUE)
		if (inherits(theme, 'try-error')) return(FALSE)
		## Try changing the tk2theme according to this value
		res <- try(tk2theme(theme), silent = TRUE)
		return(!inherits(res, "try-error"))
	} else return(FALSE)
}

.isUbuntu <- function () {
	## Note: take care not to call 'cat' on Windows: it is usually *not* there!
	if (.Platform$OS.type == "windows" || grepl("^mac", .Platform$pkgType))
		return(FALSE)	# This is either Windows or Mac OS X!
	grepl("^Ubuntu", suppressWarnings(try(system("cat /etc/issue",
		intern = TRUE, ignore.stderr = TRUE), silent = TRUE))[1])	
}

.mergeList <- function (l1, l2)
{
	## For named lists, overwrite items of l1 present in l2
	nms <- names(l2)
	## Deal with named items
	if (length(nms)) {
		named <- nms != ""
		if (any(named)) {
			l2n <- l2[named]
			nmsn <- nms[named]
			for (i in 1:length(nmsn)) l1[[nmsn[i]]] <- l2n[[nmsn[i]]]
		}
		## Keep only non named items in l2
		l2 <- l2[!named]
	}
	## Deal with non named items in l2
	if (length(l2)) { # Unnamed list
		n1 <- length(l1)
		n2 <- length(l2)
		for (i in 1:n2) l1[[n1 + i]] <- l2[[i]]
	}
	return(l1)
}

.configStd <- function (x, lstval)
{
	## These config parameters are considered as data
	## Image
	if (!is.null(lstval$image)) {
		tkconfigure(x, image = lstval$image)
		lstval$image <- NULL
	}
	## Text
	if (!is.null(lstval$text)) {
		tkconfigure(x, text = lstval$text)
		lstval$text <- NULL
	}
	## Textvariable
	if (!is.null(lstval$textvariable)) {
		tkconfigure(x, textvariable = lstval$textvariable)
		lstval$textvariable <- NULL
	}
	## Values
	if (!is.null(lstval$values)) {
		values(x) <- lstval$values
		lstval$values <- NULL
	}
	## Value
	if (!is.null(lstval$value)) {
		value(x) <- lstval$value
		lstval$value <- NULL
	}
	## Selection
	if (!is.null(lstval$selection)) {
		selection(x) <- lstval$selection
		lstval$selection <- NULL
	}
	## Label (not a Tk attribute)
	if (!is.null(lstval$label)) {
		label(x) <- lstval$label
		lstval$label <- NULL
	}
	## Tag (not a Tk attribute)
	if (!is.null(lstval$name)) {
		tag(x) <- lstval$tag
		lstval$tag <- NULL
	}
	## Tooltip
	if (!is.null(lstval$tip)) {
		tip(x) <- lstval$tip
		lstval$tip <- NULL
	}
	## Disabled (is tk 'state' parameter indeed)
	if (!is.null(lstval$disabled)) {
		disabled(x) <- lstval$disabled
		lstval$disabled <- NULL
	}
	## Return modified value list
	return(lstval)
}

.wraplength <- function (w, width)
{
	## Calculate wraplength required for tk2label widgets
	## width is expressed in characters, but wraplength must be given in pixels
	## This is stupid and requires additional computation to calculate the
	## width in pixel of an average character, like "0" to do the conversion!
	## Get the average size of one character in the current font used
	
	## If width is not set, just return a large value for wraplength
	if (!length(width)) return(1000)
	
	## Get the font and measure it	
	font <- tclvalue(tkcget(w, "-font"))
	if (font == "") font <- tk2style("tk2label", "font")
	if (font == "") {
		charsize <- 8 # Use an everage value
	} else charsize <- as.numeric(tkfont.measure(tkfont.actual(font), "0"))
	
	## Optimal wraplength is width * charsize
	return(width * charsize)
}

.TempEnv <- function ()
{
    pos <-  match("SciViews:TempEnv", search())
    if (is.na(pos)) {  # Must create it
        `SciViews:TempEnv` <- list()
        Attach <- function (...) get("attach", mode = "function")(...)
        Attach(`SciViews:TempEnv`, pos = length(search()) - 1)
        rm(`SciViews:TempEnv`)
        pos <- match("SciViews:TempEnv", search())
    }
    return(pos.to.env(pos))
}

.assignTemp <- function (x, value, replace.existing = TRUE)
    if (replace.existing || !exists(x, envir = .TempEnv(), mode = "any",
		inherits = FALSE))
        assign(x, value, envir = .TempEnv())

.getTemp <- function (x, default = NULL, mode = "any", item = NULL)
{
    if (is.null(item)) Mode <- mode else Mode <- "any"
    if  (exists(x, envir = .TempEnv(), mode = Mode, inherits = FALSE)) {
        dat <- get(x, envir = .TempEnv(), mode = Mode, inherits = FALSE)
        if (is.null(item)) return(dat) else {
            item <- as.character(item)[1]
            if (inherits(dat, "list") && item %in% names(dat)) {
                dat <- dat[[item]]
                if (mode != "any" && mode(dat) != mode) dat <- default
                return(dat)
            } else {
                return(default)
            }
        }
    } else {  # Variable not found, return the default value
        return(default)
    }
}
