### tk2commands.R - Additional tk commands to manipulate tk2 widgets
### Copyright (c), Philippe Grosjean (phgrosjean@sciviews.org)
### Licensed under LGPL 3 or above
###
### Changes:
### - 2007-01-01: fisrt version (for tcltk2_1.0-0)
###
### To do:
### - Rework all this...
### - Style option of Tile widgets?
### - Implement style element options ...

tk2column <- function (widget, action = c("add", "configure", "delete", "names",
"cget", "nearest"), ...)
{
    Action <- action[1]
    tcl(widget, "column", Action, ...)
}

tk2list.set <- function (widget, items)
{
	## Set a list of values for a widget (e.g., combobox)
	if (inherits(widget, "ttk2combobox")) {
        ## ttk::combobox uses -values parameter
        tkconfigure(widget, values = as.character(items))
    } else {
        ## Try to use the defaul method
		## First, clear the list
		tcl(widget, "list", "delete", 0, "end")
		## Then, insert all its elements
		items <- as.character(items)
		for (item in items) tcl(widget, "list", "insert", "end", item)
    }
}

tk2list.insert <- function (widget, index = "end", ...)
{
	## Insert one or more items in a list
	if (inherits(widget, "ttk2combobox")) {
        ## ttk::combobox uses -values parameter
		Items <- as.character(unlist(list(...)))
		if (length(Items) < 1) return()	# Nothing to insert
		List <- as.character(tcl(widget, "cget", "-values"))
		if (length(List) < 2 && List == "") {
			## The list in empty, simply add these items
			List <- Items
		} else if (index == "end" || index > length(List) - 1) {
			List <- c(List, Items)
		} else if (index == 0){
			## Insert items at the beginning of the list
			List <- c(Items, List)
		} else {
			## Insert items inside the list
			List <- c(List[1:index], Items, List[(index + 1):length(List)])
		}
		## Reassign this modified list to the combobox
		tkconfigure(widget, values = List)
	} else {
		tcl(widget, "list", "insert", index, ...)
	}
}

tk2list.delete <- function (widget, first, last = first)
{
	## Delete one or more items from a list
	if (inherits(widget, "ttk2combobox")) {
        ## ttk::combobox uses -values parameter
		List <- as.character(tcl(widget, "cget", "-values"))
		if (length(List) < 2 && List == "") return(List)	# The list in empty
		if (last == "end") last <- length(List) else last <- last + 1
		List <- List[-((first + 1):last)]
		## Reassign this modified list to the combobox
		tkconfigure(widget, values = List)
	} else {
		tcl(widget, "list", "delete", first, last)
	}
}

tk2list.get <- function (widget, first = 0, last = "end")
{
	## Get the list of elements in a widget (e.g., combobox)
	if (inherits(widget, "tk2combobox")) {
        ## ttk::combobox uses -values parameter
		List <- as.character(tcl(widget, "cget", "-values"))
		if (length(List) < 2 && List == "") return(List)
		if (last == "end") last <- length(List) else last <- last + 1
		return(List[(first +1):last])
	} else {
		as.character(tcl(widget, "list", "get", first, last))
	}
}

tk2list.size <- function (widget)
{
	## Get the length of the list of elements in a widget (e.g., combobox)
	if (inherits(widget, "tk2combobox")) {
        ## ttk::combobox uses -values parameter
		List <- as.character(tcl(widget, "cget", "-values"))
		return(length(List))
	} else {
		as.numeric(tcl(widget, "list", "size"))
	}
}

tk2state.set <- function (widget, state = c("normal", "disabled", "readonly"))
{
	## Change the state of a widget
	state <- as.character(state[1])
	tkconfigure(widget, state = state)
}

tk2insert.multi <- function (widget, where = "end", items)
{
    ## We insert one or several lines in a multicolumn widget
    items <- as.matrix(items)
    ## A vector is coerced into a column matrix and we want a row matrix here
    if (ncol(items) == 1) items <- t(items)
    ## Convert the matrix into [list {el1} {el2} {el3}] [list {el4}, {el5}, {el6}], ...
    makeTclList <- function(x) paste("[list {", paste(x, collapse = "} {"), "}]", sep = "")
    TclList <- paste(apply(items, 1, makeTclList), collapse = "\\\n")
    .Tcl(paste(widget, "insert", where, TclList))
}

tk2notetraverse <- function (nb)
{
	res <- tcl("ttk::notebook::enableTraversal", nb)
	return(invisible(res))
}

tk2notetab <- function (nb, tab)
{
    if (inherits(nb, "tk2notebook")) {
        ## We need the tab index, so, look for it
		ntab <- as.numeric(tcl(nb, "index", "end"))
		if (ntab < 1) return(NULL)
		tabidx <- -1
		for (i in 0:(ntab - 1))
			if (tclvalue(tcl(nb, "tab", i, "-text")) == tab) {
				tabidx <- i
				break
			}
		if (tabidx > -1) {
			tabid <- paste(nb$ID, tabidx + 1, sep = ".")
			## Create a simili tkwin object referring to this page
            w <- list()
            w$ID <- tabid
			w$env <- new.env()
			w$env$num.subwin <- 0
			w$env$parent <- nb
            class(w) <- c("ttk2notetab", "tk2container", "tkwin")
            return(w)
		} else return(NULL)  # Tab not found!
    } else stop ("'nb' must be a 'tk2notebook' object")
}

tk2notetab.select <- function (nb, tab)
{
    ## Select a tab in a notebook
    if (inherits(nb, "tk2notebook")) {
        ## Tile notebook
		## We need the tab index, so, look for it
		ntab <- as.numeric(tcl(nb, "index", "end"))
		if (ntab < 1) return(invisible(FALSE))
		tabidx <- -1
		for (i in 0:(ntab - 1))
			if (tclvalue(tcl(nb, "tab", i, "-text")) == tab) {
				tabidx <- i
				break
			}
		if (tabidx > -1) {
			tkselect(nb, tabidx)
			return(invisible(TRUE))
		} else return(invisible(FALSE))
    } else stop ("'nb' must be a 'tk2notebook' object")
}

tk2notetab.text <- function (nb)
{
    ## Select a tab in a notebook
    if (inherits(nb, "tk2notebook")) {
		return(tclvalue(tcl(nb, "tab", "current", "-text")))
    } else stop ("'nb' must be a 'tk2notebook' object")
}

## Themes management
tk2theme.elements <- function ()
	return(as.character(.Tcl("ttk::style element names")))

tk2theme.list <- function ()
	return(as.character(.Tcl("ttk::style theme names")))

tk2theme <- function (theme = NULL)
{
    if (is.null(theme)) {  # Get it
        res <- getOption("tk2theme")
    } else {  # Set it to theme
        ## First, check if the theme is already loaded... or try loading it
		loadedThemes <- tk2theme.list()
		if (!theme %in% loadedThemes) {
			## Could be plastik, keramik, keramik_alt, clearlooks, radiance 
			res <- try(tclRequire(paste0("ttk::theme::", theme)), silent = TRUE)
			if (inherits(res, "try-error"))
				stop("Ttk theme ", theme, " is not found")
		}
		## Themes (like radiance) change TkDefaultFont => reset it for the others
		if (theme == "radiance") {
			tkfont.configure("TkDefaultFont", family = "Ubuntu", size = 11)
		} else tk2font.set("TkDefaultFont", tk2font.get("TkSysDefaultFont"))
		## Change theme
		.Tcl(paste("ttk::style theme use", theme))
        ## And save current theme in option "tk2theme"
        options(tk2theme = theme)
		## Make sure to homogenize background for old tk widgets (suggested by Milan Bouchet-Valat)
		## Note: foreground not defined for plastik and keramik => workaround
		fg <- tclvalue(.Tcl("ttk::style lookup TLabel -foreground"))
		if (fg == "") fg <- "#000000"
		afg <- tclvalue(.Tcl("ttk::style lookup TLabel -foreground active"))
		if (afg == "") afg <- "#000000"
		ffg <- tclvalue(.Tcl("ttk::style lookup TLabel -foreground focus"))
		if (ffg == "") ffg <- "#000000"
		hfg <- tclvalue(.Tcl("ttk::style lookup TLabel -foreground hover"))
		if (hfg == "") hfg <- "#000000"
		.Tcl(paste("tk_setPalette",
			"background",
				tclvalue(.Tcl("ttk::style lookup TLabel -background")),
			"foreground", fg,
			"activeBackground",
				tclvalue(.Tcl("ttk::style lookup TLabel -background active")),
			"activeForeground", afg,
			"disabledForeground",
				tclvalue(.Tcl("ttk::style lookup TLabel -foreground disabled")),
			"highlightBackground", "white",
				#tclvalue(.Tcl("ttk::style lookup TLabel -background focus")),
			"highlightColor", ffg,
			"insertBackground",
				tclvalue(.Tcl("ttk::style lookup TLabel -foreground active")),
			"selectBackground",
				tclvalue(.Tcl("ttk::style lookup TText -selectbackground")),
			"selectForeground",
				tclvalue(.Tcl("ttk::style lookup TText -selectforeground")),
			"selectColor",
				tclvalue(.Tcl("ttk::style lookup TText -selectforeground")),
			"throughColor", hfg),
			"fieldBackground",
				tclvalue(.Tcl("ttk::style lookup TEntry -fieldbackground")))
		
		## Set menu font the same as label font
		font <- tclvalue(.Tcl("ttk::style lookup TLabel -font"))
		if (!length(font) || font == "") font <- "TkDefaultFont"
		tk2font.set("TkMenuFont", tk2font.get(font))
		
		## Return the theme
		res <- theme
    }
    res
}
### Note: to change a style element: .Tcl('ttk::style configure TButton -font "helvetica 24"')
### Create a derived style: ttk::style configure Emergency.TButton -font "helvetica 24" -foreground red -padding 10
### Changing different states:
###ttk::style map TButton \ 
###	-background [list disabled #d9d9d9  active #ececec] \ 
###	-foreground [list disabled #a3a3a3] \ 
###	-relief [list {pressed !disabled} sunken] \ 
###	;

## Function to look for a ttk style
tk2style <- function (class, style, state = c("default", "active",
"disabled", "focus", "!focus", "pressed", "selected", "background", "readonly",
"alternate", "invalid", "hover", "all"), default = NULL)
{
	## Get a ttk style in the current theme
	## Class is either the TTk class, or the tk2 function name
	## TODO: add tk2toolbutton and tk2sizegrip!
	class <- switch(class,
		tk2button = "TButton",
		tk2label = "TLabel",
		tk2toolbutton = "Toolbutton",
		tk2menubutton = "TMenubutton",
		tk2checkbutton = "TCheckbutton",
		tk2radiobutton = "TRadiobutton",
		tk2entry = "TEntry",
		tk2combobox = "TCombobox",
		tk2notebook = "TNotebook",
		tk2labelframe = "TLabelframe",
		tk2scrollbar = "TScrollbar",
		tk2scale = "TScale",
		tk2progress = "TProgressbar",
		#tk2spinbox = "TSpinbox",
		tk2tree = "Treeview",
		tk2frame = "TFrame",
		tk2panedwindow = "TPanedwindow",
		tk2separator = "TSeparator",
		#"TSizegrip",
		as.character(class)[1] # Supposed to be the ttk class
		## Not ttk widgets: tk2canvas, tk2ctext, tk2edit, tk2listbox,
		## tk2mclistbox, tk2menu, tk2menuentry, tk2spinbox, tk2table
	)
	style = paste("-", as.character(style)[1], sep = "")
	state = match.arg(state)
	if (is.null(default)) default <- ""
	
	## styles creates a named vector (items in even elements, labels = odd)
	styles <- function (x) {
		st <- as.character(x)
		l <- length(st)
		if (l == 0) return(character(0))
		if (l == 1) return(c(default = st))
		if (l %% 2 > 0) stop("Didn't get an even number of items: ", st)
		stnames <- st[seq(1, l - 1, by = 2)]
		st <- st[seq(2, l, by = 2)]
		names(st) <- stnames
		return(st)
	}
	
	## First look at the map for this class
	res <-  styles(tcl("ttk::style", "map", class, style))
	res2 <-  styles(tcl("ttk::style", "map", ".", style))
	res <- c(res, res2[!names(res2) %in% names(res)])
	res2 <-  styles(tcl("ttk::style", "configure", class, style))
	res <- c(res, res2[!names(res2) %in% names(res)])
	res2 <-  styles(tcl("ttk::style", "configure", ".", style))
	res <- c(res, res2[!names(res2) %in% names(res)])
	if (length(res) == 0) res <- c(default = default)
	
	## If state != "all", try to resolve the right state
	if (state != "all") {
		## If the given state is there, use it
		if (state %in% names(res)) {
			return(res[state])
		} else if ("default" %in% names(res)) {
			return(res["default"])
		} else {
			return(c(default = as.character(default)[1]))
		}
	} else return(res)
}

tk2dataList <- function (x)
{
	## List data parameters for a given tk2widget
	## Data manage the content of the widgets
	## Common items are label, tag, and tip
	## image: widgets that can display images
	## text, textvariable: display text
	## values, value and selection
	## command: the command to run
	## validate, validatecommand, invalidcommand: validation mechanism
	## variable: varaible associated with value
	## postcommand: specific to comboboxes, to fill them!
	## onvalue & offvalue: specific to checkbutton
	## default: specific for button (default button in a dialog box)
	## show: specific to entry for password... clash with treeview show => ???
	## mode, maximum, value: for progressbars
	## from, to, increment, : for spinbox & scale + format
	## Look in text widget what we keep!
	if (is.tk2widget(x)) cl <- class(x)[1] else cl <- as.character(x)[1]
	res <- switch(cl,
		tk2button = c("image", "text", "textvariable", "command", "default"),
		tk2canvas = character(0),
		tk2checkbutton = c("image", "text", "textvariable", "variable",
			"command", "onvalue", "offvalue"),
		tk2combobox = c("postcommand", "textvariable", "values"),
		tk2ctext = c("values", "value", "selection", "maxundo", "undo",
			"spacing1", "spacing2", "spacing3", "tabs", "tabstyle"), # language
		tk2entry = c("invalidcommand", "textvariable", "validate",
			"validatecommand", "show"),
		tk2label = c("image", "text", "textvariable"),
		tk2labelframe = c("text"),
		tk2listbox = c("values", "value", "selection"),
		tk2mclistbox = c("values", "value", "selection"),
		tk2notebook = character(0),
		tk2notetab = c("image", "text"),
		tk2panedwindow = character(0),
		tk2progress = c("mode", "maximum", "value", "variable"),
		tk2radiobutton = c("image", "text", "textvariable",
			"command", "value", "variable"),
		tk2scale = c("command", "from", "to", "value", "variable"),
		tk2scrollbar = c("command"),
		tk2separator = character(0),
		#tk2sizegrip = character(0),
		tk2spinbox = c("validate", "validatecommand", "from", "to", "increment",
			"values", "format", "command"),
		tk2table = c("values", "value", "selection"),
		tk2text = c("values", "value", "selection", "maxundo", "undo",
			"spacing1", "spacing2", "spacing3", "tabs", "tabstyle"),
		tk2tree = c("values", "value", "selection"),
		stop("Unknown tk2widget, provide a tk2widget object or its class")	
	)
	## Add label, tag & tip for all
	res <- c(res, "label", "tag", "tip")
	return(res)
}

tk2configList <- function (x)
{
	## List config parameters for a given tk2widget
	## Note: most of the appearance is controlled by the theme, we keep here
	## only a subset of items that are most useful considering themed widgets:
	## height, width or length: the size of the widget
	## compound: how image and text are composed
	## justify and wrap: control of text flow
	## orient: for widgets that can be horizontal or vertical
	## selectmode: for widgets that allow for multiple selections
	## show: tree and/or headings for the treeview widget
	if (is.tk2widget(x)) cl <- class(x)[1] else cl <- as.character(x)[1]
	res <- switch(cl,
		tk2button = c("compound", "width"),
		tk2canvas = c("height", "width"),
		tk2checkbutton = c("compound", "width"),
		tk2combobox = c("justify", "height", "width"),
		tk2ctext = c("height", "width"),
		tk2entry = c("justify", "width"),
		tk2label = c("compound", "justify", "width", "wraplength"), # Use wrap!
		tk2labelframe = c("height", "width"),
		tk2listbox = c("height", "width", "selectmode"),
		tk2mclistbox = c("height", "width", "selectmode"),
		tk2notebook = c("height", "width"),
		tk2notetab = c("compound"),
		tk2panedwindow = c("orient", "height", "width"),
		tk2progress = c("length", "orient"),
		tk2radiobutton = c("compound", "width"),
		tk2scale = c("length", "orient"),
		tk2scrollbar = c("orient"),
		tk2separator = character(0),
		#tk2sizegrip = character(0),
		tk2spinbox = c("wrap"),
		tk2table = c("height", "width"),
		tk2text = c("height", "width"),
		tk2tree = c("height", "selectmode", "show"), # show tree and/or headings
		stop("Unknown tk2widget, provide a tk2widget object or its class")
	)
	## Add cursor and takefocus that are common to all
	# Should we really add these?
	#res <- c(res, "cursor", "takefocus")
	return(res)
}

setLanguage <- function (lang)
{
	## Change locale for both R and Tcl/Tk
	Sys.setenv(language = lang)
	try(Sys.setlocale("LC_MESSAGES", lang), silent = TRUE)  # Fails on Windows!
	res <- tclRequire("msgcat")
	if (inherits(res, "tclObj")) {
		tcl("::msgcat::mclocale", lang)
		return(TRUE)
	} else {
		return(FALSE)
	}
}

getLanguage <- function ()
{
	## Try to recover current language used for messages and GUI stuff in R
	lang <- Sys.getenv("language")
	if (lang == "") lang <- Sys.getlocale("LC_MESSAGES")
	## This is a bad hack that probably does not work all the time, but at least,
	## it works under Windows for getting "fr" for French language
	if (lang == "") lang <- tolower(substr(Sys.getlocale("LC_TIME"), 1, 2))
	return(lang)
}

is.tk <- function ()
	return(tclvalue(.Tcl("catch { package present Tk }")) == "0")

is.ttk <- function ()
{
	res <- is.tk() && as.numeric(tcl("set", "::tk_version")) >= 8.5
	return(res)
}
