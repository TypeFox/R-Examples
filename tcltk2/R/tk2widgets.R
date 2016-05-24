### tk2widgets.R - Support for the ttk widgets
### Copyright (c), Philippe Grosjean (phgrosjean@sciviews.org)
### Licensed under LGPL 3 or above
###
### Changes:
### - 2012-01-09: reworked tk2label to provide "full options"
###
### - 2012-01-07: listbox now behaves like a ttk widget, although it is not!
###
### - 2009-06-30: only use ttk (no tile or plain widgets)
###
### - 2007-01-01: first version (for tcltk2_1.0-0)
###
### To do:
### - Rework all this and add new widgets like sizegrip, tkplot, ...

tk2button <- function (parent, tip = "", ...)
{
	if (!is.ttk()) stop("Tcl/Tk >= 8.5 is required")
	w <- tkwidget(parent, "ttk::button", ...)
	if (tip != "") tk2tip(w, tip)
	class(w) <- c("tk2button", "tk2widget", class(w))
	return(w)
}

tk2canvas <- function (parent, tip = "", ...)
{
	if (!is.ttk()) stop("Tcl/Tk >= 8.5 is required")
### TODO: use autoscroll here!
	## Default background to fieldbackground
	if (any(names(list(...)) == "background")) {
        w <- tkwidget(parent, "canvas", ...)
    } else {
        background <- tclvalue(.Tcl("ttk::style lookup TEntry -fieldbackground"))
		if (background == "") background <- "white"
		w <- tkwidget(parent, "canvas",
            background = background, ...)
    }
	if (tip != "") tk2tip(w, tip)
	class(w) <- c("tk2canvas", "tk2widget", class(w))
	return(w)
}

tk2checkbutton <- function (parent, tip = "", ...)
{
### TODO: associate with a variable and set both states values
	if (!is.ttk()) stop("Tcl/Tk >= 8.5 is required")
	w <-tkwidget(parent, "ttk::checkbutton", ...)
	if (tip != "") tk2tip(w, tip)
	class(w) <- c("tk2checkbutton", "tk2widget", class(w))
	return(w)
}

tk2combobox <- function (parent, tip = "", ...)
{
### TODO: associate the list and results with a variable and intialize the widget
	if (!is.ttk()) stop("Tcl/Tk >= 8.5 is required")
	w <- tkwidget(parent, "ttk::combobox", ...)
	if (tip != "") tk2tip(w, tip)
	class(w) <- c("tk2combobox", "tk2widget", class(w))
	return(w)
}

### TODO: a centry widget
tk2entry <- function (parent, tip = "", ...)
{
### TODO: add cut/copy/paste/select all/clear context menu
	if (!is.ttk()) stop("Tcl/Tk >= 8.5 is required")
	w <- tkwidget(parent, "ttk::entry", cursor = "xterm", ...)
	if (tip != "") tk2tip(w, tip)
	class(w) <- c("tk2entry", "tk2widget", class(w))
	return(w)
}

tk2frame <- function (parent, ...)
{
	if (!is.ttk()) stop("Tcl/Tk >= 8.5 is required")
	w <- tkwidget(parent, "ttk::frame", ...)
	class(w) <- c("tk2frame", "tk2widget", class(w))
	return(w)
}

tk2label <- function (parent, tip, label, tag, cfglist, wrap = FALSE, ...)
{
	## Also image, text, textvariable, label & tag for data
	## width = -10, compound = "left", justify = "left", wrap = FALSE for config
	## Special treatment from wrap (wraplength)
	if (!is.ttk()) stop("Tcl/Tk >= 8.5 is required")
	## If a config list is provided, rework arguments
	if (!missing(cfglist)) {
		args <- .mergeList(cfglist, list(...))
		args$parent <- parent
		if (!missing(tip)) args$tip <- tip
		if (!missing(label)) args$label <- label
		if (!missing(tag)) args$tag <- tag
		return(do.call(tk2label, args))
	}
	## Create the widget and the corresponding R object
	w <- tkwidget(parent, "ttk::label", ...)
	class(w) <- c("tk2label", "tk2widget", class(w))
	if (!missing(tip)) tip(w) <- tip
	if (!missing(label)) label(w) <- label
	if (!missing(tag)) tag(w) <- tag
	## Special treatment for 'wrap' argument that does not exists in ttk::label
	if (!missing(wrap)) {
		## We need width
		width <- abs(as.integer(tclvalue(tkcget(w, "-width"))))
		if (isTRUE(wrap)) wraplength <- .wraplength(w, width) else wraplength <- 0
		tkconfigure(w, wraplength = wraplength)
		## If width is not reapplied after wraplength, the text is not always
		## wrapped in the widget (is this a bug?)
		if (wraplength > 0 && length(width)) tkconfigure(w, width = width)
	}	
	return(w)
}

tk2labelframe <- function (parent, ...)
{
	if (!is.ttk()) stop("Tcl/Tk >= 8.5 is required")
	w <- tkwidget(parent, "ttk::labelframe", ...)
	class(w) <- c("ttk2labelframe", "tk2widget", class(w))
	return(w)
}

## TODO: I need to rework this on the same scheme as tk2label
tk2listbox <- function (parent, values, value, selection,
selectmode = c("extended", "single", "browse", "multiple"), height = 5,
tip = "", scroll = "both", autoscroll = "x", enabled = TRUE, ...)
{
	## Check conditions and arguments
	if (!is.ttk()) stop("Tcl/Tk >= 8.5 is required")
	selectmode <- match.arg(selectmode)
	scrolls <- c("x", "y", "both", "none")
	scroll <- match.arg(scroll, scrolls)
	autoscroll <- match.arg(autoscroll, scrolls)
	if (missing(values) || length(values) == 0) {
		values <- NULL
	} else {
		values <- as.character(values)
	}
	## We provide either value or selection... translate value into selection
	if (missing(selection)) selection <- NULL
	if (!missing(value)) {
		if (is.null(values) || length(values) < 1) selection <- NULL else
			selection <- (1:length(values))[values %in% as.character(value)]
	}
	if (length(selection) == 0 || is.null(values) || length(values) < 1) {
		selection <- NULL
	} else {
		selection <- sort(as.integer(round(selection)))
		if (selection[1] < 1)
			stop("Numerical selections must be indices > 0")
		if (selection[1] > length(values)) selection <- NULL
		if (selectmode == "single" && length(selection) > 1) {
			warning("Only lowest selection used in single selection mode")
			selection <- selection[1]
		}
	}
	
	## Location of the widget depends if we add scrollbars or not
	background <- tclvalue(.Tcl("ttk::style lookup TEntry -fieldbackground"))
	if (background == "") background <- "white"
	if (scroll == "none") {
		w <- tkwidget(parent, "listbox", font = "TkDefaultFont",
			borderwidth = 1, relief = "sunken", activestyle = "dotbox",
			selectmode = selectmode, height = height, exportselection = 0,
			background = background, ...)
	} else { # We need to create a tk2frame as parent of the listbox
		wf <- tk2frame(parent)
		w <- tkwidget(wf, "listbox", font = "TkDefaultFont",
			borderwidth = 1, relief = "sunken", activestyle = "dotbox",
			selectmode = selectmode, height = height, exportselection = 0,
			background = background, ...)
	}
	## Make it react to tk2theme changes, and integrate the listbox as much
	## as possible with current ttk theme
	#restyleListbox <- function (W) {
	#	## Restyle the listbox according to current ttk style
	#	## Note: font is set to TkDefaultFont => already managed there!
	#	tkconfigure(W,
	#		foreground = tk2style("tk2entry", "foreground",
	#			default = "#000000"),
	#		borderwidth = tk2style("", "borderwidth", default = 0),
	#		disabledforeground = tk2style("tk2entry", "foreground",
	#			"disabled", default = "#a3a3a3"),
	#		highlightbackground = tk2style("tk2entry", "selectbackground",
	#			default = "#c3c3c3"),
	#		highlightcolor = tk2style("tk2entry", "selectbackground",
	#			default = "#c3c3c3"),
	#		selectbackground = tk2style("tk2entry", "selectbackground",
	#			default = "#c3c3c3"),
	#		selectforeground = tk2style("tk2entry", "selectforeground",
	#			default = "#ffffff")
	#	)
	#}
	## Restyle it now
	#restyleListbox(w)
	restyleListbox <- function (W) {
		background <- tclvalue(.Tcl("ttk::style lookup TEntry -fieldbackground"))
		if (background == "") background <- "white"
		tkconfigure(W, background = background)
	}
	
	## If there are values and/or selections, populate the list now
	for (item in values)
		tkinsert(w, "end", item)
	if (!is.null(selection)) {
		for (sel in selection)
			tkselection.set(w, sel - 1) # Because Tcl uses 0-based indexing!
		tksee(w, selection[1]) # Ensure that the first selected item is visible	
	}
	
	## Possibly add a tooltip
	if (tip != "") tk2tip(w, tip)
	
	## Do we add scrollbars?
	if (scroll == "none") {
		## Apply bindings to original listbox
		tkbind(w, "<<ThemeChanged>>", restyleListbox)
		tkbind(w, "<1>", function (W) tkfocus(W)) # Needed for mouseweel action
		## Do we disable it?
		if (!isTRUE(enabled)) tkconfigure(w, state = "disabled")
		## Done... just return the widget
		class(w) <- c("tk2listbox", "tk2widget", class(w))
		return(w)
	} else {
		## Add (autohide) scrollbar(s)
		tcl("scrolledWidget", w, wf, scroll, autoscroll)
		## Apply bindings to frame container
		tkbind(wf, "<<ThemeChanged>>", restyleListbox)
		tkbind(wf, "<1>", function (W) tkfocus(W)) # Needed for mouseweel action
		class(wf) <- c("tk2listbox", "tk2widget", class(w))
		## Do we disable it?
		if (!isTRUE(enabled)) tkconfigure(wf, state = "disabled")
		return(wf)
	}
}

tk2mclistbox <- function (parent, tip ="", ...)
{
### TODO: a tile equivalent of this widget
### TODO: or adjust the header: font, color and frame
	if (!is.ttk()) stop("Tcl/Tk >= 8.5 is required")
	res <- tclRequire("mclistbox")
	if (!inherits(res, "tclObj"))
		stop("Impossible to load the Tcl mclistbox package; check your Tcl/Tk installation")
	
	background <- tclvalue(.Tcl("ttk::style lookup TEntry -fieldbackground"))
	if (background == "") background <- "white"
	w <- tkwidget(parent, "mclistbox::mclistbox", font = "TkDefaultFont",
		background = background, ...)
	tkconfigure(w, relief = "flat")
	if (tip != "") tk2tip(w, tip)
	class(w) <- c("tk2mclistbox", "tk2widget", class(w))
	return(w)
}

#tk2menu <- function (parent, ...)
#{
#	if (!is.ttk()) stop("Tcl/Tk >= 8.5 is required")
#	w <- tkwidget(parent, "menu", ...)
#	class(w) <- c("tk2menu", "tk2widget", class(w))
#	return(w)
#}

tk2menu <- function (parent, activebackground, activeforeground, ...) 
{
    if (!is.ttk()) stop("Tcl/Tk >= 8.5 is required")
    w <- tkwidget(parent, "menu", ...)
    if (missing(activebackground))
		activebackground <- tk2style("tk2button", "selectbackground")
    if (activebackground == "") activebackground = "darkblue" # Default value
    if (missing(activeforeground))
		activeforeground <- tk2style("tk2button", "selectforeground")
    if (activeforeground == "") activeforeground = "white" # Default value
    tkconfigure(w, activebackground = activebackground,
		activeforeground = activeforeground)
    class(w) <- c("tk2menu", "tk2widget", class(w))
    return(w)
}


tk2menubutton <- function (parent, tip = "", ...)
{
	if (!is.ttk()) stop("Tcl/Tk >= 8.5 is required")
	w <- tkwidget(parent, "ttk::menubutton", ...)
	if (tip != "") tk2tip(w, tip)
	class(w) <- c("tk2menubutton", "tk2widget", class(w))
	return(w)
}

tk2message <- function (parent, text = "", justify = c("left", "center", "right"),
width = -1, aspect = 150, tip = "", ...)
{
	if (!is.ttk()) stop("Tcl/Tk >= 8.5 is required")
	justify <- as.character(justify[1])
	w <- tkwidget(parent, "message", text = text, justify = justify,
		width = width, aspect = aspect, font = "TkDefaultFont", ...)
	if (tip != "") tk2tip(w, tip)
	class(w) <- c("tk2message", "tk2widget", class(w))
	return(w)
}

tk2notebook <- function (parent, tabs, ...)
{
	if (!is.ttk()) stop("Tcl/Tk >= 8.5 is required")
	w <- tkwidget(parent, "ttk::notebook", ...)
	## Add pages
	tabs <- as.character(tabs)
	for (tab in tabs) {
		tframe <- tk2frame(w)
		tkadd(w, tframe, text = tab, sticky = "nsew")
	}
	tk2notetraverse(w)  # Enable keyboard traversal for this notebook
	class(w) <- c("tk2notebook", "tk2widget", class(w))
	return(w)
}

tk2panedwindow <- function (parent, orientation = c("horizontal", "vertical"), ...)
{
	if (!is.ttk()) stop("Tcl/Tk >= 8.5 is required")
	orientation <- as.character(orientation[1])
	w <- tkwidget(parent, "ttk::panedwindow", orient = orientation, ...)
	class(w) <- c("tk2panedwindow", "tk2widget", class(w))
	return(w)
}

tk2progress <- function (parent, orientation = c("horizontal", "vertical"),
tip = "", ...)
{
	if (!is.ttk()) stop("Tcl/Tk >= 8.5 is required")
	orientation <- as.character(orientation[1])
	w <- tkwidget(parent, "ttk::progressbar", ...)
	if (tip != "") tk2tip(w, tip)
	class(w) <- c("tk2progress", "tk2widget", class(w))
	return(w)
}

tk2radiobutton <- function (parent, tip = "", ...)
{
### TODO: associate with a variable and set both states values
	if (!is.ttk()) stop("Tcl/Tk >= 8.5 is required")
	w <-tkwidget(parent, "ttk::radiobutton", ...)
	if (tip != "") tk2tip(w, tip)
	class(w) <- c("tk2radiobutton", "tk2widget", class(w))
	return(w)
}

tk2scale <- function (parent, orientation = c("horizontal", "vertical"),
tip = "", ...)
{
	if (!is.ttk()) stop("Tcl/Tk >= 8.5 is required")
	orientation <- as.character(orientation[1])
	w <- tkwidget(parent, "ttk::scale", orient = orientation, ...)
	if (tip != "") tk2tip(w, tip)
	class(w) <- c("tk2scale", "tk2widget", class(w))
	return(w)
}

tk2scrollbar <- function (parent, orientation = c("horizontal", "vertical"), ...)
{
	if (!is.ttk()) stop("Tcl/Tk >= 8.5 is required")
	orientation <- as.character(orientation[1])
	w <- tkwidget(parent, "ttk::scrollbar", orient = orientation, ...)
	class(w) <- c("tk2scrollbar", "tk2widget", class(w))
	return(w)
}

tk2separator <- function (parent, orientation = c("horizontal", "vertical"), ...)
{
	if (!is.ttk()) stop("Tcl/Tk >= 8.5 is required")
	orientation <- as.character(orientation[1])
	w <-tkwidget(parent, "ttk::separator", orient = orientation, ...)
	class(w) <- c("tk2separator", "tk2widget", class(w))
	return(w)
}

tk2spinbox <- function (parent, tip = "", ...)
{
	if (!is.ttk()) stop("Tcl/Tk >= 8.5 is required")
	## Default background to fieldbackground
	if (any(names(list(...)) == "background")) {
        w <- tkwidget(parent, "spinbox", font = "TkDefaultFont",
			relief = "solid", borderwidth = 1, ...)
    } else {
        background <- tclvalue(.Tcl("ttk::style lookup TEntry -fieldbackground"))
		if (background == "") background <- "white"
		w <- tkwidget(parent, "spinbox", font = "TkDefaultFont",
			relief = "solid", borderwidth = 1,
            background = background, ...)
    }
	
	if (tip != "") tk2tip(w, tip)
	class(w) <- c("tk2spinbox", "tk2widget", class(w))
	return(w)
}

tk2table <- function (parent, ...)
{
### TODO: a pure Tcl equivalent
### TODO: a tile equivalent (some customization, if possible)
	if (!is.ttk()) stop("Tcl/Tk >= 8.5 is required")
	if (inherits(tclRequire("Tktable", warn = FALSE), "tclObj")) {
		w <- tkwidget(parent, "table", font = "TkDefaultFont", ...)
		class(w) <- c("tk2table", "tk2widget", class(w))
		return(w)
	} else stop("Tcl package 'Tktable' must be installed first")
}

tk2tablelist <- function (parent, ...) 
{
    if (!is.ttk()) 
        stop("Tcl/Tk >= 8.5 is required")
    if (inherits(tclRequire("tablelist_tile", warn = FALSE), "tclObj")) {
        ## Default background to fieldbackground
		if (any(names(list(...)) == "background")) {
			w <- tkwidget(parent, "tablelist::tablelist",
				font = "TkDefaultFont", ...)
		} else {
			background <- tclvalue(.Tcl("ttk::style lookup TEntry -fieldbackground"))
			if (background == "") background <- "white"
			w <- tkwidget(parent, "tablelist::tablelist",
				font = "TkDefaultFont", background = background, ...)
		}
		
        class(w) <- c("tk2tablelist", "tk2widget", class(w))
        return(w)
    }
    else stop("Tcl package 'tablelist' must be installed first")
}

tk2text <- function (parent, tip = "", ...)
{
### TODO: autohide scrollbars
	if (!is.ttk()) stop("Tcl/Tk >= 8.5 is required")
	
	## Default background to fieldbackground
	if (any(names(list(...)) == "background")) {
		w <- tkwidget(parent, "text", font = "TkTextFont", ...)
	} else {
		background <- tclvalue(.Tcl("ttk::style lookup TEntry -fieldbackground"))
		if (background == "") background <- "white"
		w <- tkwidget(parent, "text", font = "TkTextFont",
			background = background, ...)
	}
	
	tkconfigure(w, relief = "flat")
	if (tip != "") tk2tip(w, tip)
	class(w) <- c("tk2text", "tk2widget", class(w))
	return(w)
}

tk2ctext <- function (parent, tip = "", ...)
{
### TODO: autohide scrollbars
	if (!is.ttk()) stop("Tcl/Tk >= 8.5 is required")
	tclRequire("ctext")
	
	## Default background to fieldbackground
	if (any(names(list(...)) == "background")) {
		w <- tkwidget(parent, "ctext", font = "TkFixedFont", ...)
	} else {
		background <- tclvalue(.Tcl("ttk::style lookup TEntry -fieldbackground"))
		if (background == "") background <- "white"
		w <- tkwidget(parent, "ctext", font = "TkFixedFont",
			background = background, ...)
	}
	
	tkconfigure(w, relief = "flat")
	if (tip != "") tk2tip(w, tip)
	class(w) <- c("tk2ctext", "tk2widget", class(w))
	return(w)
}

### TODO: rework this, using ttk::treeview
tk2tree <- function (parent, tip = "", ...)
{
### TODO: better icons!
	## Reasonable default icons for files and folders
	if (!is.ttk()) stop("Tcl/Tk >= 8.5 is required")
	images <- as.character(tcl("image", "names"))
	if (!"Tree:dir" %in% images)
		.Tcl("image create photo Tree:dir -data {R0lGODdhEAAQAPIAAAAAAHh4eLi4uPj4APj4+P///wAAAAAAACwAAAAAEAAQAAADPVi63P4wLkKCtTTnUsXwQqBtAfh910UU4ugGAEucpgnLNY3Gop7folwNOBOeiEYQ0acDpp6pGAFArVqthQQAO///}")
	if (!"Tree:file" %in% images)
		.Tcl("image create photo Tree:file -data {R0lGODdhEAAQAPIAAAAAAHh4eLi4uPj4+P///wAAAAAAAAAAACwAAAAAEAAQAAADPkixzPODyADrWE8qC8WN0+BZAmBq1GMOqwigXFXCrGk/cxjjr27fLtout6n9eMIYMTXsFZsogXRKJf6uP0kCADv/}")

### TODO: correct support of font
	w <- tkwidget(parent, "Tree:create")  #, font = "TkDefaultFont", ...)
	tkconfigure(w, relief = "flat")
	if (tip != "") tk2tip(w, tip)
	class(w) <- c("tk2tree", "tk2widget", class(w))
	return(w)
}
