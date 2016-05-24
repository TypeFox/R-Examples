## Various methods for tk2widgets
## See also tk2tip.R for methods tip(w) and tip(w) <- value

is.tk2widget <- function (x)
	return(inherits(x, "tk2widget"))

print.tk2widget <- function (x, ...)
{
	if (disabled(x)) txt <- " (disabled)" else txt <- ""
	cat("A tk2widget of class '", class(x)[1], "'", txt, "\n", sep = "")
	cat("State: ", state(x), "\n", sep = "")
	cursize <- size(x)
	if (cursize > 0) cat("Size: ", cursize, "\n", sep = "")
	val <- value(x)
	if (!is.null(val)) {
		cat("Value:\n")
		print(value(x))
	}
	return(invisible(x))
}

tk2cfglist <- function (...)
{
	res <- list(...)
	class(res) <- c("tk2cfglist", class(res))
	return(res)
}
	
print.tk2cfglist <- function (x, ...)
{
	if (!length(x)) {
		cat("An empty tk2widget cfglist\n")
	} else {
		cat("A tk2widget cfglist with:\n\n")
		print(unclass(x))
	}
	return(invisible(x))
}

state <- function (x, ...)
	UseMethod("state")

state.tk2widget <- function (x, ...)
	return(as.character(tcl(x, "state")))

## TODO: a state.tk2listbox, because there is no state property defined for it!

label <- function (x, ...)
	UseMethod("label")
	
label.tk2widget <- function (x, ...)
	return(x$env$label)
	
`label<-` <- function (x, value)
	UseMethod("label<-")
	
`label<-.tk2widget` <- function (x, value)
{
	x$env$label <- as.character(value)[1]
	return(x)
}

tag <- function (x, ...)
	UseMethod("tag")
	
tag.tk2widget <- function (x, ...)
	return(x$env$tag)
	
`tag<-` <- function (x, value)
	UseMethod("tag<-")
	
`tag<-.tk2widget` <- function (x, value)
{
	x$env$tag <- value
	return(x)
}

disabled <- function (x, ...)
	UseMethod("disabled")

disabled.tk2widget <- function (x, ...)
	return(tclvalue(tkcget(x, "-state")) == "disabled")
	
`disabled<-` <- function (x, value)
	UseMethod("disabled<-")
	
`disabled<-.tk2widget` <- function (x, value)
{
	if (isTRUE(value)) state <- "disabled" else state <- "normal"
	tkconfigure(x, state = state)
	return(x)
}
	
values <- function (x, ...)
	UseMethod("values")
	
values.tk2widget <- function (x, ...)
	return(NULL) # Default value, for widgets that do not support this!

values.tk2listbox <- function (x, ...)
	return(as.character(tkget(x, 0, "end")))

`values<-` <- function (x, value)
	UseMethod("values<-")

`values<-.tk2widget` <- function (x, value)
	stop("This tk2widget does not seem to support values")

`values<-.tk2listbox` <- function (x, value)
{
	## Save current selection
	cursel <- selection(x)
	tclServiceMode(FALSE)
	on.exit(tclServiceMode(TRUE))
	isDisabled <- disabled(x)
	on.exit(disabled(x) <- isDisabled, add = TRUE)
	if (isDisabled) disabled(x) <- FALSE
	## Change items (no attempt to match them -possible future improvement!-)
	tkdelete(x, 0, "end")
	for (item in as.character(value))
		tkinsert(x, "end", item)
	## Try to reapply selection
	for (sel in cursel)
		tkselection.set(x, sel - 1)
	return(x)
}

value <- function (x, ...)
	UseMethod("value")

value.tk2widget <- function (x, ...)
	return(NULL) # Default value is NULL for tk2widgets

value.tk2listbox <- function (x, ...)
	return(values(x)[selection(x)])

`value<-` <- function (x, value)
	UseMethod("value<-")

`value<-.tk2widget` <- function (x, value)
	stop("This tk2widget does not seem to support setting its value")

`value<-.tk2listbox` <- function (x, value)
{
	items <- items(x)
	if (length(items) > 0)
		selection(x) <- (1:length(items))[items %in% value] 
	return(x)
}

selection <- function (x, ...)
	UseMethod("selection")

selection.tk2widget <- function (x, ...)
	return(NULL) # For tk2widgets that do not support selection

selection.tk2listbox <- function (x, ...)
	return(as.integer(tkcurselection(x)) + 1)

`selection<-` <- function (x, value)
	UseMethod("selection<-")

`selection<-.tk2widget` <- function (x, value)
	stop("This tk2widget does not seem to support setting its selection")
	
`selection<-.tk2listbox` <- function (x, value)
{
	## Prepare
	tclServiceMode(FALSE)
	on.exit(tclServiceMode(TRUE))
	isDisabled <- disabled(x)
	on.exit(disabled(x) <- isDisabled, add = TRUE)
	if (isDisabled) disabled(x) <- FALSE
	## Clear selection only
	if (is.null(value) || length(value) < 1) {
		tkselection.clear(x, 0, "end")
		return(x)
	}
	## Check data
	value <- sort(as.integer(round(value)))
	if (value[1] < 1)
		stop("Selections must be indices > 0")
	if (value[length(value)] > size(x)) return(x)
	## Change selection
	tkselection.clear(x, 0, "end")
	if (tclvalue(tkcget(x, "-selectmode")) == "single" && length(value) > 1) {
		warning("Single selection mode... only lowest selection used")
		tkselection.set(x, value[1] - 1)
	} else {
		for (sel in value)
			tkselection.set(x, sel - 1)
	}
	if (!isDisabled) tksee(x, value[1] - 1)
	return(x)
}

visibleItem <- function (x, index, ...)
	UseMethod("visibleItem")

visibleItem.tk2widget <- function (x, index, ...)
	stop("This tk2widget does not seems to support the visibleItem method")

visibleItem.tk2listbox <- function (x, index, ...)
{
	## Index must be a positive integer
	index <- as.integer(round(index))
	if (is.null(index) || length(index) < 1 || index[1] < 1)
		stop("index must be a postive integer")
	tksee(x, index[1] - 1) # Because Tcl uses 0-based indices
	return()
}

size <- function (x, ...)
	UseMethod("size")

size.tk2widget <- function (x, ...)
	return(0L) # By default, a tk2widget has values of zero size (NULL)
	
size.tk2listbox <- function (x, ...)
	return(as.integer(tksize(x)))

config <- function (x, ...)
	UseMethod("config")

config.tk2widget <- function (x, cfglist, ...)
{
	## Compile a list of arguments
	args <- list(...)
	if (!missing(cfglist)) args <- .mergeList(as.list(cfglist, args))
	
	
	## Prepare an empty object
	res <- list()
	class(res) <- c("tk2cfglist", class(res))
	
	## No arguments provided... query a sublist of parameters
	if (length(args) == 0) {
		## Return the complete config (but not the data!)
		params <- tk2configList(x)
		if (!length(params)) return(res)
	} else {
		## Separate named (set) from unnamed (query only) arguments
		params <- names(args)
		if (is.null(params)) { # No named arguments, only queries
			params <- as.character(args)
		} else {
			## For those named arguments, change the config
			res <- (config(x) <- args[params != ""])
			## ... and query the others
			params <- as.character(args[params == ""])
		}
	}
	
	## Retrieve values for the queries
	if (length(params)) {
		for (i in 1:length(params)) {
			p <- tclvalue(tkcget(x, paste("-", params[i], sep = "")))
			if (!is.null(p) && p != "")
				res[[params[i]]] <- p
		}
	}
	return(res)
}
	
config.tk2label <- function (x, cfglist, ...)
{
	## wrap is special here... => how to deal with it???
	## TODO...
	return(config.tk2widget(x, cfglist, ...))
}

## TODO: config.tk2listbox()

`config<-` <- function (x, value)
	UseMethod("config<-")

`config<-.tk2widget` <- function (x, value)
{
	## The default function deleguates to tkconfigure, except for a few things
	value <- .configStd(x, value)
	value$widget <- x
	do.call(tkconfigure, value)
	return(x)
}

`config<-.tk2label` <- function (x, value)
{
	## Standard treatment
	value <- .configStd(x, value)
	
	if (!is.null(value$wrap)) {
		## wrap is not a ttk option but we use it here for convenience
		wrap <- value$wrap
		value$wrap <- NULL
	} else wrap <- NULL

	## For the other parameters, apply tkconfigure() with them
	value$widget <- x
	do.call(tkconfigure, value)
	
	## Do we still have to apply wrap?
	width <- abs(as.integer(tclvalue(tkcget(x, "-width"))))
	if (!is.null(wrap)) {
		if (isTRUE(wrap)) wraplength <- .wraplength(x, width) else wraplength <- 0
		tkconfigure(x, wraplength = wraplength)
	}
	## Reapply width to get correct text wrapping (bug in ttk::label?)
	if (length(width)) tkconfigure(x, width = width)
	return(x)
}

## TODO: `config<-.tk2listbox`