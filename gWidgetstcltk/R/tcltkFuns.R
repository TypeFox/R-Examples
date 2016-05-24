##' Coerce tclObj object to logical value
##'
##' @param x should be a "0" or "1" value
##' @param ... ignored
##' @return a logical or NA
as.logical.tclObj <- function(x, ...) as.logical(as.numeric(x))


##' Does object exists as tcl variable
##'
##' @param x character string with name of variable
##' @return logical
tclObj_exists <- function(x) as.logical(.Tcl(sprintf("info exists %s", x)))

##' create a tcl image from the file
##'
##' @param basename basename of image. We add some bit to avoid filename collisions
##' @param file file of image
##' @return image name
make_tcl_image <- function(basename, file) {
  already_defined <- function(nm) any(nm == as.character(tcl("image","names")))
  
  nm <- sprintf("::tcl::%s", basename)
  if(!already_defined(nm)) {
    tcl("image","create","photo", nm ,file=file)
  }
  return(nm)
}

##' Heuristic to determine if widget is a ttk widget
##'
##' @param x tk object or its id
##' @return logical indicating is  ttk widget or not
isTtkWidget <- function(x) {
  cl <- as.character(tkwinfo("class",x))
  (cl %in% c("Treeview")) ||  grepl("^[A-Z]{2,}", cl)
}

##' what windowing system?
##'
##' @return one of c("x11", "win32", "aqua")
windowingsystem <- function() {
  ## one of x11 (X11-based), win32 (MS Windows), or aqua (Mac OS X Aqu
  as.character(.Tcl("tk windowingsystem"))
}

## return tk widget from obj
## ** This should be a method **
getWidget = function(obj) {
  if(is(obj,"tkwin")) return(obj)

  if(is(obj,"gWidgettcltk"))
    return(getWidget(obj@widget))
  else if(is(obj,"guiWidget"))
    return(getWidget(obj@widget))
  else
    return(NA)
}

getBlock = function(obj) {

  if(is(obj,"tkwin")) return(obj)
  if(is(obj,"gWidgettcltk"))
    return(getBlock(obj@block))
  else if(is(obj,"guiWidget"))
    return(getBlock(obj@widget))
  else
    return(NA)
}


getTopParent = function(tkobj) {
  ## in env is parent variable if present
  ans <- NULL
  
  while(is.null(ans)) {
    e <- tkobj$env$parent
    if(is.list(e)  &&
       e[['ID']] =="")
      ans <- tkobj
    else tkobj <- tkobj$env$parent
  }
  return(ans)
}

getTopLevel <- function(obj) {
  if(is(obj, "guiWidget")) {
    return(getTopLevel(obj@widget))
  } else if(!is.null(obj@e$parentContainer)) {
    return(getTopLevel(obj@e$parentContainer)) 
  } else {
    return(obj)
  }
}

#######################################################
## methods to interact with underlying toolkit object
setMethod(".getToolkitWidget",
          signature(obj="gWidgettcltk", toolkit="guiWidgetsToolkittcltk"),
          function(obj, toolkit) getWidget(obj))



setMethod(".callToolkitMethod",
          signature(x="gWidgettcltk", toolkit="guiWidgetsToolkittcltk"),
          function(x, toolkit, meth_name) {
            widget <- getWidget(x)
            f <- function(...) {
              get(meth_name, parent.frame())(widget, ...)
            }
            f                           # return a function
          })

setMethod(".getToolkitProperty",
          signature(x="gWidgettcltk", toolkit="guiWidgetsToolkittcltk"),
          function(x, toolkit, property) {
            widget <- getWidget(x)
            tkcget(widget, sprintf("-%s", property))
          })

setMethod(".setToolkitProperty",
          signature(x="gWidgettcltk", toolkit="guiWidgetsToolkittcltk"),
          function(x, toolkit, property, value) {
            widget <- getWidget(x)
            l <- list(widget); l[[property]] <- value
            do.call(tkconfigure, l) 

            x
          })




## Does the top level window exists
windowExists = function(obj) {
  win = getTopParent(getWidget(obj))
  as.logical(tkwinfo("exists", win))
}

findTkIcon <- function(i) {
  if(is.null(i)) return("")
  
  stock <- getStockIcons()
  gwi <- system.file(paste("image/",i,".png",sep=""), package="gWidgetstcltk")
  if(is.null(i) || is.na(i) || i == "")
    val <- ""
  else if(file.exists(i))
    val <- i
  else if (file.exists(gwi))
    val <- gwi
  else
    val <- stock[[i, exact=TRUE]]
  ## what to return
  if(is.null(val))
    return("")
  else
    return(val)
}


##################################################
## function to add scrollbars to widget and pack into grid
addScrollbarsToWidget <- function(widget, parent) {
  xscr <- ttkscrollbar(parent, orient="horizontal",
                       command=function(...) tkxview(widget, ...))
  yscr <- ttkscrollbar(parent, orient="vertical",
                       command=function(...) tkxview(widget, ...))

  tkconfigure(widget,
              xscrollcommand=function(...) tkset(xscr,...),
              yscrollcommand=function(...) tkset(yscr,...))

  tkgrid(widget, row=0, column=0, sticky="news")
  tkgrid(yscr,row=0,column=1, sticky="ns")
  tkgrid(xscr, row=1, column=0, sticky="ew")
  tkgrid.columnconfigure(parent, 0, weight=1)
  tkgrid.rowconfigure(parent, 0, weight=1)
}

##' helper to make a treeview and populate from m
##'
##' @param parent parent container. Must use pack manager.
##' @param m character matrix. May have 1 or more columns
##' @return list with components \code{frame} (the enclosing frame for
##' size management, managed by pack); \code{tr}, a treeview widget.
populate_rectangular_treeview <- function(parent, m) {
  enc_frame <- ttkframe(parent)
  frame <- ttkframe(enc_frame)
  tkpack(frame, expand=TRUE, fill="both")
  tr <- ttktreeview(frame,
                    columns=seq_len(ncol(m)),
                    show="headings",
                    selectmode="browse"
                    )
  addScrollbarsToWidget(tr, frame)
  tkpack.propagate(enc_frame, FALSE)


  ## headings,widths
  charWidth <- as.integer(tclvalue(tcl("font","measure","TkTextFont","0")))
  sapply(seq_len(ncol(m)), function(i) {
    tcl(tr, "heading", i, text=colnames(m)[i])
    tcl(tr, "column", i, width=10 + charWidth*max(apply(m, 2, nchar)))
  })
  tcl(tr, "column", ncol(m), stretch=TRUE) # stretch last
  tcl(tr, "column", "#0", stretch=FALSE)   # no strecth on icons
  ## values
  apply(m, 1, function(vals) {
    if(length(vals) == 1) vals <- paste("{", vals, "}", sep="")
    tcl(tr, "insert", "", "end", values=vals)
  })
  return(list(tr=tr, frame=enc_frame))
}

   



### tk2tip.R - Tooltips for Tk widgets
### Copyright (c), Philippe Grosjean (phgrosjean@sciviews.org)
### Licensed under LGPL 3 or above
###
### Changes:
### - 2007-01-01: first version (for tcltk2_1.0-0)
###
### To do:
### - add and check catch instructions here

## JV: Rather than load in tcltk2 dependency, we borrow Philippe's work here
tk2tip <- function (widget, message)
{
##	if (!is.tk()) stop("Package Tk is required but not loaded")
	if (is.null(message)) message <- ""
	res <- tclRequire("tooltip")
	if (inherits(res, "tclObj")) {
		res <- tcl("tooltip::tooltip", widget, message)
		## Store tip text in the object (use NULL instead of "" for no tip)
		if (message == "") message <- NULL
		widget$env$tip <- message
	} else stop("cannot find tcl package 'tooltip'")
	return(invisible(res))
}

tk2killtip <- function ()
{
	## if (!is.tk()) stop("Package Tk is required but not loaded")
	return(invisible(tcl("tooltip::hide")))
}

## Get tip method
tip <- function (x, ...)
	UseMethod("tip")

tip.tk2widget <- function (x, ...)
	return(x$env$tip)

## Chenge tip method
`tip<-` <- function (x, value)
	UseMethod("tip<-")

`tip<-.tk2widget` <- function (x, value)
{
	tk2tip(x, value)
	return(x)
}
