# wrappers for RGtk1 functions to RGtk2

# containers

gtkChildren <- gtkContainerGetChildren
gtkParent <- gtkWidgetGetParent

# gdk manuals

gdkGetWindowSize <- function(w) {
	geom <- gdkWindowGetGeometry(w)
	c(geom["width"], geom["height"])
}

# signals

gtkObjectSignalEmit <- function(obj, signal, ...) gSignalEmit(obj, signal, ...)
gtkObjectGetSignals <- gObjectGetSignals
gtkTypeGetSignals <- gTypeGetSignals
.GtkClasses <- "GtkWidget"
getSignalInfo <- function(classes = .GtkClasses, load = TRUE) {
  sapply(classes, function(type) sapply(gtkTypeGetSignals(type), gtkSignalGetInfo))
}

gtkSignalGetInfo <- gSignalGetInfo


# args

gtkObjectGetArgInfo <- gObjectGetPropInfo
gtkObjectGetArgs <- function(obj, argNames) obj$get(argNames)
gtkObjectGetArg <- 
function(obj, argName) 
{
	props <- obj$get(argName)
	if (!is.null(props))
		props <- props[[1]]
	props
}
gtkObjectSetArgs <- function(obj, ..., .vals) obj$set(..., .vals)
names.GtkObject <- names.GObject

# dnd

gtkTargetEntry <- gtkTargetEntryNew <- 
function(target, flags, info) as.GtkTargetEntry(list(target, flags, info))

# gtk text stuff ---- not supported!
gtkTextGetText <- function(w) gtkEditableGetChars(w, 0, -1)
gtkTextClearText <- function(w, start = 0, end = -1) gtkEditableDeleteText(w, start, end)
gtkTextSetText <- function(w, contents="", append=FALSE) {
	if (append)
		pos <- nchar(gtkTextGetText(w))
	else {
		gtkTextClearText(w)
		pos <- 0
	}
	
	gtkEditableInsertText(w, contents, pos)
}

# flags

gtkWidgetGetFlags <- function(w) gtkObjectFlags(w)

# this function should not be necessary in well-structured code

findWidgetByType <-
  #
  # Recursively search a widget tree for the first occurence of
  # a widget of the specified type.
  #
function(win, gtkType = "GtkMenuBar", verbose = FALSE)
{
 if(verbose)
   print(class(win))

 if(is.function(gtkType)) {
   if(gtkType(win))
     return(win)
 } else if(as.character(gtkType) %in% class(win)) {
     return(win)
 } 

 if("GtkContainer" %in% class(win)) {
  for(i in win$GetChildren()) {
   tmp <- findWidgetByType(i, gtkType, verbose = verbose)
   if(!is.null(tmp))
     return(tmp)
  }
 }
 
 return(NULL)
}

# CList convenience access - note that this widget is deprecated 

gtkCListGetText <-
function(w, row, cols, zeroBased = TRUE)
{
  checkPtrType(w, "GtkCList")

  if(missing(cols) && is.matrix(row) && ncol(row) == 2)
     which = row
  else
     which = cbind(row, cols)

  if(!zeroBased)
    which <- which - 1
  
  storage.mode(which) <- "integer"

  .Call("R_gtkCListGetText", w, as.integer(t(which)), PACKAGE = "RGtk2")
}


gtkCListSetText <-
function(w, row, cols, values, zeroBased = TRUE)
{
  checkPtrType(w, "GtkCList")  

  if(missing(cols) && is.matrix(row) && ncol(row) == 2)
     which = row
  else
     which = cbind(row, cols)

  if(!zeroBased)
    which <- which - 1
  
  storage.mode(which) <- "integer"

  values <- rep(as.character(values), length = nrow(which))
  
  invisible(.Call("R_gtkCListSetText", w, which, values, PACKAGE = "RGtk2"))
}

# types

gtkObjectGetTypeName <- function(w)
{
 checkPtrType(w, "GtkObject")
 class(w)[1]
}
gtkObjectGetClasses <- function(w, check = TRUE)
{
 if(check)
     checkPtrType(w, "GtkObject")
 class(w)
}
gtkTypeGetClasses <- gTypeGetAncestors
gtkObjectGetType <- function(w, check = TRUE)
{
 if(check)
    checkPtrType(w, "GtkObject")
 as.GType(class(w)[1])
}
gtkGetType <- function(name) as.GType(name)

# widgets

gtkTopWindow <- 
function(title="My Window", show = TRUE)
{
	window <- gtkWindowNew("toplevel", show = show)
	window$setTitle(title)
	window
}
gtkAdd <-
function(parent, ...)
{
  widgets <- list(...)
  if(length(widgets) == 0)
    stop("No widgets to add to parent")

  if(!all(sapply(widgets, checkPtrType, "GtkWidget"))) {
    stop("Non widget objects passed to gtkAdd()")
  }

  sapply(widgets, function(widget) parent$add(widget))
}
gtkShow <- function(..., all=T)
{
	widgets <- list(...)
	func <- gtkWidgetShow
	if (all)
		func <- gtkWidgetShowAll
	sapply(widgets, func)
}

gtkAddCallback <- gtkObjectAddCallback <- 
function(w, signal, f, data = NULL, object = TRUE, after = TRUE)
  gSignalConnect(w, signal, f, data, after, object)

gtkObjectRemoveCallback <- gtkObjectDisconnectCallback <- gSignalHandlerDisconnect
gtkObjectBlockCallback <- gSignalHandlerBlock
gtkObjectUnblockCallback <- gSignalHandlerUnblock
gtkAddTimeout <- gTimeoutAdd
gtkRemoveTimeout <- gtkRemoveIdle <- gSourceRemove
gtkAddIdle <- gIdleAdd

