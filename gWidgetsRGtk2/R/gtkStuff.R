## try without stack smashing
gtktry = function(expr, silent=TRUE) {
  tryCatch(expr, error = function(e) {
    gwCat(sprintf("Error: %s\n",e))
    msg = conditionMessage(e)
    invisible(structure(msg, class = "try-error"))
  })
}


## set alignment
#Sets the alignment of the child. This property has no effect unless the child is a GtkMisc or a GtkAligment.
# xalign : the horizontal position of the child, 0.0 is left aligned, 1.0 is right aligned
# yalign : the vertical position of the child, 0.0 is top aligned, 1.0 is bottom aligned

setXYalign <- function(child, childWidget, anchor) {
  if(is(child,"GtkMisc") || is(child,"GtkAlignment")) {
    child['xalign'] <- anchor[1]
    child['yalign'] <- anchor[2]
  } else if(!is.null(childWidget)) {
    if(is(childWidget,"GtkMisc") || is(childWidget,"GtkAlignment")) {
      childWidget['xalign'] <- anchor[1]
      childWidget['yalign'] <- anchor[2]
    }
  }
}


         

## return gtk objects from others
getBlock <- function(widget) {
  if(inherits(widget,"<invalid>")) return(NULL)
  if(is(widget,"RGtkObject")) return(widget)
  if(is(widget,"gWidgetRGtk")) return(getBlock(widget@block))
  if(is(widget,"guiWidget")) return(getBlock(widget@widget))
  gwCat(gettext("Can't get block"))
  return(NULL)
}

## return NA or widget
getWidget <- function(widget) {
  if(inherits(widget,"<invalid>")) return(NULL)
  while(!is(widget,"RGtkObject")) {
    if(inherits(widget,"<invalid>")) return(NULL)
    widget = widget@widget
  }
  widget
}

## return GtkWindow if possible
getGtkWindow = function(widget) {
  if(inherits(widget,"guiContainer") || inherits(widget,"guiComponent"))
    widget = getToolkitWidget(widget)

  while(!is(widget,"GtkWindow")) {
    widget = widget$GetParent()
    if(inherits(widget,"<invalid>")) return(NULL)
  }
  return(widget)
}

## Method to interact with toolkit objects
setMethod(".getToolkitWidget",
          signature(obj="gWidgetRGtk", toolkit="guiWidgetsToolkitRGtk2"),
          function(obj, toolkit) getWidget(obj))

## setMethod(".callToolkitMethod",
##           signature(x="gWidgetRGtk", toolkit="guiWidgetsToolkitRGtk2"),
##           function(x, toolkit, meth_name) {
##             widget <- getWidget(x)
##             RGtk2:::.getAutoMethodByName(widget, meth_name, parent.frame())
##           })

setMethod(".getToolkitProperty",
          signature(x="gWidgetRGtk", toolkit="guiWidgetsToolkitRGtk2"),
          function(x, toolkit, property) {
            widget <- getWidget(x)
            RGtk2::gObjectGet(widget,property)
          })

setMethod(".setToolkitProperty",
          signature(x="gWidgetRGtk", toolkit="guiWidgetsToolkitRGtk2"),
          function(x, toolkit, property, value) {
            widget <- getWidget(x)
            widget[property] <- value
            x
          })



RtoGObjectConversion = function(obj) {
  if(inherits(obj,"gComponent")) return("GObject")
  if(is.list(obj)) return("GObject")
  
  Klasse = class(obj)[1]                # silly name?
  switch(Klasse,
         "integer"="gint",
         "numeric"="gdouble",
         "gtk"="GObject",
         "logical" = "gboolean",
         "gchararray"
         )
}

##################################################
##
## gtkTreeViewColumn stuff

setMethod("svalue",signature(obj="GtkTreeViewColumn"),
          function(obj, index=NULL, drop=NULL, ...) {
            theArgs =   list(...)
            index = ifelse(is.null(index),FALSE, as.logical(index))
            drop =  ifelse(is.null(drop), TRUE, as.logical(drop))

            ## is this a treeviewCOlumn that ggrid made?
            col.no = gtktry(tag(obj,"column.number"), silent=TRUE)
            if(inherits(col.no,"try-error"))
              return(NA)

            ## return index if requested
            if(index) return(col.no)
            
            ## else return the values
            
            gridObj = tag(obj,"gridObj")
            vals = gridObj[,col.no, visible=TRUE, drop=drop] # only show visible
            return(vals) 
          })

setMethod("id",signature(obj="GtkTreeViewColumn"),
          function(obj,  ...) {
            curname = tag(obj,"name")
            if(is.null(curname) || length(curname) == 0) {
#              gwCat(gettext("No name for this view column\n"))
              return(NA)
            } else {
              return(curname)
            }
          })


setReplaceMethod("id",signature(obj="GtkTreeViewColumn"),
          function(obj, ..., value) {
            curname = tag(obj,"name")
            if(is.null(curname) || length(curname) == 0) {
              ## not there, set it
              label = glabel(value)
              tag(obj,"widget") <- label

              ## set in view col
              widget <- getBlock(label)  ## block is event box
              obj$setWidget(widget)
#              print(class(widget))
#              tag(obj,"header") <- widget$getParent()$getParent()$getParent()
            } else {
              ## store in widget
              svalue(tag(obj,"widget"))<-value
            }
            tag(obj,"name") <- value
            return(obj)
          })

setMethod("addhandlerchanged",signature(obj="GtkTreeViewColumn"),
          function(obj, handler=NULL, action=NULL, ...) {
            lst = list()                # store ids for handlers
            lst[["cellrenderer"]] = addhandler(obj$GetCellRenderers()[[1]],
                 signal = "edited",
                 handler = handler,
                 action = action
                 )
            ## If view column comes from gdf.R then subsetBy is stored in object
            ## so changes there will propogate adding change to underlying model
            ## proved too slow as it seems to get called repeatedly, and
            ## wouldn't stop by setting return value
            gridObj = tag(obj,"gridObj")
            if(!is.null(gridObj)) {
              doSubsetBy = tag(gridObj,"doSubsetBy")
              if(!is.null(doSubsetBy) && as.logical(doSubsetBy) == TRUE) {
                subsetBy = tag(gridObj,"subsetBy")
                lst[["subsetBy"]] = addhandlerchanged(subsetBy, handler, action)
              }
            }
            return(lst)
          })


setMethod("addHandlerChanged",signature(obj="GtkTreeViewColumn"),
          function(obj, handler=NULL, action=NULL, ...) {
            addhandlerchanged(obj,handler=handler,action=action,...)
          })

setMethod("removehandler",signature(obj="GtkTreeViewColumn"),
          function(obj, ID=NULL,...) {
            removehandler(obj$GetCellRenderers()[[1]],ID,...)
          })
setMethod("removeHandler",signature(obj="GtkTreeViewColumn"),
          function(obj, ID=NULL,...) {
            removehandler(obj,ID=ID,...)
          })



## fix up [ for RGtkDataFrame
## is this needed?
setMethod("[",
          signature(x="RGtkDataFrame"),
          function(x, i, j, ..., drop=TRUE) {
            .leftBracket(x, guiToolkit("RGtk2"), i, j, ..., drop=drop)
          })
setMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="RGtkDataFrame"),
          function(x, toolkit, i, j, ..., drop=TRUE) {

            frame <- as.data.frame(x)
                                        #if (!missing(i) && length(i) > 0 && inherits(i[[1]], "GtkTreePath"))
                                        #	i <- .RGtkCall("R_gtk_tree_paths_to_indices", i)+1
            if(missing(i) && missing(j))
              frame[, , drop=drop]
            else if(missing(i))
              frame[,j, drop=drop]
            else if(missing(j))
              frame[i,, drop=drop]
            else
              frame[i,j,drop=drop]
          })



## which versino of RGtk2
getRGtk2Version = function() {
  m = installed.packages()
  ver = m["RGtk2","Version"]
  ver = unlist(strsplit(ver,"\\."))
  names(ver) <- c("major","minor","version")
  return(ver)
}

## Determin which OS
##' is windows the OS?
##'
##' @return TRUE or FALSE
isWindows <- function() {

}

##' is windows the OS?
##'
##' @return TRUE or FALSE
isMac <- function() {

}

##' is windows the OS?
##'
##' @return TRUE or FALSE
isUNIX <- function() {

}


## mouse click processing

##' Return TRUE if first mouse click
##'
##' To be called from key-press|release-event
##' @param e event for mouse press
##' @return TRUE or FALSE
isFirstMouseClick <- function(e) {
  if(!is(e, "GdkEvent"))
    stop("Must pass in an event")
  e$getButton() == 1
}

##' Return TRUE/FALSE if right mouse click
##'
##' To be called from key-press|release-event
##' @param e event for mouse press
##' @return TRUE or FALSE
isRightMouseClick <- function(e) {
  if(!is(e, "GdkEvent"))
    stop("Must pass in an event")
  
  e$GetButton() == 3 ||
  (Sys.info()["sysname"] == "Darwin" && e$GetState() == GdkModifierType['control-mask'] && e$GetButton() == 1) 
}



