## expander group, like a group, only expands, contracts if requested
## inherits from ggroup, see ggroup's arguments: horizontal, spacing, container
setClass("gExpandgroupRGtk",
         contains="gGroupRGtk",
         prototype=prototype(new("gGroupRGtk"))
         )

setMethod(".gexpandgroup",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   text="", markup=FALSE, horizontal=TRUE,
                   handler=NULL, action=NULL,
                   container = NULL, ...){

            force(toolkit)
            
            expander = gtkExpanderNew()
            if(markup)
              expander$SetUseMarkup(TRUE)
            if(text != "")
              expander$SetLabel(text)

            obj <- as.gWidgetsRGtk2(expander, horizontal=horizontal)


            theArgs = list(...)

            if(!is.null(container)) {
              if(is.logical(container) && container == TRUE)
                container = gwindow(visible=TRUE)
            
              if(!is.null(theArgs$expand) && theArgs$expand)
                add(container,obj,expand=TRUE)
              else
                add(container,obj)
            }
            
            if(!is.null(handler)) 
              tag(obj, "handler.id") <- addhandlerchanged(obj, handler, action)

            invisible(obj)
          })

as.gWidgetsRGtk2.GtkExpander <- function(widget,...) {
  ## coverting from gWidget?
  if(!is.null(tag(widget,"group"))) {
    group <- tag(widget,"group")
  } else {
    theArgs <- list(...)

    horizontal <- if(is.null(theArgs$horizontal)) TRUE else theArgs$horizontal
    spacing <- if(is.null(theArgs$spacing)) 5 else theArgs$spacing

    group = ggroup(horizontal=horizontal, spacing=spacing)
    widget$Add(getBlock(group)) # down from guiWidget to gWidgetRGtk
  }
  ## we put widget=group here to get gGroup methods, but
  ## must be careful below to use "block" when referring to expander
  obj = new("gExpandgroupRGtk",block=widget,widget=getWidget(group),
    toolkit=guiToolkit("RGtk2"))

  tag(obj,"group") <- group
  return(obj)
}

## methods

## value refers to border width
## but it used to refer to the label, we keep this here but suggest
## names be used instead
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gExpandgroupRGtk"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {
            gwCat("Use names() to access label")
            obj@block$GetLabel()        # not @widget@
          })

## if numeric -- set padding to match ggroup
## else set as a label
setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gExpandgroupRGtk",
                           value = "numeric"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   ## set as padding
                   getWidget(obj)$SetBorderWidth(value)                   
                   return(obj)
                 })
## set label, but deprecated
setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gExpandgroupRGtk"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   .Deprecated("names<-",
                               msg = "Use the names<- method to the label")
                   obj@block$SetLabel(value)
                   return(obj)
                 })

## ## names refers to label
setMethod(".names",
          signature(toolkit="guiWidgetsToolkitRGtk2", x="gExpandgroupRGtk"),
          function(x,toolkit) {
            x@block$GetLabel()
          })

setReplaceMethod(".names",
                 signature(toolkit="guiWidgetsToolkitRGtk2",x="gExpandgroupRGtk"),
                 function(x,toolkit,value) {
                   obj@block$SetLabel(value)
                   return(x)
                 })

## Is widget expanded?
setMethod(".visible",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gExpandgroupRGtk"),
                 function(obj, toolkit, ...) {
                   obj@block$getExpanded()
                 })

## control expand/close with logical
setReplaceMethod(".visible",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gExpandgroupRGtk"),
                 function(obj, toolkit, ..., value) {
                   obj@block$SetExpanded(as.logical(value))
                   return(obj)
                 })

## names refers to label
setMethod(".names",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gExpandgroupRGtk"),
          function(x, toolkit) {
            x@block$GetLabel()        # not @widget@
          })

setReplaceMethod(".names",
                 signature(toolkit="guiWidgetsToolkitRGtk2",x="gExpandgroupRGtk"),
                 function(x, toolkit, value) {
                   x@block$SetLabel(value)
                   return(x)
                 })

## set font
setReplaceMethod(".font",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gExpandgroupRGtk"),
                 function(obj, toolkit, value) {
                   label <- obj@block[[2]]
                   label <- gWidgetsRGtk2:::as.gWidgetsRGtk2(label)
                   font(label) <- value
                   return(obj)
                 })


## handlers
## putonto expander in @block
setMethod(".addhandlerchanged",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gExpandgroupRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            addhandler(obj@block, "activate",handler, action,...)
          })
