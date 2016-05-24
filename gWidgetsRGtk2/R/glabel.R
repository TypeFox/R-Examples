setClass("gLabelRGtk",
         contains="gComponentRGtk",
         prototype=prototype(new("gComponentRGtk"))
         )

## constructor
setMethod(".glabel",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   text= "", markup = FALSE, editable = FALSE, handler = NULL, 
                   action = NULL, container = NULL, 
                   ...
                   ) {

            force(toolkit)

            
            label <- gtkLabelNew()


            obj <- as.gWidgetsRGtk2(label)
#            obj = new("gLabelRGtk",block=evb, widget=label,toolkit=toolkit)

            
            if(markup) 
              tag(obj,"markup")<-TRUE
            else
              tag(obj, "markup") <- FALSE

            svalue(obj) <- text

            if(editable) {
              tag(obj,"editable") <- TRUE
              edit <- gedit()
              tag(obj, "edit") <- edit

              editWidget <- getWidget(edit)
              evb <- getBlock(obj)
              
              ## this are almost identical, as we just swap edit and label
              addHandlerChanged(edit,  handler =
                                function(h,...) {
                                  svalue(obj) <- svalue(edit)
                                  evb$Remove(evb[[1]])
                                  evb$Add(label)
                                })
              ##This is for connecting to the third mosue
              id <- addHandlerClicked(obj,
                                      handler=function(h,...) {
                                        svalue(edit) <- svalue(obj)
                                        evb$Remove(evb[[1]])                                        
                                        evb$Add(editWidget)
                                        editWidget$GrabFocus()
                                      })
              tag(obj, "handler.id") <- id

              handler <- NULL           # no editable with handler
            }
            
            if(!is.null(handler)) {
              tag(obj,"handler.id") <- addHandlerClicked(obj, handler=handler,action=action)
            }
            
            ## attach?
            if (!is.null(container)) {
              if(is.logical(container) && container == TRUE)
                container = gwindow()
              add(container, obj,...)
            }
            
            invisible(obj)
          })

## coerce gtk object
as.gWidgetsRGtk2.GtkLabel <- function(widget,...) {
  label = widget
  ## pack into an event box so that we can get signals
  ## doesn't work if there is already a parent!
  evb <- gtkEventBoxNew()
  ## Issue with labels and notebooks can be fixed here, but
  ## may mask editable event response
  ## Thanks to Felix A. for this
  evb$SetVisibleWindow(FALSE)
  
  if(is.null(label$GetParent()))
    evb$Add(label)
#  else
#    cat("Can't add gwidget handlers to this label\n")

  obj <- new("gLabelRGtk",
    block=evb, widget=label, toolkit=guiToolkit("RGtk2"))

  ## tag values -- may already be set (asgWidget(getToolkitWidget(widget)))
  vals <- c("markup"=FALSE)
  for(i in names(vals)) {
    if(!is.null(tag(label,i)))
      tag(obj,i) <- tag(label,i)
    else
      tag(obj,i) <- vals[i]
  }
  return(obj)
}


## methods
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gLabelRGtk"),
          function(obj, toolkit, index=NULL, drop=NULL, ..) {
            markup = tag(obj, "markup")
            if(is.null(markup)) markup = FALSE

            val = obj@widget$GetText()
            if(!is.empty(markup) && markup==TRUE)
              val = gsub("<[^>]*>","",val)    # strip off
            return(val)
          })

## svalue<-
setReplaceMethod(".svalue",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gLabelRGtk"),
          function(obj, toolkit, index=NULL, ..., value) {
            widget <- getWidget(obj)
            ## set the text
            markup <- tag(obj, "markup")
            if(is.null(markup))
              markup = FALSE

            ## if multiline, collapse with \n
            value <- paste(value, collapse="\n")
            
            if(as.logical(markup)==TRUE)
              widget$SetMarkup(value)
            else
              widget$SetText(value)

            return(obj)
          })


## special GTK method for rotation
setGeneric(".rotatelabel",function(obj, angle, ...) standardGeneric(".rotatelabel"))
setMethod(".rotatelabel",
          signature("gLabelRGtk"),
          function(obj, angle, ...) {  
            obj@widget$SetAngle(angle)
          }
        )

##################################################
## handlers
## need to put handler on evb -- not widget
setMethod(".addHandler",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gLabelRGtk"),
          function(obj, toolkit, signal, handler, action=NULL, ...) {
            f <- function(h,...) {
              if(h$obj@widget['sensitive'])
                handler(h,...)
            }
            ID = .addHandler(obj@block, toolkit, signal, f, action, actualobj=obj,...)
            invisible(ID)
          })

setMethod(".addhandlerclicked",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gLabelRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            .addHandler(obj, toolkit, signal="button-press-event",
                        handler=handler, action=action, actualobj=obj,...)
          })

setMethod(".addhandlerchanged",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gLabelRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            edit = tag(obj, "edit")
            if(!is.null(edit)) {
              ## we use unrealize here, the addhandlerchanged on edit wasn't
              ## working for some strage reason
              return(addhandlerunrealize(edit, handler, action))
            } else {
              ## use addhandlerclicked
              return(.addhandlerclicked(obj, toolkit, handler, action, ...))
            }
          })

### need to fuss with evb vs. label
setMethod(".adddroptarget",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gLabelRGtk"),
          function(obj, toolkit, targetType="text", handler=NULL, action=NULL, ...) {
            ## problem -- we want to add drop target to obj@block evb,
            ## but have handler refer to obj@widgeg=label. 
            addDropTarget(obj@block, toolkit, targetType, handler, action, actualobj=obj)
            
          })

setMethod(".adddropsource",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gLabelRGtk"),
          function(obj, toolkit, targetType="text", handler=NULL, action=NULL, ...) {
            ## problem -- we want to add drop target to obj@block evb,
            ## but have handler refer to obj@widgeg=label. 
            addDropSource(obj@block, toolkit, targetType, handler, action, actualobj=obj)
            
          })


## Put onto block
setMethod(".addpopupmenu",signature(toolkit="guiWidgetsToolkitRGtk2", obj="gLabelRGtk"),
          function(obj, toolkit, menulist, action=NULL, ...) {
            addPopupMenuWithSignal(obj@block, toolkit , menulist, action, actualobj=obj,...)
          })
setMethod(".add3rdmousepopupmenu",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gLabelRGtk"),
          function(obj, toolkit, menulist,action=NULL, ...) {
            add3rdMousePopupMenuWithSignal(obj@block, toolkit,
                                           menulist, action, actualobj=obj,...)
          })
##################################################
## internal function -- used by gvariables in  gcommandline
setGeneric("gaddlabel", function(obj, text="", markup=FALSE, pos=1, container=NULL, ...) standardGeneric("gaddlabel"))

setMethod("gaddlabel",
          signature("guiWidget"),
          function(obj, text="", markup=FALSE, pos=1, container=NULL, ...)
          gaddlabel(obj@widget, text, markup, pos, container, ...)
        )

setMethod("gaddlabel",
          signature("gWidgetRGtk"),
          function(obj, text="", markup=FALSE, pos=1, container=NULL, ...) {
            ## wrap widget into a new package with label
            if(pos ==2 || pos == 4) {
              group = ggroup(horizontal=TRUE,container=container,
                toolkit=obj@toolkit)
            } else {
              group = ggroup(horizontal=FALSE,container=container,
                toolkit=obj@toolkit)
            }
            
            
            if(pos ==2 || pos == 3) {
              glabel(text, markup=markup, container=group, toolkit=obj@toolkit)
              add(group, obj,expand=TRUE)
            } else {
              add(group, obj,expand=TRUE)
              glabel(text, markup=markup, container=group, toolkit=obj@toolkit)
            }
            ## group is returned. No methods added here, just a new package
            return(group)
          })
          
