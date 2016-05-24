setClass("gButtonRGtk",
         contains="gComponentRGtk",
         prototype=prototype(new("gComponentRGtk"))
         )


setMethod(".gbutton",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   text="", border=TRUE, handler=NULL, action=NULL, container=NULL,...
                   ) {

            force(toolkit)

            iconname <- getstockiconname(tolower(text))
            if(!is.na(iconname)) {
              button <- gtkButtonNewFromStock(iconname)
              button$Show()
            } else {
              button <- gtkButtonNewWithLabel(text)
            }

            ## look for border request
            if(border == FALSE) 
              button$SetRelief(2L)

            
            
            obj <- as.gWidgetsRGtk2(button)


            ##            obj <- new("gButtonRGtk",
            ##              block=button, widget=button, toolkit=toolkit)

            tag(obj, "default_fill") <- "x"
            ## add to container
            if (!is.null(container)) {
              if(is.logical(container) && container == TRUE)
                container <- gwindow(visible=TRUE, toolkit=toolkit)
              add(container, obj,...)
            }

            ## add handler
            if (!is.null(handler)) {
              tag(obj,"handler.id") <-  addhandlerclicked(obj,handler,action)
            }

            invisible(obj)
          })

## coerce gtk object
as.gWidgetsRGtk2.GtkButton <- function(widget,...) {
  parent <- widget$parent
  if(is.null(parent)) {
    parent <- gtkAlignmentNew(xscale=1, yscale=0)
    parent$add(widget)
  }

  obj <- new("gButtonRGtk",
    block=parent, widget=widget, toolkit=guiToolkit("RGtk2"))
  return(obj)
}

## constructor for actions
## proper call is gbutton(action = gaction_instnace, cont = ...)
setMethod(".gbutton",
          signature(toolkit="guiWidgetsToolkitRGtk2",
                    action = "guiComponent"),
          function(toolkit,
                   text="", border=TRUE, handler=NULL, action=NULL, container=NULL,...
                   ) {
            if(is(action@widget, "gActionRGtk")) {
              .gbutton(toolkit,  "", border, handler, action@widget, container, ...)
            } else {
              callNextMethod(toolkit, text, border, handler, action, container, ...)
            }
          })

setMethod(".gbutton",
          signature(toolkit="guiWidgetsToolkitRGtk2",
                    action = "gActionRGtk"),
          function(toolkit,
                   text="", border=TRUE, handler=NULL, action=NULL, container=NULL,...
                   ) {
            force(toolkit)

            gtkaction <- getWidget(action)
            
            button <- gtkButton()
            obj <- new("gButtonRGtk",
                       block=button, widget=button, toolkit=guiToolkit("RGtk2"))


            action@e$buttons <- c(action@e$buttons, obj)
            
            #gtkaction$connectProxy(button)
            button$setRelatedAction(gtkaction)
            button['use-action-appearance'] <- TRUE
            ## icon
            icon <- gtkaction['stock-id']
            if(!is.null(icon)) {
              image <- gtkaction$createIcon(GtkIconSize[4])
              button$setImage(image)
            }
            ## tooltip
            tip <- gtkaction['tooltip']
            if(!is.null(tip))
              tooltip(obj) <- tip

            if(!is.null(container)) {
              if(is.logical(container) && container) {
                container <- gwindow()
                add(container, obj)
              } else {
                add(container, obj, ...)
              }
            }
                
            
            return(obj)
          })

### methods
##' return button text
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gButtonRGtk"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {
            return(obj@widget$GetLabel())
          })

##' set button text
setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gButtonRGtk"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   button = obj@widget
                   image = gtkImageNew()
                   if(is.gImage(value)) {
                     filename = tag(value,"filename")
                     if(!is.na(getstockiconname(filename))) {
                       ## stock
                       image$SetFromStock(getstockiconname(filename),size=obj$size)
                     } else {
                       image$SetFromFile(filename)
                     }
                     button$SetImage(image)
                   } else if(is(value,"gLabelRGtk")) {
                     button$SetLabel(value@widget)
                   } else {
                     button$SetLabel(value)
                   }
                   return(obj)
                 })


## font -- push down to label
setReplaceMethod("font",signature(obj="gButtonRGtk"),
          function(obj, ..., value) {
            .font(obj, obj@toolkit,...) <- .fixFontMessUp(value)
            return(obj)
          })

setReplaceMethod(".font",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gButtonRGtk"),
                 function(obj, toolkit, ..., value) {
                   widget <- getWidget(obj)[[1]] # label is first child or something
                   if(is(widget, "GtkAlignment"))
                     widget <- widget[[1]][[2]] # a real hacke
                   .font(widget, toolkit, ...) <- value
                   invisible(obj)
                 })
### handlers
setMethod(".addhandlerchanged",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gButtonRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            addhandlerclicked(obj, handler, action,...)
          })

setMethod(".addhandlerclicked",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gButtonRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            addhandler(obj,"clicked",handler, action,...)
          })

## for popup menu
setMethod(".addpopupmenu",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gButtonRGtk"),
          function(obj, toolkit, menulist, action=NULL, ...) {
            addPopupMenuWithSignal(obj, toolkit, menulist, signal="clicked",...)
})
