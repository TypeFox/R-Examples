setClass("gFrameRGtk",
         contains="gGroupRGtk",
         prototype=prototype(new("gGroupRGtk"))
         )

## add a frame for packing. subclass of gGroup
setMethod(".gframe",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   text = "", markup=FALSE,
                   pos = 0, ## pos in [0,1] 0 for left, 1 for right
                   horizontal=TRUE,
                   container=NULL,
                   ...) {
            
            force(toolkit)

            frame = gtkFrameNew()


            obj <- as.gWidgetsRGtk2(frame, horizontal=horizontal,...)

            ##            group = ggroup(horizontal=horizontal, ...) # for horizontal, spacing etc.
##             frame$Add(getBlock(group))

##             ## add label to group
##             obj = new("gFrameRGtk",
##               block=frame, widget=group@widget, toolkit=toolkit)

            tag(obj,"markup") <- markup
            names(obj) <- text

            frame$SetLabelAlign(pos,0.5) # was 0, suggested value by felix andrews


            if (!is.null(container)) {
              if(is.logical(container) && container == TRUE)
                container = gwindow(visible=TRUE)
              add(container, obj, ...)
            }
            return(obj)
          })

as.gWidgetsRGtk2.GtkFrame <- function(widget,...) {

  if(is.null(tag(widget,"group"))) {
    theArgs <- list(...)
    horizontal <- if(is.null(theArgs$horizontal)) TRUE else theArgs$horizontal
    spacing <- if(is.null(theArgs$spacing)) 5 else theArgs$spacing
    group <- ggroup(horizontal=horizontal, spacing=spacing) # for horizontal, spacing etc.
    widget$Add(getBlock(group))
  } else {
    group <- tag(widget,"group")
  }

  ## add label to group
  obj <- new("gFrameRGtk",block=widget, widget=getWidget(group),
    toolkit=guiToolkit("RGtk2"))

  tag(obj,"group") <- group
  if(is.null(tag(obj,"markup"))) tag(obj,"markup") <- FALSE

  return(obj)
}
    

### methods -- inherited from ggroup

## return label name
setMethod(".names",signature(toolkit="guiWidgetsToolkitRGtk2",
                             x="gFrameRGtk"),
          function(x, toolkit) {
            if(tag(x,"markup")) {
              label <- getBlock(x)$GetLabelWidget()
              label$GetLabel()
            } else {
              getBlock(x)$GetLabel()
            }
          })


setReplaceMethod(".names",
                 signature(toolkit="guiWidgetsToolkitRGtk2",x = "gFrameRGtk"),
                 function(x,toolkit,value) {

                   frame <- getBlock(x)
                   if(is.null(tag(x,"markup")) || !tag(x,"markup")) {
                     frame$SetLabel(value)
                   } else {
                     label <- gtkLabelNew()
                     label$SetMarkup(value)
                     frame$SetLabelWidget(label)
                   }

                   return(x)
                 })

