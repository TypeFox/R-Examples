## cairo graphics device
## would like to get size from par("fin"), but this isn't so easy as it
## seems to pop up a new plot container

### Trouble when adding to a notebook. Currently when a notebook page is closed the signal to close the widget is not propogated.


setClass("gGraphicstcltk",
         contains="gComponenttcltk",
         prototype=prototype(new("gComponenttcltk"))
         )

setMethod(".ggraphics",
          signature(toolkit="guiWidgetsToolkittcltk"),
          function(toolkit,
                   width=dpi*6, height=dpi*6,
                   dpi=75, ps=12,
                   container=NULL,...) {

            force(toolkit)

            msg <- paste("There is no embeddable graphics device available for",
                         "gWidgetstcltk. However, the device created by the tkrplot",
                         "package may be embedded into gWidgets by creating a group container",
                         "with ggroup, and using the result of getToolkitWidget(group_container)",
                         "as the parent for tkrplot.",
                         sep="\n")

            ## take@widget to get glabel instance after going through gWidgets
            out <- glabel(msg, container=container)@widget
            return(out)
          })


### methods

## ## adding to a group is funny, we intercept here
## setMethod(".add",
##           signature(toolkit="guiWidgetsToolkittcltk",obj="gContainertcltk", value="gGraphicstcltk"),
##           function(obj, toolkit, value, ...) {
##             cat("can't add a ggraphics() object to a container in gWidgetsrjava")
##           })



## raise this device
setReplaceMethod(".visible",
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gGraphicstcltk"),
                 function(obj, toolkit, ..., value) {
                   if(is.logical(value) == TRUE) {
                     dev.set(tag(obj,"device"))
                   }
                   return(obj)
                 })

## save Current Page
## This uses GTK -- not R to save.
## need to have window fully shown
setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gGraphicstcltk"),
                 function(obj, toolkit, index=NULL,  ..., value) {
                   gwCat("svalue not implemented\n")
                   return(obj)
                 })


### handlers
## add this expose event for graph
setMethod(".addhandlerexpose",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gGraphicstcltk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            addhandler(obj,"expose-event",handler,action)
          })

## applies a handler to the mouse click. The handler gets extra
## argument h$x, h$y passed into it. These are in [0,1] coordinates
setMethod(".addhandlerclicked",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gGraphicstcltk"),
          function(obj, toolkit, handler, action=NULL, ...) {
          })

