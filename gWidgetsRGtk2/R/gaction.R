## reusuabel chunk of code
setClass("gActionRGtk",
         contains="gComponentRGtk",
         representation(e = "environment"),
         prototype=prototype(e=new.env())
         )


setMethod(".gaction",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   label,
                   tooltip = NULL,
                   icon = NULL,
                   key.accel = NULL,
                   handler = NULL, action = NULL,
                   parent = NULL,
                   ...) {
            
            force(toolkit)

            if(!is.null(icon))
              icon <- getstockiconname(icon)
            
            act <- gtkAction(name = make.names(label),
                             label = label,
                             tooltip = tooltip,
                             stock.id = icon)


            obj = new("gActionRGtk", block=act, widget=act, toolkit=toolkit)

            ## add for later use
            ## should be defined when used in a menu bar.
            tag(obj,"key.accel") <- key.accel
            obj@e$buttons <- list()     # for svalue<- with buttons, menu items work


            ## accel buttons
            if(!is.null(key.accel) && !is.null(parent)) {
              toplevel <- getBlock(parent)$toplevel
              ## mask Shift-1, Control-4 alt-8
              ## key sprintf("GDK_%s",key)
              ## flag GtkAccelFlags -- 1
              if(grepl("^Control", key.accel) ||
                 grepl("^Alt", key.accel) ||
                 grepl("^Shift", key.accel)) {
                tmp <- strsplit(key.accel, "-")[[1]]
                modifier <- c(Shift="shift-mask", "Control"="control-mask", Alt="mod1-mask")[tmp[1]]
                key <- sprintf("GDK_%s", tmp[2])
              } else {
                modifier <- "modifier-mask"
                key <- sprintf("GDK_%s", key.accel)
              }
              a <- gtkAccelGroup()
              toplevel$addAccelGroup(a)
              a$connect(get(key), modifier, "visible", function(...) {
                h <- list(action=action)
                handler(h, ...)
                TRUE
              })
            }

            
            if(!is.null(handler))
              addHandlerChanged(obj, handler, action)
            
            return(obj)
          })

## svalue -- get label
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gActionRGtk"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {
            widget <- getWidget(obj)
            return(widget['label'])
          })



## svalue<- set label
setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gActionRGtk"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   gtkaction <- getWidget(obj)

                   ## for menu, toolbar est label propoerty
                   gtkaction['label'] <- value

                   ## for buttons, we work harder
                   buttons <- obj@e$buttons
                   if(length(buttons) > 0)
                     sapply(buttons, function(i) {
                       if(isExtant(i))
                         svalue(i) <- value
                     })

                   return(obj)
                 })

## enabled -- inherited
setMethod(".addhandlerchanged",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gActionRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            widget <- getWidget(obj)

            ID <- gSignalConnect(widget, signal="activate",
                           f = handler,
                           data = list(action = action),
                           user.data.first = TRUE)

            invisible(ID)
          })

                             
## helper functions

.isgAction <- function(lst) {
  is(lst,"guiComponent") && is(lst@widget, "gActionRGtk") ||
  is(lst,"gActionRGtk")
}
