## reusuabel chunk of code
setClass("gActiontcltk",
         representation(widget="list",e = "environment"),
         prototype(widget=list(), e = new.env())
         )

setMethod(".tag", signature(toolkit="guiWidgetsToolkittcltk",obj="gActiontcltk"),
          function(obj, toolkit, i, drop=TRUE, ...) {
            if(missing(i)) i = NULL
            if(missing(drop)) drop <- TRUE                                    

            if(is.null(i))
              return(as.list(obj@e))
            else
              return(obj@e[[i]])
            
          })

setReplaceMethod(".tag", signature(toolkit="guiWidgetsToolkittcltk",obj="gActiontcltk"),
          function(obj, toolkit, i, replace=TRUE, ..., value) {
            if(missing(i)) i = NULL
            

            obj@e[[i]] <- value
            return(obj)

          })

setMethod(".gaction",
          signature(toolkit="guiWidgetsToolkittcltk"),
          function(toolkit,
                   label,
                   tooltip = NULL,
                   icon = NULL,
                   key.accel = NULL,
                   handler = NULL, action = NULL,
                   parent=NULL,
                   ...) {
            
            force(toolkit)

            lst <- list(label = label,
                        tooltip = tooltip,
                        icon = icon,
                        key.accel = key.accel,
                        handler = handler,
                        action = action)
            e <- new.env(); e$state <- TRUE; e$buttons <- e$menuitems <- e$toolbaritems <- list()
            e$label <- label
            obj <- new("gActiontcltk", widget = lst, e =e)

            if(!is.null(key.accel) && !is.null(parent)) {
              toplevel <- tkwinfo("toplevel", getWidget(parent))
              tkbind(toplevel, sprintf("<%s>",key.accel), function() {
                if(obj@e$state) {
                  h <- list(action=action)
                  handler(h)
                }
              })
            }
              
            return(obj)
          })

setMethod(".getToolkitWidget",
          signature(obj="gActiontcltk", toolkit="guiWidgetsToolkittcltk"),
          function(obj, toolkit) obj@widget)

## is this a gaction
.isgAction <- function(obj) {
  is(obj,"guiComponent") && is(obj@widget,"gActiontcltk")
}

## methods need to be disabled
setReplaceMethod(".enabled",
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gActiontcltk"),
                 function(obj, toolkit, ..., value) {
                   e <- obj@e
                   e$state <- as.logical(value)
                   
                   if(length(e$buttons) > 0)
                     lapply(e$buttons, function(i) enabled(i) <- as.logical(value))

                   if(length(e$toolbaritems) > 0)
                     lapply(e$toolbaritems, function(i) {
                     if(as.logical(value))
                       tkconfigure(i,state="normal")
                     else
                       tkconfigure(i, state = "disabled")
                   })
                   if(length(e$menuitems) > 0)
                     lapply(e$menuitems, function(i) {
                       if(as.logical(value))
                         tcl(i,"entryconfigure",e$label,state="normal")
                       else
                         tcl(i,"entryconfigure",e$label,state="disabled")
                     })
                   return(obj)
                 })


setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gActiontcltk"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {
            val <- obj@widget$label
            return(val)
          })

setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gActiontcltk"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   e <- obj@e
                   if(length(e$buttons) > 0)
                     lapply(e$buttons, function(i) svalue(i) <- as.character(value))

                   if(length(e$toolbaritems) > 0)
                     lapply(e$toolbaritems, function(i) {
                       tkconfigure(i, text=value)
                     })
                   if(length(e$menuitems) > 0)
                     lapply(e$menuitems, function(i) {
                       tcl(i,"entryconfigure",e$label,label=value)
                     })
                   return(obj)
                 })
