setClass("gCheckboxtcltk",
         contains="gComponenttcltk",
         prototype=prototype(new("gComponenttcltk"))
         )

## constructor
setMethod(".gcheckbox",
          signature(toolkit="guiWidgetsToolkittcltk"),
          function(toolkit,
                   text, checked=FALSE,
                   use.togglebutton=FALSE,
                   handler=NULL, action=NULL,
                   container=NULL,...) {
            
            force(toolkit)
            
            if(missing(text)) text = ""

            if(is(container,"logical") && container)
              container <- gwindow()
            if(!is(container,"guiWidget")) {
              warning("Container is not correct. No NULL containers possible\n" )
              return()
            }
            
            theArgs <- list(...)


            tt <- getWidget(container)
 #           gp = ttkframe(tt)


            
            ## widget use toolbutton or not?
            ## http://wiki.tcl.tk/17899
            if(use.togglebutton) {
              check <- ttkcheckbutton(tt, text=as.character(text), style="Toolbutton")
            } else {
              check <- ttkcheckbutton(tt, text=as.character(text))
            }
#            theLabel = ttklabel(gp, text=text)
            ## configure
            tclVar <- tclVar(as.numeric(checked))
            tkconfigure(check,variable=tclVar)

            obj <- new("gCheckboxtcltk",block=check, widget=check,
              toolkit=toolkit, ID=getNewID(), e = new.env())

            tag(obj,"tclVar") <- tclVar
            
            ## add to container
            add(container, obj,...)
            
            if (!is.null(handler)) 
              obj@e$handlerID <- addhandlerchanged(obj, handler, action=action)
  
            invisible(obj)
          })

### methods
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gCheckboxtcltk"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {
            cbVal <- as.logical(as.numeric(tclvalue(tag(obj,"tclVar"))))
            return(cbVal)
          })

setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gCheckboxtcltk"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   tclvalue(tag(obj,"tclVar")) <-
                     as.character(as.numeric(value))
                   return(obj)
                 })

## [
setMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkittcltk",x="gCheckboxtcltk"),
          function(x, toolkit, i, j, ..., drop=TRUE) {
            ## theLabel <- tag(x,"labelText")
            widget <- getWidget(x)
            val <- tclvalue(tkcget(widget, "-text"))
            return(val)
          })
            
setMethod("[",
          signature(x="gCheckboxtcltk"),
          function(x, i, j, ..., drop=TRUE) {
            .leftBracket(x, x@toolkit, i, j, ..., drop=drop)
          })

setReplaceMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkittcltk",x="gCheckboxtcltk"),
          function(x, toolkit, i, j, ..., value) {
            widget <- getWidget(x)
            tkconfigure(widget, text=paste(value, collapse="\n"))
            return(x)
          })

setReplaceMethod("[",
                 signature(x="gCheckboxtcltk"),
                 function(x, i, j,..., value) {
                   .leftBracket(x, x@toolkit, i, j, ...) <- value
                   return(x)
                 })

## inherited enabled isn't workgin                
setReplaceMethod(".enabled",
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gCheckboxtcltk"),
                 function(obj, toolkit, ..., value) {

                   ## Odd this isn't needed anymore
                   widget <- getWidget(obj)
                   if(as.logical(value))
                     tcl(widget,"state","!disabled")
                   else
                     tcl(widget,"state","disabled")
                   return(obj)
                   
                 ##   ## change both widget and label
                 ##   lapply(list(tag(obj,"check"), tag(obj,"label")), function(i) {
                 ##   if(as.logical(value))
                 ##     tcl(i,"state","!disabled")
                 ##   else
                 ##     tcl(i,"state","disabled")
                 ## })
                   
                 ##   return(obj)
                 })

### no method to change the value of text???

### handlers
setMethod(".addhandlerchanged",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gCheckboxtcltk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            changeHandler <- handler

            theArgs <- list(...); actualobj <- theArgs$actualobj
            if(is.null(actualobj))
              actualobj <- obj

            ## bind to command, not ButtonRelease-1. That binding requires a pause
            addhandler(obj,toolkit, signal="command",
                        action=action, actualobj=actualobj,
                        handler = function(h,...) {
                              changeHandler(h,...)
                          })
          })

setMethod(".addhandlerclicked",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gCheckboxtcltk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            .addhandlerchanged(obj, toolkit, handler, action)
          })
