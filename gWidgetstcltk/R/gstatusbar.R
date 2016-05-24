## StatusBar. Use value to push message, value to pop
setClass("gStatusbartcltk",
         contains="gComponenttcltk",
         prototype=prototype(new("gComponenttcltk"))
         )
## constructor
setMethod(".gstatusbar",
          signature(toolkit="guiWidgetsToolkittcltk"),
          function(toolkit,
                   text="", container=NULL, ...) {

            
            force(toolkit)            


            if(is(container,"logical") && container)
              container = gwindow()

            ## container must be a gwindow unless we pass in argument not.toplevel=TRUE
            theArgs <- list(...)
            if(!is.null(theArgs$not.toplevel) && as.logical(theArgs$not.toplevel)) {
              tt <- getBlock(container)
            } else {
              if(!(is(container,"gWindowtcltk") || is(container@widget,"gWindowtcltk"))) {
                message(gettext("gstatusbar: container must be gwindow instance\n"))
              }
              tt <- tag(container,"sb")
            }
            gp <- ttkframe(tt)
            
            sb <- ttklabel(gp, text=text)
            tkpack(sb, side="left",anchor="w", expand=TRUE, fill="x")
            
            obj = new("gStatusbartcltk",block=gp, widget=sb,
              toolkit=toolkit, ID=getNewID(), e = new.env())

            stack <- c(text)
            tag(obj,"stack") <- stack
            
            ## add to container
            add(container, obj,...)

            
            invisible(obj)
          })

### methods

## This pops label
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gStatusbartcltk"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {
            ## pop the stack
            stack <- tag(obj,"stack")
            val <- stack[1]
            if(length(stack))
              stack <- stack[-1]
            tag(obj,"stack") <- stack

            if(length(stack))
              value <- stack[1]
            else
              value <- ""
            tkconfigure(obj@widget, text=as.character(value))

            return(val)
          })


## This pushes to label
setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gStatusbartcltk"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   stack <- tag(obj,"stack")
                   stack <- c(value, stack)
                   tag(obj,"stack") <- stack
                   
                   tkconfigure(obj@widget, text=as.character(value))
                   return(obj)
                 })

