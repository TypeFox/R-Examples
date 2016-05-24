## class defined in aaaClasses for inheritance
library(tcltk)

     

## constructor
setClass("gEdittcltk",
         representation = representation("gComponentR5tcltk",
           coercewith="NULLorFunction"),
         contains="gComponenttcltk",
         prototype=prototype(new("gComponenttcltk"))
         )


setMethod(".gedit",
          signature(toolkit="guiWidgetsToolkittcltk"),
          function(toolkit,
                   text="", width=25,
                   coerce.with = NULL,
                   initial.msg = "", 
                   handler=NULL, action=NULL,
                   container=NULL,
                   ...
                   ) {

           force(toolkit)
            
            if(is(container,"logical") && container)
              container = gwindow()
            if(!is(container,"guiWidget")) {
              warning("Container is not correct. No NULL containers possible\n" )
              return()
            }


           if (is.null(text)) text<-""

            ## check that coerce.with is a function
            if(is.null(coerce.with) || is.function(coerce.with)) {
              ## okay
            } else {
              if(is.character(coerce.with)) {
                coerce.with = get(coerce.with)
              }
            }

           tt <- getWidget(container)

           e <- getRefClass("Entry")$new(tt)
           obj <- new("gEdittcltk", block=e$get_widget(), widget = e$get_widget(),
                      R5widget=e,
                      toolkit=toolkit,ID=getNewID(), e = new.env(),
                      coercewith=coerce.with)

           if(nchar(text)) {
            
             svalue(obj) <- text
           }

           
           ## initial message
           if(nchar(initial.msg) > 0 && nchar(text) == 0) {
             e$set_init_msg(initial.msg)
             e$show_init_msg()
           }

           ## width
           if(!is.null(width))
             tkconfigure(obj@widget,width=as.integer(width)) ## character count, not pixels
           
           ## Drag and drop
           ## addDropSource(obj)
           ## addDropTarget(obj)
           
           add(container, obj,...)
           
           if (!is.null(handler)) 
             tag(obj, "handler.id") <- addhandlerchanged(obj,handler,action)
           
           
           invisible(obj)
            
            
          })

## methods
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gEdittcltk"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {
            widget <- obj@R5widget
            val <- widget$get_value()
            ## val = tclvalue(tag(obj,"tclVar"))
            if(val == "<NA>")
              val <- NA
            
            coercewith <- obj@coercewith
            
            if(is.null(coercewith))
              return(val)

            if(is.character(coercewith))
              coercewith <- get(coercewith)
            
            if(!is.function(coercewith))
              return(val)

            
            return(coercewith(val))
          })

## svalue<-
setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gEdittcltk"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   if(is.na(value))
                     value <- "<NA>"

                   widget <- obj@R5widget
                   widget$set_value(value)
                   ## tclvalue(tag(obj, "tclVar")) <- value
                   return(obj)
          })


## left bracket implement completion
setMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkittcltk",x="gEdittcltk"),
          function(x, toolkit, i, j, ..., drop=TRUE) {
            widget <- x@R5widget
            vals <- widget$words
            if(missing(i))
              vals
            else
              vals[i]
          })
            
setMethod("[",
          signature(x="gEdittcltk"),
          function(x, i, j, ..., drop=TRUE) {
            if(missing(i))
              .leftBracket(x,x@toolkit, ...)
            else
              .leftBracket(x,x@toolkit, i, ...)
          })

setReplaceMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkittcltk",x="gEdittcltk"),
          function(x, toolkit, i, j, ..., value) {
            widget <- x@R5widget
            vals <- widget$words
            # vals <- tag(x, "typeAhead")
            
            if(missing(i))
              vals <- value
            else
              vals[i] <- value
            widget$set_words(vals)
            ## tag(x, "typeAhead") <- vals
            return(x)
          })

setReplaceMethod("[",
                 signature(x="gEdittcltk"),
                 function(x, i, j,..., value) {
                   .leftBracket(x, x@toolkit, i, j, ...) <- value
                   return(x)
                 })


setReplaceMethod(".size", 
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gEdittcltk"),
                 function(obj, toolkit, ..., value) {
                   if(is.numeric(value))
                     tkconfigure(obj@widget,width=ceiling(value[1]/widthOfChar)) 
                                        
                   else
                     message(gettext("size needs a numeric vector c(width,...)\n"))
                   return(obj)
                 })


##' visible<- if FALSE, for password usage
setReplaceMethod(".visible",signature(toolkit="guiWidgetsToolkittcltk", obj="gEdittcltk"),
          function(obj, toolkit, ..., value) {
            widget <- getWidget(obj)
            if(as.logical(value))
              tkconfigure(widget, show="")
            else
              tkconfigure(widget, show="*")
            return(obj)
          })


##################################################
## handlers

## changed is called after a commit (svalue, Return key in widget -- not drop down menu)
## keystroke is called when widget display changes

## Use Virtual Event for KeyRelease, as other one is used by class above

## use "R5classes" handlers here
setMethod(".addhandlerchanged",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gEdittcltk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            widget <- obj@R5widget
            FUN <- function() {
              h <- list(obj=obj, action=action)
              handler(h)
            }
            widget$add_handler("<<Changed>>", FUN)
##            .addHandler(widget$get_widget(), toolkit, signal="<<Changed>>", handler, action, actualobj=obj)
          })

setMethod(".addhandlerkeystroke",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gEdittcltk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            FUN = function(d) {
              h = list(obj = obj, action = action, key=d) # add in key
              handler(h)
            }
            widget <- obj@R5widget            

            widget$add_handler("<<KeyRelease>>", FUN)
            ##            .addHandler(widget$get_widget(), toolkit, signal="<<KeyRelease>>", FUN, action, actualobj=obj)
          })

## This just dispatches down to R5 classes. It should be put into a class

