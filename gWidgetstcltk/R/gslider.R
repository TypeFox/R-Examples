## FIX up for non-integer values

setClass("gSlidertcltk",
         contains="gComponenttcltk",
         prototype=prototype(new("gComponenttcltk"))
         )

setMethod(".gslider",
          signature(toolkit="guiWidgetsToolkittcltk"),
          function(toolkit,
                   from=0, to=100, by = 1,
                   value=from,
                   horizontal=TRUE,
                   handler=NULL, action=NULL,
                   container=NULL, ...) {
            force(toolkit)

            ## if from a single value, then from, to ,by specify sequence
            if(length(from) == 1)
              x <- seq(from, to, by)
            else
              x <- from
            
            ## x needs sorting, make unique
            x <- sort(unique(x))                  # do I need to do for different types
            ind <- seq_along(x)
            value <- which(as.character(value) == as.character(x))
            
  
            if(is(container,"logical") && container)
              container = gwindow()
            if(!is(container,"guiWidget")) {
              warning("Container is not correct. No NULL containers possible\n" )
              return()
            }
            
  
            if(horizontal)
              orientation <- "horizontal"
            else
              orientation <- "vertical"
            
            tt <- getWidget(container)
            SliderValue <- tclVar(as.character(value))
            
            ## ## use old school. ttk:::scale doesn't allow steps, using other values.
            ## slider <- tkscale(tt, from=1L, to=length(x),
            ##                   showvalue=FALSE, variable=SliderValue,
            ##                   resolution=1L, orient=orientation)

            slider <- tkwidget(tt, "ttk::scale", from=1L, to=length(x), variable=SliderValue,
                               orient=orientation)
            
            obj <- new("gSlidertcltk",block=slider, widget=slider,
                       toolkit=toolkit, ID=getNewID(), e = new.env())
            tag(obj,"..tclVar") <- SliderValue
            tag(obj, "..byIndexValues") <- x

            ## ## modify label
            ## modifyLabel <- function() {
            ##   tkconfigure(slider, label=format(svalue(obj), digts=3))
            ## }
            ## modifyLabel()
            ## tkbind(slider, "<Motion>", modifyLabel)
            
            
            add(container, obj,...)
            
            if (!is.null(handler))  {
              id <- addhandlerchanged(obj, handler, action)
            }
            
            return(obj)
          })


### methods
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gSlidertcltk"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {
            rbValue = tag(obj,"..tclVar")
            val <- as.numeric(tclvalue(rbValue))
            if(is.null(index) || !index) {
              x <- tag(obj, "..byIndexValues")
              val <- x[val]
            }
            return(val)
          })

setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gSlidertcltk"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   ## can set by index or match
                   if(is.null(index) || index==FALSE) {
                     value <- match(value, tag(obj, "..byIndexValues"))
                   } else {
                     value <- value
                   }
                                         
                   n <- length(tag(obj, "..byIndexValues"))
                   if(!is.na(value) &&
                      value >= 1 &&
                      value <= n) 
                     tclvalue(tag(obj,"..tclVar")) <- value
                   ## ## update label
                   ## tkconfigure(getWidget(obj), label=format(svalue(obj), digts=3))
                   return(obj)
               })

##' return values
setMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkittcltk",x="gSlidertcltk"),
          function(x, toolkit, i, j, ..., drop=TRUE) {
            tag(x, "..byIndexValues")
          })

## Method to replace values of spin button
setReplaceMethod("[",
                 signature(x="gSlidertcltk"),
                 function(x, i, j,..., value) {
                   .leftBracket(x, x@toolkit, i, j, ...) <- value
                   return(x)
                 })



## Method to replace values of spin button
setReplaceMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkittcltk",x="gSlidertcltk"),
          function(x, toolkit, i, j, ..., value) {
            obj <- x
            widget <- getWidget(obj)
            curVal <- svalue(obj)

            value <- sort(unique(value))
            tag(obj, "..byIndexValues") <- value
            tkconfigure(widget, from=1, to=length(value))

            svalue(obj) <- curVal
            
            return(obj)
          })





### handlers
setMethod(".addhandlerchanged",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gSlidertcltk"),
          function(obj, toolkit, handler, action=NULL, ...) {
#            .addHandler(obj,toolkit, signal="<ButtonRelease-1>",handler,action)
            .addHandler(obj,toolkit, signal="command",handler,action)
          })




