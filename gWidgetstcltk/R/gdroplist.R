## editable has entry widget that can be edited
setClass("gDroplisttcltk",
         contains="gComponenttcltk",
         prototype=prototype(new("gComponenttcltk"))
         )

setMethod(".gdroplist",
          signature(toolkit="guiWidgetsToolkittcltk"),
          function(toolkit,
                   items, selected = 1, # use 0 for blank
                   editable=FALSE,
                   coerce.with = NULL,
                   handler=NULL, action=NULL,
                   container=NULL,
                   ...               # do.quote = TRUE for quote of answer
                   ) {

            force(toolkit)

            if(is(container,"logical") && container)
              container = gwindow()
            if(!is(container,"guiWidget")) {
              warning("Container is not correct. No NULL containers possible\n" )
              return()
            }

            ## items can be vector of items or data frame with
            ## one col: items
            ## two cols: items, icons
            ## three cols: items, icons, tooltip
            ## four or more cols: toolkit specific

            if(inherits(items,"data.frame")) {
              items <- items[,1, drop=TRUE]
            }
            ## no icons, tooltip in tcltk
            
            
            ## items must be a vector here
            items = as.vector(items)              # undoes factor
            items = unique(items)                 # unique
            
            theArgs = list(...)
            ## keep this, but don't advertise
            if(!is.null(theArgs$do.quote)) {
              coerce.with = function(x) paste("'",x,"'",sep="") # no space
            }
            
            if(editable)
              state <- "normal"
            else
              state <- "readonly"


            if(!is.null(theArgs$width))
              width <- theArgs$width
            else if(length(items))
              width <- max(sapply(items,nchar))  + 5
            else
              width <- NULL
            
            tt <- getWidget(container)
            gp <- ttkframe(tt)
            cbVar <- tclVar()
            if(length(items) == 1)
              values <- as.tclObj(as.character(items))
            else
              values <- as.character(items)
            cb <- ttkcombobox(gp,
                              values = values,
                              textvariable = cbVar,
                              state = state)
            if(!is.null(width))
              tkconfigure(cb, width=width)

            tkgrid(cb,row=0, column=0, sticky="we") # stretch horizontally. Use news for both
            tkgrid.columnconfigure(gp,0, weight=1)
            
            obj = new("gDroplisttcltk",block=gp,widget=cb,
              toolkit=toolkit,ID=getNewID(), e = new.env())

            tag(obj,"coerce.with") <- coerce.with
            tag(obj,"editable") <- editable
            tag(obj,"tclVar") <- cbVar
            tag(obj,"items") <- items

            addDropTarget(obj, handler = function(h,...)
                           svalue(obj) <- h$dropdata)

            
            add(container, obj, ...)

            if(!is.null(theArgs$width))
              size(obj) <- c(theArgs$width,0)
            

            ## can add numeric or for Richie, a value
            if(is.numeric(selected)) { 
              svalue(obj, index = TRUE) <- selected 
            } else { 
              svalue(obj) <- as.character(selected) 
            } 
            

            
            if (!is.null(handler)) {
              id <- addhandlerchanged(obj, handler, action)
              tag(obj, "handler.id") <- id
            }
            
            invisible(obj)
          })
          
### methods
## value is for getting/setting the selected value
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gDroplisttcltk"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {
            widget <- getWidget(obj)

            ind <-  as.numeric(tclvalue(tcl(widget, "current"))) + 1 # 0-based

                        
            ## if index
            if(!is.null(index) && index) {
              return(ind)
            }


            editable <- as.character(tkcget(widget, "-state")) != "readonly"
            if(editable) {
              val <- tclvalue(tcl(widget,"get"))
            } else {
              if(ind == 0) {
                ## no selection
                return(NA)
              }

              ## else get values from items -- not get to avoid conversion
              items <- tag(obj,"items")
              val <- items[ind]
            }

            
            ## add in an as.numeric flag, getwidget when editable
            theArgs = list(...)         # deprecated

            coerce.with = tag(obj, "coerce.with")

            ## do we coerce return value?
            if(is.null(coerce.with))
              return(val)
            else if(is.function(coerce.with))
              return(coerce.with(val))
            else if(is.character(coerce.with))
              return(do.call(coerce.with,list(val)))
            else
              return(val)               # what else?
            
          })

## set the displayed value to value
## if index=TRUE and value=0, seet to no state
setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gDroplisttcltk"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   ## we can only handle vectors for value -- not data frame
                   ## with value, label, icon info


                   theArgs = list(...)

                   widget <- getWidget(obj)
                   tclVar <- tag(obj, "tclVar")
                   
                   n <- length(obj)
                   if(n < 1) return(obj)
                   
                   if(is.null(index))
                     index <- FALSE
                   index <- as.logical(index)

                   ##  if editable do differently
                   ## editable not implented
                   editable <- as.character(tkcget(widget, "-state")) != "readonly"
                   ## if index, set
                   if(index) {
                     if(value > 0 && value <= n)
                       tcl(widget,"current", as.numeric(value) - 1)
                     else               # set to no state
                       tclvalue(tclVar) = ""
                     ##tcl(widget,"set", "") # aka -1 for get
                   } else {
                     if(!is.null(editable) && editable) {
                       ## editable
                       tclvalue(tcl(widget,"set",as.character(value)))
                     } else {
                       ## not editable, check its there
                       vals <- tag(obj,"items")
                       if(value %in% vals) {
                         tclvalue(tcl(widget,"set",as.character(value)))
                       } else {
                         message(sprintf("%s is not a valid item",value),"\n")
                       }
                     }
                   }
                   ## notify event handlers unless set to 0
                   if(!(index && value < 0))
                     tkevent.generate(getWidget(obj), "<<ComboboxSelected>>")
                   
                   return(obj)
                 })




setMethod("length",
          signature(x="gDroplisttcltk"),
          function(x) {
            .length(x, x@toolkit)
          })
setMethod(".length",
          signature(toolkit="guiWidgetsToolkittcltk",x="gDroplisttcltk"),
          function(x, toolkit) {
            return(length(tag(x,"items")))
          })


## the methods [ and [<- refer to the pre-defined values in the drop list.
## [
setMethod("[",
          signature(x="gDroplisttcltk"),
          function(x, i, j, ..., drop=TRUE) {
            .leftBracket(x, x@toolkit, i, j, ..., drop=drop)
          })
setMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkittcltk",x="gDroplisttcltk"),
          function(x, toolkit, i, j, ..., drop=TRUE) {

            n = length(x)               # no. items
            if(n == 0)
              character(0)              # Thanks Yves
            
            items = tag(x,"items")
            
            if(missing(i))
              return(items)
            else
              return(items[i])
          })


## replaces the values in droplist
## values is a vector of values -- not a dataframe
#set.values.gDropList = function(obj, values, ...) {
setReplaceMethod("[",
                 signature(x="gDroplisttcltk"),
                 function(x, i, j,..., value) {
                   .leftBracket(x, x@toolkit, i, j, ...) <- value
                   return(x)
                 })

setReplaceMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkittcltk",x="gDroplisttcltk"),
          function(x, toolkit, i, j, ..., value) {
            if(is.data.frame(value))
              value <- value[,1,drop=TRUE]

            widget <- getWidget(x)
            ind <- svalue(x, index=TRUE)
            
            if(missing(i)) {
              if(length(value) == 1)
                 tcl(widget,"configure",values=as.tclObj(value))
              else
                tcl(widget,"configure",values=value)
              tag(x,"items") <- value
              if(ind > 0)
                svalue(x, index=TRUE) <- ind
            } else {
              items = x[]
              items[i] <- value
              x[] <- items ## recurse
            }

            return(x)
          })



###################################################
  
### handlers
setMethod(".addhandlerchanged",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gDroplisttcltk"),
          function(obj, toolkit, handler, action=NULL, ...) {

            .addHandler(obj,toolkit,"<<ComboboxSelected>>",handler,action,...)

            editable <- as.character(tkcget(getWidget(obj), "-state")) != "readonly"
            if(editable) ## tag(obj,"editable"))
              .addHandler(obj, toolkit, signal="<Return>", handler, action)
            
          })

setMethod(".addhandlerclicked",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gDroplisttcltk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            .addhandlerchanged(obj,toolkit, handler,action)
          })

