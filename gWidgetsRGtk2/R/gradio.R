## Use a reference class


##################################################
## Radio widget stuff

RadioWidgetGtk <- setRefClass("RadioWidgetGtk",
                     contains="GWidgetGtk",
                     fields=list(
                       inner_block="ANY",    # replaceble box container
                       items="ANY",          # store the items
                       horizontal="logical", # layout direction
                       obj = "ANY"           # gradio object for callbacks
                       ),
                     methods=list(
                       initialize=function(items, selected=1, horizontal=TRUE) {
                         horizontal <<- horizontal
                         block <<- gtkHBox()
                         inner_block <<- gtkHBox(); block$packStart(inner_block)
                         if(!missing(items)) {
                           set_items(items, selected)
                         }
                         .self
                       },
                       get_items = function() {
                         "Return items"
                         items
                       },
                       set_items=function(items, selected=NULL) {
                         if(length(items) == 0) return()
                         if(is.null(selected))
                           selected <- get_index()
                         items <<- items

                         block$remove(inner_block)
                         inner_block <<- if(horizontal) gtkHBox() else gtkVBox()
                         block$packStart(inner_block)

                         widget <<- gtkRadioButton(label=items[1])
                         ## Keep rbs around until after sapply statement
                         rbs <- lapply(items[-1], gtkRadioButtonNewWithLabelFromWidget, group = widget)
                         sapply(rev(widget$getGroup()), gtkBoxPackStart, object = inner_block)

                         ## add handlers
                         lapply(widget$getGroup(), gSignalConnect, signal="toggled", f=function(self, w, ...) {
                           if(w$getActive())
                             self$notify_observers(...)
                         }, data=.self, user.data.first=TRUE)

                         set_index(selected)
                         
                         invisible()
                       },
                       get_index = function() {
                         "Return index of selected"
                         which(sapply(rev(widget$getGroup()), gtkToggleButtonGetActive))
                       },
                       set_index = function(i) {
                         "Set index of selection"
                         i <- as.integer(i)
                         l <- rev(widget$getGroup())
                         if(1 <= i && i <= length(l))
                           l[[i]]$setActive(TRUE)

                         invisible()
                       }
                       ))


setClass("gRadioRGtk",
         contains="gComponentWithRefClassWithItemsRGtk",
         prototype=prototype(new("gComponentWithRefClassWithItemsRGtk"))
         )






## constructor
setMethod(".gradio",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   items, selected=1, horizontal=FALSE,
                   handler=NULL, action=NULL,
                   container=NULL,       
                   ...
                   ) {

            force(toolkit)

            ref_widget <- RadioWidgetGtk$new(items, selected, horizontal)
            obj <- new("gRadioRGtk",block=ref_widget$block, widget=ref_widget$block,
                       ref_widget=ref_widget,
                       toolkit=guiToolkit("RGtk2"))


            if(is.data.frame(items))
              items <- items[,1, drop=TRUE] # first column

            obj[] <- items
            svalue(obj,index=TRUE) <- selected

            ## do we add to the container?
            if (!is.null(container)) {
              if(is.logical(container) && container == TRUE)
                container = gwindow(visible=TRUE, toolkit=obj@toolkit)
              add(container,  obj,...)
            }
  
            ## add handler
            if(!is.null(handler))
              addHandlerChanged(obj,  handler, action)

            
            invisible(obj)
          })

##' NOOP
as.gWidgetsRGtk2.GtkRadioButton <- function(widget,...) {}

## methods
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gRadioRGtk"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {

            index = ifelse(is.null(index),FALSE,as.logical(index))

            ind <- obj@ref_widget$get_index()

            if(index)
              return(ind)
            else
              return(obj[ind])
          })

## svalue<-
setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gRadioRGtk"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   
                   if(is.data.frame(value))
                     value <- value[,1, drop=TRUE]
                   
                   index <- getWithDefault(index, is.logical(value))
                   if(!index) 
                     value <- match(value, obj[])

                   obj@ref_widget$set_index(value[1])

                   return(obj)
          })


setMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gRadioRGtk"),
          function(x, toolkit, i, j, ..., drop=TRUE) {
            items <- x@ref_widget$get_items()
            
            if(missing(i))
              items
            else
              items[i]
          })
            
setMethod("[",
          signature(x="gRadioRGtk"),
          function(x, i, j, ..., drop=TRUE) {
            .leftBracket(x, x@toolkit, i, j, ..., drop=drop)
          })

setReplaceMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gRadioRGtk"),
          function(x, toolkit, i, j, ..., value) {

            items <- x[]

            if(!missing(i))
              items[i] <- value
            else
              items <- value
            x@ref_widget$set_items(items)
            
            ## all done
            return(x)
          })

setReplaceMethod("[",
                 signature(x="gRadioRGtk"),
                 function(x, i, j,..., value) {
                   .leftBracket(x, x@toolkit, i, j, ...) <- value
                   return(x)
                 })

## length
setMethod(".length",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gRadioRGtk"),
          function(x,toolkit) {
            length(x[])
          })


##################################################
## handlers

## need to deal with changing buttons via [<-
## added a handlers cache that we can manipulate
setMethod(".addhandlerclicked",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gRadioRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            o <- Observer$new(o=handler, obj=obj, action=action)
            obj@ref_widget$add_observer(o)
            o
          })


setMethod(".addhandlerchanged",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gRadioRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            .addhandlerclicked(obj, toolkit, handler, action, ...)
          })
 

## ## There is an issue here. When we set values via [<- the handlers are gone!
## setMethod(".addhandlerclicked",
##           signature(toolkit="guiWidgetsToolkitRGtk2",obj="gRadioRGtk"),
##           function(obj, toolkit, handler, action=NULL, ...) {

##             radiogp <- getWidget(obj)
##             btns <- rev(radiogp$GetGroup())
            
##             IDs = sapply(btns, function(x) {
##               gtktry(connectSignal(x,
##                                    signal="toggled",
##                                    f=function(h,w,...) {
##                                      ## only call handler for change to active
##                                      ## not just toggle
##                                      if(w$GetActive())
##                                        handler(h,w,...)
##                                    },
##                                    data=list(obj=obj, action=action,...),
##                                    user.data.first = TRUE,
##                                    after = FALSE), silent=FALSE)
##             })
            
##             handler.ID = tag(obj, "handler.id")
##             if(is.null(handler.ID))
##               handler.ID =list()
##             for(i in 1:length(IDs))
##               handler.ID[[length(handler.ID)+1]] = IDs[[i]]
##             tag(obj, "handler.id", replace=FALSE) <- handler.ID
            
##             invisible(IDs)
##           })
 

##################################################
## ## constructor
## setMethod(".gradio",
##           signature(toolkit="guiWidgetsToolkitRGtk2"),
##           function(toolkit,
##                    items, selected=1, horizontal=FALSE,
##                    handler=NULL, action=NULL,
##                    container=NULL,       
##                    ...
##                    ) {

##             force(toolkit)

##             if(horizontal)
##               g <- gtkHBox()
##             else
##               g <- gtkVBox()

##             radiogp <- gtkRadioButton(group=NULL, label=items[1]) # initial

##             obj <- as.gWidgetsRGtk2(radiogp, block=g)

##             if(is.data.frame(items))
##               items <- items[,1, drop=TRUE] # first column

##             obj[] <- items
##             svalue(obj,index=TRUE) <- selected

##             tag(obj, ".handlers") <- list() # list of handlers keyed by ID
            
##             ## do we add to the container?
##             if (!is.null(container)) {
##               if(is.logical(container) && container == TRUE)
##                 container = gwindow(visible=TRUE, toolkit=obj@toolkit)
##               add(container,  obj,...)
##             }
  
##             ## add handler
##             if(!is.null(handler))
##               addHandlerChanged(obj,  handler, action)

            
##             invisible(obj)
##           })

## ##' coercion method from a gtkRadioButton widget. Pass in container via bloc
## as.gWidgetsRGtk2.GtkRadioButton <- function(widget,...) {
##   theArgs <- list(...)
##   if(!is.null(theArgs$block))
##     block <- theArgs$block
##   else
##     block <- gtkHBox()                  # or vbox!

##   obj <- new("gRadioRGtk",block=block, widget=widget,
##     toolkit=guiToolkit("RGtk2"))
##   return(obj)
## }

## ## methods
## setMethod(".svalue",
##           signature(toolkit="guiWidgetsToolkitRGtk2",obj="gRadioRGtk"),
##           function(obj, toolkit, index=NULL, drop=NULL, ...) {

##             index = ifelse(is.null(index),FALSE,as.logical(index))

##             radiogp <- getWidget(obj)
##             btns <- rev(radiogp$GetGroup())
##             ind <- sapply(btns, function(i) i$GetActive())

##             if(index)
##               return(which(ind))
##             else
##               return(obj[ind])
##           })

## ## svalue<-
## setReplaceMethod(".svalue",
##                  signature(toolkit="guiWidgetsToolkitRGtk2",obj="gRadioRGtk"),
##                  function(obj, toolkit, index=NULL, ..., value) {

##                    if(is.data.frame(value))
##                      value <- value[,1, drop=TRUE]
                   
##                    radiogp <- getWidget(obj)
##                    btns <- rev(radiogp$GetGroup())
##                    items <- obj[]
                   
##                    if(!is.null(index) && index==TRUE) {
##                      if(value %in% 1:length(obj))
##                        btns[[as.numeric(value)]]$SetActive(TRUE)
##                      else
##                        cat(sprintf("index outside of range\n"))
##                    } else {
##                      if(value %in% items) {
##                        whichIndex = min(which(value == items))
##                        btns[[whichIndex]]$SetActive(TRUE)
##                      } else {
##                        cat(sprintf("Value %s is not among the items\n",value))
##                      }
##                    }
                   
                   
##                    return(obj)
##           })


## setMethod(".leftBracket",
##           signature(toolkit="guiWidgetsToolkitRGtk2",x="gRadioRGtk"),
##           function(x, toolkit, i, j, ..., drop=TRUE) {
##             radiogp <- getWidget(x)
##             btns <- rev(radiogp$GetGroup())
##             btns <- btns[1:tag(x,".n")]
##             items <- sapply(btns, function(i) i$GetLabel())
            
##             if(missing(i))
##               items
##             else
##               items[i]
##           })
            
## setMethod("[",
##           signature(x="gRadioRGtk"),
##           function(x, i, j, ..., drop=TRUE) {
##             .leftBracket(x, x@toolkit, i, j, ..., drop=drop)
##           })

## setReplaceMethod(".leftBracket",
##           signature(toolkit="guiWidgetsToolkitRGtk2",x="gRadioRGtk"),
##           function(x, toolkit, i, j, ..., value) {

##             radiogp <- getWidget(x)
##             gp <- getBlock(x)
##             btns <- rev(radiogp$GetGroup())
##             n = length(x)
            
##             ## set items
##             if(missing(i)) {

##               ## The radio group doesn't like to reduce the size. We trick it by
##               ## keeping track of a variable length in ".n" tag
##               n <- length(value)
##               if(n < 2) {
##                 cat(sprintf("Length of items must be 2 or more\n"))
##                 return(x)
##               }

##               ## we store the length of items in the .n value.
##               ## When shortening a lenght by setting the group, the GTK
##               ## widget does not truncate. We use this to do so. (leaving some
##               ## possible orphans in the radioButtonGroup object.)
##               tag(x, ".n") <- n

##               ## clear old
##               sapply(gp$getChildren(), function(i) gp$remove(i))
##               ## make new
##               radiogp1 <- gtkRadioButton(group=NULL, label=value[1])
##               sapply(value[-1], function(i) {
##                 radiogp1$newWithLabelFromWidget(i)
##               })
##               ## replace -- doesn't clear, just replaces first n (even if more than n)
##               radiogp$setGroup(radiogp1$getGroup())
              
##               ## now add to container
##               btns <- rev(radiogp$getGroup())[1:n] # no more than n of them
##               sapply(btns, function(i) gp$PackStart(i))

##               ## need to add in the handlers
##               ## Always call to see if a handler exists
##               sapply(btns, function(i) {
##                 gSignalConnect(i, "toggled", f=function(obj, w, ...) {
##                   if(w$getActive()) {
##                     ## call handlers from h
##                     handlers <- tag(obj, ".handlers")
##                     if(length(handlers)) {
##                       ## handler is list with blocked, handler, action component
##                       sapply(handlers, function(handler) {
##                         if(!handler$blocked)
##                           handler$handler(list(obj=obj, action=handler$action), ...)
##                       })
##                     }
##                   }
##                 },
##                                data=x,
##                                user.data.first=TRUE,
##                                after=FALSE
##                                )
##               })
##             } else {
##               ## update just the i values
##               i <- i[i <= n]
##               for(j in 1:length(i)) 
##                 btns[[j]]$SetLabel(value[j])
##             }
            
##             ## all done
##             return(x)
##           })

## setReplaceMethod("[",
##                  signature(x="gRadioRGtk"),
##                  function(x, i, j,..., value) {
##                    .leftBracket(x, x@toolkit, i, j, ...) <- value
##                    return(x)
##                  })

## ## length
## setMethod(".length",
##           signature(toolkit="guiWidgetsToolkitRGtk2",x="gRadioRGtk"),
##           function(x,toolkit) {
##             tag(x, ".n")
## #            radiogp <- getWidget(x)
## #            btns <- rev(radiogp$GetGroup())
## #            length(btns)
##           })


## ## enabled must go on each button
## ## enabled <-
## setReplaceMethod(".enabled",
##                  signature(toolkit="guiWidgetsToolkitRGtk2",obj="gRadioRGtk"),
##                  function(obj, toolkit, ..., value) {
##                    radiogp <- getWidget(obj)
##                    btns <- rev(radiogp$GetGroup())
##                    sapply(btns, function(i) {
##                      i$SetSensitive(as.logical(value))
##                      })
##                    return(obj)
##                  })

## setReplaceMethod(".visible",
##                  signature(toolkit="guiWidgetsToolkitRGtk2",obj="gRadioRGtk"),
##                  function(obj, toolkit, ..., value) {
##                    radiogp <- getWidget(obj)
##                    btns <- rev(radiogp$GetGroup())
##                    sapply(btns, function(i) {
##                      if(value)
##                        i$show()
##                      else
##                        i$hide()
## #                     i$SetSensitive(as.logical(value))
##                    })
##                    return(obj)
##                  })

## ##################################################
## ## handlers

## ## need to deal with changing buttons via [<-
## ## added a handlers cache that we can manipulate
## setMethod(".addhandlerclicked",
##           signature(toolkit="guiWidgetsToolkitRGtk2",obj="gRadioRGtk"),
##           function(obj, toolkit, handler, action=NULL, ...) {

##             handlers <- tag(obj, ".handlers")
##             if(length(handlers))
##               nhandlers <- max(as.numeric(names(handlers)))
##             else
##               nhandlers <- 0

##             newhandler <- list(blocked=FALSE,
##                                handler=handler,
##                                action=action)
##             ID <- as.character(nhandlers + 1)
##             handlers[[ID]] <- newhandler
##             tag(obj, ".handlers") <- handlers

##             invisible(ID)
            
 
##           })
 



## setMethod(".addhandlerchanged",
##           signature(toolkit="guiWidgetsToolkitRGtk2",obj="gRadioRGtk"),
##           function(obj, toolkit, handler, action=NULL, ...) {
##             .addhandlerclicked(obj,toolkit,handler,action,...)
##           })



## setMethod(".removehandler",
##           signature(toolkit="guiWidgetsToolkitRGtk2",obj="gRadioRGtk"),
##           function(obj, toolkit, ID=NULL, ...) {
##             handlers <- tag(obj, ".handlers")
            
##             if(is.null(ID)) {
##               handlers <- list()        # remove all
##             } else {
##               sapply(ID, function(id) {
##                 handlers[[id]] <<- NULL
##               })
##             }
##             tag(obj, ".handlers") <- handlers          
##           })

## setMethod(".blockhandler",
##           signature(toolkit="guiWidgetsToolkitRGtk2",obj="gRadioRGtk"),
##           function(obj, toolkit, ID=NULL, ...) {
##             handlers <- tag(obj, ".handlers")
##             if(is.null(ID)) {
##               ID <- names(handlers)
##             }
##             sapply(ID, function(id) {
##               handlers[[id]]$blocked <<- TRUE
##             })
##             tag(obj, ".handlers") <- handlers          
##           })

## setMethod(".unblockhandler",
##           signature(toolkit="guiWidgetsToolkitRGtk2",obj="gRadioRGtk"),
##           function(obj, toolkit, ID=NULL, ...) {
##             handlers <- tag(obj, ".handlers")
##             if(is.null(ID)) {
##               ID <- names(handlers)
##             }
##             sapply(ID, function(id) {
##               handlers[[id]]$blocked <<- FALSE
##             })
##             tag(obj, ".handlers") <- handlers          
##           })

## ## ## There is an issue here. When we set values via [<- the handlers are gone!
## ## setMethod(".addhandlerclicked",
## ##           signature(toolkit="guiWidgetsToolkitRGtk2",obj="gRadioRGtk"),
## ##           function(obj, toolkit, handler, action=NULL, ...) {

## ##             radiogp <- getWidget(obj)
## ##             btns <- rev(radiogp$GetGroup())
            
## ##             IDs = sapply(btns, function(x) {
## ##               gtktry(connectSignal(x,
## ##                                    signal="toggled",
## ##                                    f=function(h,w,...) {
## ##                                      ## only call handler for change to active
## ##                                      ## not just toggle
## ##                                      if(w$GetActive())
## ##                                        handler(h,w,...)
## ##                                    },
## ##                                    data=list(obj=obj, action=action,...),
## ##                                    user.data.first = TRUE,
## ##                                    after = FALSE), silent=FALSE)
## ##             })
            
## ##             handler.ID = tag(obj, "handler.id")
## ##             if(is.null(handler.ID))
## ##               handler.ID =list()
## ##             for(i in 1:length(IDs))
## ##               handler.ID[[length(handler.ID)+1]] = IDs[[i]]
## ##             tag(obj, "handler.id", replace=FALSE) <- handler.ID
            
## ##             invisible(IDs)
## ##           })
 

