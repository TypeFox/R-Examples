## Use reference class, like gradio

CbgWidgetGtk <- setRefClass("CbgWidgetGtk",
                     contains="GWidgetGtk",
                     fields=list(
                       inner_block="ANY",    # replaceble box container
                       items="ANY",          # store the items
                       horizontal="logical", # layout direction
                       obj = "ANY"           # gradio object for callbacks
                       ),
                     methods=list(
                       initialize=function(items, checked=FALSE, horizontal=TRUE) {
                         horizontal <<- horizontal
                         block <<- gtkHBox()
                         inner_block <<- gtkHBox(); block$packStart(inner_block)
                         if(!missing(items)) {
                           set_items(items)
                           set_index(which(rep(checked, length.out=length(items))))
                         }
                         .self
                       },
                       get_items = function() {
                         "Return items"
                         items
                       },
                       set_items=function(items) {
                         if(length(items) == 0) return()
                         items <<- items

                         block$remove(inner_block)
                         inner_block <<- if(horizontal) gtkHBox() else gtkVBox()
                         block$packStart(inner_block)

                         widget <<- lapply(items, gtkCheckButton)
                         lapply(widget, gtkBoxPackStart, object = inner_block)
                         
                         ## add handlers
                         lapply(widget, gSignalConnect, signal="toggled", f=function(self, w, ...) {
                           self$notify_observers(...)
                         }, data=.self, user.data.first=TRUE)

                         invisible()
                       },
                       get_index = function() {
                         "Return indices of selected"
                         which(sapply(widget, gtkToggleButtonGetActive))
                       },
                       set_index = function(i) {
                         "Set selection indices"
                         ind <- rep(FALSE, length(items))
                         ind[as.integer(i)] <- TRUE
                         sapply(seq_along(widget), function(j) widget[[j]]$setActive(ind[j]))

                         invisible()
                       }
                       ))


setClass("gCheckboxgroupRGtk",
         contains="gComponentWithRefClassWithItemsRGtk",
         prototype=prototype(new("gComponentWithRefClassWithItemsRGtk"))
         )


setMethod(".gcheckboxgroup",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   items, checked = FALSE,
                   horizontal=FALSE, use.table=FALSE,
                   handler = NULL, action = NULL, container = NULL, ...) {

            force(toolkit)

            if(as.logical(use.table)) {
              obj <- .gcheckboxgrouptable(toolkit,
                                          items, checked=checked,
                                          handler=handler, action=action,
                                          container=container, ...)
              return(obj)
            }

            
            if(missing(items))
              stop(gettext("Need items to be defined"))

            if(is.data.frame(items))
              items <- items[, 1, drop=TRUE]

            
            checked = rep(checked, length(items))

            ref_widget <- CbgWidgetGtk$new(items, checked, horizontal)

            obj = new("gCheckboxgroupRGtk",
              block=ref_widget$block,
              widget=ref_widget$block,
              ref_widget=ref_widget,
              toolkit=toolkit)


            ## do we add to the container?
            if (!is.null(container)) {
              if(is.logical(container) && container == TRUE)
                container = gwindow(visible=TRUE, toolkit=obj@toolkit)
              add(container,  obj,...)
            }
  
            ## add handler
            if(!is.null(handler))
              ID = addhandlerchanged(obj, handler=handler, action=action, ...)
            
            return(obj)
          })


### methods
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gCheckboxgroupRGtk"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {

            index <- getWithDefault(index, FALSE)
            vals <- obj@ref_widget$get_index()
            
            if(index) {
              vals
            } else {
              obj[vals]
            }
          })

## toggles state to be T or F
setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gCheckboxgroupRGtk"),
                 function(obj, toolkit, index=NULL, ..., value) {

                   if(is.data.frame(value))
                     value <- value[,1,drop=TRUE]

                   index <- getWithDefault(index, FALSE)

                   ## compute values -- logical vector with length n
                   if(!index) {
                     if(!is.logical(value)) {
                       ## characters
                       value <- match(value, obj[])
                     } else {
                       value = rep(value, length.out=length(obj)) ## recycle
                       value <- which(value)
                     }
                   }
                   obj@ref_widget$set_index(value)

                   return(obj)
                 })

## [ and [<- refer to the names -- not the TF values

setMethod("[",
          signature(x="gCheckboxgroupRGtk"),
          function(x, i, j, ..., drop=TRUE) {
            .leftBracket(x, x@toolkit, i, j, ..., drop=drop)
          })
setMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gCheckboxgroupRGtk"),
          function(x, toolkit, i, j, ..., drop=TRUE) {
            items <- x@ref_widget$get_items()
            if(missing(i))
              return(items)
            else
              return(items[i])
          })

## assigns names
setReplaceMethod("[",
                 signature(x="gCheckboxgroupRGtk"),
                 function(x, i, j,..., value) {
                   .leftBracket(x, x@toolkit, i, j, ...) <- value
                   return(x)
                 })

setReplaceMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gCheckboxgroupRGtk"),
          function(x, toolkit, i, j, ..., value) {
            items = x[]
            if(!missing(i))
              items[i] <- value
            else
              items <- value

            x@ref_widget$set_items(items)
  
             return(x)
          })

## length
setMethod(".length",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gCheckboxgroupRGtk"),
          function(x, toolkit) {
            length(x[])
          })


## handlers
setMethod(".addhandlerclicked",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gCheckboxgroupRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            o <- Observer$new(o=handler, obj=obj, action=action)
            obj@ref_widget$add_observer(o)
            o
          })

setMethod(".addhandlerchanged",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gCheckboxgroupRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            .addhandlerclicked(obj, toolkit, handler, action, ...)
          })



##################################################
##################################################

### Checkbox group in a table
setClass("gCheckboxgroupTableRGtk",
         contains="gComponentRGtk",
         prototype=prototype(new("gComponentRGtk"))
         )

setGeneric(".gcheckboxgrouptable", function(toolkit, items, checked=FALSE,
                                            handler=NULL, action=NULL,
                                            container=NULL, ...)
           standardGeneric(".gcheckboxgrouptable"))

setMethod(".gcheckboxgrouptable",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   items, checked = FALSE,
                   handler = NULL, action = NULL, container = NULL, ...) {


            force(toolkit)

            
            tbl <- gtkTreeViewNew(TRUE)
            tbl$SetRulesHint(TRUE)      # shade
            
            store <- rGtkDataFrame(.makeItems())
            tbl$setModel(store)
            tbl$setHeadersVisible(FALSE)

            sw <- gtkScrolledWindowNew()
            sw$SetPolicy("GTK_POLICY_AUTOMATIC","GTK_POLICY_AUTOMATIC")
            sw$Add(tbl)


            ## set up the view columns
            vc <- gtkTreeViewColumnNew()
            tbl$insertColumn(vc, 0)
            cr <- gtkCellRendererToggle()
            vc$PackStart(cr, TRUE)
            cr['activatable'] <- TRUE                  # needed
            vc$addAttribute(cr, "active", 1)            
            item.toggled <- function(tbl, cell, path, data) {
              store <- tbl$getModel()
              row <- as.numeric(path) + 1
              store[row,2] <- !store[row, 2]
            }
            gSignalConnect(cr, "toggled", item.toggled, data=tbl, user.data.first=TRUE)

            
            
            cr <- gtkCellRendererTextNew()
            vc <- gtkTreeViewColumnNew()
            vc$PackStart(cr, TRUE)
            vc$addAttribute(cr, "text", 0)            
            tbl$insertColumn(vc, 1)

            ## how to add icons, tooltips!
            
            ## make combination widget with all the values
            obj = new("gCheckboxgroupTableRGtk", block=sw, widget=tbl,
              toolkit=toolkit)

            obj[] <- items       
            svalue(obj) <- checked
            
            
            if(!is.null(handler))
              tag(obj, "handler.id") <- addhandlerchanged(obj,handler,action)

            if(!is.null(container)) {
              if(is.logical(container)) {
                if(container) {
                  container <- gwindow()
                } else {
                  return(obj)
                }
              }
              add(container, obj, ...)
            }
            
            return(obj)
          })

## helper
.makeItems <- function(items, icons, tooltips, checked=rep(FALSE, length(items))) {
  if(missing(items) ||
     (is.data.frame(items) && nrow(items) == 0) ||
     (length(items) == 0)
     ) {
    out <- data.frame(items=character(0),
                      checked=logical(0),
                      icons=character(0),
                      tooltips=character(0),
                      stringsAsFactors=FALSE)
  } else if(is.data.frame(items)) {
    ## check
    out <- items
    if(ncol(out) == 1) 
      out$checked <- as.logical(rep(checked, length=nrow(items)))
    if(ncol(out) == 2)
      out$icons <- rep("", nrow(items))
    if(ncol(out) == 3)
      out$tooltip <- rep("", nrow(items))
  } else {
    ## piece together
    items <- as.character(items)
    
    if(missing(icons))
      icons <- ""
    icons <- rep(icons, length=length(items))

    if(missing(tooltips))
      tooltips <- ""
    icons <- rep(tooltips, length=length(items))

    checked <- rep(checked, length=length(items))

    out <- data.frame(items=items, checked=checked, icons=icons, tooltips=tooltips,
                      stringsAsFactors=FALSE)
  }
  return(out)
}
    


### methods
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gCheckboxgroupTableRGtk"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {
            n <- length(obj)
            if(n == 0)
              return(logical(0))

            tbl <- getWidget(obj)
            store <- tbl$getModel()
            vals <- store[,2, drop=TRUE]
            index <- getWithDefault(index, FALSE)

            if(index) {
              return(which(vals))       # return indices
            } else {
              obj[vals]                 # vals is logical
            }
          })

## toggles state to be T or F
setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gCheckboxgroupTableRGtk"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   
                   n <- length(obj)
                   if(n == 0)
                     return(obj)

                   tbl <- getWidget(obj)
                   store <- tbl$getModel()

                   index <- getWithDefault(index, FALSE)
                   if(!index) {
                     if(is.logical(value)) {
                       value <- rep(value, length.out=n)
                       value <- which(value)
                     } else {
                       value <- match(value, obj[])
                     }
                   }

                   ## value is index, we want logical
                   ind <- rep(FALSE, n)
                   ind[value] <- TRUE
                   store[,2] <- ind
                   
                   return(obj)
                 })

## [ and [<- refer to the names -- not the TF values
## Here we can have a vector of names -- or a data frame
## 1st column names, 2nd icon, third tooltip -- like gcombobox
setMethod("[",
          signature(x="gCheckboxgroupTableRGtk"),
          function(x, i, j, ..., drop=TRUE) {
            .leftBracket(x, x@toolkit, i, j, ..., drop=drop)
          })
setMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gCheckboxgroupTableRGtk"),
          function(x, toolkit, i, j, ..., drop=TRUE) {
            if(length(x) == 0)
              return(character(0))

            tbl <- getWidget(x)
            store <- tbl$getModel()
            
            items <- store[,1, drop=TRUE]
            
            if(missing(i))
              return(items)
            else
              return(items[i])
          })

## assigns names
setReplaceMethod("[",
                 signature(x="gCheckboxgroupTableRGtk"),
                 function(x, i, j,..., value) {
                   .leftBracket(x, x@toolkit, i, j, ...) <- value
                   return(x)
                 })

setReplaceMethod(".leftBracket",
                 signature(toolkit="guiWidgetsToolkitRGtk2",x="gCheckboxgroupTableRGtk"),
                 function(x, toolkit, i, j, ..., value) {
                   ## value can be a vector or data frame
                   ## if a data.frame we have
                   ## text, stockicon, tooltip
                   items <- .makeItems(value)
                   
                   tbl <- getWidget(x)
                   store <- tbl$getModel()
                   
                   if(missing(i)) {
                     ## replace the store
                     newStore <- rGtkDataFrame(items)
                     tbl$setModel(newStore)
                   } else {
                     if(is.logical(i))
                       i = which(i)
                     
                     ## set items
                     m <- nrow(items)
                     if(m == 0)
                       return(x)
                     
                     store[i,] <- items
                   }
                 
                   return(x)
                 })


setMethod(".length",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gCheckboxgroupTableRGtk"),
          function(x,toolkit) {
            tbl <- getWidget(x)
            store <- tbl$getModel()
            dim(store)[1]
          })



## Handlers must just pass down to each item in the list.
setMethod(".addhandlerchanged",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gCheckboxgroupTableRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            .addhandlerclicked(obj, toolkit, handler=handler,action=action,...)
          })

## clicked is changed
setMethod(".addhandlerclicked",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gCheckboxgroupTableRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            ## push down to cr
            tbl <- getWidget(obj)
            vc <- tbl$getColumn(0)
            cr <- vc$getCellRenderers()[[1]]
            ID <- gSignalConnect(cr, "toggled", function(h,...) handler(h),
                                 user.data.first=TRUE,
                                 data=list(obj=obj, action=action))
            invisible(ID)
          })

setMethod(".removehandler",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gCheckboxgroupTableRGtk"),
          function(obj, toolkit, ID=NULL, ...) {
            tbl <- getWidget(obj)
            vc <- tbl$getColumn(0)
            cr <- vc$getCellRenderers()[[1]]
            gSignalHandlerDisconnect(cr, ID)
          })

setMethod(".blockhandler",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gCheckboxgroupTableRGtk"),
          function(obj, toolkit, ID=NULL, ...) {
            tbl <- getWidget(obj)
            vc <- tbl$getColumn(0)
            cr <- vc$getCellRenderers()[[1]]
            gSignalHandlerBlock(cr, ID)
          })

setMethod(".unblockhandler",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gCheckboxgroupTableRGtk"),
          function(obj, toolkit, ID=NULL, ...) {
            tbl <- getWidget(obj)
            vc <- tbl$getColumn(0)
            cr <- vc$getCellRenderers()[[1]]
            gSignalHandlerUnblock(cr, ID)
          })

##################################################
## build widget based on gcheckbox
## setMethod(".gcheckboxgroup",
##           signature(toolkit="guiWidgetsToolkitRGtk2"),
##           function(toolkit,
##                    items, checked = FALSE,
##                    horizontal=FALSE, use.table=FALSE,
##                    handler = NULL, action = NULL, container = NULL, ...) {

##             force(toolkit)

##             if(as.logical(use.table)) {
##               obj <- .gcheckboxgrouptable(toolkit,
##                                           items, checked=checked,
##                                           handler=handler, action=action,
##                                           container=container, ...)
##               return(obj)
##             }

            
##             if(missing(items))
##               stop(gettext("Need items to be defined"))

##             if(is.data.frame(items))
##               items <- items[, 1, drop=TRUE]

            
##             checked = rep(checked, length(items))

##             group = ggroup(horizontal = horizontal, container=container, ...)
            
##             lst = list()
##             n = length(items)
##             for(i in 1:n) {
##               newItem = gcheckbox(items[i], checked=checked[i])
##               lst[[ as.character(items[i]) ]] = newItem
##               add(group, newItem)
##             }
  

##             ## make combination widget with all the values
##             obj = new("gCheckboxgroupRGtk",block=group, widget=group, toolkit=toolkit)
  
##             tag(obj, "items") <- items
##             tag(obj, "itemlist") <- lst
##             tag(obj, "handlerList") <- list()
##             tag(obj, "handlerCount") <- 0

##             ## add handler
##             if(!is.null(handler))
##               ID = addhandlerchanged(obj, handler=handler, action=action, ...)
            
##             return(obj)
##           })


## ### methods
## setMethod(".svalue",
##           signature(toolkit="guiWidgetsToolkitRGtk2",obj="gCheckboxgroupRGtk"),
##           function(obj, toolkit, index=NULL, drop=NULL, ...) {
##             theArgs = list(...)
            
##             lst = tag(obj, "itemlist")
##             vals = sapply(lst, svalue)         # logicals
            
##             if(!is.null(index) && index == TRUE) {
##               return(which(vals))
##             } else {
##               return(tag(obj,"items")[vals])
##             }
##           })

## ## toggles state to be T or F
## setReplaceMethod(".svalue",
##                  signature(toolkit="guiWidgetsToolkitRGtk2",obj="gCheckboxgroupRGtk"),
##                  function(obj, toolkit, index=NULL, ..., value) {
##                    if(is.data.frame(value))
##                      value <- value[,1,drop=TRUE]
                   
##                    lst = tag(obj,"itemlist")
##                    n <- length(obj)
##                    ## compute values -- logical vector with length n
##                    if(!is.null(index) && index) {
##                      ## indices
##                      values <- rep(FALSE, n)
##                      values[value] <- TRUE
##                    } else if(!is.logical(value)) {
##                      ## characters
##                     ind <- match(value, obj[])
##                     ind <- ind[!is.na(ind)]
##                     values <- rep(FALSE,length=n)
##                     values[ind] <- TRUE
##                    } else {
##                      ## logical vector, we recycle
##                      values = rep(value, length.out=n) ## recycle
##                    }
##                    ## apply to each checkbox
##                    sapply(1:n, function(i) svalue(lst[[i]]) <- values[i])

##                    return(obj)
##                  })

## ## [ and [<- refer to the names -- not the TF values

## setMethod("[",
##           signature(x="gCheckboxgroupRGtk"),
##           function(x, i, j, ..., drop=TRUE) {
##             .leftBracket(x, x@toolkit, i, j, ..., drop=drop)
##           })
## setMethod(".leftBracket",
##           signature(toolkit="guiWidgetsToolkitRGtk2",x="gCheckboxgroupRGtk"),
##           function(x, toolkit, i, j, ..., drop=TRUE) {
##             items = tag(x,"items")
##             if(missing(i))
##               return(items)
##             else
##               return(items[i])
##           })

## ## assigns names
## setReplaceMethod("[",
##                  signature(x="gCheckboxgroupRGtk"),
##                  function(x, i, j,..., value) {
##                    .leftBracket(x, x@toolkit, i, j, ...) <- value
##                    return(x)
##                  })

## setReplaceMethod(".leftBracket",
##           signature(toolkit="guiWidgetsToolkitRGtk2",x="gCheckboxgroupRGtk"),
##           function(x, toolkit, i, j, ..., value) {
##             items = tag(x,"items")
##             lst = tag(x,"itemlist")
##             n = length(items)

##             ## if i is missing, we can relabel if length the same
##             ## otherwise we delete and start again
##             ## We will need to add the handlers back
            
##             if(missing(i)) {
##               if(length(value) != n) {
##                 group <- x@widget
##                 ## delete
##                 sapply(rev(lst), function(child)
##                        delete(group, child))
##                 ## add
##                 lst <- list()
##                 for(i in 1:length(value)) {
##                   newItem = gcheckbox(value[i], checked=FALSE)
##                   lst[[ as.character(value[i]) ]] = newItem
##                   add(group, newItem)
##                 }
##                 tag(x, "items") <- value
##                 tag(x, "itemlist") <- lst

##                 ## addhandlers
##                 handlerList <- tag(x,"handlerList")
##                 if(length(handlerList) > 0) {
##                   for(j in handlerList) {
##                     sapply(lst, function(i)
##                            addhandlerchanged(i,
##                                              handler=j$handler, action=j$action,
##                                              actualobj=x, ...))
##                   }
##                 }
##                 ## return
##                 return(x)
##               } else {
##                 ## back to our regularly scheduled programming
##                 i = 1:n
##               }
##             }
  
##             if(is.logical(i))
##               i = which(i)

##             items[i] = value
##             sapply(1:n, function(i) 
##                    lst[[i]][] <- items[i]
##                    )
##             tag(x,"items") <- items
##             tag(x,"itemlist") <- lst
  
##              return(x)
##           })

## ## length
## setMethod(".length",
##           signature(toolkit="guiWidgetsToolkitRGtk2",x="gCheckboxgroupRGtk"),
##           function(x,toolkit) {
##             length(tag(x,"items"))
##           })


## ## handlers
## setMethod(".addhandlerchanged",
##           signature(toolkit="guiWidgetsToolkitRGtk2",obj="gCheckboxgroupRGtk"),
##           function(obj, toolkit, handler, action=NULL, ...) {
##             handlerList <- tag(obj,"handlerList")
##             ct <- tag(obj,"handlerCount")
##             ID <- as.character(ct+1)
##             handlerList[[ID]] <- list(
##                                       handler=handler,
##                                       action=action
##                                       )
##             tag(obj,"handlerList") <- handlerList
##             ## now call on each
##             lst = tag(obj,"itemlist")
##             IDs <- lapply(lst, function(i)
##                    addhandlerchanged(i,handler=handler, action=action, actualobj=obj, ...))
##             invisible(IDs)
##           })
          

## setMethod(".removehandler",
##           signature(toolkit="guiWidgetsToolkitRGtk2",obj="gCheckboxgroupRGtk"),
##           function(obj, toolkit, ID=NULL, ...) {
##             tag(obj,"handlerList") <- NULL
##             lst <- tag(obj,"itemlist")
##             sapply(1:length(lst), function(i)
##                    removehandler(lst[[i]], ID[[i]])
##                  )
##           })

## setMethod(".blockhandler",
##           signature(toolkit="guiWidgetsToolkitRGtk2",obj="gCheckboxgroupRGtk"),
##           function(obj, toolkit, ID=NULL, ...) {

##             lst <- tag(obj,"itemlist")
##             sapply(1:length(lst), function(i)
##                    blockhandler(lst[[i]], ID[[i]])
##                    )
##           })

## setMethod(".unblockhandler",
##           signature(toolkit="guiWidgetsToolkitRGtk2",obj="gCheckboxgroupRGtk"),
##           function(obj, toolkit, ID=NULL, ...) {

##             lst <- tag(obj,"itemlist")
##             sapply(1:length(lst), function(i)
##               unblockhandler(lst[[i]], ID[[i]])
##             )
##           })



## ##################################################
## ##################################################

## ### Checkbox group in a table
## setClass("gCheckboxgroupTableRGtk",
##          contains="gComponentRGtk",
##          prototype=prototype(new("gComponentRGtk"))
##          )

## setGeneric(".gcheckboxgrouptable", function(toolkit, items, checked=FALSE,
##                                             handler=NULL, action=NULL,
##                                             container=NULL, ...)
##            standardGeneric(".gcheckboxgrouptable"))

## setMethod(".gcheckboxgrouptable",
##           signature(toolkit="guiWidgetsToolkitRGtk2"),
##           function(toolkit,
##                    items, checked = FALSE,
##                    handler = NULL, action = NULL, container = NULL, ...) {


##             force(toolkit)

            
##             tbl <- gtkTreeViewNew(TRUE)
##             tbl$SetRulesHint(TRUE)      # shade
            
##             store <- rGtkDataFrame(.makeItems())
##             tbl$setModel(store)
##             tbl$setHeadersVisible(FALSE)

##             sw <- gtkScrolledWindowNew()
##             sw$SetPolicy("GTK_POLICY_AUTOMATIC","GTK_POLICY_AUTOMATIC")
##             sw$Add(tbl)


##             ## set up the view columns
##             vc <- gtkTreeViewColumnNew()
##             tbl$insertColumn(vc, 0)
##             cr <- gtkCellRendererToggle()
##             vc$PackStart(cr, TRUE)
##             cr['activatable'] <- TRUE                  # needed
##             vc$addAttribute(cr, "active", 1)            
##             item.toggled <- function(tbl, cell, path, data) {
##               store <- tbl$getModel()
##               row <- as.numeric(path) + 1
##               store[row,2] <- !store[row, 2]
##             }
##             gSignalConnect(cr, "toggled", item.toggled, data=tbl, user.data.first=TRUE)

            
            
##             cr <- gtkCellRendererTextNew()
##             vc <- gtkTreeViewColumnNew()
##             vc$PackStart(cr, TRUE)
##             vc$addAttribute(cr, "text", 0)            
##             tbl$insertColumn(vc, 1)

##             ## how to add icons, tooltips!
            
##             ## make combination widget with all the values
##             obj = new("gCheckboxgroupTableRGtk", block=sw, widget=tbl,
##               toolkit=toolkit)

##             obj[] <- items       
##             svalue(obj) <- checked
            
            
##             if(!is.null(handler))
##               tag(obj, "handler.id") <- addhandlerchanged(obj,handler,action)

##             if(!is.null(container)) {
##               if(is.logical(container)) {
##                 if(container) {
##                   container <- gwindow()
##                 } else {
##                   return(obj)
##                 }
##               }
##               add(container, obj, ...)
##             }
            
##             return(obj)
##           })

## ## helper
## .makeItems <- function(items, icons, tooltips, checked=rep(FALSE, length(items))) {
##   if(missing(items) ||
##      (is.data.frame(items) && nrow(items) == 0) ||
##      (length(items) == 0)
##      ) {
##     out <- data.frame(items=character(0),
##                       checked=logical(0),
##                       icons=character(0),
##                       tooltips=character(0),
##                       stringsAsFactors=FALSE)
##   } else if(is.data.frame(items)) {
##     ## check
##     out <- items
##     if(ncol(out) == 1) 
##       out$checked <- as.logical(rep(checked, length=nrow(items)))
##     if(ncol(out) == 2)
##       out$icons <- rep("", nrow(items))
##     if(ncol(out) == 3)
##       out$tooltip <- rep("", nrow(items))
##   } else {
##     ## piece together
##     items <- as.character(items)
    
##     if(missing(icons))
##       icons <- ""
##     icons <- rep(icons, length=length(items))

##     if(missing(tooltips))
##       tooltips <- ""
##     icons <- rep(tooltips, length=length(items))

##     checked <- rep(checked, length=length(items))

##     out <- data.frame(items=items, checked=checked, icons=icons, tooltips=tooltips,
##                       stringsAsFactors=FALSE)
##   }
##   return(out)
## }
    


## ### methods
## setMethod(".svalue",
##           signature(toolkit="guiWidgetsToolkitRGtk2",obj="gCheckboxgroupTableRGtk"),
##           function(obj, toolkit, index=NULL, drop=NULL, ...) {
##             n <- length(obj)
##             if(n == 0)
##               return(logical(0))

##             tbl <- getWidget(obj)
##             store <- tbl$getModel()
##             vals <- store[,2, drop=TRUE]
##             index <- getWithDefault(index, FALSE)

##             if(index) {
##               return(which(vals))       # return indices
##             } else {
##               obj[vals]                 # vals is logical
##             }
##           })

## ## toggles state to be T or F
## setReplaceMethod(".svalue",
##                  signature(toolkit="guiWidgetsToolkitRGtk2",obj="gCheckboxgroupTableRGtk"),
##                  function(obj, toolkit, index=NULL, ..., value) {
                   
##                    n <- length(obj)
##                    if(n == 0)
##                      return(obj)

##                    tbl <- getWidget(obj)
##                    store <- tbl$getModel()
                   
##                    index <- getWithDefault(index, is.numeric(value))
##                    if(index) {
##                      tmp <- rep(FALSE, n)
##                      tmp[value] <- TRUE
##                      value <- tmp
##                    }
##                    ## recycle
##                    value <- as.logical(rep(value, length=n))
##                    store[,2] <- value
                   
##                    return(obj)
##                  })

## ## [ and [<- refer to the names -- not the TF values
## ## Here we can have a vector of names -- or a data frame
## ## 1st column names, 2nd icon, third tooltip -- like gcombobox
## setMethod("[",
##           signature(x="gCheckboxgroupTableRGtk"),
##           function(x, i, j, ..., drop=TRUE) {
##             .leftBracket(x, x@toolkit, i, j, ..., drop=drop)
##           })
## setMethod(".leftBracket",
##           signature(toolkit="guiWidgetsToolkitRGtk2",x="gCheckboxgroupTableRGtk"),
##           function(x, toolkit, i, j, ..., drop=TRUE) {
##             if(length(x) == 0)
##               return(character(0))

##             tbl <- getWidget(x)
##             store <- tbl$getModel()
            
##             items <- store[,1, drop=TRUE]
            
##             if(missing(i))
##               return(items)
##             else
##               return(items[i])
##           })

## ## assigns names
## setReplaceMethod("[",
##                  signature(x="gCheckboxgroupTableRGtk"),
##                  function(x, i, j,..., value) {
##                    .leftBracket(x, x@toolkit, i, j, ...) <- value
##                    return(x)
##                  })

## setReplaceMethod(".leftBracket",
##                  signature(toolkit="guiWidgetsToolkitRGtk2",x="gCheckboxgroupTableRGtk"),
##                  function(x, toolkit, i, j, ..., value) {
##                    ## value can be a vector or data frame
##                    ## if a data.frame we have
##                    ## text, stockicon, tooltip
##                    items <- .makeItems(value)
                   
##                    tbl <- getWidget(x)
##                    store <- tbl$getModel()
                   
##                    if(missing(i)) {
##                      ## replace the store
##                      newStore <- rGtkDataFrame(items)
##                      tbl$setModel(newStore)
##                    } else {
##                      if(is.logical(i))
##                        i = which(i)
                     
##                      ## set items
##                      m <- nrow(items)
##                      if(m == 0)
##                        return(x)
                     
##                      store[i,] <- items
##                    }
                 
##                    return(x)
##                  })


## setMethod(".length",
##           signature(toolkit="guiWidgetsToolkitRGtk2",x="gCheckboxgroupTableRGtk"),
##           function(x,toolkit) {
##             tbl <- getWidget(x)
##             store <- tbl$getModel()
##             dim(store)[1]
##           })



## ## Handlers must just pass down to each item in the list.
## setMethod(".addhandlerchanged",
##           signature(toolkit="guiWidgetsToolkitRGtk2",obj="gCheckboxgroupTableRGtk"),
##           function(obj, toolkit, handler, action=NULL, ...) {
##             .addhandlerclicked(obj, toolkit, handler=handler,action=action,...)
##           })

## ## clicked is changed
## setMethod(".addhandlerclicked",
##           signature(toolkit="guiWidgetsToolkitRGtk2",obj="gCheckboxgroupTableRGtk"),
##           function(obj, toolkit, handler, action=NULL, ...) {
##             ## push down to cr
##             tbl <- getWidget(obj)
##             vc <- tbl$getColumn(0)
##             cr <- vc$getCellRenderers()[[1]]
##             ID <- gSignalConnect(cr, "toggled", function(h,...) handler(h),
##                                  user.data.first=TRUE,
##                                  data=list(obj=obj, action=action))
##             invisible(ID)
##           })

## setMethod(".removehandler",
##           signature(toolkit="guiWidgetsToolkitRGtk2",obj="gCheckboxgroupTableRGtk"),
##           function(obj, toolkit, ID=NULL, ...) {
##             tbl <- getWidget(obj)
##             vc <- tbl$getColumn(0)
##             cr <- vc$getCellRenderers()[[1]]
##             gSignalHandlerDisconnect(cr, ID)
##           })

## setMethod(".blockhandler",
##           signature(toolkit="guiWidgetsToolkitRGtk2",obj="gCheckboxgroupTableRGtk"),
##           function(obj, toolkit, ID=NULL, ...) {
##             tbl <- getWidget(obj)
##             vc <- tbl$getColumn(0)
##             cr <- vc$getCellRenderers()[[1]]
##             gSignalHandlerBlock(cr, ID)
##           })

## setMethod(".unblockhandler",
##           signature(toolkit="guiWidgetsToolkitRGtk2",obj="gCheckboxgroupTableRGtk"),
##           function(obj, toolkit, ID=NULL, ...) {
##             tbl <- getWidget(obj)
##             vc <- tbl$getColumn(0)
##             cr <- vc$getCellRenderers()[[1]]
##             gSignalHandlerUnblock(cr, ID)
##           })
