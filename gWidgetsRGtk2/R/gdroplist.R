## editable has entry widget that can be edited
setClass("gDroplistRGtk",
         contains="gComponentRGtk",
         prototype=prototype(new("gComponentRGtk"))
         )

##' Combobox widget
##' @param items vector of names; 1-column data.frame of names; 2-column names, icons; 3-column names, icons, tooltip
##' @param selected index of initial, 0 if blank
##' @param editible -- are we editable?
setMethod(".gdroplist",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   items, selected = 1, # use 0 for blank
                   editable=FALSE,
                   coerce.with = NULL,
                   handler=NULL, action=NULL,
                   container=NULL,
                   ...               # do.quote = TRUE for quote of answer
                   ) {

            force(toolkit)

            ## Changed this. Make objects a data.frame if
            ## two columns, second is a stock-icon name.
            ## ideally third column would specify tooltip
            
            ## items must be a vector or data frame
            if(!inherits(items,"data.frame")) {
              items = as.vector(items)              # undoes factor
              items = unique(items)                 # unique
              items = data.frame(items, stringsAsFactors=FALSE)
            }
            
            doIcons = ifelse(ncol(items) >= 2, TRUE, FALSE)
            if(ncol(items) == 3) {
              types <- c(data="gchararray",icons="gchararray", tooltip="gchararray")
            } else if(ncol(items) == 2) {
              types <- c(data="gchararray",icons="gchararray")
            } else {
              types <- c(dataonly="gchararray")
            }
            
            theArgs = list(...)
            
            ## keep this, but don't advertise
            if(!is.null(theArgs$do.quote)) {
              coerce.with = function(x) paste("'",x,"'",sep="",collapse="")
            }
            
            
            ## droplist is not happy with datastore class
            ## droplist was not happy with numeric vectors! seems strange
            
            if(editable) {
              store = gtkListStoreNew(types)
              combo <- gtkComboBoxEntryNewWithModel(store, 0)
             
              
              ## now add icon if there
              if(ncol(items)  >= 2) {
                ## icon renderer
                cellrenderer = gtkCellRendererPixbufNew()
                combo$PackStart(cellrenderer, expand=FALSE)
                combo$AddAttribute(cellrenderer, "stock-id", 1)
              }
              
              entry = combo$GetChild()
              entry$SetEditable(TRUE)
              ## add in drop target to entry
              ## we can't pass in obj here, so we find via scoping
              dropHandler =   function(h,...) {
                theName = id(h$dropdata)
                ## override value -- in case it is a widget
                tag(obj, "value") <- h$dropdata # find obj via scoping
                svalue(obj) <- "" 
                return(TRUE)
              }
              .adddroptarget(entry, toolkit, targetType="object",handler=dropHandler)
#              .adddroptarget(entry, toolkit, targetType="object")
              
            } else {
              store = gtkTreeStoreNew(types)
              combo <- gtkComboBoxNewWithModel(store)
              if(doIcons) {
                ## icon renderer
                cellrenderer = gtkCellRendererPixbufNew()
                combo$PackStart(cellrenderer, expand=FALSE)
                combo$AddAttribute(cellrenderer, "stock-id", 1)
              }
              ## pack in text -- not done if no entry
              cellrenderer = gtkCellRendererTextNew()
              combo$PackStart(cellrenderer, expand=TRUE)
              combo$AddAttribute(cellrenderer,"text", 0)
            }

            ## add tooltip if there
            if(ncol(items) >= 3) {
              ## no easy way. Thought that we could do the following, but
              ## it doesn't work. Tooltip is on combobox when not expanded for searching
              ## combo['has-tooltip'] <-TRUE
              ## gSignalConnect(combo, "query-tooltip", function(w, x, y, bool, tool, ...) {
              ##   ## look up text from w, x, y
              ##   tool$setText("text")
              ##   TRUE
              ## })
            }

            
            obj <- as.gWidgetsRGtk2(combo)
            
##             obj = new("gDroplistRGtk",block=combo,widget=combo, toolkit=toolkit)

##             tag(obj,"store") <- store
##             tag(obj,"combo") <- combo
##             tag(obj,"editable") <- editable

            tag(obj, "items") <- items
            tag(obj, "doIcons") <- doIcons
            tag(obj, "coerce.with") = coerce.with
            tag(obj, "default_fill") <- "x"
            
            ## load up the store
            if(length(items) > 0) {
              obj[] <- items
            }
            ## should I have actiirst be blank? Use 0 (to make -1) for this
            combo$Show()
            combo$SetActive(selected-1)

            ## set size if really small under windows
            if(.Platform$OS == "windows") {
              if(dim(items)[1] > 0) {
                colChars <- max(sapply(items[,1,drop=TRUE],nchar))
                if(colChars < 3)
                  combo['width-request'] <- 15*(4 + colChars)
              }
            }
            
            ## add drophandler -- switch if drop matches
            adddroptarget(obj, handler = function(h,...) {
              name = id(h$dropdata)
              theValues = obj[]
              if(!is.na(name) && !is.null(name) && name %in% theValues) {
                svalue(obj) <- name
              }
            })
            
            
            if (!is.null(container)) {
              if(is.logical(container) && container == TRUE)
                container = gwindow(visible=TRUE)
              add(container, obj,...)
            }

            ## pass in a size via width=
            if(!is.null(theArgs$width)) 
              size(obj) <- c(width=theArgs$width, height=0)
            
            if (!is.null(handler)) {
              id <- addhandlerchanged(obj, handler, action)
              tag(obj, "handler.id") <- id
            }
            
            invisible(obj)
          })

as.gWidgetsRGtk2.GtkComboBoxEntry <- function(widget,...) {
  obj <- .as.gWidgetsRGtk2.gdroplist(widget,...)
  tag(obj,"editable") <- TRUE
  return(obj)
}

as.gWidgetsRGtk2.GtkComboBox <- function(widget,...) {

  obj <- .as.gWidgetsRGtk2.gdroplist(widget,...)
  tag(obj,"editable") <- FALSE
  return(obj)

}

.as.gWidgetsRGtk2.gdroplist <- function(widget) {
  parent <- widget$parent
  if(is.null(parent)) {
    parent <- gtkAlignmentNew(xscale=1, yscale=0)
    parent$add(widget)
  }
  
  obj <- new("gDroplistRGtk",block=parent,widget=widget,
    toolkit=guiToolkit("RGtk2"))

  store <- widget$GetModel()
  tag(obj,"store") <- store
  tag(obj,"combo") <- widget

  ## get items then store. This is only useful for coerced values
  ## as otherwise set in constructor
  items <- c()

  iter <- store$GetIterFirst()
  if(is.logical(iter$retval) && iter$retval) {
    items <- store$GetValue(iter$iter,0)$value
    ret <- store$IterNext(iter$iter)
    while(ret) {
      items <- c(items,store$GetValue(iter$iter,0)$value)
      ret <- store$IterNext(iter$iter)
    }
  }
  
  tag(obj,"items") <- data.frame(items=items, stringsAsFactors=FALSE)
  
  return(obj)
}  

### methods
## value is for getting/setting the selected value
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gDroplistRGtk"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {
            ## add in an as.numeric flag, getwidget when editable
            theArgs = list(...)         # deprecated
            coerce.with = tag(obj, "coerce.with")
  
            ## do things depending on whether there is an entry or not
            ## if editable, then entry is widget and combo may be found by tag("combo")
            if(tag(obj,"editable")) {
              if(is.null(index) || index==FALSE) {
##                entry = obj@widget      # entry is widget
                entry = obj@widget$GetChild()
                if(!is.null(theArgs$getwidget)) {
                  gwCat("DEBUG: getwidget is deprecated\n")
                }
                if(!is.null(theArgs$as.numeric)) {
                  gwCat("DEBUG: as.numeric as an argument is deprected. Use coerce.with\n")
                }
                
                ## else we return text
                val = entry$GetText()

                coerce.with<-tag(obj,"coerce.with")
                if(is.null(coerce.with))
                  return(val)
                else if(is.function(coerce.with))
                  return(coerce.with(val))
                else if(is.character(coerce.with))
                  return(do.call(coerce.with,list(val)))
                else
                  warning("Error: coerce.with is a function or character")
              } else {
                ## return the index or NA
                combobox = tag(obj,"combo") # of obj@widget$GetParent()
                active = combobox$GetActive()
                if(active < 0)
                  return(NA)
                else
                  return(active+1)
              }
            } else {
              ## from pygtk manual
              combobox = obj@widget
              model = combobox$GetModel()
              selected = combobox$GetActive()
              items = obj[]

              ## selected is the index. It is 0 based
              if(selected < 0) {
                return(NULL)                      # none selected
              } else {
                ## do we return the index?
                if(!is.null(index) && index==TRUE) {
                  return(selected + 1)
                } else {
                  val = items[selected+1]
                  coerce.with<-tag(obj,"coerce.with")
                  if(is.null(coerce.with))
                    return(val)
                  else if(is.function(coerce.with))
                    return(coerce.with(val))
                  else if(is.character(coerce.with))
                    return(do.call(coerce.with,list(val)))
                  else
                    warning("Error: coerce.with is a function or character")
                }
              }
            }
          })

## set the displayed value to value
setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gDroplistRGtk"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   theArgs = list(...)

                   ##  if editable do differently
                   if(tag(obj,"editable")) {
                     if(is.null(index) || index == FALSE)  {
##                       entry = obj@widget
                       entry = obj@widget$GetChild()
                       entry$SetText(value)              # gtk Call
                     } else {
                       ## set the index
                       combobox = tag(obj,"combo") # or obj@widget$GetParent()
                       combobox$SetActive(value-1)
                     }
                   } else {
                     combobox = obj@widget
                     items = obj[]         # drops icons
                     if(!is.null(index) && as.logical(index)) { # either value or index is non-null
                       combobox$SetActive(value-1)
                     } else {
                       if(any(value == items)) {
                         combobox$SetActive(min(which(value==items)) - 1)
                       } else {
                         combobox$AppendText(value)
                         combobox$SetActive(length(items))
                       }
                     }
                   }
                   return(obj)
                 })

## the methods [ and [<- refer to the pre-defined values in the drop list.
## [
setMethod("[",
          signature(x="gDroplistRGtk"),
          function(x, i, j, ..., drop=TRUE) {
            .leftBracket(x, x@toolkit, i, j, ..., drop=drop)
          })
setMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gDroplistRGtk"),
          function(x, toolkit, i, j, ..., drop=TRUE) {

            items = tag(x,"items")
            if(missing(i))
              return(items[,1,drop=TRUE])
            else
              return(items[i,1,drop=TRUE])
          })


## replaces the values in droplist
## values is a vector of values -- not a dataframe
#set.values.gDropList = function(obj, values, ...) {
setReplaceMethod("[",
                 signature(x="gDroplistRGtk"),
                 function(x, i, j,..., value) {
                   .leftBracket(x, x@toolkit, i, j, ...) <- value
                   return(x)
                 })

setReplaceMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gDroplistRGtk"),
          function(x, toolkit, i, j, ..., value) {
            ## items is a data frame
            ## ncol=2 if doIcons, else 1
            olditems = tag(x,"items")

            old_value <- svalue(x)

            
            ## coerce value to a data frame, if not one
            if(!inherits(value,"data.frame"))
              value = as.data.frame(value, stringsAsFactors=FALSE)
            
            if(missing(i)) {
              items = value
              n = nrow(items)
              i = 1:n
            } else {
              items = olditems
              items[i,] = value
            }
            ## update items
            tag(x,"items") <- items
            
            ## now update widget
            store = tag(x,"store")
            store$Clear()
            n = nrow(items)

            doIcons = tag(x,"doIcons")
            allIcons = getStockIcons()
            if(n  > 0)  {
              for(j in 1:n) {
#                iter = store$Append(parent=NULL)
                iter = store$Append()
                store$SetValue(iter$iter, column = 0, items[j,1])
                if(doIcons)
                  store$SetValue(iter$iter, column = 1,
                                 allIcons[[as.character(items[j,2])]]) # convert to name
                if(ncol(items) >= 3) {
                  store$setValue(iter$iter, column =2, items[j,3])
                }
              }
            }

            ## set value if we can
            if(!is.null(old_value) &&
               old_value %in% items[,1, drop=TRUE])
              svalue(x) <- old_value
            
            return(x)
          })

setMethod("length",
          signature(x="gDroplistRGtk"),
          function(x) {
            .length(x, x@toolkit)
          })
setMethod(".length",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gDroplistRGtk"),
          function(x,toolkit) {
            return(length(x[]))
          })

###################################################
  
### handlers
setMethod(".addhandlerclicked",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gDroplistRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            addhandler(obj,"changed",handler,action,...)
          })

## want changed by activate -- or arrow for editable -- not keystroke
setMethod(".addhandlerchanged",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gDroplistRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            if(tag(obj,"editable")) {
              id = addhandler(obj,"changed",handler = function(h,...) {
                if(obj@widget$GetActive() != -1) {
                  handler(h,...)
                }
              },action)  # clicked -- not keystroke
              ## put handler on entry too
              gtktry(connectSignal(obj@widget$GetChild(),
                            signal="activate",
                            f=handler,
                            data=list(obj=obj, action=action,...),
                            user.data.first = TRUE,
                            after = FALSE),
                  silent=TRUE)
              invisible(id)
            } else {
              addhandler(obj,"changed",handler,action)
            }
          })

setMethod(".addhandlerkeystroke",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gDroplistRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            ## put handler on entry 
            gtktry(connectSignal(obj@widget$GetChild(),
                          signal="changed",
                          f=handler,
                          data=list(obj=obj, action=action,...),
                          user.data.first = TRUE,
                          after = FALSE),
                silent=TRUE)
          })
