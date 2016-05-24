## Some reference classes for clearing up programming

##' Observable class sets up objects that can be observed. Inherited by Model
Observable <- setRefClass("Observable",
                          fields=list(
                            ..observers="list",
                            ..blocked_observers = "list",
                            ..blocked="logical"
                            ),
                          methods=list(
                            add_observer=function(o, signal="DEFAULT") {
                              "Add an observer. Return id for block/remove/..."
                              if(!is(o, "Observer"))
                                stop("Not an observer")
                              l <- ..observers
                              if(is.null(l[[signal]]))
                                l[[signal]] <- list(o)
                              else
                                l[[signal]] <- c(l[[signal]], o)
                              ..observers <<- l
                              list(signal=signal, o=o)
                            },
                            remove_observer=function(id) {
                              "Remove observer"
                              if(!is(id$o, "Observer"))
                                stop("Call with an observer id")
                              
                              signal <- id$signal
                              ind <- lapply(..observers[[signal]], function(i) identical(i, id$o))
                              if(any(unlist(ind)) )
                                ..observers[[signal]][[which(ind)]] <<- NULL
                              
                            },
                            block_observer=function(id) {
                              "Block observers. If o missing, block all"
                              if(missing(id) || is.null(id)) {
                                ..blocked <<- TRUE
                              } else {
                                if(is.null(..blocked_observers[[id$signal]]))
                                  ..blocked_observers[[id$signal]] <<- list(id$o)
                                else
                                  ..blocked_observers[[id$signal]] <<-
                                    c(..blocked_observers[[id$signal]], o)
                              }
                            },
                            unblock_observer=function(id) {
                              "Unblock observer. If id missing, unblock global block"
                              if(missing(id) || is.null(id)) {
                                ..blocked <<- FALSE
                              } else {
                                signal <- id$signal
                                ind <- lapply(..blocked_observers[[signal]], function(i) identical(i, id$o))
                                if(any(unlist(ind))) 
                                  ..blocked_observers[[signal]][[which(ind)]] <<- NULL
                              }
                            },
                            notify_observers=function(..., signal="DEFAULT") {
                              "Call each non-blocked observer"
                              if(length(..blocked) && ..blocked)
                                return()
                              lapply(..observers[[signal]], function(o) {
                                ind <- lapply(..blocked_observers[[signal]], function(i) identical(i, o))
                                if(!any(unlist(ind))) 
                                  o$update(...)
                              })
                            }
                            )
                          )

##' Observer class is used to observe an observable
Observer <- setRefClass("Observer",
                        fields=list(
                          o = "ANY",    # want "function", but doesn't work with proto objects
                          obj="ANY",
                          action="ANY"
                          ),
                        methods=list(
                          initialize=function(o,  obj, action=NULL) {
                            o <<- o
                            obj <<- obj
                            action <<- action
                            .self
                          },
                          update=function(...) {
                            "Call self."
                            h <- list(obj=obj, action=action)
                            o(h, ...)
                          }
                          )
                        )

##' Base class for widgets -- just a widget and block, but could
##' put much more here later
GWidgetGtk <- setRefClass("GWidgetGtk",
                          contains="Observable",
                          fields=list(
                            widget="ANY",
                            block="ANY"
                            ),
                          methods=list(
                            getWidget=function() widget,
                            getBlock=function() block
                            )
                          )


## Base class for widgets using a reference class, as gradio
## This should be moved elsewhere -- once we have more than one
setClass("gComponentWithRefClassRGtk",
         representation=representation(ref_widget="Observable"),
         contains="gComponentRGtk",
         prototype=prototype(new("gComponentRGtk"))
         )


## Some methods
setMethod(".removehandler",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gComponentWithRefClassRGtk"),
          function(obj, toolkit, ID=NULL, ...) {
            obj@ref_widget$remove_observer(ID)
          })

setMethod(".blockhandler",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gComponentWithRefClassRGtk"),
          function(obj, toolkit, ID=NULL, ...) {
            obj@ref_widget$block_observer(ID)
          })

setMethod(".unblockhandler",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gComponentWithRefClassRGtk"),
          function(obj, toolkit, ID=NULL, ...) {
            obj@ref_widget$unblock_observer(ID)
          })


##################################################
## Base class for gradio, gcheckboxgroup
setClass("gComponentWithRefClassWithItemsRGtk",
         contains="gComponentWithRefClassRGtk",
         prototype=prototype(new("gComponentWithRefClassRGtk"))
         )

##' enabled <- enables/disable the block
##' 
##' @param value a logical
setReplaceMethod(".enabled",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gComponentWithRefClassWithItemsRGtk"),
                 function(obj, toolkit, ..., value) {
                   block <- obj@ref_widget$block
                   block['sensitive'] <- as.logical(value)
                   return(obj)
                 })

##' visible<- hides or shows the block
##'
##' @param value a logical
setReplaceMethod(".visible",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gComponentWithRefClassWithItemsRGtk"),
                 function(obj, toolkit, ..., value) {
                   block <- obj@ref_widget$block
                   if(as.logical(value))
                     block$show()
                   else
                     block$hide()
                   
                   return(obj)
                 })

