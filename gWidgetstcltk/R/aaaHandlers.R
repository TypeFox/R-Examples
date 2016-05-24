## Handler code
## redid 7/2010

##' run handlers for this signal
##'
##' @param obj is gWidgets object
##' @param signal the signal (handler list keyed by signal)
##' @param h list with proper components from call
runHandlers <- function(obj, signal, h, ...) {
  ## check if enabled
  W <- getWidget(obj)
  if(isTtkWidget(W))
    enabled <- enabled_ttkwidget(W)
  else
    enabled <- enabled_tkwidget(W)

  if(enabled) {
    l <- tag(obj, "..handlers")
    signalList <- l[[signal]]      # first run last?
    lapply(signalList, function(i) {
      if(!i$blocked) {
        i$handler(h, ...)
      }
    })
  }
}

##' add a handler to an object
##'
##' The basic idea is that a list of handlers (keyed by the signal) is kept along with a flag
##' indicating whether the handler is blocked or not
##' The binding is done to call the runHandlers function so that this flag can be intercepted
##' For signal="command" we use the command option of the widget, otherwise we bind with tkbind
setMethod(".addHandler",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
          function(obj, toolkit,
                   signal, handler, action=NULL, ...) {
            ## use tkbind

            l <- tag(obj, "..handlers")
            if(is.null(l))
              l <- list()

            signalList <- l[[signal]]
            if(is.null(signalList))
              signalList <- list()


            ## each component of signalList is a list with ID, blocked, handler, action=action
            id <- digest(Sys.time())    # some unique key
            hList <- list(ID=id,
                          blocked=FALSE,
                          handler=handler,
                          action=action)
            
            signalList[[length(signalList) + 1]] <- hList # append
            l[[signal]] <- signalList
            tag(obj, "..handlers") <- l

            id <- list(id=id, signal=signal) # need this to block/remove/unblock
            
            theArgs = list(...)

            actualobj <- getWithDefault(theArgs$actualobj, obj)

            ## theArgs may have an extra with name=key, value
            FUN <- theArgs$FUN
            handler <- force(handler)

            if(is.null(FUN)) {
              FUN <- function(...) {
                h = list(
                  obj=actualobj,
                  action=action)
                
                runHandlers(obj, signal, h, ...)
              }
            }

            if(signal == "command")
              tkconfigure(getWidget(obj), command=FUN)
            else
              tkbind(getWidget(obj), signal, FUN)

            ## return id
            invisible(id)
          })

## for tcltk objects
setMethod(".addHandler",
          signature(toolkit="guiWidgetsToolkittcltk",obj="tcltkObject"),
          function(obj, toolkit, signal, handler, action=NULL, ...) {
            
            theArgs = list(...)
            theobj = 

            ## tcltk object
            tcltk_obj <- obj
            obj <- theArgs$actualobj

            ## copied from above
            l <- tag(obj, "..handlers")
            if(is.null(l))
              l <- list()

            signalList <- l[[signal]]
            if(is.null(signalList))
              signalList <- list()


            ## each component of signalList is a list with ID, blocked, handler, action=action
            id <- digest(Sys.time())    # some unique key
            hList <- list(ID=id,
                          blocked=FALSE,
                          handler=handler,
                          action=action)
            
            signalList[[length(signalList) + 1]] <- hList # append
            l[[signal]] <- signalList
            tag(obj, "..handlers") <- l

            id <- list(id=id, signal=signal) # need this to block/remove/unblock

            ## add handler
            handler <- force(handler)
            FUN <- theArgs$FUN
            if(is.null(FUN)) {
              FUN <- function(...) {
                ## check if enabled
                if(isTtkWidget(tcltk_obj))
                  enabled <- enabled_ttkwidget(tcltk_obj)
                else
                  enabled <- enabled_tkwidget(tcltk_obj)
                
                if(enabled) {
                  h = list(obj=obj, action=action)
                  handler(h,...)
                }
              }
            }
          

            if(signal == "command")
              tkconfigure(tcltk_obj, command=FUN)
            else
              tkbind(tcltk_obj, signal, FUN)

              ## return
            invisible(id)


          })

##' idle handler is different == have this hack to keep calling "after"
setMethod(".addhandleridle",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
          function(obj, toolkit,
                   handler=NULL, action=NULL, interval=1000, ...) {

            signal <- "idle"
            
            l <- tag(obj, "..handlers")
            if(is.null(l))
              l <- list()
            
            signalList <- l[[signal]]
            if(is.null(signalList))
              signalList <- list()
            
            
            ## each component of signalList is a list with ID, blocked, handler, action=action
            id <- digest(Sys.time())    # some unique key
            hList <- list(ID=id,
                          blocked=FALSE,
                          handler=handler,
                          action=action)
            
            signalList[[length(signalList) + 1]] <- hList # append
            l[[signal]] <- signalList
            tag(obj, "..handlers") <- l
            
            
            ## use tcl("apply",time, function) in a while loop here
            h = list()
            h$obj=obj
            h$action=action
            
            f <- function() {
              if(!windowExists(obj))
                return() # otherwise, issue when destroyed

              l <- tag(obj, "..handlers")
              sigList <- l[['idle']]
              ind <- sapply(sigList, function(i) i$ID == id)
              if(any(ind)) {
                if(!sigList[[which(ind)]]$blocked)
                  sigList[[which(ind)]]$handler(h)
                tcl("after", interval, f)
              }
            }
            ## start it off
            f()
            ## ID
            id <- list(id=id, signal="idle")
          })


##' Function to call to update the "blocked" flag on a handler. runHandlers consults this
##' before making the call
.blockUnblock <- function(obj, ID, block=TRUE, ...) {
  l <- tag(obj, "..handlers")

  if(is.null(ID)) {
    ## do all IDS
    lapply(names(l), function(signal) {
      sigList <- l[[signal]]
      if(length(sigList)) {
        for(i in sigList)
          .blockUnblock(obj, list(id=i$ID, signal=signal), block, ...)
      }
    })
    return()
  } else if(is.null(ID$id) && !is.null(ID[[1]]$obj)) {
    ## might be a list of IDs (gradio, gcheckboxgroup), we check here
    lapply(ID, function(i) {
      .blockUnblock(i$obj, i$id, block)
    })
    return()
  } else {
    ## single ID
    id <- ID$id
    signal <- ID$signal
    if(is.null(id) || is.null(signal))
      return()
    
    if(is.null(l[[signal]]))
      return()                 # no signal list
    ind <- sapply(l[[signal]], function(i) {
      i$ID == id
    })
    
    if(!any(ind))
      return()   
    
    for(i in which(ind)) {
      l[[signal]][[i]]$blocked <- block
    }
  }
  
  tag(obj, "..handlers") <- l
}

##' call to block a handler by ID. If ID=NULL, all handlers are blocked
setMethod(".blockhandler",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
          function(obj, toolkit, ID=NULL, ...) {
            .blockUnblock(obj, ID, block=TRUE)
            invisible()
          })

setMethod(".blockhandler",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gComponentR5tcltk"),
          function(obj, toolkit, ID=NULL, ...) {
            #widget <- tag(obj, "widget")
            widget <- obj@R5widget            
            widget$block_handler(ID)
          })

##' call to unblock a handler by ID. If ID=NULL, all handlers are unblocked
setMethod(".unblockhandler",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
          function(obj, toolkit, ID=NULL, ...) {
            .blockUnblock(obj, ID, block=FALSE)
            invisible()
          })


setMethod(".unblockhandler",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gComponentR5tcltk"),
          function(obj, toolkit, ID=NULL, ...) {
#            widget <- tag(obj, "widget")
            widget <- obj@R5widget            
            widget$unblock_handler(ID)
          })

##' method to remove a handler
setMethod(".removehandler",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
          function(obj, toolkit, ID=NULL, ...) {

            if(is.null(ID)) {
              ## remove all handlers. Get id, signal then call this recursively
              l <- tag(obj, "..handlers")
              lapply(names(l), function(signal) {
                sigList <- l[[signal]]
                if(length(sigList)) {
                  for(i in sigList)
                    .removehandler(obj, toolkit, ID=list(id=i$ID, signal=signal))
                }
              })
            } else if(is.null(ID$id) && !is.null(ID[[1]]$obj)) {
              ## might be a list of IDs (gradio, gcheckboxgroup), we check here
              lapply(ID, function(i) {
                removehandler(i$obj, ID=i)
              })
              return()
            } else {
              ## single ID
              ## ID here has two components: id, signal
              id <- ID$id
              signal <- ID$signal
              
              if(is.null(id) || is.null(signal))
                return()
              
              l <- tag(obj, "..handlers")
              if(is.null(l[[signal]]))
                return()                 # no signal list
              ind <- sapply(l[[signal]], function(i) {
                i$ID == id
              })
              
              if(!any(ind))
                return()                  # no match on id

              ## remove list that stores the handler
              for(i in which(ind))
                l[[signal]][[i]] <- NULL
              ## now save
              tag(obj, "..handlers") <- l
            }
          })
          
setMethod(".removehandler",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gComponentR5tcltk"),
          function(obj, toolkit, ID=NULL, ...) {
            ##widget <- tag(obj, "widget")
            widget <- obj@R5widget
            widget$remove_handler(ID)
          })


