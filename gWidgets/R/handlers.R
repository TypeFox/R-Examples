##' @imports guiComponents
##' 

############### removeHandler ###################################
##' generic to remove a handler
##' @alias removeHandler
setGeneric("removehandler",function(obj, ID=NULL, ...)
           standardGeneric("removehandler"))

##' base method to remove a handler
##' @alias removeHandler
setMethod("removehandler", signature("guiWidget"),
          function(obj, ID=NULL, ...) {
            .removehandler(obj@widget, obj@toolkit, ID, ...)
          })

##' dispatch to toolkit
##' @alias removeHandler
setGeneric(".removehandler",function(obj, toolkit, ID=NULL, ...)
           standardGeneric(".removehandler"))

##' generic to define method to remove a handler by ID
setGeneric("removeHandler",function(obj, ID=NULL, ...)
           standardGeneric("removeHandler"))

##' base method to remove a handler
##'
##' @param obj object with event handler set for it
##' @param ID ID of handler returned when handler is set. If NULL, all handlers are removed
##' @return NULL, called for side effect
##' @export
setMethod("removeHandler", signature("guiWidget"),
          function(obj, ID=NULL, ...) {
            .removehandler(obj@widget, obj@toolkit, ID, ...)
          })

############### blockHandler ###################################
##' Generic to block a handler from being called until block is unblocked
##' @alias blockHandler
setGeneric("blockhandler",function(obj, ID=NULL, ...)
           standardGeneric("blockhandler"))
##' base method for blocking a handler by ID
##' @alias blockHandler
setMethod("blockhandler", signature("guiWidget"),
          function(obj, ID=NULL, ...) {
            .blockhandler(obj@widget, obj@toolkit, ID, ...)
          })
##' method for toolkit dispatch
##' @alias blockHandler
setGeneric(".blockhandler",function(obj, toolkit, ID=NULL, ...)
           standardGeneric(".blockhandler"))

##' Generic to define method to block a handler from being called
##' 
##' @param obj object with event handler set for it
##' @param ID ID of handler returned when handler is set. If NULL, all handlers are blocked
##' @return NULL, called for side effect
##' @export
setGeneric("blockHandler",function(obj, ID=NULL, ...)
           standardGeneric("blockHandler"))

##' base method to block a handler from being called.
setMethod("blockHandler", signature("guiWidget"),
          function(obj, ID=NULL, ...) {
            .blockhandler(obj@widget, obj@toolkit, ID, ...)
          })

##################################################
##' Generic to define method to unblock a blocked handler
setGeneric("unblockhandler",function(obj, ID=NULL, ...)
           standardGeneric("unblockhandler"))

##' base method to unblock a blocked handler
setMethod("unblockhandler", signature("guiWidget"),
          function(obj, ID=NULL, ...) {
            .unblockhandler(obj@widget, obj@toolkit, ID, ...)
          })

##' method for toolkit dispatch
##' @alias unblockHandler
setGeneric(".unblockhandler",function(obj, toolkit, ID=NULL, ...)
           standardGeneric(".unblockhandler"))

##' Generic for unblocking a blocked handler
##'
##' @param obj object with event handler
##' @param ID id of handler set through addHandlerXXX call
##' @param ... ignored
##' @return NULL, called for side effect
##' @export
setGeneric("unblockHandler",function(obj, ID=NULL, ...)
           standardGeneric("unblockHandler"))
setMethod("unblockHandler", signature("guiWidget"),
          function(obj, ID=NULL, ...) {
            .unblockhandler(obj@widget, obj@toolkit, ID, ...)
          })


##################################################
## addhandler is now exported
setGeneric("addhandler",function(obj, signal, handler, action=NULL, ...) standardGeneric("addhandler"))
setMethod("addhandler",signature(obj="guiWidget"),
          function(obj, signal, handler, action=NULL, ...) {
            toolkit = obj@toolkit
            .addhandler(obj@widget, toolkit, handler, action, ...)
          })
## dispatch with toolkit
setGeneric(".addhandler",function(obj,  toolkit, signal, handler, action=NULL,...) standardGeneric(".addhandler"))

##' Add a handler using toolkit-specific signal
##'
##' This is not portable across toolkits, as the signal passed in a toolkit specific.
##' @param obj gWidgets object to get event handler
##' @param signal toolkit signal defining event
##' @param handler function to call when event occurs. Uses standard signature (h,...)
##' @param action used to parameterize handler call, passed in as \code{h$action}
##' @param ... ignored
##' @return an ID of the handler for blocking or removal
##' @export
setGeneric("addHandler",function(obj, signal, handler, action=NULL, ...) standardGeneric("addHandler"))

##
setMethod("addHandler",signature(obj="guiWidget"),
          function(obj, signal, handler, action=NULL, ...) {
            toolkit = obj@toolkit
            .addhandler(obj@widget, toolkit, signal, handler, action, ...)
          })

           

## addhandlerchanged
setGeneric("addhandlerchanged",function(obj, handler=NULL, action=NULL, ...) standardGeneric("addhandlerchanged"))
setMethod("addhandlerchanged",signature(obj="guiWidget"),
          function(obj, handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .addhandlerchanged(obj@widget, toolkit, handler, action, ...)
          })
## dispatch with toolkit
setGeneric(".addhandlerchanged",function(obj, toolkit,...) standardGeneric(".addhandlerchanged"))


##' Set handler for most typical event
##'
##' The "Changed" handler is set to be called for the most typical
##' event (arbitrarily defined of course). This event varies from
##' widget to widget. It is the same event that the widget
##' constructor's \code{handler} argument is called for.
setGeneric("addHandlerChanged",function(obj, handler=NULL, action=NULL, ...) standardGeneric("addHandlerChanged"))
setMethod("addHandlerChanged",signature(obj="guiWidget"),
          function(obj, handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .addhandlerchanged(obj@widget, toolkit, handler, action, ...)
          })



## addhandlerkeystroke
setGeneric("addhandlerkeystroke",function(obj, handler=NULL, action=NULL, ...) standardGeneric("addhandlerkeystroke"))
setMethod("addhandlerkeystroke",signature(obj="guiWidget"),
          function(obj, handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .addhandlerkeystroke(obj@widget, toolkit,handler, action, ...)
          })
## dispatch with toolkit
setGeneric(".addhandlerkeystroke",function(obj, toolkit,...) standardGeneric(".addhandlerkeystroke"))
#caps
setGeneric("addHandlerKeystroke",function(obj, handler=NULL, action=NULL, ...) standardGeneric("addHandlerKeystroke"))
setMethod("addHandlerKeystroke",signature(obj="guiWidget"),
          function(obj, handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .addhandlerkeystroke(obj@widget, toolkit,handler, action, ...)
          })



## addhandlerclicked
setGeneric("addhandlerclicked",function(obj, handler=NULL, action=NULL, ...) standardGeneric("addhandlerclicked"))
setMethod("addhandlerclicked",signature(obj="guiWidget"),
          function(obj, handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .addhandlerclicked(obj@widget, toolkit,handler, action, ...)
          })
## dispatch with toolkit
setGeneric(".addhandlerclicked",function(obj, toolkit,...) standardGeneric(".addhandlerclicked"))
## caps
setGeneric("addHandlerClicked",function(obj, handler=NULL, action=NULL, ...) standardGeneric("addHandlerClicked"))
setMethod("addHandlerClicked",signature(obj="guiWidget"),
          function(obj, handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .addhandlerclicked(obj@widget, toolkit,handler, action, ...)
          })



## addhandlerdoubleclick
setGeneric("addhandlerdoubleclick",function(obj, handler=NULL, action=NULL, ...) standardGeneric("addhandlerdoubleclick"))
setMethod("addhandlerdoubleclick",signature(obj="guiWidget"),
          function(obj, handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .addhandlerdoubleclick(obj@widget,toolkit,handler, action, ...)
          })
## dispatch with toolkit
setGeneric(".addhandlerdoubleclick",function(obj, toolkit,...) standardGeneric(".addhandlerdoubleclick"))
## caps
setGeneric("addHandlerDoubleclick",function(obj, handler=NULL, action=NULL, ...) standardGeneric("addHandlerDoubleclick"))
setMethod("addHandlerDoubleclick",signature(obj="guiWidget"),
          function(obj, handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .addhandlerdoubleclick(obj@widget,toolkit,handler, action, ...)
          })



## addhandlerrightclick
setGeneric("addhandlerrightclick",function(obj, handler=NULL, action=NULL, ...) standardGeneric("addhandlerrightclick"))
setMethod("addhandlerrightclick",signature(obj="guiWidget"),
          function(obj, handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .addhandlerrightclick(obj@widget,toolkit,handler, action, ...)
          })
## dispatch with toolkit
setGeneric(".addhandlerrightclick",function(obj,toolkit,...) standardGeneric(".addhandlerrightclick"))
## caps
setGeneric("addHandlerRightclick",function(obj, handler=NULL, action=NULL, ...) standardGeneric("addHandlerRightclick"))
setMethod("addHandlerRightclick",signature(obj="guiWidget"),
          function(obj, handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .addhandlerrightclick(obj@widget,toolkit,handler, action, ...)
          })


###
## Column clicks
## clicked
setGeneric("addhandlercolumnclicked",function(obj, handler=NULL, action=NULL, ...) standardGeneric("addhandlercolumnclicked"))
setMethod("addhandlercolumnclicked",signature(obj="guiWidget"),
          function(obj, handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .addhandlercolumnclicked(obj@widget, toolkit,handler, action, ...)
          })
## dispatch with toolkit
setGeneric(".addhandlercolumnclicked",function(obj, toolkit,...) standardGeneric(".addhandlercolumnclicked"))
## caps
setGeneric("addHandlerColumnClicked",function(obj, handler=NULL, action=NULL, ...) standardGeneric("addHandlerColumnClicked"))
setMethod("addHandlerColumnClicked",signature(obj="guiWidget"),
          function(obj, handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .addhandlercolumnclicked(obj@widget, toolkit,handler, action, ...)
          })



## addhandlerCOLUMNdoubleclick
setGeneric("addhandlercolumndoubleclick",function(obj, handler=NULL, action=NULL, ...) standardGeneric("addhandlercolumndoubleclick"))
setMethod("addhandlercolumndoubleclick",signature(obj="guiWidget"),
          function(obj, handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .addhandlercolumndoubleclick(obj@widget,toolkit,handler, action, ...)
          })
## dispatch with toolkit
setGeneric(".addhandlercolumndoubleclick",function(obj, toolkit,...) standardGeneric(".addhandlercolumndoubleclick"))
## caps
setGeneric("addHandlerColumnDoubleclick",function(obj, handler=NULL, action=NULL, ...) standardGeneric("addHandlerColumnDoubleclick"))
setMethod("addHandlerColumnDoubleclick",signature(obj="guiWidget"),
          function(obj, handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .addhandlercolumndoubleclick(obj@widget,toolkit,handler, action, ...)
          })

## columnrightclick


## addhandlerCOLUMNdoubleclick
setGeneric("addhandlercolumnrightclick",function(obj, handler=NULL, action=NULL, ...) standardGeneric("addhandlercolumnrightclick"))
setMethod("addhandlercolumnrightclick",signature(obj="guiWidget"),
          function(obj, handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .addhandlercolumnrightclick(obj@widget,toolkit,handler, action, ...)
          })
## dispatch with toolkit
setGeneric(".addhandlercolumnrightclick",function(obj, toolkit,...) standardGeneric(".addhandlercolumnrightclick"))
## caps
setGeneric("addHandlerColumnRightclick",function(obj, handler=NULL, action=NULL, ...) standardGeneric("addHandlerColumnRightclick"))
setMethod("addHandlerColumnRightclick",signature(obj="guiWidget"),
          function(obj, handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .addhandlercolumnrightclick(obj@widget,toolkit,handler, action, ...)
          })



## Selections
setGeneric("addhandlerselect",function(obj, handler=NULL, action=NULL, ...) standardGeneric("addhandlerselect"))
setMethod("addhandlerselect",signature(obj="guiWidget"),
          function(obj, handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .addhandlerselect(obj@widget,toolkit,handler, action, ...)
          })
## dispatch with toolkit
setGeneric(".addhandlerselect",function(obj,toolkit,...) standardGeneric(".addhandlerselect"))
## caps
setGeneric("addHandlerSelect",function(obj, handler=NULL, action=NULL, ...) standardGeneric("addHandlerSelect"))
setMethod("addHandlerSelect",signature(obj="guiWidget"),
          function(obj, handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .addhandlerselect(obj@widget, toolkit, handler, action, ...)
          })





## addhandlerFocus
setGeneric("addhandlerfocus",function(obj, handler=NULL, action=NULL, ...) standardGeneric("addhandlerfocus"))
setMethod("addhandlerfocus",signature(obj="guiWidget"),
          function(obj, handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .addhandlerfocus(obj@widget,toolkit,handler, action, ...)
          })
## dispatch with toolkit
setGeneric(".addhandlerfocus",function(obj,toolkit,...) standardGeneric(".addhandlerfocus"))
## caps
setGeneric("addHandlerFocus",function(obj, handler=NULL, action=NULL, ...) standardGeneric("addHandlerFocus"))
setMethod("addHandlerFocus",signature(obj="guiWidget"),
          function(obj, handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .addhandlerfocus(obj@widget,toolkit,handler, action, ...)
          })



## addhandlerblur
setGeneric("addhandlerblur",function(obj, handler=NULL, action=NULL, ...) standardGeneric("addhandlerblur"))
setMethod("addhandlerblur",signature(obj="guiWidget"),
          function(obj, handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .addhandlerblur(obj@widget,toolkit,handler, action, ...)
          })
## dispatch with toolkit
setGeneric(".addhandlerblur",function(obj,toolkit,...) standardGeneric(".addhandlerblur"))
## caps
setGeneric("addHandlerBlur",function(obj, handler=NULL, action=NULL, ...) standardGeneric("addHandlerBlur"))
setMethod("addHandlerBlur",signature(obj="guiWidget"),
          function(obj, handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .addhandlerblur(obj@widget,toolkit,handler, action, ...)
          })




## addhandlerdestroy
setGeneric("addhandlerdestroy",function(obj, handler=NULL, action=NULL, ...) standardGeneric("addhandlerdestroy"))
setMethod("addhandlerdestroy",signature(obj="guiWidget"),
          function(obj, handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .addhandlerdestroy(obj@widget, toolkit,handler, action, ...)
          })
## dispatch with toolkit
setGeneric(".addhandlerdestroy",function(obj,toolkit,...) standardGeneric(".addhandlerdestroy"))
##caps
setGeneric("addHandlerDestroy",function(obj, handler=NULL, action=NULL, ...) standardGeneric("addHandlerDestroy"))
setMethod("addHandlerDestroy",signature(obj="guiWidget"),
          function(obj, handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .addhandlerdestroy(obj@widget, toolkit,handler, action, ...)
          })


# addhandlerexpose
setGeneric("addhandlerexpose",function(obj, handler=NULL, action=NULL, ...) standardGeneric("addhandlerexpose"))
setMethod("addhandlerexpose",signature(obj="guiWidget"),
          function(obj, handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .addhandlerexpose(obj@widget,toolkit,handler, action, ...)
          })
## dispatch with toolkit
setGeneric(".addhandlerexpose",function(obj, toolkit,...) standardGeneric(".addhandlerexpose"))
## caps
setGeneric("addHandlerExpose",function(obj, handler=NULL, action=NULL, ...) standardGeneric("addHandlerExpose"))
setMethod("addHandlerExpose",signature(obj="guiWidget"),
          function(obj, handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .addhandlerexpose(obj@widget,toolkit,handler, action, ...)
          })

# addhandlerunrealize
##' handler when window is unrealized
##'
##' For gwindow objects this handler is called before the window is closed. If this handler
##' returns \code{TRUE} the window will be close, if \code{FALSE} the window will not be closed.
##' @param obj gWidget object
##' @param handler function with signature (h,...) to call when event occurs
##' @param action value passed to handler function in component \code{h$action}
##' @param ... ignored
##' @return ID of handler, used for blocking or removing
##' @export
##' @examples
##' \dontrun{
##' w <- gwindow("Really close?")
##' addHandlerUnrealize(w, handler=function(h,...) {
##' gconfirm("Really close?", parent=w)
##' })
##' }
setGeneric("addhandlerunrealize",function(obj, handler=NULL, action=NULL, ...) standardGeneric("addhandlerunrealize"))
setMethod("addhandlerunrealize",signature(obj="guiWidget"),
          function(obj, handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .addhandlerunrealize(obj@widget,toolkit,handler, action, ...)
          })
## dispatch with toolkit
setGeneric(".addhandlerunrealize",function(obj, toolkit,...) standardGeneric(".addhandlerunrealize"))
## caps
setGeneric("addHandlerUnrealize",function(obj, handler=NULL, action=NULL, ...) standardGeneric("addHandlerUnrealize"))
setMethod("addHandlerUnrealize",signature(obj="guiWidget"),
          function(obj, handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .addhandlerunrealize(obj@widget,toolkit,handler, action, ...)
          })




## mousemotion
##' Handler to respond to mouse motion
##'
##' Movement of mouse over widget triggers this handler
##' @param obj gWidget object
##' @param handler function with signature (h,...) to call when event occurs
##' @param action value passed to handler function in component \code{h$action}
##' @param ... ignored
##' @return ID of handler, used for blocking or removing
##' @export
setGeneric("addhandlermousemotion",function(obj, handler=NULL, action=NULL, ...) standardGeneric("addhandlermousemotion"))
setMethod("addhandlermousemotion",signature(obj="guiWidget"),
          function(obj, handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .addhandlermousemotion(obj@widget,toolkit,handler, action, ...)
          })
## dispatch with toolkit
setGeneric(".addhandlermousemotion",function(obj, toolkit,...) standardGeneric(".addhandlermousemotion"))
## caps
setGeneric("addHandlerMouseMotion",function(obj, handler=NULL, action=NULL, ...) standardGeneric("addHandlerMouseMotion"))
setMethod("addHandlerMouseMotion",signature(obj="guiWidget"),
          function(obj, handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .addhandlermousemotion(obj@widget,toolkit,handler, action, ...)
          })




# addhandleridle
setGeneric("addhandleridle",function(obj, handler=NULL, action=NULL, interval=1000, ...) standardGeneric("addhandleridle"))
setMethod("addhandleridle",signature(obj="guiWidget"),
          function(obj, handler=NULL, action=NULL, interval=1000, ...) {
            toolkit = obj@toolkit
            .addhandleridle(obj@widget, toolkit, handler=handler, action=action, interval=interval, ...)
          })
## dispatch with toolkit
setGeneric(".addhandleridle",function(obj, toolkit,handler=NULL,action=NULL,
                                      interval=1000,...) standardGeneric(".addhandleridle"))
## caps
setGeneric("addHandlerIdle",function(obj, handler=NULL, action=NULL, interval=1000, ...) standardGeneric("addHandlerIdle"))
setMethod("addHandlerIdle",signature(obj="guiWidget"),
          function(obj, handler=NULL, action=NULL, interval=1000, ...) {
            toolkit = obj@toolkit
            .addhandleridle(obj@widget, toolkit, handler=handler, action=action, interval=interval, ...)
          })

## addpopupmenu
setGeneric("addpopupmenu",function(obj, menulist, action=NULL, ...) standardGeneric("addpopupmenu"))
setMethod("addpopupmenu",signature(obj="guiWidget"),
          function(obj, menulist, action=NULL, ...) {
            toolkit = obj@toolkit
            .addpopupmenu(obj@widget, toolkit,menulist, action, ...)
          })
## dispatch with toolkit
setGeneric(".addpopupmenu",function(obj, toolkit, menulist, action=NULL, ...) standardGeneric(".addpopupmenu"))
## caps
setGeneric("addPopupmenu",function(obj, menulist, action=NULL, ...) standardGeneric("addPopupmenu"))
setMethod("addPopupmenu",signature(obj="guiWidget"),
          function(obj, menulist, action=NULL, ...) {
            toolkit = obj@toolkit
            .addpopupmenu(obj@widget, toolkit,menulist, action, ...)
          })

## add3rdmousepopupmenu
setGeneric("add3rdmousepopupmenu",function(obj, menulist, action=NULL,  ...) standardGeneric("add3rdmousepopupmenu"))
setMethod("add3rdmousepopupmenu",signature(obj="guiWidget"),
          function(obj, menulist, action=NULL,  ...) {
            .add3rdmousepopupmenu(obj@widget, obj@toolkit, menulist,
                                  action=action,  ...)
          })
## dispatch with toolkit
setGeneric(".add3rdmousepopupmenu",function(obj, toolkit,menulist, action=NULL, ...) standardGeneric(".add3rdmousepopupmenu"))
## caps
setGeneric("add3rdMousePopupmenu",function(obj, menulist, action=NULL,  ...) standardGeneric("add3rdMousePopupmenu"))
setMethod("add3rdMousePopupmenu",signature(obj="guiWidget"),
          function(obj, menulist, action=NULL,  ...) {
            .add3rdmousepopupmenu(obj@widget, obj@toolkit, menulist,
                                  action=action,  ...)
          })

## adddropsource
setGeneric("adddropsource",function(obj, targetType="text",
                                    handler=NULL, action=NULL, ...) standardGeneric("adddropsource"))
setMethod("adddropsource",signature(obj="guiWidget"),
          function(obj, targetType="text", handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .adddropsource(obj@widget, toolkit,targetType=targetType,
                           handler=handler, action=action, ...)
          })
## dispatch with toolkit
setGeneric(".adddropsource",function(obj, toolkit,targetType="text", handler=NULL, action=NULL, ...) standardGeneric(".adddropsource"))
## caps
setGeneric("addDropSource",function(obj, targetType="text",
                                    handler=NULL, action=NULL, ...) standardGeneric("addDropSource"))
setMethod("addDropSource",signature(obj="guiWidget"),
          function(obj, targetType="text", handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .adddropsource(obj@widget, toolkit,targetType=targetType,
                           handler=handler, action=action, ...)
          })


## adddropmotion
setGeneric("adddropmotion",function(obj, handler=NULL, action=NULL, ...)
           standardGeneric("adddropmotion"))
setMethod("adddropmotion",signature(obj="guiWidget"),
          function(obj,  handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .adddropmotion(obj@widget, toolkit,
                           handler=handler, action=action, ...)
          })
## dispatch with toolkit
setGeneric(".adddropmotion",function(obj, toolkit, handler=NULL, action=NULL, ...) standardGeneric(".adddropmotion"))
## caps
setGeneric("addDropMotion",function(obj, handler=NULL, action=NULL, ...)
           standardGeneric("addDropMotion"))

setMethod("addDropMotion",signature(obj="guiWidget"),
          function(obj,  handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .adddropmotion(obj@widget, toolkit,
                           handler=handler, action=action, ...)
          })

## adddroptarget
setGeneric("adddroptarget",function(obj, targetType="text",
                                    handler=NULL, action=NULL, ...) standardGeneric("adddroptarget"))
setMethod("adddroptarget",signature(obj="guiWidget"),
          function(obj, targetType="text", handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .adddroptarget(obj@widget, toolkit,targetType=targetType,
                           handler=handler, action=action, ...)
          })
## dispatch with toolkit
setGeneric(".adddroptarget",function(obj, toolkit,targetType="text", handler=NULL, action=NULL, ...) standardGeneric(".adddroptarget"))
## caps
setGeneric("addDropTarget",function(obj, targetType="text",
                                    handler=NULL, action=NULL, ...) standardGeneric("addDropTarget"))
setMethod("addDropTarget",signature(obj="guiWidget"),
          function(obj, targetType="text", handler=NULL, action=NULL, ...) {
            toolkit = obj@toolkit
            .adddroptarget(obj@widget, toolkit,targetType=targetType,
                           handler=handler, action=action, ...)
          })
