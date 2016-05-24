##' @include guiWidget.R
##'
##' 
##' A widget subclass for basic components
setClass("guiComponent",
         contains=c("guiWidget"),
         prototype=prototype(new("guiWidget"))
         )

##' define a subclass
setClass("gComponentANY",
         contains=c("gWidgetANY"),
         prototype=prototype(new("gWidgetANY"))
         )


############### methods ###################################

################ add ##################################
##' method for adding to ghelp, etc
setMethod("add",signature(obj="guiComponent"),
          function(obj, value, ...) {
            .add(obj@widget, obj@toolkit, value, ...)
          })


############### focus ###################################
##' generic to check if a widget has the focus
setGeneric("focus",function(obj, ...) standardGeneric("focus"))

##' method to check is a component has the focus 
setMethod("focus",signature(obj="guiWidget"),
          function(obj, ...) {
            toolkit = obj@toolkit
            .focus(obj@widget, toolkit,...)
          })
##' dispatch with toolkit
##' @alias focus
setGeneric(".focus",function(obj, toolkit,...) standardGeneric(".focus"))


############### focus<- ###################################

##' generic to set focus on a widget
setGeneric("focus<-",function(obj, ..., value) standardGeneric("focus<-"))

##' method for setting focus on a component
setMethod("focus<-",signature(obj="guiWidget"),
          function(obj, ..., value) {
            toolkit <- obj@toolkit
            .focus(obj@widget, toolkit,...) <- value
            return(obj)
          })
##' dispatch with toolkit
##' @alias focus<-
setGeneric(".focus<-",function(obj, toolkit,...,value)
           standardGeneric(".focus<-"))



############### tooltip<- ###################################
##' generic to set a tooltip for a component
setGeneric("tooltip<-",function(obj, ..., value) standardGeneric("tooltip<-"))

##' base method for setting a toolktip
setMethod("tooltip<-",signature(obj="guiComponent"),
          function(obj, ..., value) {
            toolkit <- obj@toolkit
            .tooltip(obj@widget, toolkit,...) <- value
            return(obj)
          })
##' dispatch with toolkit
##' @alias tooltip<-
setGeneric(".tooltip<-",function(obj, toolkit,...,value)
           standardGeneric(".tooltip<-"))


############## undo ####################################
##' Generic to implement undo feature in widgets
setGeneric("undo",function(obj, ...) standardGeneric("undo"))

##' base method for undo
setMethod("undo",signature(obj="guiComponent"),
          function(obj, ...) {
            toolkit = obj@toolkit
            .undo(obj@widget, toolkit,...)
          })
##' dispatch with toolkit
##' @alias undo
setGeneric(".undo",function(obj,toolkit,...) standardGeneric(".undo"))

############## redo ####################################
##' Generic to implement redo feature in widgets
setGeneric("redo",function(obj, ...) standardGeneric("redo"))

##' Base method for redo
setMethod("redo",signature(obj="guiComponent"),
          function(obj, ...) {
            toolkit = obj@toolkit
            .redo(obj@widget, toolkit,...)
          })
##' dispatch with toolkit
##' @alias redot
setGeneric(".redo",function(obj,toolkit,...) standardGeneric(".redo"))


