## Command objects represent a high-level operation

setClass("Command", representation(active = "logical"), prototype(active = TRUE))

### possibly supported methods
## rev(x)
## setGeneric("active", function(object, ...) standardGeneric("active"))
## setMethod("active", "Command", function(object) object@active)

## setGeneric("active<-", function(object, ... ,value) standardGeneric("active<-"))
## setReplaceMethod("active", "Command", function(object, value){
##   object@active <- value
##   object
## })
## active(), active<-()

## returns a widget for controlling and viewing this object
## setGeneric("widget", function(object, ...) standardGeneric("widget"))

## returns a widget containing an interactive visualization of the
## specified data in the context of the specified protocol
setGeneric("explore", function(object, protocol, ...)
           standardGeneric("explore"))

## how an object is identified in a user interface
setGeneric("displayName", function(object, ...) standardGeneric("displayName"))
## by default, the name of the object
setMethod("displayName", "ANY", function(object) class(object))
