##' @include guiWidget.R
##'

##' A subclass of guiWidget to hold containers
setClass("guiContainer",
         contains="guiWidget",
         prototype=prototype(new("guiWidget"))
         )


## define a subclass
setClass("gContainerANY",
         contains="gWidgetANY",
         prototype=prototype(new("gWidgetANY"))
         )

################# methods #################################


################## add ################################

##' Generic for adding child widget to parent container
setGeneric("add",function(obj,value, ...) standardGeneric("add"))

##' method for adding child to a container
setMethod("add",signature(obj="guiContainer"),
          function(obj, value, ...) {
            toolkit = obj@toolkit
            ladd <- function(obj, value, ...,do.newline, font.attr, where)
              .add(obj@widget, obj@toolkit, value, ...)
            ladd(obj, value,...)
          })

## deprecated
##' 
## setMethod("add",signature(obj="guiWidget"),
##           function(obj, value, ...) {
##             toolkit = obj@toolkit
##             .add(obj@widget, toolkit,value,...)
##           })


## setMethod("add",signature(obj="guiComponent"),
##           function(obj, value, ...) {
##             toolkit = obj@toolkit
##             ladd <- function(obj,value, ..., label, name, override.closebuttons, index, pageno, markup, anchor, expand) .add(obj@widget,obj@toolkit,value,...)
##              ladd(obj,value,...)
##           })


##' dispatch with toolkit
##' @alias add
setGeneric(".add",function(obj, toolkit,value,...) standardGeneric(".add"))

##' Generic to delete child widget from parent
setGeneric("delete",function(obj,widget, ...) standardGeneric("delete"))
##' Delete method for containers
setMethod("delete",signature(obj="guiContainer"),
          function(obj, widget, ...) {
            toolkit = obj@toolkit
            .delete(obj@widget,toolkit,widget,...)
          })
##' dispatch with toolkit
##' @alias delete
setGeneric(".delete",function(obj, toolkit,widget,...) standardGeneric(".delete"))




## dispose
setGeneric("dispose",function(obj, ...) standardGeneric("dispose"))
## add generic for Containers and sometimes widgets
setMethod("dispose",signature(obj="guiWidget"),
          function(obj, ...) {
            toolkit = obj@toolkit
            .dispose(obj@widget, toolkit,...)
          })
## dispatch with toolkit
setGeneric(".dispose",function(obj,toolkit,...) standardGeneric(".dispose"))

