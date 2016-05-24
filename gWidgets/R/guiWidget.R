##' @include toolkitClasses.R
##' 

##' Basic class for gWidgets objects
setClass("guiWidget",
         representation(
                        toolkit="guiWidgetsToolkit",
                        widget="ANY"  # could be RGtkObject, TclObject,....
                        ),
         prototype(
                   toolkit =guiToolkit(),
                   widget=NULL
                   )
         ) 
           
##' Subclass of gWidgets with NULL possibility
setClassUnion("guiWidgetOrNULL",
              c("NULL","guiWidget"))

##' for ANY toolkit
##'
##' @export
setClass("gWidgetANY",
         representation(
                        toolkit="guiWidgetsToolkit",
                        widget="ANY",  # could be RGtkObject, TclObject,....
                        block="ANY",
                        ID = "numeric"
                        ),
         prototype(
                   toolkit =guiToolkit(),
                   widget=NULL,
                   block = NULL,
                   ID = getNewID()
                   )
         ) 
           


##################################################
## methods

################ show ##################################
##' show method for guiWidget objects
##'
##' @param object guiWidget object to show
##' @return NULL, called to print summary of object
##' @export
setMethod("show",signature(object="guiWidget"),
          function(object) {
            cat("guiWidget of type:",
                class(object@widget),
                "for toolkit:",
                class(object@toolkit),
                "\n")
          })

################## svalue ################################
##' generic for selected value, svalue
##'
##' The \code{svalue} method returns the main property of a widget. This of course is context specific.
##' @param obj guiWidget object
##' @param index for selection widgets this is set to decide if value or index should be returned
##' @param drop for selection widgets do we return data frame or
##' vector (if possible). For text widgets, do we return entire text
##' or just selected text
##' @return main value of object
##' @export
setGeneric("svalue",function(obj, index=NULL, drop=NULL,...) standardGeneric("svalue"))

##' base svalue method
##'
##' @inheritParams svalue
##' @export
setMethod("svalue",signature(obj="guiWidget"),
          function(obj, index=NULL, drop=NULL, ... ) {
            toolkit = obj@toolkit
            val <- .svalue(obj@widget, toolkit, ...,index=index, drop=drop)

            if(is.logical(index) && index)
              return(as.integer(val))
            
            ## do we have coercewith?
            ## we check lots of ways, these should just be slots, someday...
            coercewith <- NULL
            if(methods::.hasSlot(obj@widget, "coercewith"))
              coercewith <- obj@widget@coercewith

            if(is.null(coercewith))
              coercewith <- tag(obj, "coercewith")

            if(is.null(coercewith))
              coercewith <- tag(obj, "coerce.with")

            ## now return
            if(is.null(coercewith))
              return(val)

            if(is.character(coercewith))
              coercewith <- get(coercewith, inherits=TRUE)
            
            if(!is.function(coercewith))
              return(val)

            ## return coerced value
            coercewith(val)
          })

##' svalue method for numeric variables
setMethod("svalue",signature(obj="numeric"),
          function(obj, index=NULL, drop=NULL, ... ) {
            return(obj)
          })

##' svalue method for character data
setMethod("svalue",signature(obj="character"),
          function(obj, index=NULL, drop=NULL, ... ) {
            if(length(obj) == 1)
              return(get(obj, envir=.GlobalEnv))
            else
              return(obj)
          })

##' package generic has toolkit, object to dispatch on
##' @alias svalue
setGeneric(".svalue",function(obj, toolkit, index=NULL, drop=NULL,  ...)
           standardGeneric(".svalue"))


################# svalue<- #################################
##' svalue<- generic
##'
##' Method for setting main property of a widget
##' @param obj guiWidget object
##' @param index For selection widgets one can set by index or by value
##' @param ... generally ignored, but some widgets use this to pass in extra values. See toolkit documentation
##' @param value new value for main property of widget
##' @return called for its side effect. Whether the "change" handler is called is toolkit and widget dependent.
##' @export
setGeneric("svalue<-",function(obj, index=NULL, ...,value) standardGeneric("svalue<-"))

##' svalue method for any widget
setReplaceMethod("svalue",signature(obj="guiWidget"),
          function(obj, index=NULL,  ...,value) {
            toolkit = obj@toolkit
            .svalue(obj@widget, toolkit, index=index,  ...) <- value
            obj
          })

##' dispatch to toolkit
##' @alias svalue<-
setGeneric(".svalue<-",function(obj, toolkit, index=NULL, ..., value)
           standardGeneric(".svalue<-"))

############### [ ###################################
##' [ method for selection widgetgs in guiWidget. 
##'
##' Method to refer to objects to select from
##' @param x guiWidget selection widget object
##' @param i index of vector, or row
##' @param j index of column, if possible
##' @param drop logical indicating if values should be dropped, when sensible
##' @return returns vector or data frame of values
##' @export
setMethod("[",
          signature(x="guiWidget"),
          function(x,i,j,...,drop=TRUE) {
            if(missing(drop)) drop <- TRUE
            return(.leftBracket(x@widget, x@toolkit,i,j,...,drop=drop))
          })

##' generic for implementing [ at toolkit level
setGeneric(".leftBracket",function(x, toolkit, i,j, ..., drop=TRUE)
           standardGeneric(".leftBracket"))

################ [<- ##################################
##' base [<- method for guiWidget.
##'
##' For selction widgets this method allows the possible selected values to be set
##' @param x guiWidget selection widget instance
##' @param i index of vector or row, if appropriate
##' @param j index of coumn, if appropriate
##' @param ... mostly ignored, but some toolkits may implement hidden arguments
##' @param value new value
setReplaceMethod("[",signature(x="guiWidget"),
          function(x,i,j,...,value) {
            toolkit = x@toolkit
            if(missing(i) && missing(j))
              .leftBracket(x@widget, toolkit,...) <- value
            else if(missing(j))
              .leftBracket(x@widget, toolkit,i,...) <- value
            else 
              .leftBracket(x@widget, toolkit,i,j,...) <- value
            return(x)
          })

##' generic for .leftBracket<- to implement [<- at toolkit level
setGeneric(".leftBracket<-",function(x, toolkit, i,j, ..., value)
           standardGeneric(".leftBracket<-"))


################## visible ################################
##' Generic to check visibility of widget
##'
##' @param obj guiWidget instance
##' @param set logical. Should one set the value. One should use
##' \code{visible<-} in most all cases, but for \code{gbasicdialog}
##' this is necessary
##' @return logical indicating if object is visible (shown)
##' @export
setGeneric("visible",function(obj, set=NULL, ...) standardGeneric("visible"))

##' visible method for widgets
setMethod("visible",signature(obj="guiWidget"),
          function(obj, set=NULL, ...) {
            toolkit = obj@toolkit
            .visible(obj@widget,toolkit, set=set, ...)
          })
##' dispatch with toolkit
##' @alias visible
setGeneric(".visible",function(obj, toolkit, set=NULL, ...) standardGeneric(".visible"))

############### visible<- ###################################
##' Generic method to adjust visibility of object
##'
##' A widget is visible if it is drawn (shown). This method is used to toggle that state.
##' @param obj guiWidget instance
##' @param ... ignored
##' @param value logical. Should object be visible or hidden
##' @note Many widgets override this method. If that is the case and
##' you wish to hide the widget, then place in a box container.
##' @return Called for side effect
##' @export
setGeneric("visible<-",function(obj, ..., value) standardGeneric("visible<-"))

##' visible<- method for widgets
setReplaceMethod("visible",signature(obj="guiWidget"),
          function(obj, ..., value) {
            toolkit = obj@toolkit
            .visible(obj@widget, toolkit, ...) <- value
            return(obj)
          })
##' dispatch with toolkit
##' @alias visible<-
setGeneric(".visible<-",function(obj, toolkit,...,value)
           standardGeneric(".visible<-"))

############### enabled ###################################
##' Generic to check if widget is sensitive to user input
##'
##' @param obj guiWiget object
##' @param ... ignored
##' @return logical indicating if widget is sensitive to user input
##' @export
setGeneric("enabled",function(obj, ...) standardGeneric("enabled"))

##' Method to check if widget is sensitive to user input
setMethod("enabled",signature(obj="guiWidget"),
          function(obj, ...) {
            toolkit = obj@toolkit
            .enabled(obj@widget, toolkit,...)
          })
##' dispatch with toolkit
##' @alias enabled
setGeneric(".enabled",function(obj, toolkit,...) standardGeneric(".enabled"))

################ enabled<- ##################################
## Generic to set whether widget is sensitive to user input
setGeneric("enabled<-",function(obj, ..., value) standardGeneric("enabled<-"))

## method to adjust whether widget is sensitive to user input
setReplaceMethod("enabled",signature(obj="guiWidget"),
          function(obj, ..., value) {
            toolkit = obj@toolkit
            .enabled(obj@widget, toolkit,...) <- value
            return(obj)
          })
##' dispatch with toolkit
##' @alias enabled<-
setGeneric(".enabled<-",function(obj, toolkit,...,value)
           standardGeneric(".enabled<-"))

############### editable ###################################
##' Generic to check if widget can be edited
setGeneric("editable",function(obj, ...) standardGeneric("editable"))

##' Method to check if widget can be edited
setMethod("editable",signature(obj="guiWidget"),
          function(obj, ...) {
            toolkit = obj@toolkit
            .editable(obj@widget, toolkit,...)
          })
##' dispatch with toolkit
##' @alias editable
setGeneric(".editable",function(obj, toolkit,...) standardGeneric(".editable"))

################ editable<- ##################################
## Generic to set whether widget can be edited
setGeneric("editable<-",function(obj, ..., value) standardGeneric("editable<-"))

## method to adjust whether widget can be edited
setReplaceMethod("editable",signature(obj="guiWidget"),
          function(obj, ..., value) {
            toolkit = obj@toolkit
            .editable(obj@widget, toolkit,...) <- value
            return(obj)
          })
##' dispatch with toolkit
##' @alias editable<-
setGeneric(".editable<-",function(obj, toolkit,...,value)
           standardGeneric(".editable<-"))

############## size ####################################
##' Generic for size method
setGeneric("size",function(obj, ...) standardGeneric("size"))

##' Method to return preferred size of widget
setMethod("size",signature(obj="guiWidget"),
          function(obj, ...) {
            toolkit = obj@toolkit
            .size(obj@widget, toolkit,...)
          })
##' dispatch with toolkit
##' @alias size
setGeneric(".size",function(obj, toolkit,...) standardGeneric(".size"))

############### sizes<- ###################################
## Generic for method to adjust size of widget
setGeneric("size<-",function(obj, ..., value) standardGeneric("size<-"))

##' Method to adjust size of widget
##'
##' 
setMethod("size<-",signature(obj="guiWidget"),
          function(obj, ..., value) {
            toolkit = obj@toolkit
            .size(obj@widget, toolkit,...) <- value
            return(obj)
          })

##' dispatch with toolkit
##' @alias size<-
setGeneric(".size<-",function(obj, toolkit,...,value)
           standardGeneric(".size<-"))



################# font #################################
##' Generic to get the font information
setGeneric("font",function(obj, ...) standardGeneric("font"))
## base method for font generic
setMethod("font",signature(obj="guiWidget"),
          function(obj, ...) {
            toolkit <- obj@toolkit
            .font(obj@widget, toolkit,...)
          })

##' dispatch with toolkit
##' @alias font
setGeneric(".font",function(obj,toolkit,...) standardGeneric(".font"))

################# font<- #################################
##' Generic to set font properties of widget
setGeneric("font<-",function(obj, ..., value) standardGeneric("font<-"))

##' base method for setting font properties of a widget
setMethod("font<-",signature(obj="guiWidget"),
          function(obj, ..., value) {
            toolkit <- obj@toolkit
            .font(obj@widget, toolkit,...) <- value ## DEPRECATED.fixFontMessUp(value)
            return(obj)
          })

##' dispatch with toolkit
##' @alias font<-
setGeneric(".font<-",function(obj, toolkit,...,value)
           standardGeneric(".font<-"))

############## tag ####################################
##' Generic for retrieving data from an object
setGeneric("tag",function(obj,i, drop=TRUE, ...) standardGeneric("tag"))

## base method for retrieving data from a widget
setMethod("tag",signature(obj="guiWidget"),
          function(obj,i,drop=TRUE, ...) {
            toolkit = obj@toolkit
            .tag(obj@widget, toolkit,i, drop=drop,...)
          })
##' dispatch with toolkit
##' @alias tag
setGeneric(".tag",function(obj, toolkit,i, drop=TRUE,...) standardGeneric(".tag"))


################ tag<- ##################################
##' Generic to set arbitrary data in an object
setGeneric("tag<-",function(obj, i, replace=TRUE, ..., value) standardGeneric("tag<-"))

##' base method for assigning data to an object
setMethod("tag<-",signature(obj="guiWidget"),
          function(obj, i, replace=TRUE, ..., value) {
            toolkit <- obj@toolkit
            .tag(obj@widget, toolkit,i, replace, ...) <- value
            return(obj)
          })

##' dispatch with toolkit
##' @alias tag<-
setGeneric(".tag<-",function(obj, toolkit,i, replace=TRUE,...,value)
           standardGeneric(".tag<-"))

################ id ##################################
##' Generic to retrieve id of object
setGeneric("id",function(obj, ...) standardGeneric("id"))
##' Base method to retrieve id from an object
setMethod("id",signature(obj="guiWidget"),
          function(obj, ...) {
            toolkit <- obj@toolkit
            .id(obj@widget, toolkit,...)
          })
##' dispatch with toolkit
##' @alias id
setGeneric(".id",function(obj, toolkit,...) standardGeneric(".id"))


################ id<- ##################################
##' Generic to set id for an object
setGeneric("id<-",function(obj, ..., value) standardGeneric("id<-"))

##' base method to set id for an object
setMethod("id<-",signature(obj="guiWidget"),
          function(obj, ..., value) {
            toolkit = obj@toolkit
            .id(obj@widget,toolkit,...) <- value
            return(obj)
          })
##' dispatch with toolkit
##' @alias id<-
setGeneric(".id<-",function(obj, toolkit,...,value)
           standardGeneric(".id<-"))

############## isExtant ####################################
##' generic to check if a widget still exists
setGeneric("isExtant",function(obj, ...) standardGeneric("isExtant"))

##' base method to check if a widget still exists
setMethod("isExtant",signature(obj="guiWidget"),
          function(obj, ...) {
            toolkit = obj@toolkit
            .isExtant(obj@widget,toolkit, ...)
          })

##' dispatch with toolkit
##' @alias isExtant
setGeneric(".isExtant",function(obj, toolkit, ...) standardGeneric(".isExtant"))


################## toolkitProvidesWidget ################################
##' function to check if a toolkit provides a widget
##'
##' @export
##' @param widgetname Character. Name of widget
##' @param toolkit the toolkit to check
##' @returns a logical indicating if the widget is implemented in the toolkit
toolkitProvidesWidget <- function(
                                  widgetname,
                                  toolkit=guiToolkit()){
  .toolkitProvidesWidget (toolkit, widgetname)
}

##' generic for toolkit dispatch
##' @alias toolkitProvidesWidget
setGeneric( '.toolkitProvidesWidget' ,
           function(toolkit,
                    widgetname)
           standardGeneric( '.toolkitProvidesWidget' ))


##' implementation for ANY toolkit
##' @alias toolkitProvidesWidget
setMethod(".toolkitProvidesWidget",
          signature(toolkit="ANY"),
          function(toolkit,
                   widgetname) {
            notThere <- list(guiWidgetsToolkitQt=c("ggraphics","ggraphicsnotebook"),
                             guiWidgestToolkitRGtk2=c("gsvg", "ghtml"),
                             guiWidgetsrToolkitJava=c("gsvg", "ghtml", "ggraphics", "ggraphicsnotebook"),
                             guiWidgetsToolkittcltk=c("gsvg", "ghtml", "ggraphics", "ggraphicsnotebook",
                               "gdfnotebook")
                             )

            notThere <- notThere[[class(toolkit)]]
            return(!widgetname %in% notThere)
          })



##################################################
## Usual R methods made into methods for dispatch

################### update ###############################
##' generic for update
setGeneric("update")

##' base method for update to dispatch on guiWidget instances
setMethod("update",signature(object="guiWidget"),
          function(object, ...) {
            .update(object@widget, object@toolkit, ...)
          })

##' generic for toolkit dispatch
##' @alias update
setGeneric(".update",function(object, toolkit,  ...)
           standardGeneric(".update"))

##' Base method to find length of guiWidget instances
setMethod("length",signature(x="guiWidget"),
          function(x) {
            .length(x@widget, x@toolkit)
          })

##' generic for toolkit dispatch
##' @alias length
setGeneric(".length",function(x, toolkit)
           standardGeneric(".length"))

############## dim ####################################
##' base method to find dim(ension) of guiWidget instances
setMethod("dim",signature(x="guiWidget"),
          function(x) {
            .dim(x@widget, x@toolkit)
          })

##' generic for toolkit dispatch
##' @alias dim
setGeneric(".dim",function(x, toolkit)
           standardGeneric(".dim"))

##' base method to find dimnames attribute of guiWidget instance
#setGeneric("dimnames")
setMethod("dimnames",signature(x="guiWidget"),
          function(x) {
            .dimnames(x@widget, x@toolkit)
          })

##' generic for toolkit dispatch
##' @alias dimnames
setGeneric(".dimnames",function(x, toolkit)
           standardGeneric(".dimnames"))

##' base method to set dimnames for guiWidget instance
#setGeneric("dimnames<-")
setReplaceMethod("dimnames",signature(x="guiWidget"),
                 function(x,value) {
                   .dimnames(x@widget, x@toolkit) <- value
                   return(x)
                 })
##' generic for toolkit dispatch
##' @alias dimnames<-
setGeneric(".dimnames<-",function(x, toolkit, value) {
  standardGeneric(".dimnames<-")
})


## names
## as of 2.5.0 this became primiive
if(as.numeric(R.Version()$major) <= 2 &
   as.numeric(R.Version()$minor) <= 4.1) {
  setGeneric("names")
  setGeneric("names<-")
}

##' Base method for getting names of guiWidget instances
setMethod("names",signature(x="guiWidget"),
          function(x) {
            .names(x@widget, x@toolkit)
          })

##' generic for toolkit dispatch
##' @alias names
setGeneric(".names",function(x, toolkit)
           standardGeneric(".names"))

##' base method to set names of guiWidget instance
setReplaceMethod("names",signature(x="guiWidget"),
                 function(x,value) {
                   .names(x@widget, x@toolkit) <- value
                   return(x)
                 })
##' generic for toolkit dispatch
##' @alias names<-
setGeneric(".names<-",function(x, toolkit, value) {
  standardGeneric(".names<-")
})

##################################################
## Work with underlying toolkit objects

##################################################
##' Generic for method to return toolkit widget from guiWidget instance
setGeneric("getToolkitWidget",function(obj) standardGeneric("getToolkitWidget"))

##' base method to return toolkit widget from guiWidget instance
setMethod("getToolkitWidget",signature(obj="guiWidget"),
          function(obj) {
            .getToolkitWidget(obj@widget, obj@toolkit)
          })

##' generic for toolkit dispatch
##' @alias getToolkitWidget
setGeneric(".getToolkitWidget",function(obj, toolkit)
           standardGeneric(".getToolkitWidget"))


################## access via $ and [[ #############################
##' Invoke a method on the underlying toolkit object
##'
##' @param x guiwidget
##' @param meth_name name of method
##' @return function that one can call, as in \code{x$method(a=1)}
##' @export
"$.guiWidget" <- function(x, meth_name)
            .callToolkitMethod(x@widget, x@toolkit, meth_name)

setGeneric(".callToolkitMethod",function(x, toolkit, meth_name)
           standardGeneric(".callToolkitMethod"))

##' Get an underlying property of the toolkit object
##' @param x guiWidget
##' @param ... first argument has property
##' @return value
##' @export
##'
"[[.guiWidget" <- function(x, ...) {
  property <- list(...)[[1]]
  .getToolkitProperty(x@widget, toolkit=x@toolkit, property=property)
}

setGeneric(".getToolkitProperty",function(x, toolkit, property)
           standardGeneric(".getToolkitProperty"))

##' Set underlying toolkit property
##' @param x guiWidget
##' @param i property to set
##' @param j ignored
##' @param value new value of property
##' @export
"[[<-.guiWidget" <- function(x, i, j, value) {
  .setToolkitProperty(x@widget, toolkit=x@toolkit, property=i, value=value)
  x
}

setGeneric(".setToolkitProperty",function(x, toolkit, property, value)
           standardGeneric(".setToolkitProperty"))

