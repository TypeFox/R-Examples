##' @includes guiComponent.R

##' Button class
setClass("gButton",
         contains="guiComponent",
         prototype=prototype(new("guiComponent"))
         )

##' Button widget constructor
##' 
##' @param text Button text. If text matches a stock icon, that icon is also displayed.
##' @param border Logical. Is button drawn with a border
##' @param handler handler to call when button is invoked
##' @param action passed to handler when called. If a \code{gaction}
##' instance, then the action will define the text and the handler
##' @param container parent container
##' @param ... passed to \code{add} method of container
##' @param toolkit toolkit
##' @export
##' @example ~/pmg/r-forge/gwidgets/pkg/gWidgets/inst/tests/ex-gbutton.R
##' @return An object of class \code{gButton} for which the following methods are overridden:
##'
##' \enumerate{
##' \item \code{svalue} gets button text
##' 
##' \item \code{svalue<-} sets button text
##' }
##' The "change" handler is called when the button is clicked.
gbutton =function(
  text = "", border=TRUE, handler = NULL, action = NULL, container = NULL,      ...,
  toolkit=guiToolkit()){
  force(toolkit)                        # load package


  text <- paste(as.character(text), collapse="\n")

  
  widget =  .gbutton (toolkit,
    text, border, handler, action, container,...)

  obj <- new( 'gButton',widget=widget,toolkit=toolkit) 

  return(obj)
}

##' gbutton generic for toolkit dispatch
##' @alias gbutton
##' @export
setGeneric( '.gbutton' , function(toolkit,
                                  text = "", border=TRUE, handler = NULL, action = NULL, container = NULL,...)
           standardGeneric( '.gbutton' ))



################ methods ##################################

##' svalue method for button
##'
##' Main property of a button is the label
##' @param obj object
##' @param index ignored
##' @param drop ignored
##' @return character. String on button
##' @exports
setMethod("svalue", signature(obj="gButton"),
          function(obj, index=NULL, drop=NULL, ... ) {
            callNextMethod()
          })




##' set label in a button
##'
##' @param obj
##' @param index ignored
##' @param ... ignored
##' @param value character. If string matches stock icon name, then an icon will be added to button.
##' @return void
##' @exports
setReplaceMethod("svalue", signature(obj="gButton"),
          function(obj, index=NULL, ...,value) {
            value <- paste(value, collapse="\n")
            callNextMethod()
          })




################## defaultWidget ################################
##' Generic to check if a button is the default button
setGeneric("defaultWidget",function(obj, ...) standardGeneric("defaultWidget"))

##' base method to check if a button is the default
setMethod("defaultWidget",signature(obj="gButton"),
          function(obj, ...) {
            toolkit <- obj@toolkit
            .defaultWidget(obj@widget, toolkit,...)
          })
##' dispatch with toolkit
##' @alias defaultWidget
setGeneric(".defaultWidget",function(obj, toolkit,...) standardGeneric(".defaultWidget"))

################### defaultWidget<- ###############################
##' generic to set a button as the default one
setGeneric("defaultWidget<-",function(obj, ..., value) standardGeneric("defaultWidget<-"))
##' base method to set a button as default
setMethod("defaultWidget<-",signature(obj="gButton"),
          function(obj, ..., value) {
            toolkit = obj@toolkit
            .defaultWidget(obj@widget, toolkit,...) <- value
            return(obj)
          })
##' dispatch with toolkit
##' @alias defaultWidget<-
setGeneric(".defaultWidget<-",function(obj, toolkit,...,value)
           standardGeneric(".defaultWidget<-"))





