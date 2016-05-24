##' @include guiComponents.R

##' Label class
setClass("gLabel",
         contains="guiComponent",
         prototype=prototype(new("guiComponent"))
         )

##' constructor for label widget
##'
##' @param text character. Text for label. Coerced to character and pasted together with newlines
##' @param markup logical. For some toolkits, one can specify marked up text.
##' @param editable logical. For some toolkits, label can be edited
##' when this is \code{TRUE}. Generally found to be an unexpected
##' interface for users, so use is discouraged.
##' @param handler function. For some toolkits, this handler will be
##' called when the label is clicked on. In general, labels are
##' expected to be static objects so this use is discouraged.
##' @param action passed to \code{handler}, when called
##' @param container parent container
##' @param ... generally ignored
##' @param toolkit underlying toolkit. Usually not specified
##' @return \code{gLabel} object to manipulate and creates widget on screen
##' @export
##' @examples
##' w <- gwindow()
##' g <- ggroup(container=w, horizontal=FALSE)
##' l1 <- glabel("Some label", container=g)
##' l2 <- glabel(c("pasted with", "new lines"), container=g)
##' svalue(l1) <- "New text for some label")
##' svalue(l1)
glabel = function(
  text= "", markup = FALSE, editable = FALSE, handler = NULL, 
    action = NULL, container = NULL, 
  ..., toolkit=guiToolkit()) {

  ## collapse if more than one line
  text <- paste(text, collapse="\n")
  
  widget = .glabel(toolkit,
    text= text, markup = markup, editable = editable, handler = handler, 
    action = action, container = container, 
    ...)
  obj = new("gLabel",widget=widget,toolkit=toolkit)
  return(obj)
}

##' glabel generic for toolkit
setGeneric(".glabel",function(toolkit,
                              text= "", markup = FALSE, editable = FALSE, handler = NULL, 
                              action = NULL, container = NULL, 
                              ...) standardGeneric(".glabel"))

##' svalue<- generic
##'
##' Ensure value is character vector. Pastes values together by collapsing with a new line.
##' @use_svalue_otherwise
setReplaceMethod("svalue",signature(obj="gLabel"),
                 function(obj, index=NULL,  ...,value) {
                   ## enforce that value is character
                   value <- paste(as.character(value), collapse="\n")
                   callNextMethod()
                 })
