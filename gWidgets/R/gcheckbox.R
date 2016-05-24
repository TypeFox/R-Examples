##' @include guiComponents.R

##' Checkbox class
setClass("gCheckbox",
         contains="guiComponent",
         prototype=prototype(new("guiComponent"))
         )

##' constructor for checkbox widget
##'
##' A checkbox widget has a checkbox or toggle button to indicate selection or not
##' @param text label text
##' @param checked is button selected
##' @param use.togglebutton Use a toggle button (shows depressed) not a check box
##' @param handler Callback called when toggle is changed.
##' @param action passed to handler
##' @param container parent container
##' @param ... passed to \code{add} method of container
##' @param toolkit toolkit
##' @example ~/pmg/r-forge/gwidgets/pkg/gWidgets/inst/tests/ex-gcheckbox.R
##' @export
##' @return Returns an object of class \code{gCheckbox} for which the
##' following methods are overridden:
##' %
##' \enumerate{
##' \item \code{svalue} gets state by boolean value
##' 
##' \item \code{svalue<-} sets state by boolean value
##' 
##' \item \code{[} Returns label
##' 
##' \item \code{[<-} sets label
##' }
gcheckbox =function(
  text, checked = FALSE, use.togglebutton=FALSE, handler = NULL, action = NULL, container = NULL, ... ,
  toolkit=guiToolkit()){

  ## ensure text is given
  if(missing(text))
    text <- ""
  text <- as.character(text)[1]
  ## checked is logical
  checked <- as.logical(checked)[1]
  
  widget =  .gcheckbox (toolkit,
    text=text, checked=checked,
    use.togglebutton=use.togglebutton,
    handler=handler, action=action, container=container, ...
    )
  obj = new( 'gCheckbox',widget=widget,toolkit=toolkit) 
  return(obj)
}


##' Generic for toolkit dispatch
##' @alias gcheckbox
setGeneric( '.gcheckbox' , function(toolkit,
                                    text, checked = FALSE, use.togglebutton=FALSE, handler = NULL, action = NULL,
                                    container = NULL, ... ) standardGeneric( '.gcheckbox' ))



##' svalue method
##'
##' Ensure value is logical
setReplaceMethod("svalue",signature(obj="gCheckbox"),
          function(obj, index=NULL,  ...,value) {
            value <- as.logical(value)[1]
            callNextMethod(obj, index, ..., value=value)
          })

            
##' ensure value is character of length 1
setReplaceMethod("[",signature(x="gCheckbox"),
          function(x,i,j,...,value) {
            value <- as.character(value)[1]
            callNextMethod(x, i, j, ..., value=value)
          })
