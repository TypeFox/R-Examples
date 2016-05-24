##' @include guiComponents.R

##' single line text edit class
setClass("gAction",
         contains="guiComponent",
         prototype=prototype(new("guiComponent"))
         )

##' An action constructor
##'
##' @param label label for action
##' @param tooltip toolktip for actin
##' @param icon icon (stock icon name) for icon
##' @param key.accel keyboard accelerator. If given, parent must be specified.
##' @param handler handler to call when action is invoked
##' @param action values passed to parameterize action
##' @param parent parent window. Needed if keyboard accelerator used.
##' @param ... 
##' @param toolkit 
##' @export
##' @return a gaction instance
gaction <- function(
                    label, tooltip=NULL, icon = NULL, key.accel = NULL,
                    handler = NULL, action = NULL, parent=NULL, ...,
                    toolkit=guiToolkit()) {
  widget <- .gaction (toolkit,
                      label, tooltip, icon, key.accel, handler, action, parent, ...
                      )
  obj <- new( 'gAction',widget=widget,toolkit=toolkit) 
  return(obj)
}


##' generic for toolkit dispatch
##' @alias gaction
setGeneric( '.gaction' ,
           function(toolkit,
                    label, tooltip = NULL, icon = NULL, key.accel=NULL,
                    handler = NULL, action = NULL, parent=NULL, 
                    ... )
           standardGeneric( '.gaction' ))


