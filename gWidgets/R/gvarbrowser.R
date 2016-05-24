##' @include guiComponents.R

##' Class for a workspace or variable browser
setClass("gVarBrowser",
         contains="guiComponent",
         prototype=prototype(new("guiComponent"))
         )

##' Constructor for workspace variable browser
##'
##' @export
gvarbrowser <- function(
                        handler = NULL,
                        action = "summary",
                        container = NULL ,...,
                        toolkit=guiToolkit()){
  widget <- .gvarbrowser (toolkit,
                          handler=handler, action=action, container = container, ...
                          )
  obj <- new( 'gVarBrowser',widget=widget,toolkit=toolkit) 
  return(obj)
}


##' generic for toolkit dispatch
##' @alias gvarbrowser
setGeneric( '.gvarbrowser' ,
           function(toolkit,
                    handler = NULL,
                    action = "summary", container = NULL,... )
           standardGeneric( '.gvarbrowser' ))


