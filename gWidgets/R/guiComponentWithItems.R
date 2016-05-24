##' @include guiComponent.R
##'
##' 
##' A widget subclass for components with items
setClass("guiComponentWithItems",
         contains="guiComponent",
         prototype=prototype(new("guiComponent"))
         )

## [, [<-

