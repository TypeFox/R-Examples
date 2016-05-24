##' @include guiComponents.R

##' class to display heiarchical data in a tree widget
setClass("gTree",
         contains="guiComponent",
         prototype=prototype(new("guiComponent"))
         )

##' constructor for widget to display heirarchical data
##'
##' @exports
gtree <- function(
                  offspring = NULL, hasOffspring = NULL, offspring.data = NULL,
                  col.types = NULL, icon.FUN = NULL, chosencol = 1, multiple = FALSE,
                  handler = NULL, action = NULL, container = NULL, ... ,
                  toolkit=guiToolkit()){
  widget <- .gtree (toolkit,
                    offspring=offspring, hasOffspring=hasOffspring,
                    offspring.data=offspring.data, col.types=col.types, icon.FUN=icon.FUN,
                    chosencol=chosencol, multiple=multiple,
                    handler=handler, action=action, container=container ,...
                    )
  obj <- new( 'gTree',widget=widget,toolkit=toolkit) 
  return(obj)
}


##' generic for toolkit dispatch
##' @alias gtree
setGeneric( '.gtree' ,
           function(toolkit,
                    offspring = NULL, hasOffspring = NULL,
                    offspring.data = NULL,
                    col.types = NULL, icon.FUN = NULL, chosencol = 1,
                    multiple = FALSE,
                    handler = NULL, action = NULL, container = NULL, ... )
           standardGeneric( '.gtree' ))


## methods

##' svalue get selected values. index=TRUE: 1-based path; index=FALSE|NULL chararacter value
##' svalue<- set path by index
##' [ get text based path (by ids)
##' update update tree at root
