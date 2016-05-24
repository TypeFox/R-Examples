##' @include guiContainer.R

##' Class for top level windows
setClass("gWindow",
         contains="guiContainer",
         prototype=prototype(new("guiContainer"))
         )

##' constructor for top-level window
##'
##' @export
gwindow <- function(title="Window" , visible=TRUE, name=title,
                    width = NULL, height = NULL, parent = NULL,
                    handler = NULL, action = NULL,
                    ...,
                    toolkit=guiToolkit()
                    ) {
  theArgs <- list(...)
  if(!is.null(theArgs$location)) {
    parent <- theArgs$location
    cat(gettext("location argument is renamed to 'parent'\n"))
  }
  
  ## THe visible=TRUE default is not the best. I'd change it if I could go back in time, but
  ## c'est la vie. Anyways, for those that it really bugs there is this check
  if(!is.null(getOption("gWidgets:gwindow-default-visible-is-false")))
    visible <- FALSE
  
  
  win <- .gwindow(toolkit,title, visible,width, height, parent, handler, action, ...)
  obj <- new("gWindow",widget=win,toolkit=toolkit)
  return(obj)
}

##' define a toolkit constructor, dispatch on toolkit
##' @alias gwindow
setGeneric(".gwindow",function(toolkit, title, visible, width, height, parent, handler, action,...) standardGeneric(".gwindow"))


## methods
## add
## svalue: title property
## size, size<-
## visible, visible<-
## update: recompute window size
