##' @include methods.R
NULL


##' constructor for a two-paned container
##'
##' A container for holding two child widgets where the space
##' allocated to each can be manipulated by the user with a pane
##' between the widgets, or programatticaly via \code{svalue<-}.  The
##' value specified to \code{svalue<-} can be a number in $[0,1]$, in
##' which case it is a proportion or an integer, in which case it is a
##' pixel size (from the left or the top). The ambiguous case \code{1}
##' or \code{1L} is determined by class. The value of \code{svalue} is in proportion units.
##' 
##' Child widgets are added in the usual way, typically through the
##' \code{container} argument of a constructor. Only two children may
##' be added. Children expand and fill the allocated space.
##' @param horizontal direction of layout
##' @inheritParams gwidget
##' @export
##' @note Setting the size is often only possible after the container
##' has been realized on the screen. In the example, this call of
##' \code{svalue<-} is done after the parent window is made visible
##' for this reason. There were arguments to specify the children at
##' construction, but these have been removed.
##' @examples
##' \dontrun{
##' w <- gwindow("gpanedgroup", visible=FALSE)
##' pg <- gpanedgroup(cont=w)
##' gbutton("left", cont=pg)
##' gbutton("right", cont=pg)
##' 
##' visible(w) <- TRUE
##' svalue(pg) <- 0.33
##' }
gpanedgroup <- function(horizontal = TRUE,  container = NULL , ...,
                        toolkit=guiToolkit()){


  deprecated_args <-
    list(widget1=c("widget1 (and 2) argument has been deprecated.","Child components are added as with other containers."),
         widget2=c("widget1 (and 2) argument has been deprecated.","Child components are added as with other containers.")
         )
  check_deprecated(deprecated_args, ...)
  
  obj <- .gpanedgroup (toolkit,
                       horizontal=horizontal, 
                       container=container, ...
                )

  check_return_class(obj, "GPanedGroup")
  obj   
  
}


##' generic for toolkit dispatch
##'
##' @export
##' @rdname gpanedgroup
.gpanedgroup <-  function(toolkit,
                          horizontal = TRUE, container = NULL, ... )
  UseMethod( '.gpanedgroup' )


