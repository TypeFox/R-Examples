##' @include methods.R
NULL


##' A widget for displaying an image file
##'
##' @param filename basename of file
##' @param dirname dirname of file
##' @param stock.id stock id of icon (if non NULL)
##' @param size size of icon when a stock id (toolkit dependent)
##' @param handler handler if image is clicked on. 
##' @param action  passed to handler
##' @param container parent container
##' @param ... passed to add method of parent
##' @param toolkit toolkit
##' @export
gimage <- function(
                   filename = "", dirname = "", stock.id=NULL, size = "", handler = NULL,
                   action = NULL, container = NULL, ... ,
                   toolkit=guiToolkit()){
  obj <- .gimage (toolkit,
                  filename=filename, dirname=dirname, stock.id=stock.id, size=size,
                  handler=handler, action=action, container=container ,...
                  )
  check_return_class(obj, "GImage")
  obj
}

##' generic for toolkit dispatch
##'
##' @export
##' @rdname gimage
.gimage <- function(toolkit,
                    filename = "", dirname = "", stock.id=NULL, size = "",
                    handler = NULL,action = NULL, container = NULL, ... )
           UseMethod( '.gimage' )


