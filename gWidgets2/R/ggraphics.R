##' @include methods.R
NULL

##' Constructor for an embeddable graphics device
##'
##' Some toolkits provide an embeddable graphics device. When this is the case, this widget provides same.
##' @param width width of device (pixels)
##' @param height hieght of widget (pixels)
##' @param dpi dots per inch
##' @param ps pointsize
##' @inheritParams gwidget
##' @export
##' @examples
##' \dontrun{
##' ## This shows how to use the device within a notebook
##' 
##' w <- gwindow("notebook example")
##' nb <- gnotebook(cont=w)
##' 
##' devs <- lapply(1:5, function(i) ggraphics(cont=nb, label=as.character(i)))
##' 
##' addHandlerChanged(nb, handler=function(h,...) {
##'   ## Tricky part is svalue(h$obj) is not the new page number -- but old
##'   ## so we use the pageno component here
##'     gg <- h$obj[h$pageno]
##'     visible(gg) <- TRUE
##' })
##' 
##'
##' }
ggraphics <- function(
                      width = dpi * 6, height = dpi * 6, dpi = 75, ps = 12,
                      handler = NULL,action = NULL, container = NULL, ... ,
                      toolkit=guiToolkit()){
  obj <- .ggraphics (toolkit,
                     width=width, height=height, dpi=dpi, ps=ps,
                     handler = NULL,action = NULL, container=container ,...
                     )
  check_return_class(obj, "GGraphics")
  return(obj)
}


##' generic for toolkit dispatch
##'
##' @export
##' @rdname ggraphics
.ggraphics <-  function(toolkit,
                        width = dpi * 6, height = dpi * 6, dpi = 75, ps = 12,
                        handler = NULL,action = NULL, container = NULL, ... )
  UseMethod( '.ggraphics' )


##' change handler for ggraphics
##'
##' The change handler for ggraphics is called when a rubber-band selection is completed
##' @inheritParams addHandler
##' @export
##' @rdname gWidgets-handlers
##' @method addHandlerChanged default
##' @S3method addHandlerChanged default
addHandlerChanged.GGraphics <- function(obj, handler, action=NULL, ...)
  obj$add_handler_changed(handler, action=action, ...)

##' click handler for ggraphics
##'
##' The click handler is called on a mouse click. The handler object should pass in value for \code{x}, \code{y}
##' @inheritParams addHandler
##' @export
##' @rdname gWidgets-handlers
##' @method addHandlerClicked default
##' @S3method addHandlerClicked default
addHandlerClicked.default <-  function(obj, handler, action=NULL, ...)
    obj$add_handler_clicked(handler, action=action, ...)




## ##' constructor for notebook to hold multiple graphics devices
## ##'
## ##' @export
## ggraphicsnotebook <- function(
##                               width = dpi * 6, height = dpi * 6, dpi = 75, container = NULL,      ... ,
##                               toolkit=guiToolkit()){
##   widget <- .ggraphicsnotebook (toolkit,
##                                 width=width, height=height, dpi=dpi, container=container ,...
##                                 )
##   obj <- new( 'gGraphicsNotebook',widget=widget,toolkit=toolkit) 
##   return(obj)
## }


## ##' generic for toolkit dispatch
## ##' @alias ggraphicsnotebook
## setGeneric( '.ggraphicsnotebook' ,
##            function(toolkit,
##                     width = dpi * 6, height = dpi * 6, dpi = 75,
##                     container = NULL,      ... )
##            standardGeneric( '.ggraphicsnotebook' ))
