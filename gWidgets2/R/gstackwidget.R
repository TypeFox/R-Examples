##' @include methods.R
NULL

##' Constructor for a stack of widgets
##'
##' This widget is like a notebook -- it holds a stack of pages, but
##' does not provide the tabs to work with. Most methods are inherited
##' from gnotebook's.
##' @inheritParams gcontainer
##' @export
##' @examples
##' \dontrun{
##' w <- gwindow("stack widget", visible=FALSE)
##' add_page <- function(cont, i) {
##'   g <- gvbox(container=cont)
##'   glabel(sprintf("page %s",i), container=g)
##'   bg <- ggroup(container=g); addSpring(bg)
##'   lb <- gbutton("Previous", container=bg, handler=function(h,...) {
##'     svalue(cont) <- max(1, i - 1)
##'   })
##'   rb <- gbutton("Next", container=bg, handler=function(h,...) {
##'     svalue(cont) <- min(i + 1, length(cont))
##'   })
##' }
##' sw <- gstackwidget(cont=w)
##' sapply(1:5, add_page, cont=sw)
##' visible(w) <- TRUE
##' }
gstackwidget <- function(container = NULL, ... ,
                      toolkit=guiToolkit()){

  obj <- .gstackwidget(toolkit,
                container=container ,...
                )
  check_return_class(obj, "GStackWidget")
  obj
  
}


##' generic for toolkit dispatch
##'
##' @export
##' @rdname gstackwidget
.gstackwidget <-  function(toolkit,
                        container = NULL, ... )
           UseMethod( '.gstackwidget' )




## toolkit class should inherit from GNotebook
## but just in case


##' Remove current page from stackwidget
##'
##' Dispose deletes the current page, not the entire notebook
##' object. To delete a specific page, a combination of
##' \code{svalue<-} and \code{dispose} may be used.
##' @export
##' @rdname gnotebook
##' @method dispose GStackWidget
##' @S3method dispose GStackWidget
dispose.GStackWidget <- function(obj, ...) {
  obj$remove_current_page()
}


