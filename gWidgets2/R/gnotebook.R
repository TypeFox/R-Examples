##' @include methods.R
NULL

##' Constructor for a tabbed notebook container
##'
##' The tabbed notebook container allows one to hold many different
##' pages with a mechanism -- tabs -- to switch between them. In
##' \code{gWidgets2} new pages are added through the \code{add}
##' method. This is usually called implicitly in the child object's
##' constructor. One passes in the tab label through the extra
##' \code{label} argument. Labels may be subsequently changed through
##' \code{names<-}.
##' 
##' @param tab.pos integer. Position of tabs, 1 on bottom, 2 left, 3
##' top, 4 right. (If supported)
##' @param container parent container
##' @param ... passed to \code{add} method for container
##' @param toolkit underlying toolkit
##' @seealso \code{\link{gstackwidget}} for a similar widget without
##' tabs.
##' @note In \pkg{gWidgets2} the button arguments of the
##' \code{gWidgets} constructor are removed. One passes the close
##' button request to the \code{add} method.
##' @export
##' @examples
##' \dontrun{
##' 
##' w <- gwindow("notebook example", visible=FALSE)
##' nb <- gnotebook(container=w)
##' gbutton("Page one", label="tab 1", container=nb) ## note label argument
##' gbutton("Page two", label="tab 2", container=nb)
##' svalue(nb) <- 1
##' addHandlerChanged(nb, handler=function(h,...) {
##'   message(sprintf("On page %s", h$page.no)) ## svalue(h$obj) not always right
##' })
##' svalue(nb) <- 2 ## or use "Page two"
##' dispose(nb)
##' length(nb)
##' 
##' }
gnotebook <- function(
                      tab.pos = 3, 
                      container = NULL, ... ,
                      toolkit=guiToolkit()){

  obj <- .gnotebook (toolkit,
              tab.pos=tab.pos, 
              container=container ,...
              )

  check_return_class(obj, "GNotebook")
  obj   
  
}


##' generic for toolkit dispatch
##'
##' @export
##' @rdname gnotebook
.gnotebook <-  function(toolkit,
                        tab.pos = 3, 
                        container = NULL, ... )
  UseMethod( '.gnotebook' )


##' add method for notebooks
##'
##' Children added to notebooks need a label, a position and
##' optionally a close button (if supported). The arguments expand,
##' fill, anchor are not specified -- children expand and fill the allocated space.
##' @param obj gnotebook object
##' @param child some child component to add
##' @inheritParams add
##' @note To keep the signature the same as the generic, several arguments are passed in via ...:
##' 
##' \describe{
##' 
##' \item{label}{ A character. Label text for tab}
##' 
##' \item{i}{An integer in \code{0} to \code{length(obj)} indicating
##' the position to insert child. The new page is inserted to the
##' right of page  number \code{i}. When \code{i=0}, the page appears
##' at the front, when \code{i} is not specified it appears at the
##' end.
##' }
##' 
##' \item{close.button}{A logical. If \code{TRUE} -- and the toolkit
##' supports it -- the page tab will include a close button.
##' }
##' }
##' @return none. called for its side effect.
##' @export
##' @rdname gnotebook
##' @method add GNotebook
##' @S3method add GNotebook
add.GNotebook <- function(obj, child, expand, fill, anchor, ...) {
  ## process passed in args
  args <- list(...)
  label <- getWithDefault(args$label, "")
  i <- getWithDefault(args$i, length(obj))
  close.button <- getWithDefault(args$close.button, FALSE)
  
  obj$add_child(child, label, i, close.button, ...)
}


##' Remove current page from notebook
##'
##' Dispose deletes the current page, not the entire notebook
##' object. To delete a specific page, a combination of
##' \code{svalue<-} and \code{dispose} may be used.
##' @export
##' @rdname gnotebook
##' @method dispose GNotebook
##' @S3method dispose GNotebook
dispose.GNotebook <- function(obj, ...) {
  obj$remove_current_page()
}


##' get tab names of notebook
##'
##' The \code{names} of a notebook are the page tab labels. These may
##' be retrieved and set through the \code{names} method.
##' @export
##' @rdname gnotebook
##' @method names GNotebook
##' @S3method names GNotebook
"names.GNotebook" <- function(x) x$get_names()


##' svalue method
##'
##' Set the currently raised tab by index (the default) or name
##' @param index  \code{TRUE} refer to tab by 1-based
##' index; \code{FALSE} allows reference by tab label.
##' @param value assignment value
##' @export
##' @usage \method{svalue}{GNotebook} (obj, index=TRUE, ...) <- value
##' @rdname gnotebook
##' @method svalue<- GNotebook
##' @S3method svalue<- GNotebook
"svalue<-.GNotebook" <- function(obj, index=TRUE,  ...,value) {
    if (!index) {
        index = TRUE
        value = match(value, names(obj))
    }
    NextMethod()
}

##' "[" method
##'
##' The notebook object contains pages referenced by index. This allows access to underlying page.
##' @param x \code{GNotebook} object
##' @param i row index. Either integer or character
##' @param j ignored
##' @param drop ignored
##' @export
##' @rdname gnotebook
##' @method [ GNotebook
##' @S3method [ GNotebook
"[.GNotebook" <- function(x, i, j, ..., drop=TRUE) {
    if (is.character(i))
        i <- match(i, names(x))
    NextMethod()
}

##' add change handler
##'
##' The change handler for the notebook is called when the page
##' changes. Tthe new page number is passed back in the \code{page.no}
##' component of 'h', which in some cases may differ from the value
##' given by \code{svalue} within the handler call.
##' @export
##' @rdname gnotebook
##' @param handler handler
##' @param action passed along to handler via \code{h[["action"]]}.
##' @method addHandlerChanged GNotebook
##' @S3method addHandlerChanged GNotebook
addHandlerChanged.GNotebook <- function(obj, handler, action=NULL, ...) {
  NextMethod()
}
