##' @include methods.R
NULL

##' constructor providing a widget for displaying a line in a GUI
##'
##' The \code{gseparator} widget provides a horizontal or vertical
##' line to visually divide child components of its parent
##' container. In addition to box containers this can be used within
##' toolbars (where one uses \code{parent} and not \code{container}). 
##' @param horizontal Logical. Is separator drawn horizontally?
##' @inheritParams gcontainer
##' @export
##' @examples
##' \dontrun{
##' 
##' w <- gwindow("Within page", visible=FALSE)
##' g <- gvbox(container=w)
##' glabel("Lorem ipsum ...", cont=g)
##' gseparator(cont=g)
##' bg <- ggroup(cont=g); addSpring(bg)
##' gbutton("dismiss", container=bg, handler=function(h,...) dispose(w))
##' visible(w) <- TRUE
##' 
##' w1 <- gwindow("within layout", visible=FALSE)
##' lyt <- glayout(container=w1)
##' lyt[1,1] <- "label"
##' lyt[2,1:2] <- gseparator(container=lyt)
##' lyt[3,2] <- "asdf"
##' visible(w1) <- TRUE
##' 
##' w2 <- gwindow("Within toolbar", visible=FALSE)
##' l <- list(file=gaction("File", parent=w2),
##'           sep=gseparator(parent=w2),
##'           quit=gaction("quit", parent=w2))
##' gtoolbar(l, cont=w2)
##' glabel("Lorem ipsum ...", container=w2)
##' visible(w2) <- TRUE
##' }
gseparator <- function(
                       horizontal = TRUE, container = NULL, ... ,
                       toolkit=guiToolkit()){
  obj <- .gseparator (toolkit,
               horizontal=horizontal, container=container ,...
               )
  check_return_class(obj, "GSeparator")
  obj
  
}


##' generic for toolkit dispatch
##'
##' @export
##' @rdname gseparator
.gseparator <- function(toolkit,
                    horizontal = TRUE, container = NULL, ... )
           UseMethod( '.gseparator' )
