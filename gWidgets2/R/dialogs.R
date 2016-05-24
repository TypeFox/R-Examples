##' @include methods.R
NULL



##' Alert dialog to display transient messages
##' 
##' @param msg character. main message. If length is 2, second component is used for detail, providing it is available.
##' @param title Title (may not be displayed)
##' @param delay length of time (in seconds) to display
##' @param parent parent object to show near
##' @param ... ignored
##' @param toolkit toolkit
##' @seealso \code{\link{gmessage}}, \code{\link{gconfirm}}, 
##' \code{\link{gbasicdialog}}, \code{\link{galert}}
##' @export
galert = function(
  msg,
  title = "message",
  delay = 3,
  parent = NULL, 
  ..., toolkit=guiToolkit()) {
  .galert(toolkit,msg, title, delay=delay, parent=parent,
            ...)
}
##' generic for toolkit dispatch
##'
##' @export
##' @rdname galert
.galert <- function(toolkit,
                    msg, title="message", delay=3, parent=NULL, ...)
  UseMethod(".galert")


##' Constructor for modal message dialog
##' 
##' @export
##' @param msg Character. message to display. 
##' @param title Character. Title
##' @param icon What icon to show
##' @param parent Hint as to where to display
##' @param ... ignored
##' @param toolkit toolkit
##' @return NULL
##' @rdname gmessage
##' @seealso \code{\link{gmessage}}, \code{\link{gconfirm}}, 
##' \code{\link{gbasicdialog}}, \code{\link{galert}}
gmessage <- function(msg,
                     title = "message",
                     icon = c("info", "warning", "error", "question"),
                     parent=NULL,
                     ..., toolkit=guiToolkit()) {
  .gmessage(toolkit, msg, title, icon=match.arg(icon), parent, ...)
  NULL
}

##' generic for toolkit dispatch
##'
##' @return NULL
##' @export
##' @rdname gmessage
.gmessage <- function(toolkit,
                      msg, title="message", icon="",
                      parent = NULL,
                      ...)
  UseMethod(".gmessage")

############## ginput ####################################

##' Constructor for modal dialog to collect a line of text
##'
##' @param msg Character. Message to display.
##' @param text Character. Initial text
##' @param title Character. Title of window
##' @param icon which icon to display
##' @param parent gives hint as to where to place dialog
##' @param ... ignored
##' @param toolkit toolkit
##' @return value typed into box or \code{character(0)}
##' @export
##' @rdname ginput
##' @seealso \code{\link{gmessage}}, \code{\link{gconfirm}}, 
##' \code{\link{gbasicdialog}}, \code{\link{galert}}
ginput <- function(
                   msg, text="",
                   title = "Input",
                   icon = c("info", "warning", "error", "question"),
                   parent = NULL,
                   ..., toolkit=guiToolkit()) {
  val <- .ginput(toolkit,
                 msg, text=text, title=title, icon=match.arg(icon), parent,
                 ...)
  if(!is.character(val))
    stop("Should return a character")
  val
}

##' generic for toolkit dispatch
##'
##' @export
##' @rdname ginput
.ginput <- function(toolkit,
                    msg, text="",
                    title = "Input",
                    icon = c("info", "warning", "error", "question"),
                    parent = NULL,
                    ...)
           UseMethod(".ginput")

################# gconfirm #################################
##' Constructor for modal dialog to get confirmation
##'
##' @inheritParams ginput
##' @return logical inidicating confirmation
##' @export
##' @rdname gconfirm
##' @seealso \code{\link{gmessage}}, \code{\link{gconfirm}}, 
##' \code{\link{gbasicdialog}}, \code{\link{galert}}
gconfirm = function(
  msg,
  title = "Confirm",
  icon = c("info", "warning", "error", "question"),
  parent=NULL,
  ..., toolkit=guiToolkit()) {
  
  out <- .gconfirm(toolkit,
                   msg=msg,
                   title=title,
                   icon=match.arg(icon),
                   parent=parent,
                   ...)
  if(!is.logical(out))
    stop("gconfirm returns a logical")
  out[1]
}

##' generic for toolkit dispatch
##'
##' @export
##' @rdname gconfirm
.gconfirm <- function(toolkit,
                      msg,
                      title = "Confirm",
                      icon = c("info", "warning", "error", "question"),
                      parent=NULL,
                      ...)
           UseMethod(".gconfirm")


##' Constructor for modal dialog that can contain an arbitrary widget
##'
##' The basic dialog is basically a modal window. To use there is a 3
##' step process: 1) Create a container by calling this constructor,
##' say \code{dlg}; 2) use \code{dlg} as a container for your
##' subsequent GUI; 3) set the dialog to be modal by calling
##' \code{visible(dlg)}. (One can't call \code{visible(dlg)
##' <- TRUE}.)
##' 
##' @param title title for window
##' @param parent parent to display by
##' @param do.buttons FALSE to suppress buttons when no parent
##' @param handler handler called when \code{Ok} button invoked
##' @param action passed to handler for OK button
##' @param ... ignored
##' @param toolkit toolkit
##' @return A \code{GBasicDialog} instance with a visible method
##' @seealso \code{\link{gmessage}}, \code{\link{gconfirm}}, 
##' \code{\link{gbasicdialog}}, \code{\link{galert}}
##' @export
##' @examples
##' \dontrun{
##' ## a modal dialog for editing a data frme 
##' fix_df <- function(DF, ...) {
##'   dfname <- deparse(substitute(DF))
##'   w <- gbasicdialog(..., handler=function(h,...) {
##'     assign(dfname, df[,], .GlobalEnv)
##'   })
##'   g <- ggroup(cont=w, horizontal=FALSE)
##'   glabel("Edit a data frame", cont=g)
##'   df <- gdf(DF, cont=g, expand=TRUE)
##'   size(w) <- c(400, 400)
##'   out <- visible(w)
##' }
##' 
##' m <- mtcars[1:3, 1:4]
##' fix_df(m)
##' }
gbasicdialog <- function(
                         title = "Dialog", 
                         parent = NULL,
                         do.buttons=TRUE,
                         handler = NULL, action = NULL,
                         ..., toolkit=guiToolkit()) {
  
  obj <- .gbasicdialog(toolkit,
                       title=title, parent=parent, do.buttons=do.buttons,
                       handler=handler, action=action, 
                       ...)

  check_return_class(obj, "GBasicDialog")
  obj
}
  
##' generic for toolkit dispatch
##'
##' @export
##' @rdname gbasicdialog
.gbasicdialog <- function(toolkit,
                          title = "Dialog", 
                          parent = NULL,
                          do.buttons=TRUE,
                          handler = NULL, action = NULL, 
                          ...)
           UseMethod(".gbasicdialog")


##' make basic dialog visible and modal
##'
##' We overrided the basic use of \code{visible} for  the
##' \code{gbasicdialog} container to have it become visible and modal
##' after this call. The better suited call \code{visible(dlg) <-
##' TRUE} does not work as wanted, for we want to capture the return
##' value.
##' @param obj dialog object
##' @return logical indicating which button was pushed (or TRUE if no buttons present)
##' @export
##' @rdname gbasicdialog
##' @method visible GBasicDialog
##' @S3method visible GBasicDialog
visible.GBasicDialog <- function(obj, ...) {
  obj$set_visible(TRUE)
}

##' dispose of dialog
##'
##' dispose method for a basic dialog
##' @return NULL
##' @export
##' @rdname gbasicdialog
##' @method dispose GBasicDialog
##' @S3method dispose GBasicDialog
dispose.GBasicDialog <- function(obj, ...) {
  obj$dispose()
}
