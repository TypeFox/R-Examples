##' @include guiComponents



##' Alert dialog to display transient messages
##' 
##' @param message main message. 
##' @param title Title (may not be displayed)
##' @param delay length of time (in seconds) to display
##' @param parent parent object to show near
##' @param ... ignored
##' @param toolkit toolkit
##' @rdname gWidgets-dialogs
galert = function(
  message,
  title = "message",
  delay = 3,
  parent = NULL, 
  ..., toolkit=guiToolkit()) {
  .galert(toolkit,message, title, delay=delay, parent=parent,
            ...)
}
##' generic for toolkit dispatch
##' @alias galert
setGeneric(".galert",
           function(toolkit,
                    message, title="message", delay=3, parent=NULL, ...)
           standardGeneric(".galert"))


##' Constructor for modal message dialog
##' 
##' @export
##' @param message Character. message to display. 
##' @param title Character. Title
##' @param icon What icon to show
##' @param parent Hint as to where to display
##' @param handler called if Ok button clicked
##' @param action passed to handler
##' @param ... ignored
##' @param toolkit toolkit
##' @rdname gWidgets-dialogs
gmessage = function(
  message,
  title = "message",
  icon = c("info", "warning", "error", "question"),
  parent=NULL,
  handler = NULL, action = NULL,
  ..., toolkit=guiToolkit()) {
  .gmessage(toolkit,message, title, icon, parent, handler, action,
            ...)
}

##' generic for toolkit dispatch
##' @alias gmessage
setGeneric(".gmessage",
           function(toolkit,
                    message="", title="", icon="",
                    parent = NULL,
                    handler=NULL, action=NULL, ...)
           standardGeneric(".gmessage"))

############## ginput ####################################

##' Constructor for modal dialog to collect a line of text
##'
##' @export
##' @param message Character. Message to display.
##' @param text Character. Initial text
##' @param title Character. Title of window
##' @param icon which icon to display
##' @param parent gives hint as to where to place dialog
##' @param handler called on \code{Ok}
##' @param action passed to handler
##' @param ... ignored
##' @param toolkit toolkit
##' @rdname gWidgets-dialog
ginput <- function(
  message, text="",
  title = "Input",
  icon = c("info", "warning", "error", "question"),
  parent = NULL,
  handler = NULL, action = NULL,
  ..., toolkit=guiToolkit()) {
  .ginput(toolkit,
          message, text=text, title=title, icon, parent, handler, action, 
          ...)
}

##' generic for toolkit dispatch
##' @alias ginput
##' @export
##' @rdname gWidgets-dialog
setGeneric(".ginput",
           function(toolkit,
                    message=message, text=text, title=title, icon=icon,
                    parent = parent,
                    handler=handler, action=action, ...)
           standardGeneric(".ginput"))

################# gconfirm #################################
##' Constructor for modal dialog to get confirmation
##'
##' @export
gconfirm = function(
  message,
  title = "Confirm",
  icon = c("info", "warning", "error", "question"),
  parent=NULL,
  handler = NULL, action = NULL,
  ..., toolkit=guiToolkit()) {
  .gconfirm(toolkit,message=message, icon=icon, parent=parent, handler=handler, action=action,
            ...)
}

##' generic for toolkit dispatch
##' @alias gconfirm
setGeneric(".gconfirm",
           function(toolkit,
                    message=message, title=title, icon=icon,
                    parent = parent,
                    handler=handler, action=action, ...)
           standardGeneric(".gconfirm"))

################# gbasicdialog #################################

## define subclass for basic dialog
setClass("guiDialog",
         contains="guiContainer",
         prototype=prototype(new("guiContainer"))
         )




##' Constructor for modal dialog that can contain an arbitrary widget
##'
##' The basic dialog is basically a modal window. To use there is a 3
##' step process: 1) Create a container by calling this constructor,
##' say \code{dlg}; 2) use \code{dlg} as a container for your
##' subsequent GUI; 3) set the dialog to be modal by calling
##' \code{visible(dlg, set=TRUE)}. (One can't call \code{visible(dlg)
##' <- TRUE}.)
##' 
##' @export
##' @param title title for window
##' @param widget widget to add (Only if toolkit supports it)
##' @param parent parent to display by
##' @param do.buttons FALSE to suppress buttons when no parent
##' @param handler handler called when \code{Ok} button invoked
##' @param action passed to handler
##' @param ... ignored
##' @param toolkit toolkit
##' @rdname gWidgets-dialogs
gbasicdialog <- function(
                         title = "Dialog", widget,
                         parent = NULL,
                         do.buttons=TRUE,
                         handler = NULL, action = NULL,
                         ..., toolkit=guiToolkit()) {
  
  if(missing(widget)) {
    obj <- .gbasicdialognoparent(toolkit, title, parent, handler, action,  ..., do.buttons=do.buttons)
    obj <- new( 'guiDialog',widget=obj,toolkit=toolkit) 
  } else {
    obj <- .gbasicdialog(toolkit,
                  title=title, widget=widget,parent=parent,
                  handler=handler, action=action, 
                         ...)
  }
  return(obj)
}

##' generic for toolkit dispatch
##' 
##' @alias gbasicdialog
setGeneric(".gbasicdialog",
           function(toolkit,
                    title = "Dialog", widget, parent,
                    handler = NULL, action = NULL, 
                    ...)
           standardGeneric(".gbasicdialog"))

##' generic for toolkit dispatch when there is no parent
##' 
##' @alias gbasicdialog
setGeneric(".gbasicdialognoparent",
           function(toolkit,
                    title = "Dialog",  parent,
                    handler = NULL, action = NULL, 
                    ...)
           standardGeneric(".gbasicdialognoparent"))



