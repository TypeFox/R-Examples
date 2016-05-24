##' @include methods.R
NULL

##' An action constructor
##'
##' A action object encapsulates an action (a callback) adding textual
##' and graphic information. Actions may be proxied in buttons, menu
##' bars or tool bars.
##' @param label label for action
##' @param tooltip toolktip for actin
##' @param icon icon (stock icon name) for icon
##' @param key.accel keyboard accelerator. If given, parent must be specified.
##' @param handler handler to call when action is invoked
##' @param action values passed to parameterize action
##' @param parent parent window. Needed if keyboard accelerator used.
##' @inheritParams gcontainer
##' @export
##' @return a gaction instance
gaction <- function(
                    label, tooltip=NULL, icon = NULL, key.accel = NULL,
                    handler = NULL, action = NULL, parent=NULL, ...,
                    toolkit=guiToolkit()) {
  .gaction (toolkit,
            label, tooltip, icon, key.accel, handler, action, parent, ...
            )
}


##' generic for toolkit dispatch
##' 
##' @export
##' @rdname gaction
.gaction <- function(toolkit, label, tooltip = NULL, icon = NULL,
                     key.accel=NULL, handler = NULL, action = NULL,
                     parent=NULL, ... ) {
  UseMethod( '.gaction' )
}


