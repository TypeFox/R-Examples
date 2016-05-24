##' @include guiToolkit.R
##' @include misc.R
NULL

## simple message function
define_me <- function(...) {
  curcall <- as.character(sys.call()[[1]])[3]
  if(!is.null(getOption("gWidgets2:debug")))
    message(sprintf("Method %s not defined for class %s\n",
                    curcall,
                    class(get(".self"))          # issue with warning
                    ))
}

##' Reference class to create an observer of an observable object
##'
##' An observer can be observed
##' @param ... passed to constructor
##' @aliases Observer 
##' @rdname S4-classes
##' @name Observer-class
Observer <- setRefClass("Observer",
                        fields=list(
                          o = "ANY",   
                          obj="ANY"
                          ),
                        methods=list(
                          initialize=function(o=NULL, obj=NULL) {
                            initFields(o=o, obj=obj)
                            callSuper()
                          },
                          update=function(...) {
                            "Call self."
                            o(obj, ...)
                          }
                          )
                        )

##' Handler is a special observer with obj and actino passed as first argument
##'
##' @aliases Handler
##' @rdname S4-classes
##' @name Handler-class
Handler <- setRefClass("Handler",
                       contains="Observer",
                       fields=list(
                         action="ANY"
                         ),
                        methods=list(
                          initialize=function(o=NULL, obj=NULL, action=NULL) {
                            initFields(action=action)
                            callSuper(o=o, obj=obj)
                          },
                          update=function( ..., extra_args) {
                            "Call self."
                            h <- list(obj=obj, action=action)
                            if(!missing(extra_args))
                              h <- merge.list(h, extra_args, overwrite=FALSE)
                            else
                              h <- merge.list(h, list(...), overwrite=FALSE)
#                            if(!missing(extra_args)) {
#                              h <- merge.list(h, extra_args, overwrite=FALSE)
#                            }
                            o(h, ...)
                          }
                          )
                        )

##' constructor for handler object
##'
##' @param receiver object receiving event
##' @param handler function to call
##' @param action used to parametrize handler call
##' not exported, call using :::
observer <- function(receiver, handler, action=NULL) 
  Handler$new(handler, receiver, action)


##' Observable class sets up objects that can be observed. Inherited by template
##'
##' @aliases Observable
##' @exportClass Observable
##' @rdname S4-classes
##' @name Observable-class
Observable <- setRefClass("Observable",
                          fields=list(
                            ..observers="list",
                            ..blocked_observers = "list",
                            ..blocked="integer"
                            ),
                          methods=list(
                            add_observer=function(o, signal="DEFAULT") {
                              "Add an observer. Return id for block/remove/..."
                              if(!is(o, "Observer"))
                                stop("Not an observer")
                              l <- ..observers
                              if(is.null(l[[signal]]))
                                l[[signal]] <- list(o)
                              else
                                l[[signal]] <- c(l[[signal]], o)
                              ..observers <<- l
                              list(signal=signal, o=o)
                            },
                            ## these rmove/block/unblock all observers
                            remove_observers=function(signal) {
                              if(missing(signal))
                                ## remove all
                                ..observers <<- list()
                              else
                                ..observers[[signal]] <<- NULL
                            },
                            block_observers=function() {
                              "Block all observers"
                              if(is("..blocked", "uninitializedField") || length(..blocked) == 0) {
                                ..blocked <<- 1L
                              } else {
                                ..blocked <<- ..blocked + 1L
                              }
                            },
                            unblock_observers=function() {
                              "Remove block of all observers. Keeps count, so may need to call again"
                              if(is("..blocked", "uninitializedField") || length(..blocked) == 0) {
                                ..blocked <<- 0L
                              } else {
                                ..blocked <<- max(..blocked - 1L, 0L)
                              }
                              invisible(..blocked)
                            },
                            ## These block/unblock one at a time
                            remove_observer=function(id) {
                              "Remove observer"
                              if(!is(id$o, "Observer"))
                                stop("Call with an observer id")
                              
                              signal <- id$signal
                              ind <- unlist(lapply(..observers[[signal]], identical, id$o))
                              if(any(ind))
                                ..observers[[signal]][[which(ind)]] <<- NULL
                              
                            },
                            block_observer=function(id) {
                              "Block observers. If o missing, block all"
                              if(missing(id) || is.null(id)) {
                                block_observers()
                              } else {
                                if(is.null(..blocked_observers[[id$signal]]))
                                  ..blocked_observers[[id$signal]] <<- list(id$o)
                                else
                                  ..blocked_observers[[id$signal]] <<-
                                    c(..blocked_observers[[id$signal]], id$o)
                              }
                            },
                            unblock_observer=function(id) {
                              "Unblock observer. If id missing, unblock global block"
                              if(missing(id) || is.null(id)) {
                                unblock_observers()
                              } else {
                                signal <- id$signal
                                ind <- unlist(lapply(..blocked_observers[[signal]], identical, id$o))
                                if(any(ind))
                                  ..blocked_observers[[signal]][[which(ind)]] <<- NULL
                              }
                            },
                            notify_observers=function(..., signal="DEFAULT") {
                              "Call each non-blocked observer"
                              if(!is("..blocked", "uninitializedField") && length(..blocked) && ..blocked > 0)
                                return()
                              QT <- lapply(..observers[[signal]], function(o) {
                                ind <- lapply(..blocked_observers[[signal]], function(i) identical(i, o))
                                if(!any(unlist(ind))) 
                                  o$update(...)
                              })
                            }
                            )
                          )

##' Basic interface for a widget. These are methods referenced by the S3 methods
##'
##' This interface is inherited by the base GComponent classes in the
##' toolkit implementations. The methods defined here are referenced
##' by the S3 methods. For exampe, \code{svalue} dispatches to
##' \code{get_value} or \code{get_index}.
##'
##' We combine both widget and container methods here. It isn't
##' perfect, but they do share quite a bit. Perhaps, we could make the
##' container class subclass the basic interface.
##' @exportClass BasicToolkitInterface
##' @aliases BasicToolkitInterface
##' @rdname S4-classes
##' @name BasicToolkitInterface-class
BasicToolkitInterface <- setRefClass("BasicToolkitInterface",
                                     contains="Observable",
                                     fields=list(
                                       toolkit="ANY",
                                       widget="ANY",
                                       block="ANY",
                                       parent="ANY", # NULL for gwindow, else parent container
                                       default_expand="LogicalCharacterOrNULL",
                                       default_fill="LogicalCharacterOrNULL",
                                       coerce_with="FunctionOrNULL"
                                       ),
                                     methods=list(
                                       ## (drop=NULL, ...), (value, drop=TRUE, ...)
                                       get_value=function(drop=NULL, ...) {
                                         "Get main value of widget. From `svalue` when index = FALSE or NULL"
                                         define_me()
                                       }, # svalue
                                       set_value=function(value, ..., drop=NULL) {
                                         "for `svalue<-` when index = FALSE or NULL"
                                         define_me()
                                       }, 
                                       get_index=function(drop=NULL, ...) {
                                         "svalue; index=TRUE"
                                         define_me()
                                       },
                                       set_index=function(value, ..., drop=NULL) {
                                         define_me()   # svalue <-; index=TRUE
                                       },
                                       ## (i, j, ..., drop=NULL)
                                       get_items=function(i,j, .., drop=TRUE) {
                                         define_me()   # [
                                       },
                                       ## (value, i, j, ...)
                                       set_items=function(value, i, j, ...) {
                                         define_me()   # [<-
                                       },
                                       ## () and (value)
                                       get_enabled=function() {
                                         "is widget sensistive to user input"
                                         define_me()
                                       },
                                       ## enabled<-                                       
                                       set_enabled=function(value, ...) {
                                         "specify with logical if widget is sensistive to user input"
                                         define_me
                                       }, 
                                       get_visible=define_me, # visible
                                       set_visible=define_me, # visible<-
                                       get_editable=define_me, # editable
                                       set_editable=define_me, # editable<-
                                       get_focus=define_me,    # foucs
                                       set_focus=define_me,    # focus<-
                                       get_font=define_me,    # font
                                       set_font=define_me,    # font<-
                                       get_length=define_me,  # length
                                       set_length=define_me,  # length<-
                                       ##
                                       get_dim=define_me,     # dim
                                       get_names=define_me,   # names
                                       set_names=define_me,   # names<-
                                       get_dimnames=define_me, # dimnames
                                       set_dimnames=define_me, # dimnames <-
                                       get_attr=define_me,    # tag
                                       set_attr=define_me,    # tag<-
                                       ##
                                       set_invalid=define_me, # for validation
                                       is_invalid=function() {}, # for validation
                                       update_widget=define_me, # update
                                       is_extant=function() TRUE,   # isExtant
                                       undo=define_me,          # undo
                                       redo=define_me,          # redo
                                       add_child=define_me,     # add child to container (if present)
                                       set_parent=function(parent) parent <<- parent,
                                       ## (signal, handler, action=NULL, decorator, emitter)
                                       add_handler=function(signal, handler, action, ...) {
                                         "Add a handler to be called for the event indicated by signal"
                                         define_me()
                                       },
                                       ## (handler, action, ...)
                                       add_handler_changed=define_me,
                                       add_handler_clicked=define_me,
                                       add_handler_double_clicked=define_me,
                                       add_handler_right_clicked=define_me,
                                         add_handler_control_clicked=define_me,
                                         add_handler_shift_clicked=define_me,
                                         
                                       add_handler_column_clicked=define_me,
                                       add_handler_column_double_clicked=define_me,
                                       add_handler_column_right_clicked=define_me,
                                       add_handler_select=function(handler, action=NULL, ...) { # selection vs. selection_changed
                                         add_handler_changed(handler, action=NULL, ...)
                                       },
                                       add_handler_selection_changed=define_me,
                                       add_handler_focus=define_me,
                                       add_handler_blur=define_me,
                                       add_handler_destroy=define_me,
                                       add_handler_unrealize=define_me,
                                       add_handler_expose=define_me,
                                       add_handler_keystroke=define_me,
                                       add_handler_mouse_motion=define_me,
                                       add_popup_menu=define_me,
                                       add_3rd_mouse_popup_menu=define_me,
                                       add_drop_source=define_me,
                                       add_drop_target=define_me,
                                       add_drag_Motion=define_me,
                                       emit_signal=define_me, # gSignalEmit, ... to invoke. Signal missing do change one
                                       ## show method
                                       show=function() {
                                         cat(sprintf("An object of class %s\n", class(.self)[1]))
                                       }
                                       ))


## needed for internal guys: gfilter, gdfnotebook, ggraphicsnotebook
## Causes warnings when installin gWidgets2XXX packages about conflicting classes
## 
## @exportClass GComponent
## @rdname S4-classes
## @name GComponent-class
##GComponent <- setRefClass("GComponent",
##                          contains="BasicToolkitInterface")


## Class for default widgets
## Used by g*notebook, and gfilter
GDefaultWidget <- setRefClass("GDefaultWidget",
                                contains="BasicToolkitInterface")
                          


## For odd reasons, we don't want to use GComponent methods here.

##' Return items
##'
##' @export
##' @rdname gWidgets2-S3methods
##' @method [ GDefaultWidget
"[.GDefaultWidget" <- function(x, i, j, ...) x$get_items(i, j, ...)


##' Set object's items
##'
##' @export
##' @usage \method{[}{GDefaultWidget} (x, i, j, ...) <- value
##' @rdname gWidgets2-S3methods
##' @method [<- GDefaultWidget
##' @S3method [<- GDefaultWidget
"[<-.GDefaultWidget" <- function(x, i, j, ..., value) x$set_items(value, i, j, ...)

##' method for getWidget
##'
##' @rdname getToolkitWidget
##' @export
##' @method getWidget GDefaultWidget
##' @S3method getWidget GDefaultWidget
getWidget.GDefaultWidget <- function(obj) getWidget(obj$widget)

