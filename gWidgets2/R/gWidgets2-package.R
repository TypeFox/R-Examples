##' The \pkg{gWidgets2} package provides a programming interface for
##' making graphical user interfaces within R. The package is a
##' rewrite of the \pkg{gWidgets} package, introducing a few external
##' changes but a significant number of internal ones. The package
##' relies on one of several underlying toolkit packages providing
##' access to the graphical libraries. These will include \pkg{RGtk2},
##' \pkg{tcltk}, \pkg{qtbase}, and a collection of browser widgets
##' provided by \code{ExtJS}. As of now, only \pkg{gWidgets2RGtk2} is
##' near completion.
##'
##' The package provides constructors to produce controls, the widgets
##' that a user interacts with; containers, GUI objects used to
##' organize the controls within a window; and dialogs, simple one-off
##' windows for gathering quick user feedback. These objects are
##' manipulated through various methods. The package provides a few
##' new generics and, as much as possible, leverages existing methods
##' for R.
##' 
##'
##' \subsection{Control constructors}{
##' 
##' Controls are created by constructors. The package API includes the
##' following. As much as possible these are implemented in the
##' toolkit packages, but there may be instances where that is not
##' possible.
##' 
##' \describe{
##' \item{\code{\link{gbutton}}}{Provides a basic button to initiate an action}
##' 
##' \item{\code{\link{gcalendar}}}{Provides a text entry area with selector for a date}
##' 
##' \item{\code{\link{gcheckbox}}}{Provides a labeled checkbox to allow a user to toggle a selection}
##' 
##' \item{\code{\link{gcheckboxgroup}}}{Provides a collection of
##' checkboxes allowing selection or 0, 1, or several from many items}
##' 
##' \item{\code{\link{gcombobox}}}{Provides a drop down list of
##' choices to select from and possible and entry area for free
##' response}
##' 
##' \item{\code{\link{gdf}}}{Provides a data frame editing widget}
##' 
##' 
##' \item{\code{\link{gedit}}}{Provides a single line text entry widget}
##' 
##' \item{\code{\link{ggraphics}}}{Provides an embeddable graphic device}
##' 
##' \item{\code{\link{gimage}}}{Provides a widget to hold images}
##' 
##' \item{\code{\link{glabel}}}{Provides a widget to hold labels for other controls}
##' 
##' \item{\code{\link{gmenu}}}{Provides menus for top-level windows and popup menus}
##' 
##' \item{\code{\link{gradio}}}{Provides a means to select one of many items}
##' 
##' \item{\code{\link{gseparator}}}{Provides a visual line to separate off parts of a window}
##' 
##' \item{\code{\link{gslider}}}{Provides a means to select one value from a (seeming) continuum of values}
##' 
##' \item{\code{\link{gspinbutton}}}{Provides means to select a value from s sequence of values}
##' 
##' \item{\code{\link{gstatusbar}}}{Provides a widget to display status messages in a top-level window}
##' 
##' \item{\code{\link{gtable}}}{Provides a widget to display tabular data for selection}
##' 
##' \item{\code{\link{gtext}}}{Provides a multiline text-editing widget}
##' 
##' \item{\code{\link{gtimer}}}{Provides a one-shot or repeatable timer}
##' 
##' \item{\code{\link{gtoolbar}}}{Provides toolbars for top-level windows}
##' 
##' \item{\code{\link{gtree}}}{Provides a display for heirarchicial data}
##' 
##' \item{\code{\link{gvarbrowser}}}{Provides a widget showing a shapshot of the current global workspace}
##' 
##' \item{\code{\link{gaction}}}{Provides a means to encapsulate actions for use with menu bars, tool bars and buttons.}
##' }
##'
##' Containers are used to organize controls with a window. The package provides the following:
##' 
##' \describe{
##' \item{\code{\link{gexpandgroup}}}{Provides a container with an option to disclose or hide its children}
##' 
##' \item{\code{\link{gframe}}}{Provides a framed box container}
##' 
##' \item{\code{\link{ggroup}}}{Provides a horizontal or vertical box container for packing in child components}
##' 
##' \item{\code{\link{glayout}}}{Provides a container to organize data by row and column cell}
##' 
##' \item{\code{\link{gnotebook}}}{Provides a notebook container}
##' 
##' \item{\code{\link{gpanedgroup}}}{Provides a divided container with adjustable divider}
##' 
##' \item{\code{\link{gstackwidget}}}{Provides a container like a notebook, but without tab labels}
##' 
##' \item{\code{\link{gwindow}}}{Provides a top-level window}
##'
##' }
##' }

##' \subsection{Dialog constructors}{
##'
##' Dialogs in \pkg{gWidgets2} are typically modal, meaning they block
##' input to any other window and the R process. They do not return
##' objects to be manipulated through methods, but rather return
##' values selected by the user.
##' 
##' \describe{
##' 
##' \item{\code{\link{gmessage}}}{Produces a simple dialog to display a message}
##'
##' \item{\code{\link{gconfirm}}}{Produces a dialog for a user to confirm an action}
##'
##' \item{\code{\link{ginput}}}{Provides a dialog to gather user in=put}
##'
##' \item{\code{\link{gbasicdialog}}}{Provides a means to produce general modal dialogs}
##'
##' \item{\code{\link{galert}}}{Provides a short transient message dialog}
##'
##' \item{\code{\link{gfile}}}{Provides a dialog to select a filename or directory name}
##' 
##' }
##' }
##' 
##' \subsection{Methods}{
##'
##' Except for dialogs, the constructors produce objects for which
##' several methods are defined that allow the programmer access to
##' getting and setting of the object state. For the most part these
##' are S3 methods. The actual returned object is a reference class
##' instance, as provided by an underlying toolkit. These may have
##' toolkit-specific methods defined as reference class methods (i.e.,
##' call them using \code{$meth_name}). Any such methods are
##' documented in the toolkit packages. 
##'
##' \describe{
##' \item{\code{\link{svalue}}, \code{\link{svalue<-}}}{
##'  The main new method. This is used to retrieve or set the main property associated with a widget
##' }
##' 
##' \item{\code{\link{enabled}}, \code{\link{enabled<-}}}{
##' A widget is enabled if it is sensitive to user input. Non-enabled widgets typically are rendered in a greyed out state. 
##' }
##' 
##' \item{\code{\link{visible}}, \code{\link{visible<-}}}{ The generic
##' idea of a visible widget is one that is drawn. However, several
##' classes override this to mean part of the widget is visible or not
##' visible.  }
##' 
##' \item{\code{\link{focus}}, \code{\link{focus<-}}}{ A widget with
##' focus receives any keyboard input.  }
##' 
##' \item{\code{\link{editable}}, \code{\link{editable<-}}}{A widget is editable if it can receive keyboard input.}
##' 
##' \item{\code{\link{font}}, \code{\link{font<-}}}{The font for an object is specifed through this method using a convention illustrated in the help page.}
##' 
##' \item{\code{\link{size}}, \code{\link{size<-}}}{The size of a widget is retrieved or requested through these methods}
##' 
##' \code{\link{tooltip}}, \code{\link{tooltip<-}}{A tooltip provides contextual information when a mouse hovers over an object}
##' 
##' \item{\code{\link{undo}}, \code{\link{redo}}}{Some widgets support an undo and redo stack}
##' 
##' \item{\code{\link{isExtant}}}{A method to check if the GUI part of
##' a widget still exists. (A constructor produces an R object and GUI
##' object through the toolkit.)}
##' 
##' \item{\code{\link{tag}}, \code{\link{tag<-}}}{A method used to set
##' attributes for an object that are stored in an environment so that
##' they are passed by reference, not copy. This allows event handlers
##' to manipulate an object's attributes outside the scope of the
##' callback.}
##' 
##' \item{\code{\link{getToolkitWidget}}}{Returns the underlying
##' toolkit object that is packaged into a \pkg{gWidgets2} object}
##' 
##' 
##' \item{\code{\link{add}}}{Method used to add a child component to a
##' parent container}
##' 
##' \item{\code{\link{delete}}}{Method used to delete a component from
##' its parent}
##' 
##' \item{\code{\link{dispose}}}{Method used to delete a component}
##' 
##' }
##' 
##' The package overloads some familar R methods.
##'
##' \describe{
##'
##' \item{\code{length}, \code{length<-}}{Returns the length of an
##' object, typically related to the number of children a container
##' has, or the length of the items that a user can selection from.}
##'
##' \item{\code{dim}}{Used to return row and column size information
##' as applicable.}
##'
##' \item{\code{names}, \code{names<-}}{Used to set the names
##' associated to an object. These may be column names in the table
##' widget, or tab names in the notebook container.}
##'
##' \item{\code{dimnames}, \code{dimnames<-}}{Used to set row and
##' column names, as applicable.}
##'
##' \item{\code{[}, \code{[<-}}{Used to update the underlying items
##' that a selection widget offers. Also used to specify layout in
##' \code{glayout}}
##'
##' \item{\code{update}}{Call to update the state of a widget, when
##' applicable.}
##'
##' } }
##' 
##'
##' \subsection{Event handlers}{
##'
##' Graphical User Interfaces are made interactive by assigning a
##' function (a callback) to be called when some event happens. In
##' \pkg{gWidgets2} the \code{addHandlerXXX} methods are defined to
##' assign this callback to a type of event specified through the
##' \code{XXX}, detailed below. The generic \code{addHandlerChanged}
##' is the most common event for a widget. This event can also have a
##' handler specified through the \code{handler} argument of the
##' widget constructor.
##'
##' In \pkg{gWidgets2} handlers are functions which when called are
##' passed a list as the first argument, and possibly toolkit-specific
##' arguments for subsequent arguments. As such the signature
##' typically looks like \code{(h,...)}, where the list \code{h} has
##' components \code{obj}, containing a reference to the widget
##' emitting the event and \code{action} passing in any information
##' specified to the \code{action} argument. Some events also pass
##' back extra information, such as \code{x} and \code{y} for
##' position, or \code{key} for key events, as appropriate.
##' 
##' \describe{
##'
##' \item{\code{\link{addHandlerChanged}}}{Assigns callback for the
##' most generic event}
##' 
##' \item{\code{\link{addHandlerClicked}}}{Assigns callback for a
##' mouse click event}
##' 
##' \item{\code{\link{addHandlerDoubleclick}}}{Assigns callback for a
##' mouse double-click event}
##' 
##' \item{\code{\link{addHandlerRightclick}}}{Assigns callback for a
##' mouse right-click event}
##' 
##' \item{\code{\link{addHandlerColumnclicked}}}{Assigns callback for
##' a column-click event}
##' 
##' \item{\code{\link{addHandlerColumnDoubleclicked}}}{Assigns
##' callback for a column-double-click event}
##' 
##' \item{\code{\link{addHandlerColumnRightclicked}}}{Assigns callback
##' for a column-right-click event}
##' 
##' \item{\code{\link{addHandlerSelect}}}{Assigns callback when the
##' underlying selection is changed}
##' 
##' \item{\code{\link{addHandlerFocus}}}{Assigns a callback for when a
##' widget receives focus}
##' 
##' \item{\code{\link{addHandlerBlur}}}{Assigns a callback for when a
##' widget loses focus}
##' 
##' \item{\code{\link{addHandlerDestroy}}}{Assigns a callback for when
##' a widget is destroyed}
##' 
##' \item{\code{\link{addHandlerUnrealize}}}{For \code{gwindow} this
##' is called before the destroy event and may stop that from
##' happening. }
##' 
##' \item{\code{\link{addHandlerExpose}}}{Assigns callback to be
##' called when a widget is exposed}
##' 
##' \item{\code{\link{addHandlerKeystroke}}}{Assigns callback to be
##' called when a key event occurs}
##' 
##' \item{\code{\link{addHandlerMouseMotion}}}{Assigns callback to be
##' called when a mouse moves over a widget}
##'
##' \item{\code{\link{addHandler}}}{Base method to add a callback
##' though rarely called, as it is toolkit specific}
##' 
##' \item{\code{\link{addHandlerIdle}}}{Assign a callback to be called
##' at periodic intervals. See also \code{\link{gtimer}}}
##' 
##' \item{\code{\link{addPopupMenu}}}{Add a popup menu}
##' 
##' \item{\code{\link{add3rdmousePopupMenu}}}{Add a popup menu for the
##' right mouse (context menu)}
##' 
##' \item{\code{\link{addDropSource}}}{Specify widget as a source
##' (drag area) for drag and drop}
##' 
##' \item{\code{\link{addDropTarget}}}{Specify widget as a target
##' (drop area) for drag and drop}
##' 
##' \item{\code{\link{addDragMotion}}}{Assign callback for event that
##' a drag event crosses a widget}
##' 
##' \item{\code{\link{blockHandlers}}, \code{\link{blockHandler}}}{Block all
##' handlers for a widget (or by single ID)}
##' 
##' \item{\code{\link{unblockHandlers}}, \code{\link{unblockHandler}}}{Unblock any
##' blocked handlers (or by single ID)}
##' 
##' \item{\code{\link{removeHandler}}}{Remove a handler by it ID}
##'
##' }
##' }
##' 
##' @author
##' John Verzani \email{jverzani@@gmail.com}
##'
##' Maintainer: John Verzani \email{jverzani@@gmail.com}
##' @aliases gWidgets2-package
##' @name gWidgets2-package
##' @docType package
##' @title gWidgets2. An API for programming GUIs
##' @keywords package
NULL
