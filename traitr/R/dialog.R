##  Copyright (C) 2010 John Verzani
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 2 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  A copy of the GNU General Public License is available at
##  http://www.r-project.org/Licenses/

##' @include itemgroup.R
roxygen()

##################################################
##' A Dialog wraps a top-level window around a collection of items which may be item groups
##' 
##' One can specify the parent when making the UI
##' Buttons are specified through the buttons property
##' @export
Dialog <- ItemGroup$proto(class=c("Dialog", ItemGroup$class),
                          ## window title
                          .doc_title=paste(
                            desc("Title of dialog window")
                            ),
                          title="Dialog title",
                          ## called by default help button handler
                          .doc_help_string=paste(
                            desc("String called by default Help button handler")
                            ),
                          help_string="",
                          ## !is.null to create a status bar
                          .doc_status_text=paste(
                            desc("If non <code>NULL</code> this text appears in bottom statusbar.",
                                 "When there is a status bar, see also <code>set_status_text</code>")
                            ),
                          status_text=NULL,
                          ## If !is.null, specifies list to pass to gmenu
                          .doc_menu_list=paste(
                            desc("If non <code>NULL</code> this list is passed to <code>gmenu</code>",
                                 "to make a menu bar")
                            ),
                          menu_list = NULL,
                          ## If !is.null, specifies list to pass to gtoolbar
                          .doc_toolbar_list=paste(
                            desc("If non <code>NULL</code> this list is passed to <code>gtoolbar</code>",
                                 "to make a toolbar.")
                            ),
                          toolbar_list=NULL,
                          ## Buttons are OK, Cancel, Help, or any other name. The button handler
                          ## is NAME_handler defined as follows
                          .doc_buttons=paste(
                            desc("A property listing the buttons to place into dialog.",
                                 "A button with name <code>name</code> will call",
                                 "a method <code>name_handler</code>. Default handlers are defined for",
                                 "OK, Cancel, Help, Undo, Redo, but usually the OK handler would be ",
                                 "redefined for a dialog instance.",
                                 "The default set of buttons is OK, Cancel, Help.",
                                 "Use button named <code>SPACE</code> to add 12px space between",
                                 "buttons. Use button named <code>SPRING</code>to add spring between buttons (",
                                 "pushes buttons to left and right.",
                                 "set the property default_button to make a button the default.",
                                 "Use character(0) for no default buttons"
                                 )
                            ),
                          buttons=c("OK", "Cancel","SPACE", "Help"),#, "Undo","Redo"),
                          doc_default_button=paste(
                            desc("Name of default button. Leave empty (<code>NULL</code>) for none.")
                            ),
                          default_button=NULL,
                          ## Default button handlers. OK_handler needs to be
                          ## overridden in the instance
                          doc_OK_handler=paste(
                            desc("Handler called when 'OK' button is clicked -- if 'OK' button among <code>buttons</code>property.")
                            ),
                          OK_handler=function(.) {
                            print(.$to_R())
                          },
                          doc_Cancel_handler=paste(
                            desc("Handler called when 'Cancel' button is clicked -- if 'Cancel' button among <code>buttons</code>property.")
                            ),                          
                          Cancel_handler = function(.) {
                            dispose(.$get_widget('toplevel'))
                          },
                          doc_Help_handler=paste(
                            desc("Handler called when 'Help' button is clicked -- if 'Help' button among <code>buttons</code>property.")
                            ),
                          Help_handler = function(.) {
                            gmessage(c("Help", .$help_string), icon="info", parent=.$get_widget('toplevel'))
                          },
                          doc_Undo_handler=paste(
                            desc("Handler called when 'Undo' button is clicked -- if 'Undo' button among <code>buttons</code>property.")
                            ),
                          Undo_handler=function(.) .$undo(),
                          doc_Redo_handler=paste(
                            desc("Handler called when 'Redo' button is clicked -- if 'Redo' button among <code>buttons</code>property.")
                            ),
                          Redo_handler=function(.) .$redo(),
                          ## on_realized is called after make_gui
                          on_realized = function(.) {
                            .$update_ui()
                          },
                          
                          ## This is a list, not named, of Items or ItemLists
                          ## names are found from items.
                          items=list(),  # Items or ItemGroup
                          widget_list=list(), # where widgets are stored
                          get_widget=function(.,key) .$widget_list[[key]],
                          ## make_gui makes the dialog (make_ui has different arguments)
                          ## make the dialog.
                          ## parent can be NULL, a gwindow instance, or a dialog instance
                          ## set visible to FALSE to suppress drawing. Use .$visible(TRUE) to show
                          .doc_make_gui=paste(
                            desc("Make gui for the dialog",
                                  "Creates top level window located near <code>parent</code>, if given.",
                                  "Creates buttons. When clicked, these call handlers of similar name. That is,",
                                  "'OK' button calls 'OK_handler'. Button name is stripped of non letters before call,",
                                  "so 'Next >>' would call 'Next_handler'.",
                                  "Call on_realized method, if defined after construction"
                                  )
                            ),
                          ### JV: XXX Add in attr for gwindow call
                          make_gui=function(., gui_layout=.$make_default_gui_layout(), parent=NULL, visible=TRUE) {
                            .$init_model() # initialize model ## also called in aDialog?

                            
                            if(!is.null(parent)) {
                              if(is.proto(parent) && exists("traitr", parent) && parent$is("Dialog"))
                                parent <- parent$get_widget("toplevel")
                            }
                            widgets <- list() # where widgets are to store
                            ## make window
                            if(.$has_slot("toplevel")) {
                              w <- .$toplevel
#                              svalue(w) <- .$title
#                              visible(w) <- FALSE
                            } else {
                              w <- gwindow(.$title, parent=parent, visible=FALSE)
                            }
                            widgets[['toplevel']] <- w

                            ## uses gWidgets -- not actionItems 
                            ## menu
                            if(!is.null(.$menu_list))
                              widgets[['menu']] <- gmenu(.$menu_list, container = w)
                            if(!is.null(.$toolbar_list))
                              widgets[['toolbar']] <- gtoolbar(.$toolbar_list, container = w)
                            
                            ## main part
                            g <- ggroup(horizontal=FALSE, cont=widgets[['toplevel']])


                            .$next_method("make_gui")(., container=g, gui_layout=gui_layout)



                            ## buttons
                            ## handlers called via buttonName_handler method defined
                            bg <- ggroup(cont=g)
                            sapply(.$buttons, function(i) {
                              if(i == "SPACE") {
                                addSpace(bg, 10)
                                return()
                              } else if(i == "SPRING") {
                                addSpring(bg)
                                return()
                              }
                              
                              widgets[[i]] <<- (b <- gbutton(i, cont=bg))
                              addHandlerClicked(b, function(h,...) {
                                . <- h$action$self
                                button_name <- h$action$button_name
                                ## strip off all but characters
                                button_name <- gsub("[^a-zA-Z0-9]","",button_name)
                                .$do_call(sprintf("%s_handler", button_name), list())
                              }, action=list(self=., button_name=i))
                              if(!is.null(.$default_button) &&
                                 .$default_button == i)
                                defaultWidget(i) <- TRUE
                            })
                            ## status bar if requested
                            if(!is.null(.$status_text)) {
                              widgets[['statusbar']] <- gstatusbar(.$status_text, cont=widgets[['toplevel']])
                              .$.doc_set_status_text <- paste(
                                desc("Method to update text in status bar"),
                                param("value","New value for status bar")
                                )
                              .$set_status_text <- function(., value) {
                                sb <- .$get_widget('statusbar')
                                svalue(sb) <- value
                              }
                            }
                            ## set visible if requested
                            if(visible)
                              visible(widgets[['toplevel']]) <- visible

                            ## set widgets
                            .$widget_list <- widgets
                            ## call hook for realized
                            .$do_call("on_realized")

                            .$assign_if_null("model_value_changed", function(.) .$update_ui())
#                            .$init() ## initialize controller, ... ## called in ItemGroup

                            ## we listen to ourselves
                            ## already in init_controller for 
#                            .$add_observer(.)

                            invisible()
                          },
                          ## close the gui
                          .doc_close_gui=paste(
                            desc("Method call to delete GUI window")
                            ),
                          close_gui = function(.) {
                            l <- .$widget_list
                            try(dispose(l$toplevel), silent=TRUE)
                          },
                          ## toggle visibility of top level window of dialog
                          .doc_visible=paste(
                            desc("Method call to toggle visibility of toplevel widget"),
                            param("value","Logical indicating if window should be visible")
                            ),
                          visible=function(., value=TRUE) {
                            widget <- .$get_widget('toplevel')
                            visible(widget) <- as.logical(value)
                          },
                          ## add in undo/redo then call for itemGroup
                          update_ui=function(.) {
                            ## undo/redo buttons
                            undo <- .$get_widget("Undo");
                            redo <- .$get_widget("Redo");
                            if(!is.null(undo) && isExtant(undo))
                              enabled(undo) <- .$undo_can_undo()
                            if(!is.null(redo) && isExtant(redo))
                              enabled(redo) <- .$undo_can_redo()

                            .$next_method("update_ui")(.)
                          } 
                          )

##' Create a Dialog instance
##'
##' A dialog is like an item group, in that it combines items
##' into a model. However, an item group is meant to be incorporated
##' into other GUIS, whereas a dialog creates its own window and
##' decorations. A dialog has default buttons, and options for adding
##' in menubars, toolbars, and statusbars. The choice of buttons can
##' be specified at construction.
##' \cr
##' 
##' Methods:
##' 
##' The main method that a dialog has is
##' its \code{OK_handler} which is a method called when the "OK"
##' button is clicked (one of the default buttons). 
##'
##' The getters and setters for the main value for an item are
##' \code{get_NAME} and \code{set_NAME}, where \code{NAME} is the
##' item name. The name is specified when the item is constructed
##' (through its \code{name} property) or more conveniently, taken
##' from the name of the component in the \code{items} list that
##' defines the items for  dialog or item group.
##' 
##' The method \code{to_R} returns the items' values as a list (useful in combination with \code{do.call}).
##'
##' The method \code{get_item_by_name} returns an item object by its
##' name. (Names should be unique.) This is useful if more properties
##' than the main one for an item are needed to be set. (The main
##' value is set via the setter.) The example shows how the validate
##' property of some items can be set.
##' 
##' The method \code{is_valid} is \code{TRUE} if all items validate and \code{FALSE} otherwise.
##'
##' The method \code{model_value_changed(.)} is called whenever an
##' item property is changed either through the GUI. A dialog observes
##' itself.
##'
##' For each item one can listen for changes with the method \code{property_NAME_value_changed(., value, old_value)}.
##' 
##' Properties that are of interest:
##' \enumerate{
##' 
##' \item{\code{status_text}}{ If non-NULL, when GUI is drawn, a status bar will be made with this text. The method \code{set_status_text} can be used to update the status}
##' 
##' \item{\code{menu_list}}{ A menu list to specify a menubar. (See \code{\link{gmenu}}.)}
##' 
##' \item{\code{toolbar_list}}{ A menu list to specify a toolbar. (See \code{\link{gtoolbar}}.)}
##' 
##' \item{\code{buttons}}{ A list of buttons names. The default is
##' \code{c("OK", "SPACE", "Cancel", "Help")}. The special names
##' \code{SPACE} and \code{SPRING} adjust their positioning, otherwise
##' the values are button names. When a button is clicked, the handler
##' \code{buttonname_handler} is called, where the buttonname is
##' stripped on non-alphanumeric characters. The basic buttons and
##' \code{Redo} and \code{Undo} have default handlers. Likely, only
##' \code{OK_handler} will need redefining. The property
##' \code{default_button} can be specified to make a button the
##' default one (so that it is activated when a user presses the enter
##' key).}
##' }
##' @param items List of item instances to create the model for the
##' dialog object. May also be an item group
##' (\code{\link{anItemGroup}}).
##' @param title Title of dialog
##' @param help_string String for default Help button
##' @param buttons Character vector of button names. "OK","Cancel","Help","Undo","Redo" are some standard ones.
##'        "SPACE" and "SPRING" adjust the layout.
##' @param ... How to pass in other properties and methods of the dialog object. For example \code{OK_handler}.
##' 
##' 
##' @return Returns a proto object. See its \code{show_help} method for details.
##' @export
##' @examples
##' ##
##' ## a Basic example
##' dlg <- aDialog(items=list(
##'                  a = numericItem(0),
##'                  b = stringItem("a")
##'                  ),
##'                title="The title",
##'                help_string="Help on this dialog"
##'                )
##' \dontrun{dlg$make_gui()}
##' ##
##' ##
##' ## example with model_value_changed
##' plotIt <- function(n, mean, sd, ...) hist(rnorm(n, mean, sd))
##' dlg <- aDialog(items=list(
##'   n = integerItem(10),
##'   mean = numericItem(0),
##'   sd = numericItem(1),
##'   out=graphicDeviceItem()
##' ),
##' buttons="Cancel",
##' model_value_changed=function(.) if(.$is_valid()) do.call("plotIt", .$to_R())
##' )
##' ##
##' ## validation for n, sd
##' n <- dlg$get_item_by_name("n")
##' n$validate <- function(., rawvalue) {
##'   if(rawvalue <= 1) stop("n must be positive integer") else rawvalue
##' }
##' sd <- dlg$get_item_by_name("sd")
##' sd$validate <- function(., rawvalue) {
##'   if(rawvalue <- 0) stop("sd must be positive") else rawvalue
##' }
##' \dontrun{dlg$make_gui()}
##' ##
##' ##
##' ## subtle point about scope. Proto methods can be defined via $<- or [[<- but there is a difference.
##' ## $<- does not have lexical scope whereas [[<- does. The $<- might be more natural to type,
##' ## but [[<- might be more natural to use. In this example,
##' ## The "b" button does not work, as it can't find the
##' ## function a -- the frame of evaluation is the environment dlg (not its enclosing frame).
##' ## Thanks to Gabor for his help with this.
##' scope_example <- function() {
##' a <- function(...) print("hi")
##' dlg <- aDialog(items=list(),
##'                 buttons=c("a","b","c"),
##'                 a_handler=function(.) a(),   ## like [[<-, not $<-
##'                 title="a, c work; b doesn't"
##'                 )
##'  dlg$b_handler <- function(.) a()  ## $<- has unbound variables found in dlg
##'  dlg[['c_handler']] <- a           ## [[<- uses lexical scope for unbound variables
##"  dlg$make_gui()
##' }
##' \dontrun{scope_example()}
##' ## See ?anItemGroup for an example of a modal dialog
  



aDialog <- function(items=list(),
                    title="",
                    help_string="",
                    buttons=c("OK","SPACE","Cancel","Help"),
                    ...
                    ) {
  dlg <- Dialog$proto(items=items,
                      title=title,
                      help_string=help_string,
                      buttons=buttons,
                      ...,
                      funEnvir=FALSE)   # Many thanks to Gabor for
                                        # pointing out this argument
                                        # and to J Hallman for finding
                                        # its necessity.
  dlg$init_model()
  dlg
}



##' Automatically create a dialog for a function
##'
##' Function must have a special markup for its argument. A named
##' argument a..b is interpreted with b determining the type of item
##' to use. We support numeric, string, choice, range, ???  Within the
##' body of the function, the variable a..b should be referred to by
##' a. The idea is that you write and debug the function as usual,
##' then simply modify the argument list to include the types. This
##' function will not work for functions whose arguments use lazy
##' evaluation referring to other argument's values.
##' 
##' All arguments should have a default
##' A choice items should have its default with a vector. The first argument is the selected one
##' A range item is specified with values c(from=., to=..., by=..., [value=from]). If value not give, then from is used.
##' The OK_handler will call f.
##' @param f function to make dialog for. Its arguments must be specified in a certain way.
##' @param title Title for dialog window
##' @param help_string String for help information
##' @param make_gui If \code{TRUE} or \code{add_graphic_device=TRUE} then call dialogs \code{make_gui} method
##' @param add_graphic_device If \code{TRUE} add an graphicDeviceItem to dialog
##' @param ... passed to \code{make_gui} when no graphic device asked for
##' 
##' @return Returns an instance of \code{aDialog}.
##' @export
##' @examples
##' f <- function(x..numeric=1, y..string="a") print(list(x,y))
##' \dontrun{dialogMaker(f)}
##' ## can have missing arguments
##' f <- function(x, y..numeric=1) print(list(x,y))
##' \dontrun{dialogMaker(f)}
##' ## a choice item. Sizing is funny for tables
##' f <- function(x..choice=letters) print(x)
##' \dontrun{dialogMaker(f)}
##' ## range items
##' f <- function(x..numeric=0, mu..numeric=0,
##'               alternative..choice=c("two.sided","less","greater"),
##'               conf.level..range=c(.80,1.00, .01, .95)) {
##'               out <- capture.output(t.test(x, alt=alternative, conf.level=conf.level))
##'               print(out)
##' }
##' \dontrun{dialogMaker(f, title="CI from t.test with summarized values")}
##' 
dialogMaker <- function(f, title="Dialog", help_string="",
                        make_gui=TRUE,
                        add_graphic_device=FALSE,
                        ...
                        ) {

  ## method to get value from a component of the formals list
  ## this is a hack. Perhaps we should use typeof, as B. Ripley suggests on
  ## the mailing list.
  ## Definitely won't handle situations like function(a=1, b=a, ...)
  
  getValFromFormalsComponent <- function(i) UseMethod("getValFromFormalsComponent")
  getValFromFormalsComponent.default <- function(i) i
  getValFromFormalsComponent.call <- function(i) eval(i)
  getValFromFormalsComponent.name <- function(i) get(as.character(i), inherits=TRUE)

  ## help
  if(missing(help_string)) {
    help_string <- paste("Specify the function arguments, then click",
                         "the OK button to call the function.",
                         collapse="\n")
  }
  
  ##
  l <- formals(f)
  itemList <- list()
  for(i in names(l)) {
    if(grepl("\\.{2,2}", i)) {
      val <- unlist(strsplit(i, "\\.\\."))
      value <- eval(getValFromFormalsComponent(l[[i]]))
      itemList[[val[1]]] <- switch(val[2],
                                   "numeric"=numericItem(value, eval_first=TRUE),
                                   "string"=stringItem(value),
                                   "choice"={
                                     if(length(value) > 10)
                                       choiceItem(value[1], value, attr=list(size=c(200,200)))
                                     else
                                       choiceItem(value[1], value)
                                   },
                                   "range"={
                                     if(length(value) == 3)
                                       value[4] <- value[1]
                                     rangeItem(value=value[4],from=value[1], to=value[2], by=value[3])
                                   }
                            )
    } else {
      itemList[[i]] <- expressionItem(as.character(l[[i]]))
    }
  }

  ## now fix formals of f
  nms <- names(l)
  nms <- gsub("\\.\\..*","",nms)
  names(l) <- nms
  formals(f) <- l
  env <- new.env()                      # we define g in an environment
  assign("g", f, envir=env)

  if(add_graphic_device)
    itemList[['graphicDevice']] <- graphicDeviceItem(show_label=FALSE)
  
  dlg <- aDialog(items=itemList,
                 title=title,
                 help_string=help_string,
                 .env=env,
                 OK_handler=function(.) {
                   do.call("g", .$to_R(), envir=.$.env)
                 })

  if(add_graphic_device) {
    ## fix layout
    tmp <- lapply(names(l), function(i) i)
    tmp$horizontal=FALSE
    view <- aGroup(do.call("aGroup", tmp ),
                   aGroup("graphicDevice"))
    dlg$make_gui(gui_layout=view)
  } else if(make_gui) {
    dlg$make_gui(...)
  }

  invisible(dlg)
}
  
