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

##' @include container.R
roxygen()

## Basic editor is a view only. The Item is the model and controller

##' Base Trait for Editor.
##'
##' An editor is a basic view for a widget, essentially a map from gedit, say, to a view. 
##' @export
Editor <- View$proto(class=c("Editor", View$class),
                     ## used in make_ui to make widget, or override make_ui ifd esird
                     .doc_editor_name=paste(
                       desc("Name of gwidget constructor for basic use of make_ui")
                       ),
                     editor_name = "gedit",
                     .doc_editor_style=paste(
                       desc("If non <code>NULL</code>, calls method <code>make_ui_STYLE</code>, if present,",
                            "instead of <code>make_ui</code> method. Example use might be to offer a compact",
                            "style.")
                       ),
                     editor_style=NULL,
                     ## widget name
                     .doc_view_widget_name=paste(
                       desc("Name of primary widget for lookup for setting value from controller")
                       ),
                     view_widget_name="editor",
                     ## attributes are used to specify arguments for constructor
                     attr = list(expand=FALSE),
                     ##
                     ## make_ui sets up widget.
                     .doc_make_ui=paste(
                       desc("Makes widget. Store in widgets property by name",
                            "The container is just a ggroup instance. The container already knows to",
                            "place in a table, notebook etc along with a label. That is handled in the",
                            "Item constructor. This just creates the widget")
                     ),
                     make_ui = function(., container, attr=.$attr, context=., ...) {
                       theSize <- attr$size; attr$size <- NULL
                       attr$container <- container

                       widget <- do.call(.$editor_name, attr)
                       visible(widget) <- TRUE
                       if(!is.null(theSize))
                         size(widget) <- theSize

                       if(!is.null(context$tooltip) && nchar(context$tooltip) > 0)
                         tooltip(widget) <- context$tooltip
                     
                       .$append("widgets", widget, key=.$view_widget_name)

                     },
                     ## this is called by container.R, then it calls make_ui
                     ## once the the object is ready to go
                     .make_ui = function(., cont, attr=.$attr, context=., ...) {
                       ## check if container is traitr object
                       if(is.proto(cont) && exists("traitr", cont)) {
                         if(cont$is("Container"))
                           cont <- cont$container # should be gWidgets object
                         else
                           stop(gettext("The container argument, if a proto object, needs to be a Container"))
                       }

                       ## XXX should add label etc. 
                       ## set up labels, etc. Here context may be item to grab label, ...
                       if(!.$is_realized()) {

                         .$append("widgets", cont, key="parent_container") # useful for enabled, ...
                         
                         

                         ## we need to know if this is a table, or not
                         ## this uses gWidgets hack
                         if(inherits(cont, "gLayout") ||
                            (inherits(cont, "guiContainer") &&grepl("^gLayout", class(cont@widget)))
                            ) {
                           row_no <- get_with_default(tag(cont, "row_no"), 1)
                           col_no <- get_with_default(tag(cont, "col_no"), 1)
                           no_cols <- get_with_default(tag(cont, "no_cols"), 1)

                           
                           if(context$show_label) {
                             ## add label and group for widget
                             cont[row_no, 2*(col_no-1) + 1, anchor=c(1,1)] <-
                               (label <- glabel(context$label, cont = cont))
                             cont[row_no, 2*(col_no-1) + 2, anchor=c(-1,1), expand=TRUE] <-
                               (g <- ggroup(cont=cont, horizontal=FALSE, expand=TRUE, anchor=c(-1,1)))
                           } else {
                             ## add group object to hold widget, spread over 2 columns
                             cont[row_no, 2*(col_no-1) + 1:2, anchor=c(-1,1), expand=TRUE] <-
                               (g <- ggroup(cont=cont, horizontal=FALSE, expand=TRUE, fill="x"))
                           }
                           tag(cont, "row_no") <- row_no + (col_no == no_cols)
                           tag(cont, "col_no") <- (col_no %% no_cols) + 1
                         } else if(inherits(cont, "gNotebook") ||
                                   (inherits(cont, "guiContainer") && grepl("^gNotebook", class(cont@widget)))
                                   ) { 
                           attr$label <- get_with_default(context$label, "")
                           g <- ggroup(cont=cont, expand=TRUE, horizontal=FALSE)
                         } else {
                           cont <- ggroup(cont=cont, horizontal=TRUE, expand=TRUE)
                           if(context$show_label)
                             glabel(context$label, cont = cont, anchor=c(1,1))
                           g <- ggroup(cont=cont, horizontal=FALSE, expand=TRUE)
                         }
                         ## now make widget
                         ## check for editor_style
                         if(context$has_slot("editor_style") &&
                            !is.null(editor_style <- context$editor_style) &&
                            .$has_slot(paste("make_ui_",editor_style, sep=""))) {
                           .$do_call(paste("make_ui_",editor_style, sep=""), list(cont=g, attr=attr, context=context))
                         } else {
                           .$make_ui(cont=g, attr=attr, context=context, ...)
                         }
                       }
                     },
                     visible = function(., value) {
                       sapply(.$get_widgets(), function(i) visible(i) <- value)
                     },
                     enabled = function(., value) {
                       sapply(.$get_widgets(), function(i) enabled(i) <- value)
                     },
                     ## These are for the view, not the model -- they shouldn't be used
                     ## by end users.
                     .doc_get_value_from_view=paste(
                       desc("Private method for getting value from view (editor) for controller")
                       ),
                     get_value_from_view = function(.) {
                       nm <- if(!is.null(.$view_widget_name))
                         .$view_widget_name
                       else
                         1
                       index <- get_with_default(.$get_slot("by_index"), FALSE)
                       if(isExtant(.$widgets[[nm]]))
                          val <- svalue(.$widgets[[nm]], index=index)
                       else
                         val <- NULL
                     },
                     ## set value in view by name. If name is NULL, gets first
                     .doc_set_value_in_view=paste(
                       desc("Private method for setting value from view (editor) for controller",
                            "When assigning through svalue, should add blockHandler/unblockHandler calls",
                            "to prevent model from being updated, which can cause loops"),
                       param("widget_name"," name of widget, default for simple case"),
                       param("value"," value to be set into editor")
                       ),

                     set_value_in_view = function(.,widget_name=.$view_widget_name, value) {
                       if(.$is_realized()) {
                         widget <- .$get_widget_by_name(widget_name)
                         cur_val <- svalue(widget)
                         index <- get_with_default(.$get_slot("by_index"), FALSE)
                         
                         if(digest(cur_val) != digest(value)) {
                           blockHandler(widget)
                           try(svalue(widget, index=index) <- value, silent=TRUE)
                           unblockHandler(widget)
                         }
                       }
                     },
                     ## instance
                     .doc_instance=paste(
                       desc("Create a new instance of the editor. Does not share widgets")
                       ),
                     instance=function(.) {
                       obj <- .$proto()
                       obj$widgets <- list()
                       obj
                     },
                     .doc_set_valid=paste(
                       desc("Method to set editor to indicate valid value")
                       ),
                     set_valid=function(.) {},
                     .doc_set_invalid=paste(
                       desc("Method to indicate invalid value")
                       ),
                     set_invalid=function(., mesg) {
                       if(!missing(mesg))
                         galert(mesg)
                     }
                     )

## Various editor traits for different case
## see also special files (eg. itemlist) for others



##' A Base Trait for an editor using the entry widget
##'
##' @rdname Editor
EntryEditor <- Editor$proto(class=c("EntryEditor", Editor$class),
                            .doc_format_fun=paste(
                              desc("A function to call to coerce the value before displaying in entry",
                                   "widget. Set to <code>NULL</code> for no coercion."),
                              param("value"," value to coerce")
                              ),
                            format_fun = function(., value) sprintf("%s",value),
                            make_ui = function(., container, attr=.$attr, context=., ...) {
                              if(!is.null(.$format_fun))
                                attr$text <- .$format_fun(context$value)
                              else
                                attr$text <- context$value
                              .$next_method("make_ui")(., container, attr, context, ...)

                            },
                            ## draw widget in valid state
                            set_valid=function(.) {
                              widget <- .$get_widget_by_name(.$view_widget_name)
                              if(!is.proto(widget)) {
                                e <- getToolkitWidget(widget)
                                if(inherits(e, "RGtkObject")) {
                                  require(RGtk2)
                                  ## set backgroup
                                  e$modifyBg(GtkStateType['normal'],NULL)
                                }
                              }
                            },
                            .doc_set_invalid=paste(
                              desc("Method call to set widget in an invalid state. The item should do",
                              "validation during the call that takes values from the view to the model."),
                              param("mesg","Optional message to display if value for entry is invalid")
                              ),
                            set_invalid=function(., mesg) {
                              widget <- .$get_widget_by_name(.$view_widget_name)
                              if(!is.proto(widget)) {
                                e <- getToolkitWidget(widget)
                                if(inherits(e, "RGtkObject")) {
                                  ## set background
                                  require(RGtk2)
                                  e$modifyBg(GtkStateType['normal'],"red")
                                } else {
                                  .$next_method("set_invalid")(., mesg)
                                }
                              }
                            }
                            )
                            

##' Base trait for editor where there are underlying values to choose from
##'
##' editor has regular or compact style
##' @rdname Editor
ObjectWithValuesEditor <- Editor$proto(class=c("ObjectWithValuesEditor", Editor$class),
                                       editor_name="gradio",
                                       by_index=FALSE, # select, set by index
                                       ## need to put in index=FALSE here
                                       set_value_in_view = function(.,widget_name=.$view_widget_name, value) {
                                         if(.$is_realized()) {
                                           widget <- .$get_widget_by_name(widget_name)
                                           if(length(value) == 0) { # clear out
                                             blockHandler(widget)
                                             svalue(widget, index=.$by_index) <- value
                                             unblockHandler(widget)
                                           } else if(is.null(value) || is.na(value) ||
                                              (is.character(value) && length(value) == 0)
                                              ) {
                                             return()
                                           }
                                           curVal <- svalue(widget, .$by_index)
                                           if(digest(curVal) != digest(value) && length(value)) {
                                             blockHandler(widget)
                                             svalue(widget, index=.$by_index) <- value
                                             unblockHandler(widget)
                                           }
                                         }
                                         invisible()
                                       },

                                       set_values_in_view = function(., values) {
                                         widget <- .$get_widget_by_name(.$view_widget_name)
                                         if(!is.null(widget) && isExtant(widget)) {
                                           curVal <- svalue(widget, index=.$by_index)
                                           widget[] <- values
                                           if(.$by_index) {
                                             blockHandler(widget)
                                             svalue(widget, index=.$by_index) <- curVal
                                             unblockHandler(widget)
                                           } else {
                                             if(length(values) > 0 &&
                                                length(curVal) > 0 && nchar(curVal) > 0 &&
                                                curVal %in% values) {
                                               blockHandler(widget)
                                               try(svalue(widget, index=.$by_index) <- curVal, silent=TRUE)
                                               unblockHandler(widget)
                                             }
                                           }
                                         }
                                         invisible()
                                       },
                                       set_invalid=function(., ...) {},
                                       make_ui=function(., container, attr=.$attr, context=., ...) {
                                         attr$items <- context$get_values()
                                         attr$selected <- 0
                                         
                                         if(length(context$get_value())) {
                                           value <- context$get_value()
                                           values <- context$get_values()
                                           if(.$by_index) {
                                             selected <- value
                                           } else if(.$has_slot("multiple") && .$multiple) {
                                             which(values %in% value)
                                           } else {
                                             if(value %in% values) {
                                               selected <- min(which(value %in% values))
                                             } else {
                                               selected <- 1
                                             }
                                             attr$selected <- selected
                                           }
                                         }
                                         .$next_method("make_ui")(., container, attr, context, ...)
                                       },
                                       ## compact style is a label and button to pop up choice
                                       make_ui_compact=function(., container, attr=.$attr, context=., ...) {
                                         .$set_values_in_view <- function(...) {} # override
                                         cont <- ggroup(cont=container, horizontal=TRUE, expand=TRUE)
                                         
                                         attr$items <- context$get_values()
                                         attr$selected = if(context$get_value() %in% context$get_values())
                                           min(which(context$get_value() == context$get_values()))
                                         else
                                           1
                                         .$widgets[[.$view_widget_name]] <-
                                           glabel(context$value, cont=cont, expand=TRUE)
                                         .$widgets[["button"]] <- (b <- gbutton("Select value", cont=cont))
                                         addHandlerClicked(b, function(h,...) {
                                           obj <- h$action$obj
                                           context <- h$action$context
                                           widget <- obj$widgets[[obj$view_widget_name]]
                                           msg <- gettext("Select an item by double-clicking")
                                           dlg <- aDialog(items=list(item=choiceItem(context$get_value(),
                                                                       by_index=FALSE,
                                                                       values=context$get_values(),
                                                                       editor_type="gtable",
                                                                       tooltip=msg,
                                                                       show_label=FALSE,
                                                                       attr=list(size=c(300,300)))),
                                                          buttons=c("Cancel"),
                                                          title=msg
                                                          )
                                           dlg$label_widget <- widget
                                           dlg$label_context <- context
                                           dlg$make_gui(parent=container)
                                           dlg$property_item_value_changed <- function(., value, old_value) {
                                             svalue(.$label_widget) <- value
                                             context <- .$label_context
                                             context$do_call(sprintf("set_%s", context$name),list(value=value))
                                             .$close_gui()
                                           }

                                         },action=list(obj=., context=context))
                                       }
                                       )


##' Trait for Editor for TRUE/FALSE selection
##'
##' Editor has regular or compact style
##' @rdname BooleanEditor
BooleanEditor <- Editor$proto(class=c("BooleanEditor", Editor$class),
                              editor_name="gcombobox",
                              make_ui=function(., container, attr=.$attr, context, ...) {
                                attr$items <- c(TRUE, FALSE)
                                .$next_method("make_ui")(., container, attr, context, ...)
                              },
                              make_ui_compact=function(., container, attr=.$attr, context, ...) {
                                widget <- gcheckbox("", cont=container, ...) # possible use.togglebutton=TRUE
                                .$append("widgets", widget, key=.$view_widget_name)
                              }
                              )

##' Trait for making a range editor (slider, spinbox)
##'
##' @rdname Editor
RangeEditor <- Editor$proto(class=c("RangeEditor", Editor$class),
                            view_widget_name="slider",
                            ## makes combo UI. Must also adjust the
                            ## set_value_in_view method
                            make_ui = function(., container, attr=.$attr, context=., ...) {
                              g <- ggroup(cont=container, horizontal=TRUE)
                              l <- list(from=context$from,
                                        to=context$to,
                                        by=context$by,
                                        value=context$value,
                                        cont=g,
                                        expand=TRUE)
                              slider <- do.call("gslider",l)
                              visible(slider) <- TRUE
                               .$append("widgets", slider, key=.$view_widget_name)
                              by <- context$by
                              if(as.integer(by) == by) {
                                l$expand <- FALSE
                                spinner <- do.call(gspinbutton,l)
                                visible(spinner) <- TRUE
                                .$append("widgets", spinner, key="spinner")
                              }
                                                
                              if(!is.null(attr$size))
                                size(g) <- attr$size
                              
                              if(!is.null(context$tooltip) && nchar(context$tooltip) > 0)
                                tooltip(slider) <- context$tooltip
                     
                            },
                            ## need to update both slider and spinner
                            set_value_in_view = function(.,widget_name=.$view_widget_name, value) {
                              if(.$is_realized()) {
                                if(is.null(value))
                                  return()
                                widget <- .$widgets[["slider"]]
                                spinner <- .$widgets[["spinner"]] # may be NULL
                                if(digest(svalue(widget)) != digest(value)) {
                                  blockHandler(widget)
                                  try(svalue(widget) <- value, silent=TRUE)
                                  unblockHandler(widget)
                                }
                                if(!is.null(spinner) && (digest(svalue(spinner)) != digest(value))) {
                                  blockHandler(widget)
                                  try(svalue(spinner) <- value, silent=TRUE)
                                  unblockHandler(widget)
                                }
                              }
                            }
                              )
                            

##' Trait for button editor
##'
##' @rdname Editor
ButtonEditor <- Editor$proto(class=c("ButtonEditor", Editor$class),
                             editor_name="gbutton",
                             attr=list(),
                             make_ui=function(., container, attr=.$attr, context, ...) {
                               attr$text=context$value
                               .$next_method("make_ui")(., container, attr, context, ...)
                             }
                             )

##' Trait for embedding an image file
##'
##' @rdname Editor
ImageEditor <- Editor$proto(class=c("ImageEditor", Editor$class),
                            editor_name="gimage")


##' Trait for embedding graphics (Qt, RGtk2 only)
##'
##'
##' @rdname Editor
GraphEditor <- Editor$proto(class=c("GraphEditor", Editor$class),
                            attr=list(),
                            editor_name="ggraphics",
                            set_value_in_view = function(.,...) {},
                            get_value_from_view= function(.,...) {}
                            )

##' Trait for making File browser editor
##'
##' @rdname Editor
FileBrowseEditor <- Editor$proto(class=c("FileBrowseEditor", Editor$class),
                                 editor_name="gfilebrowse")

##' Trait for data selection editor
##'
##' @rdname Editor
DateEditor <- Editor$proto(class=c("DateEditor", Editor$class),
                                 editor_name="gcalendar")


##' Trait for displaying a table of information 
##'
##' No selection, just display. For selection use ChoiceItem
##'
##' @rdname Editor
TableEditor <- Editor$proto(class=c("TableEditor", Editor$class),
                            editor_name="gtable",
                            attr=list(expand=TRUE),
                            set_value_in_view=function(.,widget_name=.$view_widget_name, value) {
                              ## replace the widget -- so that size issues are there
                              if(.$is_realized()) {
                                widget <- .$get_widget_by_name(widget_name)
                                group <- .$get_widget_by_name(".cont")
                                if(!is.null(group) && !is.null(widget)) {
                                  delete(group, widget)
                                  if(!(is.data.frame(value) || is.matrix(value)))
                                    value <- as.data.frame(value)
                                  .$widgets[[widget_name]] <- gtable(value, cont=group, expand=TRUE)
                                }
                              }
                            },
                            get_value_from_view=function(., ...) {},
                            get_widget=function(., ...) {
                              .$widgets[[.$view_widget_name]]
                            },
                            make_ui=function(., container, attr=.$attr, context=., ...) {
                              df <- context$do_call(sprintf("get_%s",context$name))
                              if(!is.data.frame(df))
                                df <- data.frame(V1="", V2="")
                              .$widgets[[".cont"]] <- (g <- ggroup(cont=container, expand=TRUE))
                              attr$items <- df
                              attr$container <- g
                              attr$expand=TRUE
                              .$widgets[[.$view_widget_name]] <- (widget <- do.call("gtable",attr))
                              visible(widget) <- TRUE
                              if(!is.null(attr$size))
                                size(g) <- attr$size
                            }
                            )

                              
##' Trait for a label
##' 
##' @rdname Editor
LabelEditor <- Editor$proto(class=c("LabelEditor",Editor$class),
                            editor_name="glabel"
                            )


##' Trait for making a visual separator
##'
##' @rdname Editor
SeparatorEditor <- Editor$proto(class=c("SeparatorEditor",Editor$class),
                                editor_name="gseparator",
                                attr=list(expand=TRUE),
                                set_value_in_view = function(.,...) {},
                                get_value_from_view= function(.,...) {}
                                )
