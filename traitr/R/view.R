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

#' @include model.R
roxygen()

## Use makeUI to makeUI components, updateUI to update them

##################################################
#' Trait for View objects.
#'
#' A view "displays" the values in an associated model. The association is through a controller
#' the view does not know the controller or the model
#' Views are initialized through the \code{make_ui} method.
#' @export
View <- BaseTrait$proto(class=c("View", BaseTrait$class),
                        ## properties
                        ## list of widgets in view for access by controller by name
                        ## layout
                        .doc_attr=paste(
                          desc("Property list of attributes for constructor")
                          ),
                        attr=list(),
                        ## widgets
                        .doc_widgets=paste(
                          desc("Property list of widgets in view, see <code>get_widget_by_name</code>")
                          ),
                        widgets=list(),
                        ## methods
                        ## get all widgets
                        get_widgets = function(.) .$widgets,
                        ## get widget  by name
                        .doc_get_widget_by_name=paste(
                          desc("Look up widget by name")
                          ),
                        get_widget_by_name= function(., key) .$get_widgets()[[key]],
                        
                        ## makeUI sets up widget.
                        ## store widgets into .$widgets by name
                        .doc_make_ui=paste(
                          desc("Initialize user interface provided by the view"),
                          param("cont","gWidgets container"),
                          param("attr","attributes passed constructor")
                          ),
                        make_ui= function(., cont, attr=.$attr) {},
                         ## is_realized is TRUE if widget is realized
                        .doc_is_realized=paste(
                          desc("Is view realized and still extant?")
                          ),
                         is_realized=function(.) {
                           length(.$get_widgets()) && isExtant(.$get_widgets()[[1]])
                         },
                         ## communication between view and model done through controller
                         ## these are generic functions for a view to get and set the main value
                         ## others are specific to sub instances
                         ## These should not be called by users, they do not sychronize with model, but instead
                         ## the value refers to the view. Basic call is svalue, or svalue<-
                        .doc_get_value_from_view=paste(
                          desc("For simple views, used to transport values from view to model by the controller")
                          ),
                        get_value_from_view = function(.) {},
                        .doc_set_value_in_view=paste(
                          desc("Method to set a value in the view. Called by controller"),
                          param("widget_name","Name of widget to look up"),
                          param("value","Value to set")
                          ),
                        set_value_in_view = function(., widget_name, value) {
                          if(.$is_realized()) {
                            widget <- .$get_widget_by_name(widget_name)
                            index <- get_with_default(.$get_slot("by_index"), FALSE)
                            blockHandler(widget)
                            try(svalue(widget, index=index) <- value, silent=TRUE)
                            unblockHandler(widget)
                          }
                        },
                         ## methods to change visibility, enabled of the widgets
                         ## gWidgets specific
                        .doc_enabled=paste(
                          desc("Is view enabled"),
                          param("value","logical value")
                          ),
                         enabled = function(., value) {
                           invisible(sapply(.$get_widgets(), function(i) enabled(i) <- as.logical(value)))
                         },
                        .doc_visible=paste(
                          desc("Set if view is visible"),
                          param("value","logical value")
                          ),
                         visible = function(., value) {
                           invisible(sapply(.$get_widgets(), function(i) visible(i) <- as.logical(value)))
                         }
                         )

## constructor
#' Constructor for a View proto object
#'
#' Simply provides a more typical calling interface for the View proto object
#' @param ... passed to proto method for
#' @return Returns the View object. Call \code{obj$show_help()} to view its methods and properties.
#' @export
aView <- function(...) View$new(...)
