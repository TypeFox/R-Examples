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

#' @include view.R
roxygen()

##################################################
#' Trait for Controller objects
#'
#' A controller connects a model and an associated view to synchronize changes in one with another
#' This implementation rests on the controller having some suitably named methods
Controller <- BaseTrait$proto(class=c("Controller", BaseTrait$class),
                               ## list of adapters. Each adapter specified with a list.
                               ## list(property="propname",
                               ##      view_name="viewname",
                               ##      add_handler_name=c("addHandlerChanged"), # or NULL to suppress
                               ##      handler_user_data=NULL
                               ##      )
                               adapters=list(),  
                               ## Methods
                               ## make a new controller child must override some methods

                               ## get/set model
                              .doc_model=paste(
                                desc("Model property. See <code>set/get_model</code>")
                                ),
                               model=NULL,
                              .doc_get_model=paste(
                                desc("get the model associated to the controller")
                                ),
                               get_model=function(.) .$model,
                              .doc_set_model=paste(
                                desc("set the model associated to the controller. Moves the observer"),
                                param("model","Model instance")
                                ),
                               set_model = function(., model) {
                                 if(is.proto(model) && model$is("Model")) {
                                   .$model$remove_observer(.)
                                   .$model <- model
                                   .$model$add_observer(.)
                                   sapply(.$.adapters, function(i) i$set_model(model))
                                 }
                               },
                               ## get/set view controller object
                               view=NULL,
                              .doc_get_view=paste(
                                desc("Return view associated with the controller")
                                ),
                               get_view=function(.) .$view,
                               ## set view, remove old if present
                              .doc_set_view=paste(
                                desc("Set view associated with the controller"),
                                param("view","View instance")
                                ),
                               set_view = function(., view) {
                                 if(is.proto(view) && view$is("View")) {
                                   if(is.proto(.$get_view()) && .$get_view()$is("View"))
                                     .$remove_view()
                                   .$view <- view
                                   sapply(.$.adapters, function(i) i$set_view(view))
                                 }
                               },
                               ## function call to remove view -- disconnect handlers
                               remove_view=function(.) {
                                 sapply(.$.adapters, function(i) i$remove_view())
                               },
                                
                               ## sets function to listen for changes to view
                               ## might just define .$property_NAME_value_changed method
                              .doc_update_from_model=paste(
                                desc("Method call to set up model -> view update")
                                ),
                               update_from_model=function(.) {},

                               ## transfers changes to model back from view
                              .doc_update_from_view=paste(
                                desc("Method call to update model from view.")
                                ),
                               update_from_view = function(.) {},

                               ## function to register adapters if .$adapters non empty
                              .doc_register_adapters=paste(
                                desc("Register adapters, if any")
                                ),
                               register_adapters = function(.) {
                                 if(length(.$adapters) && !length(.$.adapters)) {
                                   .$.adapters <- lapply(.$adapters, function(i) {
                                     Adapter$proto(model=.$get_model(),
                                                   view=.$get_view(),
                                                   property=i$property,
                                                   view_widget_name=i$view_widget_name,
                                                   add_handler_name=i$add_handler_name,
                                                   handler_user_data=i$handler_user_data,
                                                   by_index=.$by_index
                                                   )
                                   })
                                 }
                                 if(length(.$.adapters))
                                   sapply(.$.adapters, function(i) i$init())
                               },
                               ## intialize controller
                               ## updates model, view and any adapters
                              .doc_init =paste(
                                desc("Initialize controller. Calls update_from_model,",
                                     "update_from_view, register_adapters")
                                ),
                               init = function(.) {
                                 if(!is.null(.$get_model())) {
                                   .$update_from_model()
                                   .$get_model()$init()
                                   .$get_model()$add_observer(.)
                                 }
                                 if(!is.null(.$get_view())) {
                                   .$update_from_view()
                                 }
                                 .$register_adapters()
                                 ## call value_changed methods to update any views
                                 nms <- .$list_methods()
                                 sapply(nms[grep("property_(.*)_value_changed$", nms)],
                                        function(i) {
                                          prop <-  gsub("property_(.*)_value_changed$","\\1",i)
##                                          XXX changed so that value is coerced -- didn't work
                                          meth_name <- sprintf("get_%s",prop)
                                          if(.$has_slot(meth_name))
                                            value <- do.call(meth_name,list(.))
                                          else
                                            value <- .$get_model()$getattr(prop)
                                          .$get_slot(i)(., value, NA)
##                                          .$get_slot(i)(., .$get_model()$getattr(prop),NA)
                                        })
                                 invisible()
                               },
                              ## when a model changes, this function is called
                              ## if the controller is an observer
                              .doc_model_value_changed=paste(
                                desc("Method called when any item in the model is changed")
                                ),
                              ## model_value_changed= function(.) {},

                              ## When a property of a model is changed through set_PROPERTYNAME
                              ## then methods of this type are called
                              .doc_property_PROPERTYNAME_value_changed=paste(
                                desc("Method called when property <code>PROPERTYNAME</code>is changed"),
                                param("value","New value (raw value, not coerced)"),
                                param("old_value", "Old value of property")
                                ),
                              ## property_PROPERTYNAME_value_changed = function(., value, old_value,...) {},

                               ## Private properties
                               .adapters=list(), # actual instances
                               .handlerIDs=list()
                               ## private methods
                               )

#' Trait for Adapter object
#' An adapter is a simple controller connecting one model property with a widget in a view
#' by default the adapter synchronizes changes
#' @export
Adapter <- Controller$proto(class=c("Adapter", Controller$class),
                            property=NULL, # property in model
                            # name of widget in view, otherwise last one
                            view_widget_name=NULL,
                            add_handler_name=c("addHandlerChanged"), ## "" for no view->model
                            handler_user_data=NULL,               # passed through to handler call
                            .handlerIDs=c(),                       # store IDs of handlers
                            ## Methods. init does the work
                            update_from_model = function(.) {
                              ## set up model to notify view For example:
                              view <- .$get_view()
                              meth_name<- sprintf("property_%s_value_changed", .$property)
                              .$set_slot(meth_name,
                                         function(., value, old_value) {
                                           view$set_value_in_view(.$view_widget_name, value)
                                         })
                                     
                              ## call method
                              .$get_slot(meth_name)(., .$get_model()$getattr(.$property), NA)
                            },
                            ## this propogates changes from the view back to the model
                            update_from_view = function(.) {
                              ## here view knows about model through controller (this adapter)
                              if(!.$get_view()$is_realized()) return()
                              if(!is.null(.$view_widget_name))
                                widget <- .$get_view()$get_widget_by_name(.$view_widget_name)
                              else
                                widget <- tail(.$get_view()$get_widgets(), n=1)[[1]]
                              ## gWidgets specific call to set up control between model and
                              ## view
                              if(is.null(.$add_handler_name))
                                .$add_handler_name="addHandlerChanged"
                              sapply(.$add_handler_name, function(i) {
                                if(i != "") {
                                  lst <- list(obj=widget,
                                              handler=function(h,...) {
                                                . <- h$action$adapter
                                                
                                                ## set property in model using name
                                                index <- get_with_default(.$get_slot("by_index"), FALSE)
                                                value <- svalue(h$obj, index=index)

                                                ## XXX added check here, not sure why it is needed. 
                                                if(isExtant(h$obj) && !is.null(value)) {
                                                  .$model$setattr(.$property, value)
                                                }
                                              },
                                              action=list(adapter=.))
                                  .$append(".handlerIDs", do.call(i, lst))
                              }
                              })
                            },
                            ## remove the view -- say be removeHandler call.
                            remove_view=function(.) {
                              if(.$has_slot(".handlerIDs"))
                                sapply(.$.handlerIDs, function() removeHandler(.$get_view(), i))
                            },
                            init=function(.) {
                              ## check that we are all there
                              if(!is.null(.$property) &&
                                 (is.proto(model <- .$get_model()) && model$is("Model")) &&
                                 (is.proto(view <- .$get_view()) && view$is("View"))) {
                                .$update_from_model()
                                .$update_from_view()
                              } else {
                                warning("Adapter does not have view, model and property")
                              }
                              .$model$add_observer(.)
                            }
                            )

## constructor
#' Constructor for a Controller proto objects
#'
#' Simply provides a more typical calling interface for the Controller proto object
#' @param ... passed to proto method for
#' @return returns the Controller object
#' @export
aController <- function(...) {
  obj <- aController$new(...)
  ## adjust class
  obj
}

