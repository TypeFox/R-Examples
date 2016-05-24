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

#' @include base.R
roxygen()

## Make a MVC framework using proto


#' Trait for a model object.
#'
#' Model objects consist of properties and the methods that manipulate them
#' Models are initialized by \code{init} so that getter and setter pairs are made
#' When setter functions are called, the model's observers are notified
#' @export
Model <- BaseTrait$proto(class=c("Model", BaseTrait$class),
                          ## properties are added to children of this
                          ## eg value="some_value",
                          ##### Methods ###
                          ## call observers
                          ## if key is passed in we notify that property was changed -- if it was
                          ## we notify controllers that model_value_changed happened
                         .doc_notify_observers=paste(
                           desc("This function is used to notify observers. It calls the observer methods",
                                "<code>property_PROPERTYNAME_value_changed</code> and <code>model_value_changed</code>",
                                "if they exist in the observer.",
                                "For private use, the method <code>.property_PROPERTY_value_cahnged</code>",
                                "is also called."),
                           param("key", "if non-null is used to signal propery method"),
                           param("value","is new value"),
                           param("old_value","is old value, useful if an undo mechanism is desired."),
                           param("notify_private","There for undo mechanism"),
                           returns("No return value")
                           ),
                          notify_observers = function(., key=NULL, value=NA, old_value=NA, notify_private=TRUE) {
                            sapply(.$list_observers(), function(i) {
                              if(digest(value) != digest(old_value)) { #serialize, then compare
                                if(!is.null(key)) {
                                  ## call private ones if...
                                  if(notify_private)
                                    i$do_call(sprintf(".property_%s_value_changed",key),
                                              list(value=value, old_value=old_value))
                                  ## call public ones
                                  i$do_call(sprintf("property_%s_value_changed",key),
                                            list(value=value, old_value=old_value))
                                }
                                ## always call these for new values
                                i$do_call("model_value_changed", list())
                              }
                            })
                            invisible()
                          },
                          ## add an observer -- a Controller instance or a child such as an adapter
                          .doc_add_observer=paste(
                            desc("This method adds an observer to the model."),
                            param("observer","A Controller instance or Model instance is the observer"),
                            returns("no return value")
                            ),
                          add_observer = function(., observer) {
                            if(is.proto(observer) && observer$is(c("Controller","Model"))) {
                              ## don't double up.
                              ind <- sapply(.$.observers, function(i) i$identical(observer))
                              if(!any(ind)) {
                                id <- length(.$.observers) + 1
                                .$.observers[[id]] <- observer
                              }
                            }
                          },
                          ## remove observer
                          .doc_remove_observer=paste(
                            desc("This method removes an observer from the model."),
                            param("observer"," observer to remove")
                            ),
                          remove_observer = function(., observer) {
                            if(!missing(observer) && (is.proto(observer) && observer$is(c("Controller","Model")))) {
                              ind <- sapply(.$.observers, function(i) i$identical(observer))
                              if(any(ind))
                                sapply(which(ind), function(i) .$.observers[[i]] <- NULL)
                            }
                          },
                          ## list observers
                         .doc_list_observers=paste(
                           desc("Return all observers of the model as a list")
                           ),
                         list_observers=function(.) .$.observers,
                          ## observables have getter/setter methods defined for them
                          ## on the fly when init is called on Modal instance
                          ## these refer to properties of the model
                          ## methods
                          .doc_getattr=paste(
                            desc("Method to get an attribute (property) value from the object. That is,",
                                 "<code>.$getattr('value')</code> is same as <code>.$value</code>"),
                            param("key","Name of property")
                            ),
                          getattr = function(., key) get(key, envir=.),
                          .doc_setattr=paste(
                            desc("Method to set a value for an attribute (property) in the object.",
                                 "Unlike assignment with <code>.$value = key </code>, this calls the",
                                 "<code>notify_observers</code>method to notify observers of a property",
                                 "value being changed."),
                            param("key","Name of property"),
                            param("value","New value of property")
                            ),
                          setattr = function(., key, value, notify_private=TRUE) {
                            old_value <- .$getattr(key)
                            assign(key, value, envir=.)
                            if(digest(value) != digest(old_value)) {#serialize, then compare
                              .$notify_observers(key=key, value=value, old_value=old_value, notify_private)
                            }
                          },
                          ## init makes getter/setter variables for all properties in the model
                          .doc_init=paste(
                            desc("Initialization method for models. Creates get_PROPERTY and set_PROPERTY methods",
                            "to access model values. These simply pass to getattr and setattr methods")
                            ),
                          init = function(.) {
                            sapply(.$list_properties(), function(i) {
                              .$assign_if_null(sprintf("get_%s",i),
                                               function(.,...) .$getattr(i)
                                               )
                              .$assign_if_null(sprintf("set_%s",i),
                                               function(., value, ...) {
                                                 .$setattr(i,value)
                                               })
                            })
                            invisible()
                          },

                          ## private
                          ## properties
                          .observers = list()
                          )

#' Constructor for a Model proto objects
#'
#' Simply provides a more typical calling interface for the Model proto object
#' @param ... passed to proto method for
#' @return returns the Model object. Call \code{obj$show_help()} to view its methods and properties.
#' @export
aModel <- function(...) Model$new(...)
