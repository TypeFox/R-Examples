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

##' @include editor.R
roxygen()

##################################################
## An item group is the main thing. A Model with a means to be viewed and
## adapters defined automatically

##' Base Trait to group items together to form a model. ItemGroups may be viewed as a model, view and controller
##' bundled together in a tidy package.
##'
##' An item group is a collection of Item instances. These
##' are specified through the \code{items} property as a list.
##'
##' ItemGroups implement the observer pattern, so are models and can have observers listen for changes
##' ItemGroups observe themselves to update the user interface on model changes. One can add other observers if desired.
##' See the \code{init} method for an example.
##'
##' If a method \code{model_value_changed} is defined, then it will be called with the ItemGroup instance being
##' in whenever a property of the model has a new value. The handlers \code{property_NAME_value_changed} is called
##' when the main value in item NAME is changed. (An item can have several properties, one of which is the main one.)
##'
##' ItemGroups have a \code{make_gui} method to make a view of the model. The layout of this GUI can be
##' specified through its \code{gui_layout} argument, or by default have a simple table layout. ItemGroup instances
##' are meant to be embedded into a GUI, so the \code{cont} argument is needed to pass in the desired container.
##'
##' @export
ItemGroup <- Model$proto(class=c("ItemGroup",  Model$class),
                         ## our properties are instances of Item
                         .doc_items=paste(
                           desc("A list containing the items in ItemGroup")
                           ),
                         items=list(),
                         ## methods
                         ## Get items
                         ## items are Item instances *and* ItemGroup instances
                         .doc_get_items=paste(
                           desc("Return named list of items")
                           ),
                         get_items = function(.) { # all items
                           items <- .$items
                           nms <- sapply(items, function(i) i$name)
                           names(items) <- make.names(nms, unique=TRUE)
                           items
                         },
                         ## gets items, flattens out ItemGroup entries
                         ## does not return some items
                         .doc_get_items_only=paste(
                           desc("Method to return items that can have values (ignoring SeparatorItems, say)")
                           ),
                         get_items_only = function(.) {
                           not_these <- c("SeparatorItem", "GraphItem")
                           l <- list()
                           for(i in .$get_items()) {
                             if(i$is("Item"))
                               l <- c(l, i)
                             else
                               l <- c(l, i$get_items_only())
                           }
                           ind <- sapply(l, function(i) i$is(not_these))
                           l <- l[!ind]
                           l[!sapply(l, is.null)]
                         },
                         ## return items by key or name, flattened ober ItemGroups
                         .doc_get_item=paste(
                           desc("Method to retrieve item by a key")
                           ),
                         get_item = function(., key) .$items[[key]], # item by key
                         .doc_get_item_by_name = paste(
                           desc("Method to retrieve item by its name")
                           ),
                         get_item_by_name = function(., name) {      # item by name
                           .$get_items()[[name]]
                         },
                         .doc_icon=paste(
                           desc("An icon for identification purposes. Use an empty string for no icon.")
                           ),
                         icon="",    # XXX modify
                         ## to get values from the items
                         .doc_to_R=paste(
                           desc("Method to take values in model and return as a named list. Names",
                                "come from the item names. Each item returns its own value, likely after",
                                "coercion to the desired type."),
                           param("drop","If <code>TRUE returns values as a flattened list. Otherwise,",
                                 "nesting is provided by itemgroup objects")
                           ),
                         to_R = function(., drop=TRUE) {
                           if(drop)
                             sapply(.$get_items_only(), function(i) i$to_R())
                           else
                             lapply(.$get_items(), function(i) i$to_R())
                         },
                         .doc_to_string=paste(
                           desc("Method to return list of string representations of items"),
                           param("drop","If <code>TRUE returns values as a flattened list. Otherwise,",
                                 "nesting is provided by itemgroup objects")
                           ),
                         to_string = function(., drop=TRUE) {
                           if(drop)
                             sapply(.$get_items_only(), function(i) i$to_string())
                           else
                             lapply(.$get_items(), function(i) i$to_string())
                         },
                         ## are all values valid
                          .doc_is_valid=paste(
                            desc("Are all items in the itemgroup valid?")
                            ),
                          is_valid=function(.) {
                            all(sapply(.$get_items_only(), function(i) i$is_valid()))
                          },
                         ## make_gui constructs the GUI from a specifed layout
                         ## parent, visible ignored
                         gui_layout=NULL, # empty by default
                         .doc_make_gui=paste(
                           desc("Method to make GUI for ItemGroup. The GUI layout is specified by a container",
                                "defaulting to a simple table listing all items (aContainer)."),
                           param("cont","container to hold the GUI. See Dialog for self contained GUIs"),
                           param("gui_layout", "A layout for the GUI. The default is a table layout"),
                           param("parent", "Not used. Here for consistency with same method in Dialog class"),
                           param("visible", "Not implemented"),
                           returns("Creates a GUI, no return value")
                           ),
                         make_gui = function(., container, gui_layout, parent=NULL, visible=TRUE, ...) {
                           .$init_model() # initialize model
                           
                           if(missing(gui_layout)) {
                             if(.$has_slot("gui_layout") && is.proto(.$gui_layout) && .$gui_layout$is("Container"))
                               gui_layout <- .$gui_layout
                             else
                               gui_layout <- .$make_default_gui_layout()
                           }

                           if(!is.proto(gui_layout))
                             stop("Layout must be a Container object")
                           if(!gui_layout$is("Container"))
                               stop("Layout must be a Container object")

                           if(is.null(gui_layout$context))
                             gui_layout$context <- .

                           gui_layout$make_ui(container=container) ## set up view
                           .$gui_layout <- gui_layout

                           .$init() ## initialize controller, ...
                         },
                         ## default layout if layout not specified
                         .doc_make_default_gui_layout=paste(
                           desc("Method describing default GUI layout. Default",
                                "is a table layout with each item being one row.")
                         ),
                         make_default_gui_layout=function(.) {
                           gui_layout <- do.call("aContainer", lapply(.$get_items(), function(i) i$make_default_gui_layout()))
                           gui_layout$context <- .
                           gui_layout
                         },
                         # visible
                         .doc_visible = paste(
                           desc("Method call to set visibility of items in group"),
                           param("value","A logical")
                           ),
                         visible=function(., value) {
                            sapply(.$get_items(), function(i) i$visible(value))
                         },
                         ## enabled
                         .doc_enabled = paste(
                           desc("Method call to set sensitivity of items in group"),
                           param("value","A logical")
                           ),
                         enabled=function(., value) {
                           sapply(.$get_items(), function(i) i$enabled(value))
                         },
                         ## in case itemgroup used as an item, this just passes off to make_gui
                         make_ui=function(., container, attr=.$attr, context=., ...) {
                           .$make_gui(container=container, attr=attr, context=context, ...)
                         },
                         ## init should set up the items controllers
                         ## as of now we call
                         ## make_ui which calls init


                         
                         ## initialize model
                         ## can be called before view is set
                         .doc_init_model=paste(
                           desc("Method to initialize model part of item group.",
                                "fixes names if given as named list in items",
                                "make get/set methods as part of item group",
                                "Initializes models for reach itemm")
                           ),                                
                         init_model=function(.) {
                           ## fix names in items if not specified
                           sapply(names(.$items), function(i) {
                             item <- .$items[[i]]
                             if(is.null(item$name) || (item$name == "Anonymous"))
                               item$name <- i
                             if(is.null(item$label) || (item$label == "Anonymous"))
                               item$label <- i
                           })

                           ## call init_model for each item
                           sapply(.$items, function(j) j$do_call("init_model"))

                           ## put parent into each item
                           sapply(.$items, function(j) j$parent <- .)
                           
                           ## make get_/set_pairs
                           sapply(.$get_items_only(), function(i) {
                             if(i$has_slot("name")) { 
                               property <- i$name
                               .$assign_if_null(paste("get_", i$name, sep=""),
                                      function(.) {
                                        .$get_item_by_name(property)$getattr(property)
#                                        i$get_model()$getattr(property)
                                      })
                               .$assign_if_null(paste("set_", i$name, sep=""),
                                                function(., value) {
                                                  .$get_item_by_name(property)$setattr(property, value)
                                                }
                                                )
                            }
                           })
                           ## add in undo
                           if(.$undo_do) {
                             sapply(.$items, function(i) {
                               property <- i$name
                               meth <- sprintf(".property_%s_value_changed", property)
                               .$assign_if_null(meth, function(., value, old_value) {
                                 .$undo_add(property, value, old_value)
                               })
                             })
                           }
                           
                           invisible()
                         },
                         ## set model (the items) for the itemgroup from another
                         .doc_set_model=paste(
                           desc("Set model for the itemgroup. The set of items should match",
                                "by type and name and not include any items groups."),
                           param("itemgroup","An itemgroup or dialog instance")
                           ),
                         set_model =function(., itemgroup) {
                           old_items <- .$get_items()
                           new_items <- itemgroup$get_items()
                           if(digest(sort(names(old_items))) == digest(sort(names(new_items)))) {
                             ## names match, lets replace models
                             for(i in names(old_items)) {
                               old_items[[i]]$set_model(new_items[[i]])
                             }
                           }
                         },
                         
                         ## initialize controllers
                         ## called after view is set
                         .doc_init_controller=paste(
                           desc("Initialize controller to connect view with model")
                           ),
                         init_controller=function(.) {
                           ## make self observe itself. This way we can add handlers
                           ## without making a controller
                           .$add_observer(.)
                           ## set up controllers
                           sapply(.$get_items(), function(i) i$init())
                           ## had code to call update_ui through a controller
                           ## this code is now in the dialog code. item groups do
                           ## not automatically have there UI updated?
                         },
                         ## initialize items -- called after making GUI
                         init=function(.) {
                           .$init_model()
                           .$init_controller()
                           .$update_ui()
                         },
                         ## instance
                         .doc_instance=paste(
                           desc("Create an instance of the item group (a clone of model, not GUI elements)")
                           ),
                         instance=function(.) {
                           obj <- .$proto()
                           obj$items <- lapply(obj$items, function(i) i$instance())
                           sapply(obj$items, function(i) i$parent <- obj) # add parent
                           if(!is.null(obj$gui_layout)) {
                             obj$gui_layout <- obj$gui_layout$instance()
##                              ## gui_layout children need to be instances too
##                              obj$gui_layout$container <- NULL
##                              obj$gui_layout$children <- lapply(obj$gui_layout$children, function(i) {
##                                if(is.proto(i))
##                                  i$instance()
##                                else
##                                  i
##                              })
                           }
                           obj
                         },
                         add_observer = function(., o) {
                           ## push down to models in items to add observer
                           sapply(.$get_items_only(), function(i) {
                             i$do_call("add_observer",list(o))
                           })
                           invisible()
                         },
                         remove_observer = function(., o) {
                           ## push down to models in items to add observer
                           sapply(.$get_items_only(), function(i) {
                             ## either an item or itemgroup so we push down
                               i$do_call("remove_observer",list(o))
                           })
                         },
                         
                         ## undo/redo stack implementation
                         ## THIS ISN"T WORKING
                         ## must add observer to listen for property value changes
                         ## as we don't call set_property to an itemgroup to change an item's value.
                         .doc_undo_do=paste(
                           desc("Logical property. If <code>TRUE</code> undo feature is implemented.")
                           ),
                         undo_do=TRUE,  # logical -- do we use undo.
                         undo_stack=list(),
                         undo_ptr = 0,
                         undo_stack_length=25, ## <= 0 to have unlimited
                         undo_trim_stack=function(.) {
                           if(.$undo_stack_length > 0) {
                             stack <- .$undo_stack
                             if((n <- length(stack)) > (m <- .$undo_stack_length)) {
                               .$undo_ptr <- max(0, .$undo_ptr-(n-m))
                               stack <- stack[(n-m):n]
                               .$undo_stack <- stack
                             }
                           }
                         },
                         undo_add=function(., property, value, old_value, ...) {
                           if(.$undo_ptr > 0)
                             stack <- .$undo_stack[1:.$undo_ptr]
                           else
                             stack <- list()
                           stack[[length(stack)+1]] <- merge(list(property=property, value=value, old_value=old_value),
                                                             list(...))
                           .$undo_stack <- stack
                           .$undo_ptr <- length(stack)
                           .$undo_trim_stack()
                         },
                         undo_can_undo=function(.) .$undo_ptr > 0,
                         undo=function(.) {
                           if(.$undo_can_undo()) {
                             val <- .$undo_stack[[.$undo_ptr]]
                             .$get_item_by_name(val$property)$setattr(val$property, val$old_value, notify_private=FALSE) ## by pass set_property
                             .$undo_ptr <- .$undo_ptr-1
##                             .$update_ui()
                           }
                         },
                         undo_can_redo=function(.) .$undo_ptr < length(.$undo_stack),
                         redo=function(.) {
                           if(.$undo_can_redo()) {
                             val <- .$undo_stack[[.$undo_ptr + 1]]
                             .$get_item_by_name(val$property)$setattr(val$property, val$value, notify_private=FALSE)
                             .$undo_ptr <- .$undo_ptr + 1
##                             .$update_ui()                             
                           }
                         },
                         
                         
                         ## methods to call down into items
                         .doc_update_ui=paste(
                           desc("Method to update user interface.")
                           ),
                         update_ui=function(.) {
                           ## update items
                           sapply(.$get_items(), function(i) i$update_ui())
                           ## update gui_layout
                           if(.$has_slot("gui_layout") && !is.null(.$gui_layout))
                             .$gui_layout$update_ui()
                         }
                         )

## a constructor
##' Constructor for ItemGroup instances
##'
##' An ItemGroup creates a model with properties given by the items
##' and a default layout for its items. This can also be specified
##' when the layout is drawn through \code{make_gui}.
##'
##' An item group bundles a list of items into a model. When
##' the model is intialized, constructors to access the model values
##' are created. These getters/setters use the item names, so that
##' \code{get_name} will get the main value for the item with name
##' attribute "name".
##' 
##' An item group has the useful methods \code{to_R} to return the
##' values in the model as a named list and \code{get_item_by_name}
##' to get the item from the list of items matching the name.
##' @param items List of Item instances or ItemGroup instances
##' @param name Name of ItemGroup.
##' @param ... Passed to ItemGroup proto trait
##' @seealso \code{\link{aContainer}} for specifying a layout
#     \code{\link{aDialog}} for an extension to an ItemGroup that creates its own window.
##' @export
##' @return A \code{proto} object. Call \code{obj$show_help()} to view its methods and properties.
##' 
##' @examples
##' \dontrun{
##' ## make a simple item group, show in non-default layout
##' i <- anItemGroup(items=list(
##'                  numericItem(0,"x"),
##'                  numericItem(0,"y"),
##'                  stringItem("","z")
##'                  ))
##' lay <- aContainer("x","y", aFrame("z", label="z in a box"))
##' ## some proto methods:
##' i$make_gui(cont=gwindow("Example of itemGroup"), gui_layout=lay)
##' i$get_x()     # get x value
##' i$set_x(10)   # set x value to 10
##' i$to_R()      # get list of x,y,z values
##' }
##' 
##' ## example of using an item group and gbasicdialog to make a modal  GUI
##' ig <- anItemGroup(items=list(
##'                     x=numericItem(2)
##'                     )
##'                   )
##' 
##' ## using gbasicdialog from gWidgets
##' \dontrun{
##' w <- gbasicdialog("testing", handler=function(h,...) {
##'   . <- h$action                         # action passes in itemgroup
##'   .$output <- sin(.$get_x())
##' },
##'                   action=ig)
##' ig$make_gui(container=w)
##' visible(w, TRUE)  ## modal now
##' print(ig$output)
##' }


anItemGroup <- function(items=list(),
                      name, ...) {

  if(missing(name))
    name <- "Anonymous"
  
  obj <- ItemGroup$proto(items=items,
                         name=name,
                         ...,
                         funEnvir=FALSE)
  obj$init_model()
  return(obj)
}
