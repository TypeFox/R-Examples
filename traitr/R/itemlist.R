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

#' @include items.R
roxygen()

## TODO: fix set_model, get_model

#' An itemList is used to store a list of similar items or itemgroups with a means to edit individually
#'
#' @export
#' @param items list of similar items, may be empty list
#' @param items_name Header name on top of table displaying item list
#' @param item_factory function to call to produce a new item, e.g. \code{function(.) numericItem(1)}
#' @param name name of itemList object
#' @param label label for itemList object
#' @param help help string
#' @param tooltip tooltip
#' @param attr attributes passed to \code{make_ui} constructor
#' @param model optional model to pass in
#' @param editor optional editor to pass in
#' @param ... passed along to \code{Item\$proto()} call
#' @note  This item's model is a a list storing child items or item groups.
#'        To create new items, the \code{item_factory} method should be provided. It provides a
#'        template for a new item, the editor allows the user to modify its values
#'        When a child item is edited the "done" button is clicked to close. The method \code{post_process}
#'        is called. (The edited changes may already have been sent back to the model.)
#'        The child items \code{to_string} method is called to make the label in the table that allows
#'        the user to select the child item to edit. This should be a character vector of length 1.
#'        The table can display an icon. Simply set the \code{icon} property of the icon to a \pkg{gWidgets}
#'        stock icon name.
#'
#'        The child items can be returned via the \code{get_value} method or the \code{get_NAME} method, where
#'        \code{NAME} is that passed into the \code{name} argument of the constructor.
#'        The \code{to_R} method can be modified to manipulate the return value. The vignette has an example
#'        where the output is coerced into a data frame. The default is a list with each child items \code{to_R}
#'        method called to form the numbered components.
#'          
#' @return A \code{proto} object. Call \code{obj\$show_help()} to view its methods and properties.
#' @examples
#' \dontrun{
#' ## make icons
#' imagedir <- system.file("images",package="traitr")
#' addStockIcons(gsub("\\\\.png","", list.files(path=imagedir)),
#'               list.files(path=imagedir, full.names=TRUE))
#' ## make item
#' item <- itemList(items=list(),
#'                  items_name="Personnel",
#'                  item_factory = function(.) {
#'                    a <- anItemGroup(items=list(
#'                                       name=stringItem(""),
#'                                       rank=choiceItem("Scout",
#'                                            values=c("Scout","Captain","General")),
#'                                       serial.number = stringItem("", label="Serial number")))
#'                    a$post_process <- function(.) {
#'                     .$icon <- tolower(.$get_rank())
#'                     }
#'                    a$to_string <- function(., drop=TRUE) .$to_R()$name
#'                    return(a)
#'                  },
#'                  name="itemlist")
#' 
#' item$make_ui(container=gwindow("itemList test"))
#' }
itemList <- function(items=list(),
                     items_name="",
                     item_factory = NULL,
                     name,
                     label=name,
                     help="",
                     tooltip="",
                     attr=list(),
                     model,
                     editor,
                     ...) {

  
  if(missing(name))
    name <- "Anonymous"
  
  
  
  obj <- Item$proto(name=name, label=label, help=help,
                    tooltip=tooltip, value=items,
                    items_name=items_name,
                    item_factory=item_factory,
                    add_handler_name=paste("addHandler", c("Changed"), sep=""),
                    ...
                    )
  obj$add_class("itemList")
  obj$clear_message="Select item to edit"

  ## Implement methods
  
  ## Append and in such a way that property_value_changed is signaled
  obj$.doc_append_item=paste(
    desc("Append an item to the item list"),
    param("value","An item or itemGroup object to append")
    )
  obj$append_item <- function(., value) {
    values <- .$get_value()
    values <- c(values, value)
    ## we don't set the value, we set the name?
    .$set_value(values)
    ## listen to changes in values
    value$add_observer(.)
  }
  ## Remove an item
  obj$.doc_remove_item=paste(
    desc("Remove an item to the item list"),
    param("value","An item index, or  item or itemGroup object to remove from item list")
    )
  obj$remove_item <- function(., value) {
    values <- .$get_value()
    if(is.numeric(value)) {
      i <- as.integer(value)
      item <- values[[i]]
      item$remove_observer(.)           #remove observer
      values <- values[-i]
    } else {
      for(i in seq_along(length(values))) {
        if(values[[i]]$identical(value))
          values[[i]] <- NULL
      }
      value$remove_observer(.)
    }
    .$set_value(values)
  }

  obj$.doc_post_process=paste(
    desc("Method called when user clicks 'done' button when editing an item. Use NULL for no call")
    )
  obj$post_process <- NULL
  
  
  ## get values
  obj$.doc_to_R=paste(
    desc("Return R object representing state of object. Returns list of same method call on each object.")
    )
  obj$to_R <- function(., drop=TRUE) {
    l <- sapply(.$get_value(), function(i) i$to_R(drop=drop))
    l <- list(l)
    names(l) <- .$name
    l
  }
  obj$.doc_to_string=paste(
    desc("String representation of object")
    )
  obj$to_string <- function(., drop=TRUE) {
    paste(sapply(.$get_value(), function(i) i$to_string(drop=drop)), collapse=" ")
  }

  ## setup and go
  
  if(!missing(attr))
    obj$attr <- merge(obj$attr, attr)
  
  if(missing(editor))
    obj$editor <- ItemListEditor$proto()

  if(!missing(model)) 
    obj$set_model(model)

  obj$init_model <- function(.) {
    .$next_method("init_model")(.)
    ## map name -> value
    .$get_value <- function(.) .$do_call(sprintf("get_%s",.$name))
    .$set_value <- function(., value) .$do_call(sprintf("set_%s", .$name), list(value))
    ## do same for model
    m <- .$get_model()
    m$..name <- .$name
    m$get_value <- function(.) .$do_call(sprintf("get_%s", .$..name))
    m$set_value <- function(., value) .$do_call(sprintf("set_%s", .$..name), list(value))
#    .$set_slot(sprintf("get_%s",.$name), .$get_slot("get_value"))
#    .$set_slot(sprintf("set_%s",.$name), .$get_slot("set_value"))
    ## observe self
    .$add_observer(.)
  }

  obj$init_model()

  ## set up controller, by overriding default
  obj$init_controller <- function(.,...) {}
#  obj$property_value_value_changed <- function(.,...) {
#    if(.$editor$is_realized())
#      .$editor$update_ui(.)
#  }
  obj$update_ui <- function(., context=.) .$editor$update_ui(context=context)
  
  obj$model_value_changed <- function(.) {
    if(.$editor$is_realized())
      .$update_ui(.)
  }

  ## observe self
#  obj$add_observer(obj)
  return(obj)
  

}

#' trait for editor for itemList
#'
#' @export
ItemListEditor <- Editor$proto(class=c("ItemListEditor", Editor$class),
                               ## return from items, values for table widget
                               get_table_entries=function(., items, name) {
                                 ## figure out icons
                                 df <- data.frame(#icon=character(0),
                                                  description=character(0),
                                                  stringsAsFactors=FALSE)
                                 if(length(items)) {
                                   df <- data.frame(#icon=sapply(items, function(i) i$icon),
                                                    description= sapply(items, function(i) i$to_string()),
                                                    stringsAsFactors=FALSE)
                                 }
                                 if(!missing(name))
                                   names(df)[1] <- name
                                 return(df)
                               },

                               ## clear edit area and selection
                               clear_edit_area = function(., message="") {
                                 widgets <- .$widgets
                                 delete(widgets[['edit_area']], widgets[['edit_area_child']])
                                 widgets[['edit_area_child']] <- glabel(message, cont=widgets[['edit_area']])
                                 .$widgets <- widgets
                                 svalue(widgets[['item_table']], index=TRUE) <- 0 # clear selection
                               },
                               ## edit an item
                               edit_item = function(., context, item) {
                                 widgets <- .$widgets
                                 delete(widgets[['edit_area']], widgets[['edit_area_child']])
                                 widgets[['edit_area_child']] <- (g <- ggroup(horizontal=FALSE,
                                                                              cont=widgets[['edit_area']], expand=TRUE))
                                 .$widgets <- widgets

                                 ## Don't want label when showing edit area
                                 inst <- item$instance()
                                 inst$show_label <- FALSE
                                 if(inst$has_slot("make_gui"))
                                   inst$make_gui(container=g)
                                 else
                                   inst$make_ui(container=g)
                                 gseparator(cont=g)
                                 bg <- ggroup(cont=g)
                                 ## XXX Need to abstract out buttons -- they should be
                                 ## item type specific
                                 gbutton(gettext("done"), cont=bg, handler=function(h,...) {
                                   ## call some method for item?
                                   item$do_call("post_process")
                                   .$clear_edit_area(context$clear_message)
                                   .$update_ui(context)
                                 })
                               },
                               ## make interface. Should also have a concise one
                               make_ui = function(., container, attr=.attr, context=., ...) {
                                 ## context is item object
                                 ## context$get_model() is the model
                                 ## 3 handlers: new, edit, remove selected
                                 new_item <- function(h,...) {
                                   ## make new item
                                   context <- h$action$context
                                   item <- context$item_factory()
                                   context$append_item(item)
                                   svalue(widgets[['item_table']], index=TRUE) <- length(context$get_model()$get_value())
                                   .$edit_item(context, item)
                                   
                                 }

                                 edit_item <- function(h,...) {
                                   ## call make_Ui, replace previous
                                   context <- h$action$context
                                   ind <- svalue(h$obj, index=TRUE)
                                   if(length(ind)) {
                                     model <- context$get_model()
                                     item <- model$get_value()[[ind]]
                                     .$edit_item(context, item)
                                   }
                                 }

                                 remove_selected_item <- function(h,...) {
                                   ## remove selection from model
                                   ind <- svalue(widgets[['item_table']], index=TRUE)
                                   context <- h$action$context
                                   context$remove_item(ind)
                                   .$clear_edit_area(context$clear_message)
                                 }


                                 ## Let's go
                                 items <- context$get_model()$get_value()
                                 df <- .$get_table_entries(items, context$get_slot("items_name"))
                                 widgets <- list()
                                 attr <- merge(list(container=container, expand=TRUE), attr)
                                 widgets[['box']] <- do.call("gpanedgroup", attr)

                                 lgroup <- ggroup(horizontal=FALSE, cont=widgets[['box']])
                                 widgets[['edit_area']] <- ggroup(horizontal=FALSE, cont=widgets[['box']])
                                 widgets[['edit_area_child']] <- glabel("", cont=widgets[['edit_area']])
                                 widgets[['item_table']] <- gtable(df, cont=lgroup, expand=TRUE,
                                                                   icon.FUN=function(frame) {
                                                                     frame <- as.data.frame(frame)
                                                                     if(nrow(frame) == 0)
                                                                       character(0)
                                                                     else
                                                                       rep("", nrow(frame))
                                                                   })
                               addHandlerClicked(widgets[['item_table']], handler=edit_item, action=list(context=context))

                                 ## Buttons
                                 bgroup <- ggroup(horizontal=TRUE, cont=lgroup)

                                 ## Make add button (XXXX condition on there being a factor function
                                 if(is.function(context$item_factory)) {
                                   widgets[['add_button']] <- gbutton("add", cont=bgroup)
                                   addHandlerClicked(widgets[['add_button']], handler=new_item, action=list(context=context))
                                 }

                                 ## Minus button
                                 widgets[['minus_button']] <- gbutton("remove", cont=bgroup)

                                 ## remove only there if selection is made
                                 enabled(widgets[['minus_button']]) <- FALSE
                                 addHandlerClicked(widgets[['item_table']], handler=function(h,...) {
                                   val <- length(svalue(h$obj, index=TRUE)) == 0
                                   enabled(h$action) <- !val
                                 }, action=widgets[['minus_button']])

                                 addHandlerClicked(widgets[['minus_button']], handler=remove_selected_item, action=list(context=context))

                                 ## store widgest, initialize
                                 .$widgets <- widgets
                                 .$clear_edit_area(context$clear_message) # initialize
                                 
                               },

                               update_ui = function(., context) {
                                 if(!context$editor$is_realized())
                                   return()
###                                 items <- context$get_model()$get_value()
                                 items <- context$get_value()
                                 df <- .$get_table_entries(items)
                                 if(length(items) && nrow(df)) {
                                   .$widgets[['item_table']][] <- df
                                   ## do icons with hack using RGtk2 backend -- ugly
                                   frame <- getToolkitWidget(.$widgets[['item_table']])
                                   if(inherits(frame, "RGtkObject")) {
                                     require(RGtk2)
                                     frame <- frame$getModel()                                   
                                     icons <- gWidgetsRGtk2:::getstockiconname(sapply(items, function(i) i$icon))
                                     frame[,2] <- icons
                                   }
                                 } else {
                                   ## clear out
                                   .$widgets[['item_table']][] <- character(0)
                                 }
                               }
                               )
