# most of these functions are implemented manually for the same reason as their C counterpart

# reason: var args, paste 'em together
gtkMessageDialogNew <-
function(parent = NULL, flags, type, buttons, ..., show = TRUE)
{
    checkPtrType(parent, "GtkWindow", nullOk = T)

    w <- .RGtkCall("S_gtk_message_dialog_new", parent, flags, type, buttons, paste(...))

    if(show)
        gtkWidgetShowAll(w)

    return(w)
}
# reason: same as above
gtkMessageDialogFormatSecondaryMarkup <-
function(object, ...)
{
        checkPtrType(object, "GtkMessageDialog")

        w <- .RGtkCall("S_gtk_message_dialog_format_secondary_markup", object, paste(...))

        return(invisible(w))
}
# reason: same as above
gtkMessageDialogFormatSecondaryText <-
function(object, ...)
{
        checkPtrType(object, "GtkMessageDialog")

        w <- .RGtkCall("S_gtk_message_dialog_format_secondary_text", object, paste(...))

        return(invisible(w))
}
# reason: same as above
gtkMessageDialogNewWithMarkup <-
function(parent, flags, type, buttons, ..., show = TRUE)
{
        checkPtrType(parent, "GtkWindow")

        w <- .RGtkCall("S_gtk_message_dialog_new_with_markup", parent, flags, type, buttons, paste(...))

        if(show)
                gtkWidgetShowAll(w)

        return(w)
}
# reason: var args, make two vectors, one for labels the other for responses
gtkDialogNewWithButtons <-
function(title = NULL, parent = NULL, flags = 0, ..., show = TRUE)
{
    title <- as.character(title)
    checkPtrType(parent, "GtkWindow", nullOk = TRUE)
		
		args <- list(...)
    args <- args[!is.null(args)]
    if (length(args) %% 2 != 0)
      stop("Must have one stock ID for every response type")
    args_split <- split(args, rep(c(1,2), length(args) / 2))
    labels <- as.character(args_split[[1]])
    responses <- as.list(as.integer(args_split[[2]]))
    
    w <- .RGtkCall("S_gtk_dialog_new_with_buttons", title, parent, flags, labels, responses)

    if(show)
      gtkWidgetShowAll(w)
    
    return(w)
}
gtkDialogAddButtons <-
function(object, ...)
{
        checkPtrType(object, "GtkDialog")

		args <- list(...)
		labels <- as.character(args[seq(1,length(args),by=2)])
		responses <- args[seq(2,length(args),by=2)]
		
        w <- .RGtkCall("S_gtk_dialog_add_buttons", object, labels, responses)

        return(invisible(w))
}


# reason: var-args for the buttons, compile into vectors of labels and responses (if given)
gtkFileChooserDialogNewWithBackend <-
function(title = NULL, parent = NULL, action, backend, ..., show = TRUE)
{
        title <- as.character(title)
        checkPtrType(parent, "GtkWindow", nullOk = T)
        backend <- as.character(backend)
        
		args <- list(...)
		
		labels <- NULL
		responses <- NULL
		if (length(args) > 1) {
			labels <- as.character(args[seq(1,length(args),by=2)])
			responses <- args[seq(2,length(args),by=2)]
		}

        w <- .RGtkCall("S_gtk_file_chooser_dialog_new_with_backend", title, parent, action, backend, labels, responses)

        if(show)
                gtkWidgetShowAll(w)

        return(w)
}
gtkFileChooserDialogNew <-
function(title = NULL, parent = NULL, action, ..., show = TRUE)
{
        w <- gtkFileChooserDialogNewWithBackend(title, parent, action, NULL, ..., show=show)

        return(w)
}
gtkRecentChooserDialogNewForManager <-
function(title = NULL, parent = NULL, manager, ..., show = TRUE)
{
  title <- as.character(title)
  checkPtrType(parent, "GtkWindow", nullOk = T)
  checkPtrType(manager, "GtkRecentManager", nullOk = T)
        
  args <- list(...)
		
  labels <- NULL
  responses <- NULL
  if (length(args) > 1) {
    labels <- as.character(args[seq(1,length(args),by=2)])
    responses <- args[seq(2,length(args),by=2)]
  }

  w <- .RGtkCall("S_gtk_recent_chooser_dialog_new_for_manager", title, parent, manager, labels, responses)

  if(show)
    gtkWidgetShowAll(w)

  return(w)
}
gtkRecentChooserDialogNew <-
function(title = NULL, parent = NULL, ..., show = TRUE)
{
        w <- gtkRecentChooserDialogNewForManager(title, parent, NULL, ..., show=show)

        return(w)
}

# reason: var-args - not yet ready for primetime

#gtkDialogSetAlternativeButtonOrder <-
#function(object, ...)
#{
#        checkPtrType(object, "GtkDialog")
#
#        w <- .RGtkCall("S_gtk_dialog_set_alternative_button_order", object, ..., -1)
#
#        return(invisible(w))
#}

# reason: var-args, implemented from scratch on C side
gtkShowAboutDialog <-
function(parent, ...)
{
        checkPtrType(parent, "GtkWindow")

        w <- .RGtkCall("S_gtk_show_about_dialog", parent, list(...))

        return(invisible(w))
}

# reason: var-args, reimplemented on C side
gtkTextBufferCreateTag <-
function(object, tag.name, ...)
{
        checkPtrType(object, "GtkTextBuffer")
        tag.name <- as.character(tag.name)
        
        w <- .RGtkCall("S_gtk_text_buffer_create_tag", object, tag.name, list(...))

        return(w)
}

# reason: var args, just compile into an array and send to alternate function
gtkListStore <- gtkListStoreNew <-
function(...)
{
        types <- checkArrType(c(...), as.GType)
        w <- .RGtkCall("S_gtk_list_store_newv", types)

        return(w)
}

# reason: var args, just compile into an array and send to alternate function
gtkTreeStoreNew <-
function(...)
{
        types <- checkArrType(c(...), as.GType)
        w <- .RGtkCall("S_gtk_tree_store_newv", types)

        return(w)
}

# reason: var args, break apart the args into cols and values vectors
gtkListStoreSet <-
function(object, iter, ...)
{
        checkPtrType(object, "GtkListStore")
        checkPtrType(iter, "GtkTreeIter")
		
		args <- list(...)
		cols <- as.integer(unlist(args[seq(1,length(args),by=2)]))
		if (length(args) == 2)
			values <- args[2]
		else values <- args[seq(2,length(args),by=2)]
		
        w <- .RGtkCall("S_gtk_list_store_set", object, iter, cols, values)

        return(invisible(w))
}

# reason: var args, let's do it completely in R

gtkTreeStoreSet <-
function(object, iter, ...)
{
        checkPtrType(object, "GtkTreeStore")
        checkPtrType(iter, "GtkTreeIter")
		
		args <- list(...)
		cols <- as.integer(unlist(args[seq(1,length(args),by=2)]))
		if (length(args) == 2)
			values <- args[2]
		else values <- args[seq(2,length(args),by=2)]
		
		w <- .RGtkCall("S_gtk_tree_store_set", object, iter, cols, values)
}

# reason: var args, compile into named list and send to C
gtkCellLayoutSetAttributes <-
function(object, cell, ...)
{
        checkPtrType(object, "GtkCellLayout")
        checkPtrType(cell, "GtkCellRenderer")
		attributes <- lapply(c(...), as.integer)
		
        w <- .RGtkCall("S_gtk_cell_layout_set_attributes", object, cell, attributes)

        return(invisible(w))
}

# reason: just var-args, so concat them as a named vector
gtkContainerAddWithProperties <-
function(object, widget, ...)
{
        checkPtrType(object, "GtkContainer")
        checkPtrType(widget, "GtkWidget")

        w <- .RGtkCall("S_gtk_container_add_with_properties", object, widget, c(...))

        return(invisible(w))
}

# reason: more var-args, like above
gtkContainerChildSet <-
function(object, child, ...)
{
        checkPtrType(object, "GtkContainer")
        checkPtrType(child, "GtkWidget")

        w <- .RGtkCall("S_gtk_container_child_set", object, child, list(...))

        return(invisible(w))
}

# reason: var args, compile and coerce to strings
gtkContainerChildGet <-
function(object, child, ...)
{
        checkPtrType(object, "GtkContainer")
        checkPtrType(child, "GtkWidget")
        props <- as.character(c(...))

        w <- .RGtkCall("S_gtk_container_child_get", object, child, props)

        return(invisible(w))
}

# reason: var args, make a vector of column ids
gtkTreeModelGet <-
function(object, iter, ...)
{
        checkPtrType(object, "GtkTreeModel")
        checkPtrType(iter, "GtkTreeIter")
		cols <- as.integer(c(...))
		
        w <- .RGtkCall("S_gtk_tree_model_get", object, iter, cols)

        return(w)
}

# reason: var args, make a vector of indices

gtkTreePathNewFromIndices <-
function(...)
{
        indices <- as.integer(c(...))

        w <- .RGtkCall("S_gtk_tree_path_new_from_indices", indices)

        return(w)
}

# reason: var args, convert to named list

gtkTreeViewInsertColumnWithAttributes <-
function(object, position, title, cell, ...)
{
        checkPtrType(object, "GtkTreeView")
        position <- as.integer(position)
        title <- as.character(title)
        checkPtrType(cell, "GtkCellRenderer")
		attributes <- lapply(c(...), as.integer)
		
        w <- .RGtkCall("S_gtk_tree_view_insert_column_with_attributes", object, position, title, cell, attributes)

        return(w)
}

# reason: var-args, collect into vectors

gtkTextBufferInsertWithTags <-
function(object, iter, text, ...)
{
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(iter, "GtkTextIter")
  text <- as.character(text)
  tags <- list(...)
  checkArrType(tags, function(x) checkPtrType(x, "GtkTextTag"))
  
  w <- .RGtkCall("S_gtk_text_buffer_insert_with_tags", object, iter, text, -1L,
                 tags)

  return(invisible(w))
}

# reason: same as above

gtkTextBufferInsertWithTagsByName <-
function(object, iter, text, ...)
{
        checkPtrType(object, "GtkTextBuffer")
        checkPtrType(iter, "GtkTextIter")
        text <- as.character(text)
		tagNames <- list(...)
		tagNames <- as.character(tagNames)

        w <- .RGtkCall("S_gtk_text_buffer_insert_with_tags_by_name", object, iter, text, tagNames)

        return(invisible(w))
}

# reason: redirect to set_with_data, GObject confuses code generator
gtkClipboardSetWithOwner <-
function(object, targets, get.func, owner = NULL)
{
        checkPtrType(object, "GtkClipboard")
        targets <- checkArrType(targets, function(x) {x <- as.GtkTargetEntry(x); x })
        get.func <- as.function(get.func)
        if (!is.null( owner )) checkPtrType(owner, "GObject")

        w <- .RGtkCall("S_gtk_clipboard_set_with_data", object, targets, get.func, owner)

        return(w)
}

# reason: position arg is omitted automatically, due to array

gtkListStoreInsertWithValuesv <-
function(object, position, columns, values)
{
  checkPtrType(object, "GtkListStore")
  position <- as.integer(position)
  columns <- as.list(as.integer(columns))
  values <- as.list(values)

  w <- .RGtkCall("S_gtk_list_store_insert_with_valuesv", object, position, columns, values)

  return(invisible(w))
}
gtkTreeStoreInsertWithValuesv <-
function(object, parent, position, columns, values)
{
  checkPtrType(object, "GtkTreeStore")
  checkPtrType(parent, "GtkTreeIter")
  position <- as.integer(position)
  columns <- as.list(as.integer(columns))
  values <- as.list(values)

  w <- .RGtkCall("S_gtk_tree_store_insert_with_valuesv", object, parent, position, columns, values, PACKAGE = "RGtk2")

  return(invisible(w))
}

# reason: here are some functions where we just leave off the text length parameter for convenience
gtkTextBufferInsertInteractive <-
function(object, iter, text, default.editable)
{
        checkPtrType(object, "GtkTextBuffer")
        checkPtrType(iter, "GtkTextIter")
        text <- as.character(text)
        default.editable <- as.logical(default.editable)

        w <- .RGtkCall("S_gtk_text_buffer_insert_interactive", object, iter, text, -1, default.editable)

        return(w)
}
gtkTextBufferInsertInteractiveAtCursor <-
function(object, text, default.editable)
{
        checkPtrType(object, "GtkTextBuffer")
        text <- as.character(text)
        default.editable <- as.logical(default.editable)

        w <- .RGtkCall("S_gtk_text_buffer_insert_interactive_at_cursor", object, text, -1, default.editable)

        return(w)
}
gtkIMContextSetSurrounding <-
function(object, text, cursor.index)
{
        checkPtrType(object, "GtkIMContext")
        text <- as.character(text)
        cursor.index <- as.integer(cursor.index)

        w <- .RGtkCall("S_gtk_im_context_set_surrounding", object, text, -1, cursor.index)

        return(invisible(w))
}

# reason: this one leaves off dimensions that must be user specified
gtkIMContextSimpleAddTable <-
function(object, data, max.seq.len, n.seqs)
{
        checkPtrType(object, "GtkIMContextSimple")
        data <- as.list(as.integer(data))
		max.seq.len <- as.integer(max.seq.len)
        n.seqs <- as.integer(n.seqs)

        w <- .RGtkCall("S_gtk_im_context_simple_add_table", object, data, max.seq.len, n.seqs)

        return(w)
}

# reason: for convenience give defaults for parameters
gtkSelectionDataSet <-
function(object, type = object[["target"]], format = 8L, data)
{
        checkPtrType(object, "GtkSelectionData")
        type <- as.GdkAtom(type)
        format <- as.integer(format)
        data <- as.list(as.raw(data)) # inefficient for large data

        w <- .RGtkCall("S_gtk_selection_data_set", object, type, format, data)

        return(invisible(w))
}

# reason: need to omit the length param here - string arrays are normally NULL-terminated
gtkIconThemeSetSearchPath <-
function(object, path)
{
        checkPtrType(object, "GtkIconTheme")
        path <- as.list(as.character(path))
        
        w <- .RGtkCall("S_gtk_icon_theme_set_search_path", object, path)

        return(invisible(w))
}
# reason: var-args
gtkDialogSetAlternativeButtonOrder <-
function(object, ...)
{
        checkPtrType(object, "GtkDialog")
        new.order <- list(...)
		
		w <- gtkDialogSetAlternativeButtonOrderFromArray(object, new.order)

        return(invisible(w))
}
# reason: more var-args
gtkListStoreInsertWithValues <-
function(object, position, ...)
{
        checkPtrType(object, "GtkListStore")
        position <- as.integer(position)

		args <- list(...)
		columns <- as.integer(args[seq(1,length(args),by=2)])
		values <- args[seq(2,length(args),by=2)]
		
        w <- gtkListStoreInsertWithValuesv(object, position, columns, values)

        return(w)
}
gtkTreeStoreInsertWithValues <-
function(object, parent, position, ...)
{
  checkPtrType(object, "GtkListStore")
  checkPtrType(parent, "GtkTreeIter")
  position <- as.integer(position)

  args <- list(...)
  columns <- as.integer(args[seq(1,length(args),by=2)])
  values <- args[seq(2,length(args),by=2)]
		
  w <- gtkTreeStoreInsertWithValuesv(object, parent, position, columns, values)

  return(w)
}
# reason: user func has no userdata, so we can't support it... maybe in the future
gtkMenuAttachToWidget <-
function(object, attach.widget)
{
        checkPtrType(object, "GtkMenu")
        checkPtrType(attach.widget, "GtkWidget")

        w <- .RGtkCall("S_gtk_menu_attach_to_widget", object, attach.widget, PACKAGE = "RGtk2")

        return(invisible(w))
}

# reason: var-args, let's just go right to gObjectSet
gtkWidgetSet <- gObjectSet

# reason: add ... argument so that we can use this on a variety of signals
gtkWidgetDestroy <-
function(object, ...)
{
	checkPtrType(object, "GtkWidget")

	w <- .RGtkCall("S_gtk_widget_destroy", object, PACKAGE = "RGtk2")

	return(invisible(w))
} 

# reason: these two are var-args and we're just going to tie them into gObjectNew
gtkWidgetNew <-
function(type, ..., show = TRUE)
{
	if (!("GtkWidget" %in% gTypeGetAncestors(type)))
		stop("GType must inherit from GtkWidget")
	gObjectNew(type, ...)
}
gtkObject <- gtkObjectNew <-
function(type, ...)
{
	if (!("GtkObject" %in% gTypeGetAncestors(type)))
		stop("GType must inherit from GtkObject")
	gObjectNew(type, ...)
}

# reason: var-args, just use _style_get_property for each
gtkWidgetStyleGet <-
function(object, ...)
{
        checkPtrType(object, "GtkWidget")
        props <- c(...)
		w <- sapply(props, function(prop) { gtkWidgetStyleGetProperty(object, prop) })
        return(invisible(w))
}

# reason: var-args, just reimplement in R
gtkTreeViewColumnNewWithAttributes <-
function(title, cell, ...)
{
		title <- as.character(title)
		checkPtrType(cell, "GtkCellRenderer")
		
		column <- gtkTreeViewColumnNew()
		column$setTitle(title)
		column$packStart(cell, TRUE)
		column$setAttributes(cell, ...)
		
        return(column)
}
# reason: var-args, just implement in R using _add_attribute for each
gtkTreeViewColumnSetAttributes <-
function(object, cell.renderer, ...)
{
        checkPtrType(object, "GtkTreeViewColumn")
        checkPtrType(cell.renderer, "GtkCellRenderer")
		attributes <- lapply(c(...), as.integer)
		print(attributes)
		
		object$clearAttributes(cell.renderer)
		w <- sapply(names(attributes), function(attr) { 
			object$addAttribute(cell.renderer, attr, attributes[[attr]]) 
		})
        
        return(invisible(w))
}

## x and y are in-out parameters
gtkTreeViewGetTooltipContext <-
  function(object, x, y, keyboard.tip)
{
  checkPtrType(object, "GtkTreeView")
  x <- as.integer(x)
  y <- as.integer(y)
  keyboard.tip <- as.logical(keyboard.tip)

  w <- .RGtkCall("S_gtk_tree_view_get_tooltip_context", object, x, y, keyboard.tip, PACKAGE = "RGtk2")

  return(w)
}
gtkIconViewGetTooltipContext <-
  function(object, x, y, keyboard.tip)
{
  checkPtrType(object, "GtkIconView")
  x <- as.integer(x)
  y <- as.integer(y)
  keyboard.tip <- as.logical(keyboard.tip)

  w <- .RGtkCall("S_gtk_icon_view_get_tooltip_context", object, x, y, keyboard.tip, PACKAGE = "RGtk2")

  return(w)
} 

## varargs
gtkStyleGet <-
  function(object, widget.type, first.property.name, ...)
{
  checkPtrType(object, "GtkStyle")
  widget.type <- as.GType(widget.type)
  property.names <- as.character(c(first.property.name, ...))

  values <- lapply(property.names, object$getStyleProperty,
                   widget.type = widget.type)
  lapply(values, `[[`, "value")
}

## varargs
gtkInfoBarNewWithButtons <-
  function(first.button.text, ..., show = TRUE)
{
  w <- gtkInfoBar(show = show)
  w$addButtons(first.button.text, ...)
  w
}
gtkInfoBarAddButtons <-
  function(object, first.button.text, ...)
{
  checkPtrType(object, "GtkInfoBar")
  args <- c(first.button.text, list(...))
  labels <- args[seq(1, length(args), 2)]
  responses <- args[seq(2, length(args), 2)]
  invisible(mapply(object$addButton, labels, responses))
}

## Create an instance of RGtkBuilder, which can find the types by name

.initClasses <- function() {
  if (boundGTKVersion() >= "2.12.0") {
    gClass("RGtkBuilder", "GtkBuilder",
           GtkBuilder = list(
             get_type_from_name = function(self, name) as.GType(name)
             )
           )
  }
}

gtkBuilderNew <- function() gObject("RGtkBuilder")

# reason: unfortunately, the 'group' must be an exact pointer match so 
# we delegate to 'gtkRadioButtonNewFromWidget' with the first element.
gtkRadioButtonNew <-
function(group = NULL, show = TRUE)
gtkRadioButtonNewFromWidget(group[[1]], show)

gtkRadioButtonNewWithLabel <-
function(group = NULL, label, show = TRUE)
gtkRadioButtonNewWithLabelFromWidget(group[[1]], label, show) 

## Unlike GtkRadioButton, the GtkRadioMenuItem from_widget
## constructors do not accept a NULL value for 'widget'. Thus, we have
## to catch NULL and call the original function (where NULL becomes an
## empty list). Otherwise, we call from_widget as in the above. In the
## GTK+ source, it looks like they have written to code to be robust
## to NULL, but there is an assertion up-front that fails. Since the
## function is not documented to accept NULL, this must be the desired
## behavior.

gtkRadioMenuItemNew <- function(group = NULL, show = TRUE) {
  if (is.null(group)) {
    w <- .RGtkCall("S_gtk_radio_menu_item_new", group, PACKAGE="RGtk2")
    if (show)
      w$show()
    w
  } else gtkRadioMenuItemNewFromWidget(group[[1]], show)
}

gtkRadioMenuItemNewWithLabel <- function(group = NULL, label, show = TRUE) {
  if (is.null(group)) {
    w <- .RGtkCall("S_gtk_radio_menu_item_new_with_label", group, label,
                   PACKAGE="RGtk2")
    if (show)
      w$show()
    w
  } else gtkRadioMenuItemNewWithLabelFromWidget(group[[1]], label, show)
}

gtkRadioMenuItemNewWithMnemonic <- function(group = NULL, label, show = TRUE) {
  if (is.null(group)) {
    w <- .RGtkCall("S_gtk_radio_menu_item_new_with_mnemonic", group, label,
                   PACKAGE="RGtk2")
    if (show)
      w$show()
    w
  } else gtkRadioMenuItemNewWithMnemonicFromWidget(group[[1]], label, show) 
}

gtkRadioToolButtonNew <-
function(group = NULL, show = TRUE) {
  if (is.null(group)) {
    w <- .RGtkCall("S_gtk_radio_tool_button_new", group, PACKAGE="RGtk2")
    if (show)
      w$show()
    w
  } else gtkRadioToolButtonNewFromWidget(group[[1]], show)
}

gtkRadioToolButtonNewFromStock <- function(group = NULL, stock.id, show = TRUE)
{
  if (is.null(group)) {
    w <- .RGtkCall("S_gtk_radio_tool_button_new_from_stock", group, stock.id,
                   PACKAGE="RGtk2")
    if (show)
      w$show()
    w
  } else gtkRadioToolButtonNewWithStockFromWidget(group[[1]], stock.id, show)
}

# getting child widgets by index
"[[.GtkContainer" <-
function(x, field, where = parent.frame())
{
  if(is.numeric(field))
    return(x$getChildren()[[field]])
  else NextMethod("[[", where = where)
}

# EXPERIMENTAL TREE MODEL ACCESS
# Currently deprecated in favor of custom RGtkDataFrame model
if (FALSE) {
# Loads data into the specified rows (or paths) and columns of a list store
gtkListStoreLoad <-
function(object, data, rows, cols = 0:(length(data)-1), paths = NULL, append=F)
{
        checkPtrType(object, "GtkListStore")
		
		w <- gtkTreeModelLoad(object, data, rows, cols, paths, append, type = "list")
		
		return(invisible(w))
}
# Loads data into the specified rows (or paths) and columns of a tree store
# Here, the rows should be vectors that describe path locations within the tree
gtkTreeStoreLoad <-
function(object, data, rows, cols = 0:(length(data)-1), paths = NULL, append=F)
{
        checkPtrType(object, "GtkTreeStore")
		
		w <- gtkTreeModelLoad(object, data, rows, cols, paths, append, type = "tree") 
		
        return(invisible(w))
}

# The (private) main function for loading data into a tree or list store
# If the column names have not been set for the model, this function sets them
# to the column names of the data frame. If append is TRUE the data is added
# to the end of the model. Otherwise, if there are no paths and the rows are not specified, 
# they are assumed to be the sequence from 0 to the number of rows in the data.
# The columns can be indices or names. Sorting is temporarily disabled for this operation.
gtkTreeModelLoad <-
function(object, data, rows, cols = 0:(length(data)-1), paths = NULL, append=F, type = c("list", "tree"))
{		
		type <- match.arg(type)
		
		col.names <- object$getData("colnames")
		if (is.null(col.names)) {
			col.names <- character(object$getNColumns())
			col.names[cols+1] <- names(data)
			object$setData("colnames", col.names)
		}
		
        data <- lapply(data, as.list)
		if (is.character(cols)) {
			m <- match(cols, col.names)
			cols[!is.na(m)] <- m[!is.na(m)]-1
		}
		cols <- as.integer(cols)
		
		if (length(cols) != length(data))
			stop("The number of specified columns (", length(cols), ")",
				" does not match the number of columns (", length(data), ") in the data")
		
		sort <- object$getSortColumnId()
		sorted <- sort["sort_column_id"] %in% cols
		if (sorted)
			object$setSortColumnId(-2, sort[["order"]]) # -2 means unsorted
		
		if (missing(paths) && (!append || type == "tree")) {
			if (missing(rows) || is.null(rows))
				rows <- 0:(length(data[[1]])-1)
			if (type == "tree")
				rows <- lapply(rows, as.integer)
			else rows <- as.integer(rows)
			w <- .RGtkCall(paste("S_gtk_", type, "_store_load", sep=""), object, data, rows, cols, append)
		} else {
			paths <- lapply(paths, function(path) { checkPtrType(path, "GtkTreePath"); path })
			w <- .RGtkCall(paste("S_gtk_", type, "_store_load_paths", sep=""), object, data, paths, cols, append)
		}
		
		if (sorted)
			object$setSortColumnId(sort[["sort_column_id"]], sort[["order"]]) 
		
        return(invisible(w))
}
# Unloads a tree model. If there are no paths, the rows must must be indices or 
# vectors of indices if accessing multiple levels. Columns may be indices or names.
# If this model has column names, they are the column names of the returned structure.
# This structure is a data frame if 'frame' is TRUE. The attribute "paths" contains
# indices or vectors of indices (in the case of multiple levels) describing the location
# of each row of the data frame in the tree model.
gtkTreeModelUnload <-
function(object, rows = NULL, cols = 0:(object$getNColumns()-1), paths, frame = T)
{
        checkPtrType(object, "GtkTreeModel")
		
		col.names <- colnames(object)
		if (is.character(cols)) {
			if (is.null(col.names))
				stop("Specifying column names is not allowed - this tree model does not have any")
			m <- match(cols, col.names)
			if (any(is.na(m)))
				stop("Invalid column names: ", paste(cols[is.na(m)], collapse=", "))
			cols[!is.na(m)] <- m[!is.na(m)]-1
		}
		cols <- as.integer(cols)
		
		if(missing(paths)) {
			rows <- lapply(rows, as.integer)
			data <- .RGtkCall("S_gtk_tree_model_unload", object, rows, cols)
		} else {
			paths <- lapply(paths, function(path) { checkPtrType(path, "GtkTreePath"); path })
			data <- .RGtkCall("S_gtk_tree_model_unload_paths", object, paths, cols)
		}
		
		if (length(rows) == 0) {
			rows <- data[[2]]
			data <- data[[1]]
		}
		
		if (frame)
			data <- as.data.frame(sapply(data, unlist))
		
		attr(data, "paths") <- rows
		
		if (length(data) > 0 && length(col.names) > 0) {
			names(data) <- col.names[cols+1]
		}
		
		return(data)
}

# we can't do "GtkTreeModel" here because it's an interface, not a class
"[<-.GtkListStore" <- "[<-.GtkTreeStore" <-
function(model, rows = NULL, cols = 0:(length(data)-1), value) {
	model$load(value, rows, cols)
	model
}
"[.GtkListStore" <- "[.GtkTreeStore" <-
function(model, rows = NULL, cols = 0:(model$getNColumns()-1)) {
	model$unload(rows, cols)
}

dimnames.GtkTreeStore <- dimnames.GtkListStore  <- function(x, ...)
{
	list(x$getData("rownames"), x$getData("colnames"))
}
}

# creating a GtkTreeIter from scratch (for implementing new models)

gtkTreeIter <- function(id, stamp)
{
  .RGtkCall("S_gtk_tree_iter", as.integer(id), as.integer(stamp))
}

# setting id's and stamps on GtkTreeIter's (for implementing new models)

gtkTreeIterGetId <- function(iter)
{
  checkPtrType(iter, "GtkTreeIter")
  .RGtkCall("S_gtk_tree_iter_get_id", iter)
}
gtkTreeIterSetId <- function(iter, id)
{
  checkPtrType(iter, "GtkTreeIter")
  id <- as.integer(id)
  .RGtkCall("S_gtk_tree_iter_set_id", iter, id)
}

gtkTreeIterGetStamp <- function(iter)
{
  checkPtrType(iter, "GtkTreeIter")
  .RGtkCall("S_gtk_tree_iter_get_stamp", iter)
}
gtkTreeIterSetStamp <- function(iter, stamp)
{
  checkPtrType(iter, "GtkTreeIter")
  stamp <- as.integer(stamp)
  .RGtkCall("S_gtk_tree_iter_set_stamp", iter, stamp)
}

# aliases for error domains
GTK_ICON_THEME_ERROR <- gtkIconThemeErrorQuark
GTK_FILE_CHOOSER_ERROR <- gtkFileChooserErrorQuark

# NOT IMPLEMENTED FUNCS #

gtkAccelGroupQuery <- 
function(object, accel.key, accel.mods) 
{
	.notimplemented("returns an array of an undocumented structure, so you probably shouldn't use it") 
}

gtkCListSetCompareFunc <- 
function(object, cmp.func) {
	.notimplemented("does not provide user-data for the comparison func, so we can't wrap an R function. Besides, it's deprecated")
}

gtkCTreeSetDragCompareFunc <-
function(object, cmp.func)
{
	.notimplemented("does not provide user data for the comparison func, so we can't wrap an R function. Besides, it's deprecated")
}

gtkItemFactoryCreateMenuEntries <-
function(entries)
{
	.notimplemented("accepts as input an array of an undocumented structure, so you probably shouldn't use it. Besides, it's deprecated")
}

gtkSettingsInstallPropertyParser <-
function(pspec, parser)
{
	.notimplemented("does not have user data for the parser callback. You probably don't need to be defining new property types in GTK style files anyway.")
}

gtkWidgetClassInstallStylePropertyParser <-
function(klass, pspec, parser)
{
	.notimplemented("does not have user data for the parser callback. You probably don't need to be defining new property types in GTK style files anyway.")
}

## version checking ##

boundGTKVersion <- function() {
  as.numeric_version(paste(.RGtkCall("boundGTKVersion"), collapse="."))
}

checkGTK <- function(version) {
  .Deprecated("boundGTKVersion() >= version")
  boundGTKVersion() >= version
}
