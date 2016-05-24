if(!exists('.virtuals')) .virtuals <- new.env()
assign("AtkAction", c("do_action", "get_n_actions", "get_description", "get_name", "get_keybinding", "set_description", "get_localized_name"), .virtuals)
assign("AtkComponent", c("contains", "ref_accessible_at_point", "get_extents", "get_position", "get_size", "grab_focus", "remove_focus_handler", "set_extents", "set_position", "set_size", "get_layer", "get_mdi_zorder", "bounds_changed"), .virtuals)
assign("AtkDocument", c("get_document_type", "get_document"), .virtuals)
assign("AtkEditableText", c("set_run_attributes", "set_text_contents", "insert_text", "copy_text", "cut_text", "delete_text", "paste_text"), .virtuals)
assign("AtkHyperlink", c("get_uri", "get_object", "get_end_index", "get_start_index", "is_valid", "get_n_anchors", "link_state", "is_selected_link", "link_activated"), .virtuals)
assign("AtkHypertext", c("get_link", "get_n_links", "get_link_index", "link_selected"), .virtuals)
assign("AtkImage", c("get_image_position", "get_image_description", "get_image_size", "set_image_description"), .virtuals)
assign("AtkObjectFactory", c("invalidate"), .virtuals)
assign("AtkImplementor", c("ref_accessible"), .virtuals)
assign("AtkObject", c("get_name", "get_description", "get_parent", "get_n_children", "ref_child", "get_index_in_parent", "ref_relation_set", "get_role", "get_layer", "get_mdi_zorder", "ref_state_set", "set_name", "set_description", "set_parent", "set_role", "remove_property_change_handler", "initialize", "children_changed", "focus_event", "state_change", "visible_data_changed", "active_descendant_changed"), .virtuals)
assign("AtkSelection", c("add_selection", "clear_selection", "ref_selection", "get_selection_count", "is_child_selected", "remove_selection", "select_all_selection", "selection_changed"), .virtuals)
assign("AtkStreamableContent", c("get_n_mime_types", "get_mime_type"), .virtuals)
assign("AtkTable", c("ref_at", "get_index_at", "get_column_at_index", "get_row_at_index", "get_n_columns", "get_n_rows", "get_column_extent_at", "get_row_extent_at", "get_caption", "get_column_description", "get_column_header", "get_row_description", "get_row_header", "get_summary", "set_caption", "set_column_description", "set_column_header", "set_row_description", "set_row_header", "set_summary", "get_selected_columns", "get_selected_rows", "is_column_selected", "is_row_selected", "is_selected", "add_row_selection", "remove_row_selection", "add_column_selection", "remove_column_selection", "row_inserted", "column_inserted", "row_deleted", "column_deleted", "row_reordered", "column_reordered", "model_changed"), .virtuals)
assign("AtkText", c("get_text", "get_text_after_offset", "get_text_at_offset", "get_character_at_offset", "get_text_before_offset", "get_caret_offset", "get_run_attributes", "get_default_attributes", "get_character_extents", "get_character_count", "get_offset_at_point", "get_n_selections", "get_selection", "add_selection", "remove_selection", "set_selection", "set_caret_offset", "text_changed", "text_caret_moved", "text_selection_changed", "text_attributes_changed", "get_range_extents", "get_bounded_ranges"), .virtuals)
assign("AtkValue", c("get_current_value", "get_maximum_value", "get_minimum_value", "set_current_value", "get_minimum_increment"), .virtuals)


atkActionIfaceDoAction <-
function(object.class, object, i)
{
  checkPtrType(object.class, "AtkActionIface")
  checkPtrType(object, "AtkAction")
  i <- as.integer(i)

  w <- .RGtkCall("S_atk_action_iface_do_action", object.class, object, i, PACKAGE = "RGtk2")

  return(w)
}

atkActionIfaceGetNActions <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkActionIface")
  checkPtrType(object, "AtkAction")

  w <- .RGtkCall("S_atk_action_iface_get_n_actions", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkActionIfaceGetDescription <-
function(object.class, object, i)
{
  checkPtrType(object.class, "AtkActionIface")
  checkPtrType(object, "AtkAction")
  i <- as.integer(i)

  w <- .RGtkCall("S_atk_action_iface_get_description", object.class, object, i, PACKAGE = "RGtk2")

  return(w)
}

atkActionIfaceGetName <-
function(object.class, object, i)
{
  checkPtrType(object.class, "AtkActionIface")
  checkPtrType(object, "AtkAction")
  i <- as.integer(i)

  w <- .RGtkCall("S_atk_action_iface_get_name", object.class, object, i, PACKAGE = "RGtk2")

  return(w)
}

atkActionIfaceGetKeybinding <-
function(object.class, object, i)
{
  checkPtrType(object.class, "AtkActionIface")
  checkPtrType(object, "AtkAction")
  i <- as.integer(i)

  w <- .RGtkCall("S_atk_action_iface_get_keybinding", object.class, object, i, PACKAGE = "RGtk2")

  return(w)
}

atkActionIfaceSetDescription <-
function(object.class, object, i, desc)
{
  checkPtrType(object.class, "AtkActionIface")
  checkPtrType(object, "AtkAction")
  i <- as.integer(i)
  desc <- as.character(desc)

  w <- .RGtkCall("S_atk_action_iface_set_description", object.class, object, i, desc, PACKAGE = "RGtk2")

  return(w)
}

atkActionIfaceGetLocalizedName <-
function(object.class, object, i)
{
  checkPtrType(object.class, "AtkActionIface")
  checkPtrType(object, "AtkAction")
  i <- as.integer(i)

  w <- .RGtkCall("S_atk_action_iface_get_localized_name", object.class, object, i, PACKAGE = "RGtk2")

  return(w)
}

atkComponentIfaceContains <-
function(object.class, object, x, y, coord.type)
{
  checkPtrType(object.class, "AtkComponentIface")
  checkPtrType(object, "AtkComponent")
  x <- as.integer(x)
  y <- as.integer(y)
  

  w <- .RGtkCall("S_atk_component_iface_contains", object.class, object, x, y, coord.type, PACKAGE = "RGtk2")

  return(w)
}

atkComponentIfaceRefAccessibleAtPoint <-
function(object.class, object, x, y, coord.type)
{
  checkPtrType(object.class, "AtkComponentIface")
  checkPtrType(object, "AtkComponent")
  x <- as.integer(x)
  y <- as.integer(y)
  

  w <- .RGtkCall("S_atk_component_iface_ref_accessible_at_point", object.class, object, x, y, coord.type, PACKAGE = "RGtk2")

  return(w)
}

atkComponentIfaceGetExtents <-
function(object.class, object, coord.type)
{
  checkPtrType(object.class, "AtkComponentIface")
  checkPtrType(object, "AtkComponent")
  

  w <- .RGtkCall("S_atk_component_iface_get_extents", object.class, object, coord.type, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkComponentIfaceGetPosition <-
function(object.class, object, coord.type)
{
  checkPtrType(object.class, "AtkComponentIface")
  checkPtrType(object, "AtkComponent")
  

  w <- .RGtkCall("S_atk_component_iface_get_position", object.class, object, coord.type, PACKAGE = "RGtk2")

  return(w)
}

atkComponentIfaceGetSize <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkComponentIface")
  checkPtrType(object, "AtkComponent")

  w <- .RGtkCall("S_atk_component_iface_get_size", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkComponentIfaceGrabFocus <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkComponentIface")
  checkPtrType(object, "AtkComponent")

  w <- .RGtkCall("S_atk_component_iface_grab_focus", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkComponentIfaceRemoveFocusHandler <-
function(object.class, object, handler.id)
{
  checkPtrType(object.class, "AtkComponentIface")
  checkPtrType(object, "AtkComponent")
  handler.id <- as.numeric(handler.id)

  w <- .RGtkCall("S_atk_component_iface_remove_focus_handler", object.class, object, handler.id, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkComponentIfaceSetExtents <-
function(object.class, object, x, y, width, height, coord.type)
{
  checkPtrType(object.class, "AtkComponentIface")
  checkPtrType(object, "AtkComponent")
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)
  

  w <- .RGtkCall("S_atk_component_iface_set_extents", object.class, object, x, y, width, height, coord.type, PACKAGE = "RGtk2")

  return(w)
}

atkComponentIfaceSetPosition <-
function(object.class, object, x, y, coord.type)
{
  checkPtrType(object.class, "AtkComponentIface")
  checkPtrType(object, "AtkComponent")
  x <- as.integer(x)
  y <- as.integer(y)
  

  w <- .RGtkCall("S_atk_component_iface_set_position", object.class, object, x, y, coord.type, PACKAGE = "RGtk2")

  return(w)
}

atkComponentIfaceSetSize <-
function(object.class, object, width, height)
{
  checkPtrType(object.class, "AtkComponentIface")
  checkPtrType(object, "AtkComponent")
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_atk_component_iface_set_size", object.class, object, width, height, PACKAGE = "RGtk2")

  return(w)
}

atkComponentIfaceGetLayer <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkComponentIface")
  checkPtrType(object, "AtkComponent")

  w <- .RGtkCall("S_atk_component_iface_get_layer", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkComponentIfaceGetMdiZorder <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkComponentIface")
  checkPtrType(object, "AtkComponent")

  w <- .RGtkCall("S_atk_component_iface_get_mdi_zorder", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkComponentIfaceBoundsChanged <-
function(object.class, object, bounds)
{
  checkPtrType(object.class, "AtkComponentIface")
  checkPtrType(object, "AtkComponent")
  bounds <- as.AtkRectangle(bounds)

  w <- .RGtkCall("S_atk_component_iface_bounds_changed", object.class, object, bounds, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkDocumentIfaceGetDocumentType <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkDocumentIface")
  checkPtrType(object, "AtkDocument")

  w <- .RGtkCall("S_atk_document_iface_get_document_type", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkDocumentIfaceGetDocument <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkDocumentIface")
  checkPtrType(object, "AtkDocument")

  w <- .RGtkCall("S_atk_document_iface_get_document", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkEditableTextIfaceSetRunAttributes <-
function(object.class, object, attrib.set, start.offset, end.offset)
{
  checkPtrType(object.class, "AtkEditableTextIface")
  checkPtrType(object, "AtkEditableText")
  attrib.set <- as.AtkAttributeSet(attrib.set)
  start.offset <- as.integer(start.offset)
  end.offset <- as.integer(end.offset)

  w <- .RGtkCall("S_atk_editable_text_iface_set_run_attributes", object.class, object, attrib.set, start.offset, end.offset, PACKAGE = "RGtk2")

  return(w)
}

atkEditableTextIfaceSetTextContents <-
function(object.class, object, string)
{
  checkPtrType(object.class, "AtkEditableTextIface")
  checkPtrType(object, "AtkEditableText")
  string <- as.character(string)

  w <- .RGtkCall("S_atk_editable_text_iface_set_text_contents", object.class, object, string, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkEditableTextIfaceInsertText <-
function(object.class, object, string, position)
{
  checkPtrType(object.class, "AtkEditableTextIface")
  checkPtrType(object, "AtkEditableText")
  string <- as.character(string)
  position <- as.list(as.integer(position))

  w <- .RGtkCall("S_atk_editable_text_iface_insert_text", object.class, object, string, position, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkEditableTextIfaceCopyText <-
function(object.class, object, start.pos, end.pos)
{
  checkPtrType(object.class, "AtkEditableTextIface")
  checkPtrType(object, "AtkEditableText")
  start.pos <- as.integer(start.pos)
  end.pos <- as.integer(end.pos)

  w <- .RGtkCall("S_atk_editable_text_iface_copy_text", object.class, object, start.pos, end.pos, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkEditableTextIfaceCutText <-
function(object.class, object, start.pos, end.pos)
{
  checkPtrType(object.class, "AtkEditableTextIface")
  checkPtrType(object, "AtkEditableText")
  start.pos <- as.integer(start.pos)
  end.pos <- as.integer(end.pos)

  w <- .RGtkCall("S_atk_editable_text_iface_cut_text", object.class, object, start.pos, end.pos, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkEditableTextIfaceDeleteText <-
function(object.class, object, start.pos, end.pos)
{
  checkPtrType(object.class, "AtkEditableTextIface")
  checkPtrType(object, "AtkEditableText")
  start.pos <- as.integer(start.pos)
  end.pos <- as.integer(end.pos)

  w <- .RGtkCall("S_atk_editable_text_iface_delete_text", object.class, object, start.pos, end.pos, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkEditableTextIfacePasteText <-
function(object.class, object, position)
{
  checkPtrType(object.class, "AtkEditableTextIface")
  checkPtrType(object, "AtkEditableText")
  position <- as.integer(position)

  w <- .RGtkCall("S_atk_editable_text_iface_paste_text", object.class, object, position, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkHyperlinkClassGetUri <-
function(object.class, object, i)
{
  checkPtrType(object.class, "AtkHyperlinkClass")
  checkPtrType(object, "AtkHyperlink")
  i <- as.integer(i)

  w <- .RGtkCall("S_atk_hyperlink_class_get_uri", object.class, object, i, PACKAGE = "RGtk2")

  return(w)
}

atkHyperlinkClassGetObject <-
function(object.class, object, i)
{
  checkPtrType(object.class, "AtkHyperlinkClass")
  checkPtrType(object, "AtkHyperlink")
  i <- as.integer(i)

  w <- .RGtkCall("S_atk_hyperlink_class_get_object", object.class, object, i, PACKAGE = "RGtk2")

  return(w)
}

atkHyperlinkClassGetEndIndex <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkHyperlinkClass")
  checkPtrType(object, "AtkHyperlink")

  w <- .RGtkCall("S_atk_hyperlink_class_get_end_index", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkHyperlinkClassGetStartIndex <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkHyperlinkClass")
  checkPtrType(object, "AtkHyperlink")

  w <- .RGtkCall("S_atk_hyperlink_class_get_start_index", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkHyperlinkClassIsValid <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkHyperlinkClass")
  checkPtrType(object, "AtkHyperlink")

  w <- .RGtkCall("S_atk_hyperlink_class_is_valid", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkHyperlinkClassGetNAnchors <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkHyperlinkClass")
  checkPtrType(object, "AtkHyperlink")

  w <- .RGtkCall("S_atk_hyperlink_class_get_n_anchors", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkHyperlinkClassLinkState <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkHyperlinkClass")
  checkPtrType(object, "AtkHyperlink")

  w <- .RGtkCall("S_atk_hyperlink_class_link_state", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkHyperlinkClassIsSelectedLink <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkHyperlinkClass")
  checkPtrType(object, "AtkHyperlink")

  w <- .RGtkCall("S_atk_hyperlink_class_is_selected_link", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkHyperlinkClassLinkActivated <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkHyperlinkClass")
  checkPtrType(object, "AtkHyperlink")

  w <- .RGtkCall("S_atk_hyperlink_class_link_activated", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkHypertextIfaceGetLink <-
function(object.class, object, link.index)
{
  checkPtrType(object.class, "AtkHypertextIface")
  checkPtrType(object, "AtkHypertext")
  link.index <- as.integer(link.index)

  w <- .RGtkCall("S_atk_hypertext_iface_get_link", object.class, object, link.index, PACKAGE = "RGtk2")

  return(w)
}

atkHypertextIfaceGetNLinks <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkHypertextIface")
  checkPtrType(object, "AtkHypertext")

  w <- .RGtkCall("S_atk_hypertext_iface_get_n_links", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkHypertextIfaceGetLinkIndex <-
function(object.class, object, char.index)
{
  checkPtrType(object.class, "AtkHypertextIface")
  checkPtrType(object, "AtkHypertext")
  char.index <- as.integer(char.index)

  w <- .RGtkCall("S_atk_hypertext_iface_get_link_index", object.class, object, char.index, PACKAGE = "RGtk2")

  return(w)
}

atkHypertextIfaceLinkSelected <-
function(object.class, object, link.index)
{
  checkPtrType(object.class, "AtkHypertextIface")
  checkPtrType(object, "AtkHypertext")
  link.index <- as.integer(link.index)

  w <- .RGtkCall("S_atk_hypertext_iface_link_selected", object.class, object, link.index, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkImageIfaceGetImagePosition <-
function(object.class, object, coord.type)
{
  checkPtrType(object.class, "AtkImageIface")
  checkPtrType(object, "AtkImage")
  

  w <- .RGtkCall("S_atk_image_iface_get_image_position", object.class, object, coord.type, PACKAGE = "RGtk2")

  return(w)
}

atkImageIfaceGetImageDescription <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkImageIface")
  checkPtrType(object, "AtkImage")

  w <- .RGtkCall("S_atk_image_iface_get_image_description", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkImageIfaceGetImageSize <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkImageIface")
  checkPtrType(object, "AtkImage")

  w <- .RGtkCall("S_atk_image_iface_get_image_size", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkImageIfaceSetImageDescription <-
function(object.class, object, description)
{
  checkPtrType(object.class, "AtkImageIface")
  checkPtrType(object, "AtkImage")
  description <- as.character(description)

  w <- .RGtkCall("S_atk_image_iface_set_image_description", object.class, object, description, PACKAGE = "RGtk2")

  return(w)
}

atkObjectFactoryClassInvalidate <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkObjectFactoryClass")
  checkPtrType(object, "AtkObjectFactory")

  w <- .RGtkCall("S_atk_object_factory_class_invalidate", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkImplementorIfaceRefAccessible <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkImplementorIface")
  checkPtrType(object, "AtkImplementor")

  w <- .RGtkCall("S_atk_implementor_iface_ref_accessible", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkObjectClassGetName <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkObjectClass")
  checkPtrType(object, "AtkObject")

  w <- .RGtkCall("S_atk_object_class_get_name", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkObjectClassGetDescription <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkObjectClass")
  checkPtrType(object, "AtkObject")

  w <- .RGtkCall("S_atk_object_class_get_description", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkObjectClassGetParent <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkObjectClass")
  checkPtrType(object, "AtkObject")

  w <- .RGtkCall("S_atk_object_class_get_parent", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkObjectClassGetNChildren <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkObjectClass")
  checkPtrType(object, "AtkObject")

  w <- .RGtkCall("S_atk_object_class_get_n_children", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkObjectClassRefChild <-
function(object.class, object, i)
{
  checkPtrType(object.class, "AtkObjectClass")
  checkPtrType(object, "AtkObject")
  i <- as.integer(i)

  w <- .RGtkCall("S_atk_object_class_ref_child", object.class, object, i, PACKAGE = "RGtk2")

  return(w)
}

atkObjectClassGetIndexInParent <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkObjectClass")
  checkPtrType(object, "AtkObject")

  w <- .RGtkCall("S_atk_object_class_get_index_in_parent", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkObjectClassRefRelationSet <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkObjectClass")
  checkPtrType(object, "AtkObject")

  w <- .RGtkCall("S_atk_object_class_ref_relation_set", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkObjectClassGetRole <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkObjectClass")
  checkPtrType(object, "AtkObject")

  w <- .RGtkCall("S_atk_object_class_get_role", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkObjectClassGetLayer <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkObjectClass")
  checkPtrType(object, "AtkObject")

  w <- .RGtkCall("S_atk_object_class_get_layer", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkObjectClassGetMdiZorder <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkObjectClass")
  checkPtrType(object, "AtkObject")

  w <- .RGtkCall("S_atk_object_class_get_mdi_zorder", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkObjectClassRefStateSet <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkObjectClass")
  checkPtrType(object, "AtkObject")

  w <- .RGtkCall("S_atk_object_class_ref_state_set", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkObjectClassSetName <-
function(object.class, object, name)
{
  checkPtrType(object.class, "AtkObjectClass")
  checkPtrType(object, "AtkObject")
  name <- as.character(name)

  w <- .RGtkCall("S_atk_object_class_set_name", object.class, object, name, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkObjectClassSetDescription <-
function(object.class, object, description)
{
  checkPtrType(object.class, "AtkObjectClass")
  checkPtrType(object, "AtkObject")
  description <- as.character(description)

  w <- .RGtkCall("S_atk_object_class_set_description", object.class, object, description, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkObjectClassSetParent <-
function(object.class, object, parent)
{
  checkPtrType(object.class, "AtkObjectClass")
  checkPtrType(object, "AtkObject")
  checkPtrType(parent, "AtkObject")

  w <- .RGtkCall("S_atk_object_class_set_parent", object.class, object, parent, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkObjectClassSetRole <-
function(object.class, object, role)
{
  checkPtrType(object.class, "AtkObjectClass")
  checkPtrType(object, "AtkObject")
  

  w <- .RGtkCall("S_atk_object_class_set_role", object.class, object, role, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkObjectClassRemovePropertyChangeHandler <-
function(object.class, object, handler.id)
{
  checkPtrType(object.class, "AtkObjectClass")
  checkPtrType(object, "AtkObject")
  handler.id <- as.numeric(handler.id)

  w <- .RGtkCall("S_atk_object_class_remove_property_change_handler", object.class, object, handler.id, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkObjectClassInitialize <-
function(object.class, object, data)
{
  checkPtrType(object.class, "AtkObjectClass")
  checkPtrType(object, "AtkObject")
  

  w <- .RGtkCall("S_atk_object_class_initialize", object.class, object, data, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkObjectClassChildrenChanged <-
function(object.class, object, change.index, changed.child)
{
  checkPtrType(object.class, "AtkObjectClass")
  checkPtrType(object, "AtkObject")
  change.index <- as.numeric(change.index)
  checkPtrType(changed.child, "AtkObject")

  w <- .RGtkCall("S_atk_object_class_children_changed", object.class, object, change.index, changed.child, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkObjectClassFocusEvent <-
function(object.class, object, focus.in)
{
  checkPtrType(object.class, "AtkObjectClass")
  checkPtrType(object, "AtkObject")
  focus.in <- as.logical(focus.in)

  w <- .RGtkCall("S_atk_object_class_focus_event", object.class, object, focus.in, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkObjectClassStateChange <-
function(object.class, object, name, state.set)
{
  checkPtrType(object.class, "AtkObjectClass")
  checkPtrType(object, "AtkObject")
  name <- as.character(name)
  state.set <- as.logical(state.set)

  w <- .RGtkCall("S_atk_object_class_state_change", object.class, object, name, state.set, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkObjectClassVisibleDataChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkObjectClass")
  checkPtrType(object, "AtkObject")

  w <- .RGtkCall("S_atk_object_class_visible_data_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkObjectClassActiveDescendantChanged <-
function(object.class, object, child)
{
  checkPtrType(object.class, "AtkObjectClass")
  checkPtrType(object, "AtkObject")
  checkPtrType(child, "AtkObject")

  w <- .RGtkCall("S_atk_object_class_active_descendant_changed", object.class, object, child, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkSelectionIfaceAddSelection <-
function(object.class, object, i)
{
  checkPtrType(object.class, "AtkSelectionIface")
  checkPtrType(object, "AtkSelection")
  i <- as.integer(i)

  w <- .RGtkCall("S_atk_selection_iface_add_selection", object.class, object, i, PACKAGE = "RGtk2")

  return(w)
}

atkSelectionIfaceClearSelection <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkSelectionIface")
  checkPtrType(object, "AtkSelection")

  w <- .RGtkCall("S_atk_selection_iface_clear_selection", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkSelectionIfaceRefSelection <-
function(object.class, object, i)
{
  checkPtrType(object.class, "AtkSelectionIface")
  checkPtrType(object, "AtkSelection")
  i <- as.integer(i)

  w <- .RGtkCall("S_atk_selection_iface_ref_selection", object.class, object, i, PACKAGE = "RGtk2")

  return(w)
}

atkSelectionIfaceGetSelectionCount <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkSelectionIface")
  checkPtrType(object, "AtkSelection")

  w <- .RGtkCall("S_atk_selection_iface_get_selection_count", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkSelectionIfaceIsChildSelected <-
function(object.class, object, i)
{
  checkPtrType(object.class, "AtkSelectionIface")
  checkPtrType(object, "AtkSelection")
  i <- as.integer(i)

  w <- .RGtkCall("S_atk_selection_iface_is_child_selected", object.class, object, i, PACKAGE = "RGtk2")

  return(w)
}

atkSelectionIfaceRemoveSelection <-
function(object.class, object, i)
{
  checkPtrType(object.class, "AtkSelectionIface")
  checkPtrType(object, "AtkSelection")
  i <- as.integer(i)

  w <- .RGtkCall("S_atk_selection_iface_remove_selection", object.class, object, i, PACKAGE = "RGtk2")

  return(w)
}

atkSelectionIfaceSelectAllSelection <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkSelectionIface")
  checkPtrType(object, "AtkSelection")

  w <- .RGtkCall("S_atk_selection_iface_select_all_selection", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkSelectionIfaceSelectionChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkSelectionIface")
  checkPtrType(object, "AtkSelection")

  w <- .RGtkCall("S_atk_selection_iface_selection_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkStreamableContentIfaceGetNMimeTypes <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkStreamableContentIface")
  checkPtrType(object, "AtkStreamableContent")

  w <- .RGtkCall("S_atk_streamable_content_iface_get_n_mime_types", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkStreamableContentIfaceGetMimeType <-
function(object.class, object, i)
{
  checkPtrType(object.class, "AtkStreamableContentIface")
  checkPtrType(object, "AtkStreamableContent")
  i <- as.integer(i)

  w <- .RGtkCall("S_atk_streamable_content_iface_get_mime_type", object.class, object, i, PACKAGE = "RGtk2")

  return(w)
}

atkTableIfaceRefAt <-
function(object.class, object, row, column)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")
  row <- as.integer(row)
  column <- as.integer(column)

  w <- .RGtkCall("S_atk_table_iface_ref_at", object.class, object, row, column, PACKAGE = "RGtk2")

  return(w)
}

atkTableIfaceGetIndexAt <-
function(object.class, object, row, column)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")
  row <- as.integer(row)
  column <- as.integer(column)

  w <- .RGtkCall("S_atk_table_iface_get_index_at", object.class, object, row, column, PACKAGE = "RGtk2")

  return(w)
}

atkTableIfaceGetColumnAtIndex <-
function(object.class, object, index)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")
  index <- as.integer(index)

  w <- .RGtkCall("S_atk_table_iface_get_column_at_index", object.class, object, index, PACKAGE = "RGtk2")

  return(w)
}

atkTableIfaceGetRowAtIndex <-
function(object.class, object, index)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")
  index <- as.integer(index)

  w <- .RGtkCall("S_atk_table_iface_get_row_at_index", object.class, object, index, PACKAGE = "RGtk2")

  return(w)
}

atkTableIfaceGetNColumns <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")

  w <- .RGtkCall("S_atk_table_iface_get_n_columns", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkTableIfaceGetNRows <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")

  w <- .RGtkCall("S_atk_table_iface_get_n_rows", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkTableIfaceGetColumnExtentAt <-
function(object.class, object, row, column)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")
  row <- as.integer(row)
  column <- as.integer(column)

  w <- .RGtkCall("S_atk_table_iface_get_column_extent_at", object.class, object, row, column, PACKAGE = "RGtk2")

  return(w)
}

atkTableIfaceGetRowExtentAt <-
function(object.class, object, row, column)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")
  row <- as.integer(row)
  column <- as.integer(column)

  w <- .RGtkCall("S_atk_table_iface_get_row_extent_at", object.class, object, row, column, PACKAGE = "RGtk2")

  return(w)
}

atkTableIfaceGetCaption <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")

  w <- .RGtkCall("S_atk_table_iface_get_caption", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkTableIfaceGetColumnDescription <-
function(object.class, object, column)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")
  column <- as.integer(column)

  w <- .RGtkCall("S_atk_table_iface_get_column_description", object.class, object, column, PACKAGE = "RGtk2")

  return(w)
}

atkTableIfaceGetColumnHeader <-
function(object.class, object, column)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")
  column <- as.integer(column)

  w <- .RGtkCall("S_atk_table_iface_get_column_header", object.class, object, column, PACKAGE = "RGtk2")

  return(w)
}

atkTableIfaceGetRowDescription <-
function(object.class, object, row)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")
  row <- as.integer(row)

  w <- .RGtkCall("S_atk_table_iface_get_row_description", object.class, object, row, PACKAGE = "RGtk2")

  return(w)
}

atkTableIfaceGetRowHeader <-
function(object.class, object, row)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")
  row <- as.integer(row)

  w <- .RGtkCall("S_atk_table_iface_get_row_header", object.class, object, row, PACKAGE = "RGtk2")

  return(w)
}

atkTableIfaceGetSummary <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")

  w <- .RGtkCall("S_atk_table_iface_get_summary", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkTableIfaceSetCaption <-
function(object.class, object, caption)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")
  checkPtrType(caption, "AtkObject")

  w <- .RGtkCall("S_atk_table_iface_set_caption", object.class, object, caption, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkTableIfaceSetColumnDescription <-
function(object.class, object, column, description)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")
  column <- as.integer(column)
  description <- as.character(description)

  w <- .RGtkCall("S_atk_table_iface_set_column_description", object.class, object, column, description, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkTableIfaceSetColumnHeader <-
function(object.class, object, column, header)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")
  column <- as.integer(column)
  checkPtrType(header, "AtkObject")

  w <- .RGtkCall("S_atk_table_iface_set_column_header", object.class, object, column, header, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkTableIfaceSetRowDescription <-
function(object.class, object, row, description)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")
  row <- as.integer(row)
  description <- as.character(description)

  w <- .RGtkCall("S_atk_table_iface_set_row_description", object.class, object, row, description, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkTableIfaceSetRowHeader <-
function(object.class, object, row, header)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")
  row <- as.integer(row)
  checkPtrType(header, "AtkObject")

  w <- .RGtkCall("S_atk_table_iface_set_row_header", object.class, object, row, header, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkTableIfaceSetSummary <-
function(object.class, object, accessible)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")
  checkPtrType(accessible, "AtkObject")

  w <- .RGtkCall("S_atk_table_iface_set_summary", object.class, object, accessible, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkTableIfaceGetSelectedColumns <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")

  w <- .RGtkCall("S_atk_table_iface_get_selected_columns", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkTableIfaceGetSelectedRows <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")

  w <- .RGtkCall("S_atk_table_iface_get_selected_rows", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkTableIfaceIsColumnSelected <-
function(object.class, object, column)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")
  column <- as.integer(column)

  w <- .RGtkCall("S_atk_table_iface_is_column_selected", object.class, object, column, PACKAGE = "RGtk2")

  return(w)
}

atkTableIfaceIsRowSelected <-
function(object.class, object, row)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")
  row <- as.integer(row)

  w <- .RGtkCall("S_atk_table_iface_is_row_selected", object.class, object, row, PACKAGE = "RGtk2")

  return(w)
}

atkTableIfaceIsSelected <-
function(object.class, object, row, column)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")
  row <- as.integer(row)
  column <- as.integer(column)

  w <- .RGtkCall("S_atk_table_iface_is_selected", object.class, object, row, column, PACKAGE = "RGtk2")

  return(w)
}

atkTableIfaceAddRowSelection <-
function(object.class, object, row)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")
  row <- as.integer(row)

  w <- .RGtkCall("S_atk_table_iface_add_row_selection", object.class, object, row, PACKAGE = "RGtk2")

  return(w)
}

atkTableIfaceRemoveRowSelection <-
function(object.class, object, row)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")
  row <- as.integer(row)

  w <- .RGtkCall("S_atk_table_iface_remove_row_selection", object.class, object, row, PACKAGE = "RGtk2")

  return(w)
}

atkTableIfaceAddColumnSelection <-
function(object.class, object, column)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")
  column <- as.integer(column)

  w <- .RGtkCall("S_atk_table_iface_add_column_selection", object.class, object, column, PACKAGE = "RGtk2")

  return(w)
}

atkTableIfaceRemoveColumnSelection <-
function(object.class, object, column)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")
  column <- as.integer(column)

  w <- .RGtkCall("S_atk_table_iface_remove_column_selection", object.class, object, column, PACKAGE = "RGtk2")

  return(w)
}

atkTableIfaceRowInserted <-
function(object.class, object, row, num.inserted)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")
  row <- as.integer(row)
  num.inserted <- as.integer(num.inserted)

  w <- .RGtkCall("S_atk_table_iface_row_inserted", object.class, object, row, num.inserted, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkTableIfaceColumnInserted <-
function(object.class, object, column, num.inserted)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")
  column <- as.integer(column)
  num.inserted <- as.integer(num.inserted)

  w <- .RGtkCall("S_atk_table_iface_column_inserted", object.class, object, column, num.inserted, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkTableIfaceRowDeleted <-
function(object.class, object, row, num.deleted)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")
  row <- as.integer(row)
  num.deleted <- as.integer(num.deleted)

  w <- .RGtkCall("S_atk_table_iface_row_deleted", object.class, object, row, num.deleted, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkTableIfaceColumnDeleted <-
function(object.class, object, column, num.deleted)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")
  column <- as.integer(column)
  num.deleted <- as.integer(num.deleted)

  w <- .RGtkCall("S_atk_table_iface_column_deleted", object.class, object, column, num.deleted, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkTableIfaceRowReordered <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")

  w <- .RGtkCall("S_atk_table_iface_row_reordered", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkTableIfaceColumnReordered <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")

  w <- .RGtkCall("S_atk_table_iface_column_reordered", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkTableIfaceModelChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkTableIface")
  checkPtrType(object, "AtkTable")

  w <- .RGtkCall("S_atk_table_iface_model_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkTextIfaceGetText <-
function(object.class, object, start.offset, end.offset)
{
  checkPtrType(object.class, "AtkTextIface")
  checkPtrType(object, "AtkText")
  start.offset <- as.integer(start.offset)
  end.offset <- as.integer(end.offset)

  w <- .RGtkCall("S_atk_text_iface_get_text", object.class, object, start.offset, end.offset, PACKAGE = "RGtk2")

  return(w)
}

atkTextIfaceGetTextAfterOffset <-
function(object.class, object, offset, boundary.type)
{
  checkPtrType(object.class, "AtkTextIface")
  checkPtrType(object, "AtkText")
  offset <- as.integer(offset)
  

  w <- .RGtkCall("S_atk_text_iface_get_text_after_offset", object.class, object, offset, boundary.type, PACKAGE = "RGtk2")

  return(w)
}

atkTextIfaceGetTextAtOffset <-
function(object.class, object, offset, boundary.type)
{
  checkPtrType(object.class, "AtkTextIface")
  checkPtrType(object, "AtkText")
  offset <- as.integer(offset)
  

  w <- .RGtkCall("S_atk_text_iface_get_text_at_offset", object.class, object, offset, boundary.type, PACKAGE = "RGtk2")

  return(w)
}

atkTextIfaceGetCharacterAtOffset <-
function(object.class, object, offset)
{
  checkPtrType(object.class, "AtkTextIface")
  checkPtrType(object, "AtkText")
  offset <- as.integer(offset)

  w <- .RGtkCall("S_atk_text_iface_get_character_at_offset", object.class, object, offset, PACKAGE = "RGtk2")

  return(w)
}

atkTextIfaceGetTextBeforeOffset <-
function(object.class, object, offset, boundary.type)
{
  checkPtrType(object.class, "AtkTextIface")
  checkPtrType(object, "AtkText")
  offset <- as.integer(offset)
  

  w <- .RGtkCall("S_atk_text_iface_get_text_before_offset", object.class, object, offset, boundary.type, PACKAGE = "RGtk2")

  return(w)
}

atkTextIfaceGetCaretOffset <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkTextIface")
  checkPtrType(object, "AtkText")

  w <- .RGtkCall("S_atk_text_iface_get_caret_offset", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkTextIfaceGetRunAttributes <-
function(object.class, object, offset)
{
  checkPtrType(object.class, "AtkTextIface")
  checkPtrType(object, "AtkText")
  offset <- as.integer(offset)

  w <- .RGtkCall("S_atk_text_iface_get_run_attributes", object.class, object, offset, PACKAGE = "RGtk2")

  return(w)
}

atkTextIfaceGetDefaultAttributes <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkTextIface")
  checkPtrType(object, "AtkText")

  w <- .RGtkCall("S_atk_text_iface_get_default_attributes", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkTextIfaceGetCharacterExtents <-
function(object.class, object, offset, coords)
{
  checkPtrType(object.class, "AtkTextIface")
  checkPtrType(object, "AtkText")
  offset <- as.integer(offset)
  

  w <- .RGtkCall("S_atk_text_iface_get_character_extents", object.class, object, offset, coords, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkTextIfaceGetCharacterCount <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkTextIface")
  checkPtrType(object, "AtkText")

  w <- .RGtkCall("S_atk_text_iface_get_character_count", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkTextIfaceGetOffsetAtPoint <-
function(object.class, object, x, y, coords)
{
  checkPtrType(object.class, "AtkTextIface")
  checkPtrType(object, "AtkText")
  x <- as.integer(x)
  y <- as.integer(y)
  

  w <- .RGtkCall("S_atk_text_iface_get_offset_at_point", object.class, object, x, y, coords, PACKAGE = "RGtk2")

  return(w)
}

atkTextIfaceGetNSelections <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkTextIface")
  checkPtrType(object, "AtkText")

  w <- .RGtkCall("S_atk_text_iface_get_n_selections", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

atkTextIfaceGetSelection <-
function(object.class, object, selection.num)
{
  checkPtrType(object.class, "AtkTextIface")
  checkPtrType(object, "AtkText")
  selection.num <- as.integer(selection.num)

  w <- .RGtkCall("S_atk_text_iface_get_selection", object.class, object, selection.num, PACKAGE = "RGtk2")

  return(w)
}

atkTextIfaceAddSelection <-
function(object.class, object, start.offset, end.offset)
{
  checkPtrType(object.class, "AtkTextIface")
  checkPtrType(object, "AtkText")
  start.offset <- as.integer(start.offset)
  end.offset <- as.integer(end.offset)

  w <- .RGtkCall("S_atk_text_iface_add_selection", object.class, object, start.offset, end.offset, PACKAGE = "RGtk2")

  return(w)
}

atkTextIfaceRemoveSelection <-
function(object.class, object, selection.num)
{
  checkPtrType(object.class, "AtkTextIface")
  checkPtrType(object, "AtkText")
  selection.num <- as.integer(selection.num)

  w <- .RGtkCall("S_atk_text_iface_remove_selection", object.class, object, selection.num, PACKAGE = "RGtk2")

  return(w)
}

atkTextIfaceSetSelection <-
function(object.class, object, selection.num, start.offset, end.offset)
{
  checkPtrType(object.class, "AtkTextIface")
  checkPtrType(object, "AtkText")
  selection.num <- as.integer(selection.num)
  start.offset <- as.integer(start.offset)
  end.offset <- as.integer(end.offset)

  w <- .RGtkCall("S_atk_text_iface_set_selection", object.class, object, selection.num, start.offset, end.offset, PACKAGE = "RGtk2")

  return(w)
}

atkTextIfaceSetCaretOffset <-
function(object.class, object, offset)
{
  checkPtrType(object.class, "AtkTextIface")
  checkPtrType(object, "AtkText")
  offset <- as.integer(offset)

  w <- .RGtkCall("S_atk_text_iface_set_caret_offset", object.class, object, offset, PACKAGE = "RGtk2")

  return(w)
}

atkTextIfaceTextChanged <-
function(object.class, object, position, length)
{
  checkPtrType(object.class, "AtkTextIface")
  checkPtrType(object, "AtkText")
  position <- as.integer(position)
  length <- as.integer(length)

  w <- .RGtkCall("S_atk_text_iface_text_changed", object.class, object, position, length, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkTextIfaceTextCaretMoved <-
function(object.class, object, location)
{
  checkPtrType(object.class, "AtkTextIface")
  checkPtrType(object, "AtkText")
  location <- as.integer(location)

  w <- .RGtkCall("S_atk_text_iface_text_caret_moved", object.class, object, location, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkTextIfaceTextSelectionChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkTextIface")
  checkPtrType(object, "AtkText")

  w <- .RGtkCall("S_atk_text_iface_text_selection_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkTextIfaceTextAttributesChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkTextIface")
  checkPtrType(object, "AtkText")

  w <- .RGtkCall("S_atk_text_iface_text_attributes_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkTextIfaceGetRangeExtents <-
function(object.class, object, start.offset, end.offset, coord.type)
{
  checkPtrType(object.class, "AtkTextIface")
  checkPtrType(object, "AtkText")
  start.offset <- as.integer(start.offset)
  end.offset <- as.integer(end.offset)
  

  w <- .RGtkCall("S_atk_text_iface_get_range_extents", object.class, object, start.offset, end.offset, coord.type, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkTextIfaceGetBoundedRanges <-
function(object.class, object, rect, coord.type, x.clip.type, y.clip.type)
{
  checkPtrType(object.class, "AtkTextIface")
  checkPtrType(object, "AtkText")
  rect <- as.AtkTextRectangle(rect)
  
  
  

  w <- .RGtkCall("S_atk_text_iface_get_bounded_ranges", object.class, object, rect, coord.type, x.clip.type, y.clip.type, PACKAGE = "RGtk2")

  return(w)
}

atkValueIfaceGetCurrentValue <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkValueIface")
  checkPtrType(object, "AtkValue")

  w <- .RGtkCall("S_atk_value_iface_get_current_value", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkValueIfaceGetMaximumValue <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkValueIface")
  checkPtrType(object, "AtkValue")

  w <- .RGtkCall("S_atk_value_iface_get_maximum_value", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkValueIfaceGetMinimumValue <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkValueIface")
  checkPtrType(object, "AtkValue")

  w <- .RGtkCall("S_atk_value_iface_get_minimum_value", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

atkValueIfaceSetCurrentValue <-
function(object.class, object, value)
{
  checkPtrType(object.class, "AtkValueIface")
  checkPtrType(object, "AtkValue")
  

  w <- .RGtkCall("S_atk_value_iface_set_current_value", object.class, object, value, PACKAGE = "RGtk2")

  return(w)
}

atkValueIfaceGetMinimumIncrement <-
function(object.class, object)
{
  checkPtrType(object.class, "AtkValueIface")
  checkPtrType(object, "AtkValue")

  w <- .RGtkCall("S_atk_value_iface_get_minimum_increment", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}