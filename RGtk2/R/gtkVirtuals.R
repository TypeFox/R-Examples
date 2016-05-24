if(!exists('.virtuals')) .virtuals <- new.env()
assign("GtkAccelGroup", c("accel_changed"), .virtuals)
assign("GtkAccessible", c("connect_widget_destroyed"), .virtuals)
assign("GtkAction", c("activate", "connect_proxy", "create_menu_item", "create_tool_item", "disconnect_proxy"), .virtuals)
assign("GtkActionGroup", c("get_action"), .virtuals)
assign("GtkAdjustment", c("changed", "value_changed"), .virtuals)
assign("GtkButton", c("pressed", "released", "clicked", "enter", "leave", "activate"), .virtuals)
assign("GtkCalendar", c("month_changed", "day_selected", "day_selected_double_click", "prev_month", "next_month", "prev_year", "next_year"), .virtuals)
assign("GtkCellEditable", c("editing_done", "remove_widget", "start_editing"), .virtuals)
assign("GtkCellLayout", c("pack_start", "pack_end", "clear", "add_attribute", "set_cell_data_func", "clear_attributes", "reorder"), .virtuals)
assign("GtkCellRenderer", c("get_size", "render", "activate", "editing_canceled", "editing_started", "start_editing"), .virtuals)
assign("GtkCellRendererText", c("edited"), .virtuals)
assign("GtkCellRendererToggle", c("toggled"), .virtuals)
assign("GtkCheckButton", c("draw_indicator"), .virtuals)
assign("GtkCheckMenuItem", c("toggled", "draw_indicator"), .virtuals)
assign("GtkCList", c("set_scroll_adjustments", "refresh", "select_row", "unselect_row", "row_move", "click_column", "resize_column", "toggle_focus_row", "select_all", "unselect_all", "undo_selection", "start_selection", "end_selection", "extend_selection", "scroll_horizontal", "scroll_vertical", "toggle_add_mode", "abort_column_resize", "resync_selection", "selection_find", "draw_row", "draw_drag_highlight", "clear", "fake_unselect_all", "sort_list", "insert_row", "remove_row", "set_cell_contents", "cell_size_request"), .virtuals)
assign("GtkColorButton", c("color_set"), .virtuals)
assign("GtkColorSelection", c("color_changed"), .virtuals)
assign("GtkComboBox", c("changed", "get_active_text"), .virtuals)
assign("GtkContainer", c("add", "remove", "check_resize", "forall", "set_focus_child", "child_type", "composite_name", "set_child_property", "get_child_property"), .virtuals)
assign("GtkCTree", c("tree_select_row", "tree_unselect_row", "tree_expand", "tree_collapse", "tree_move", "change_focus_row_expansion"), .virtuals)
assign("GtkCurve", c("curve_type_changed"), .virtuals)
assign("GtkDialog", c("response", "close"), .virtuals)
assign("GtkEditable", c("insert_text", "delete_text", "changed", "do_insert_text", "do_delete_text", "get_chars", "set_selection_bounds", "get_selection_bounds", "set_position", "get_position"), .virtuals)
assign("GtkEntry", c("populate_popup", "activate", "move_cursor", "insert_at_cursor", "delete_from_cursor", "backspace", "cut_clipboard", "copy_clipboard", "paste_clipboard", "toggle_overwrite"), .virtuals)
assign("GtkEntryCompletion", c("match_selected", "action_activated", "insert_prefix"), .virtuals)
assign("GtkExpander", c("activate"), .virtuals)
assign("GtkFontButton", c("font_set"), .virtuals)
assign("GtkFrame", c("compute_child_allocation"), .virtuals)
assign("GtkHandleBox", c("child_attached", "child_detached"), .virtuals)
assign("GtkIconTheme", c("changed"), .virtuals)
assign("GtkIconView", c("set_scroll_adjustments", "item_activated", "selection_changed", "select_all", "unselect_all", "select_cursor_item", "toggle_cursor_item", "move_cursor", "activate_cursor_item"), .virtuals)
assign("GtkIMContext", c("preedit_start", "preedit_end", "preedit_changed", "commit", "retrieve_surrounding", "delete_surrounding", "set_client_window", "get_preedit_string", "filter_keypress", "focus_in", "focus_out", "reset", "set_cursor_location", "set_use_preedit", "set_surrounding", "get_surrounding"), .virtuals)
assign("GtkInputDialog", c("enable_device", "disable_device"), .virtuals)
assign("GtkItem", c("select", "deselect", "toggle"), .virtuals)
assign("GtkLabel", c("move_cursor", "copy_clipboard", "populate_popup"), .virtuals)
assign("GtkLayout", c("set_scroll_adjustments"), .virtuals)
assign("GtkList", c("selection_changed", "select_child", "unselect_child"), .virtuals)
assign("GtkListItem", c("toggle_focus_row", "select_all", "unselect_all", "undo_selection", "start_selection", "end_selection", "extend_selection", "scroll_horizontal", "scroll_vertical", "toggle_add_mode"), .virtuals)
assign("GtkMenuItem", c("activate", "activate_item", "toggle_size_request", "toggle_size_allocate"), .virtuals)
assign("GtkMenuShell", c("deactivate", "selection_done", "move_current", "activate_current", "cancel", "select_item", "insert", "get_popup_delay"), .virtuals)
assign("GtkMenuToolButton", c("show_menu"), .virtuals)
assign("GtkNotebook", c("switch_page", "select_page", "focus_tab", "change_current_page", "move_focus_out", "reorder_tab", "insert_page"), .virtuals)
assign("GtkOldEditable", c("activate", "set_editable", "move_cursor", "move_word", "move_page", "move_to_row", "move_to_column", "kill_char", "kill_word", "kill_line", "cut_clipboard", "copy_clipboard", "paste_clipboard", "update_text", "get_chars", "set_selection", "set_position"), .virtuals)
assign("GtkOptionMenu", c("changed"), .virtuals)
assign("GtkPaned", c("cycle_child_focus", "toggle_handle_focus", "move_handle", "cycle_handle_focus", "accept_position", "cancel_position"), .virtuals)
assign("GtkPlug", c("embedded"), .virtuals)
assign("GtkProgress", c("paint", "update", "act_mode_enter"), .virtuals)
assign("GtkRadioAction", c("changed"), .virtuals)
assign("GtkRadioButton", c("group_changed"), .virtuals)
assign("GtkRadioMenuItem", c("group_changed"), .virtuals)
assign("GtkRange", c("value_changed", "adjust_bounds", "move_slider", "get_range_border", "change_value"), .virtuals)
assign("GtkRcStyle", c("create_rc_style", "parse", "merge", "create_style"), .virtuals)
assign("GtkRuler", c("draw_ticks", "draw_pos"), .virtuals)
assign("GtkScale", c("format_value", "draw_value", "get_layout_offsets"), .virtuals)
assign("GtkScrolledWindow", c("scroll_child", "move_focus_out"), .virtuals)
assign("GtkSocket", c("plug_added", "plug_removed"), .virtuals)
assign("GtkSpinButton", c("input", "output", "value_changed", "change_value", "wrapped"), .virtuals)
assign("GtkStatusbar", c("text_pushed", "text_popped"), .virtuals)
assign("GtkStyle", c("realize", "unrealize", "copy", "clone", "init_from_rc", "set_background", "render_icon", "draw_hline", "draw_vline", "draw_shadow", "draw_polygon", "draw_arrow", "draw_diamond", "draw_string", "draw_box", "draw_flat_box", "draw_check", "draw_option", "draw_tab", "draw_shadow_gap", "draw_box_gap", "draw_extension", "draw_focus", "draw_slider", "draw_handle", "draw_expander", "draw_layout", "draw_resize_grip"), .virtuals)
assign("GtkText", c("set_scroll_adjustments"), .virtuals)
assign("GtkTextBuffer", c("insert_text", "insert_pixbuf", "insert_child_anchor", "delete_range", "changed", "modified_changed", "mark_set", "mark_deleted", "apply_tag", "remove_tag", "begin_user_action", "end_user_action"), .virtuals)
assign("GtkTextTag", c("event"), .virtuals)
assign("GtkTextTagTable", c("tag_changed", "tag_added", "tag_removed"), .virtuals)
assign("GtkTextView", c("set_scroll_adjustments", "populate_popup", "move_cursor", "page_horizontally", "set_anchor", "insert_at_cursor", "delete_from_cursor", "backspace", "cut_clipboard", "copy_clipboard", "paste_clipboard", "toggle_overwrite", "move_focus"), .virtuals)
assign("GtkTipsQuery", c("start_query", "stop_query", "widget_entered", "widget_selected"), .virtuals)
assign("GtkToggleAction", c("toggled"), .virtuals)
assign("GtkToggleButton", c("toggled"), .virtuals)
assign("GtkToggleToolButton", c("toggled"), .virtuals)
assign("GtkToolbar", c("orientation_changed", "style_changed", "popup_context_menu"), .virtuals)
assign("GtkToolButton", c("clicked"), .virtuals)
assign("GtkToolItem", c("create_menu_proxy", "toolbar_reconfigured", "set_tooltip"), .virtuals)
assign("GtkTreeDragSource", c("row_draggable", "drag_data_get", "drag_data_delete"), .virtuals)
assign("GtkTreeDragDest", c("drag_data_received", "row_drop_possible"), .virtuals)
assign("GtkTreeSelection", c("changed", "changed"), .virtuals)
assign("GtkTree", c("select_child", "unselect_child"), .virtuals)
assign("GtkTreeItem", c("expand", "collapse"), .virtuals)
assign("GtkTreeModel", c("row_changed", "row_inserted", "row_has_child_toggled", "row_deleted", "rows_reordered", "get_flags", "get_n_columns", "get_column_type", "get_iter", "get_path", "get_value", "iter_next", "iter_children", "iter_has_child", "iter_n_children", "iter_nth_child", "iter_parent", "ref_node", "unref_node"), .virtuals)
assign("GtkTreeSortable", c("sort_column_changed", "get_sort_column_id", "set_sort_column_id", "set_sort_func", "set_default_sort_func", "has_default_sort_func"), .virtuals)
assign("GtkTreeView", c("set_scroll_adjustments", "row_activated", "test_expand_row", "test_collapse_row", "row_expanded", "row_collapsed", "columns_changed", "cursor_changed", "move_cursor", "select_all", "unselect_all", "select_cursor_row", "toggle_cursor_row", "expand_collapse_cursor_row", "select_cursor_parent", "start_interactive_search"), .virtuals)
assign("GtkTreeViewColumn", c("clicked"), .virtuals)
assign("GtkUIManager", c("add_widget", "actions_changed", "connect_proxy", "disconnect_proxy", "pre_activate", "post_activate", "get_widget", "get_action"), .virtuals)
assign("GtkViewport", c("set_scroll_adjustments"), .virtuals)
assign("GtkWidget", c("dispatch_child_properties_changed", "show", "show_all", "hide", "hide_all", "map", "unmap", "realize", "unrealize", "size_request", "size_allocate", "state_changed", "parent_set", "hierarchy_changed", "style_set", "direction_changed", "grab_notify", "child_notify", "mnemonic_activate", "grab_focus", "focus", "event", "button_press_event", "button_release_event", "scroll_event", "motion_notify_event", "delete_event", "destroy_event", "expose_event", "key_press_event", "key_release_event", "enter_notify_event", "leave_notify_event", "configure_event", "focus_in_event", "focus_out_event", "map_event", "unmap_event", "property_notify_event", "selection_clear_event", "selection_request_event", "selection_notify_event", "proximity_in_event", "proximity_out_event", "visibility_notify_event", "client_event", "no_expose_event", "window_state_event", "selection_get", "selection_received", "drag_begin", "drag_end", "drag_data_get", "drag_data_delete", "drag_leave", "drag_motion", "drag_drop", "drag_data_received", "popup_menu", "show_help", "get_accessible", "screen_changed", "can_activate_accel", "grab_broken_event", "composited_changed"), .virtuals)
assign("GtkWindow", c("set_focus", "frame_event", "activate_focus", "activate_default", "move_focus", "keys_changed"), .virtuals)
assign("GtkAssistant", c("prepare", "apply", "close", "cancel"), .virtuals)
assign("GtkCellRendererAccel", c("accel_edited", "accel_cleared"), .virtuals)
assign("GtkPrintOperation", c("done", "begin_print", "paginate", "request_page_setup", "draw_page", "end_print", "status_changed", "create_custom_widget", "custom_widget_apply", "preview"), .virtuals)
assign("GtkPrintOperationPreview", c("ready", "got_page_size", "render_page", "is_selected", "end_preview"), .virtuals)
assign("GtkRecentChooser", c("set_current_uri", "get_current_uri", "select_uri", "unselect_uri", "select_all", "unselect_all", "get_items", "get_recent_manager", "add_filter", "remove_filter", "list_filters", "set_sort_func", "item_activated", "selection_changed"), .virtuals)
assign("GtkRecentManager", c("changed"), .virtuals)
assign("GtkStatusIcon", c("activate", "popup_menu", "size_changed"), .virtuals)
assign("GtkBuildable", c("set_name", "get_name", "add_child", "set_buildable_property", "construct_child", "custom_tag_start", "custom_tag_end", "custom_finished", "parser_finished", "get_internal_child"), .virtuals)
assign("GtkBuilder", c("get_type_from_name"), .virtuals)
assign("GtkToolShell", c("get_icon_size", "get_orientation", "get_style", "get_relief_style", "rebuild_menu"), .virtuals)
assign("GtkActivatable", c("update", "sync_action_properties"), .virtuals)


gtkAccelGroupClassAccelChanged <-
function(object.class, object, keyval, modifier, accel.closure)
{
  checkPtrType(object.class, "GtkAccelGroupClass")
  checkPtrType(object, "GtkAccelGroup")
  keyval <- as.numeric(keyval)
  
  accel.closure <- as.GClosure(accel.closure)

  w <- .RGtkCall("S_gtk_accel_group_class_accel_changed", object.class, object, keyval, modifier, accel.closure, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkAccessibleClassConnectWidgetDestroyed <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkAccessibleClass")
  checkPtrType(object, "GtkAccessible")

  w <- .RGtkCall("S_gtk_accessible_class_connect_widget_destroyed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkActionClassActivate <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkActionClass")
  checkPtrType(object, "GtkAction")

  w <- .RGtkCall("S_gtk_action_class_activate", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkActionClassConnectProxy <-
function(object.class, object, proxy)
{
  checkPtrType(object.class, "GtkActionClass")
  checkPtrType(object, "GtkAction")
  checkPtrType(proxy, "GtkWidget")

  w <- .RGtkCall("S_gtk_action_class_connect_proxy", object.class, object, proxy, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkActionClassCreateMenuItem <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkActionClass")
  checkPtrType(object, "GtkAction")

  w <- .RGtkCall("S_gtk_action_class_create_menu_item", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkActionClassCreateToolItem <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkActionClass")
  checkPtrType(object, "GtkAction")

  w <- .RGtkCall("S_gtk_action_class_create_tool_item", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkActionClassDisconnectProxy <-
function(object.class, object, proxy)
{
  checkPtrType(object.class, "GtkActionClass")
  checkPtrType(object, "GtkAction")
  checkPtrType(proxy, "GtkWidget")

  w <- .RGtkCall("S_gtk_action_class_disconnect_proxy", object.class, object, proxy, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkActionGroupClassGetAction <-
function(object.class, object, action.name)
{
  checkPtrType(object.class, "GtkActionGroupClass")
  checkPtrType(object, "GtkActionGroup")
  action.name <- as.character(action.name)

  w <- .RGtkCall("S_gtk_action_group_class_get_action", object.class, object, action.name, PACKAGE = "RGtk2")

  return(w)
}

gtkAdjustmentClassChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkAdjustmentClass")
  checkPtrType(object, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_adjustment_class_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkAdjustmentClassValueChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkAdjustmentClass")
  checkPtrType(object, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_adjustment_class_value_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkButtonClassPressed <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkButtonClass")
  checkPtrType(object, "GtkButton")

  w <- .RGtkCall("S_gtk_button_class_pressed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkButtonClassReleased <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkButtonClass")
  checkPtrType(object, "GtkButton")

  w <- .RGtkCall("S_gtk_button_class_released", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkButtonClassClicked <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkButtonClass")
  checkPtrType(object, "GtkButton")

  w <- .RGtkCall("S_gtk_button_class_clicked", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkButtonClassEnter <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkButtonClass")
  checkPtrType(object, "GtkButton")

  w <- .RGtkCall("S_gtk_button_class_enter", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkButtonClassLeave <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkButtonClass")
  checkPtrType(object, "GtkButton")

  w <- .RGtkCall("S_gtk_button_class_leave", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkButtonClassActivate <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkButtonClass")
  checkPtrType(object, "GtkButton")

  w <- .RGtkCall("S_gtk_button_class_activate", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCalendarClassMonthChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkCalendarClass")
  checkPtrType(object, "GtkCalendar")

  w <- .RGtkCall("S_gtk_calendar_class_month_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCalendarClassDaySelected <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkCalendarClass")
  checkPtrType(object, "GtkCalendar")

  w <- .RGtkCall("S_gtk_calendar_class_day_selected", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCalendarClassDaySelectedDoubleClick <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkCalendarClass")
  checkPtrType(object, "GtkCalendar")

  w <- .RGtkCall("S_gtk_calendar_class_day_selected_double_click", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCalendarClassPrevMonth <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkCalendarClass")
  checkPtrType(object, "GtkCalendar")

  w <- .RGtkCall("S_gtk_calendar_class_prev_month", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCalendarClassNextMonth <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkCalendarClass")
  checkPtrType(object, "GtkCalendar")

  w <- .RGtkCall("S_gtk_calendar_class_next_month", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCalendarClassPrevYear <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkCalendarClass")
  checkPtrType(object, "GtkCalendar")

  w <- .RGtkCall("S_gtk_calendar_class_prev_year", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCalendarClassNextYear <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkCalendarClass")
  checkPtrType(object, "GtkCalendar")

  w <- .RGtkCall("S_gtk_calendar_class_next_year", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCellEditableIfaceEditingDone <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkCellEditableIface")
  checkPtrType(object, "GtkCellEditable")

  w <- .RGtkCall("S_gtk_cell_editable_iface_editing_done", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCellEditableIfaceRemoveWidget <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkCellEditableIface")
  checkPtrType(object, "GtkCellEditable")

  w <- .RGtkCall("S_gtk_cell_editable_iface_remove_widget", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCellEditableIfaceStartEditing <-
function(object.class, object, event)
{
  checkPtrType(object.class, "GtkCellEditableIface")
  checkPtrType(object, "GtkCellEditable")
  checkPtrType(event, "GdkEvent")

  w <- .RGtkCall("S_gtk_cell_editable_iface_start_editing", object.class, object, event, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCellLayoutIfacePackStart <-
function(object.class, object, cell, expand)
{
  checkPtrType(object.class, "GtkCellLayoutIface")
  checkPtrType(object, "GtkCellLayout")
  checkPtrType(cell, "GtkCellRenderer")
  expand <- as.logical(expand)

  w <- .RGtkCall("S_gtk_cell_layout_iface_pack_start", object.class, object, cell, expand, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCellLayoutIfacePackEnd <-
function(object.class, object, cell, expand)
{
  checkPtrType(object.class, "GtkCellLayoutIface")
  checkPtrType(object, "GtkCellLayout")
  checkPtrType(cell, "GtkCellRenderer")
  expand <- as.logical(expand)

  w <- .RGtkCall("S_gtk_cell_layout_iface_pack_end", object.class, object, cell, expand, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCellLayoutIfaceClear <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkCellLayoutIface")
  checkPtrType(object, "GtkCellLayout")

  w <- .RGtkCall("S_gtk_cell_layout_iface_clear", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCellLayoutIfaceAddAttribute <-
function(object.class, object, cell, attribute, column)
{
  checkPtrType(object.class, "GtkCellLayoutIface")
  checkPtrType(object, "GtkCellLayout")
  checkPtrType(cell, "GtkCellRenderer")
  attribute <- as.character(attribute)
  column <- as.integer(column)

  w <- .RGtkCall("S_gtk_cell_layout_iface_add_attribute", object.class, object, cell, attribute, column, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCellLayoutIfaceSetCellDataFunc <-
function(object.class, object, cell, func, func.data)
{
  checkPtrType(object.class, "GtkCellLayoutIface")
  checkPtrType(object, "GtkCellLayout")
  checkPtrType(cell, "GtkCellRenderer")
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_cell_layout_iface_set_cell_data_func", object.class, object, cell, func, func.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCellLayoutIfaceClearAttributes <-
function(object.class, object, cell)
{
  checkPtrType(object.class, "GtkCellLayoutIface")
  checkPtrType(object, "GtkCellLayout")
  checkPtrType(cell, "GtkCellRenderer")

  w <- .RGtkCall("S_gtk_cell_layout_iface_clear_attributes", object.class, object, cell, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCellLayoutIfaceReorder <-
function(object.class, object, cell, position)
{
  checkPtrType(object.class, "GtkCellLayoutIface")
  checkPtrType(object, "GtkCellLayout")
  checkPtrType(cell, "GtkCellRenderer")
  position <- as.integer(position)

  w <- .RGtkCall("S_gtk_cell_layout_iface_reorder", object.class, object, cell, position, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCellRendererClassGetSize <-
function(object.class, object, widget, cell.area)
{
  checkPtrType(object.class, "GtkCellRendererClass")
  checkPtrType(object, "GtkCellRenderer")
  checkPtrType(widget, "GtkWidget")
  cell.area <- as.GdkRectangle(cell.area)

  w <- .RGtkCall("S_gtk_cell_renderer_class_get_size", object.class, object, widget, cell.area, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCellRendererClassRender <-
function(object.class, object, window, widget, background.area, cell.area, expose.area, flags)
{
  checkPtrType(object.class, "GtkCellRendererClass")
  checkPtrType(object, "GtkCellRenderer")
  checkPtrType(window, "GdkDrawable")
  checkPtrType(widget, "GtkWidget")
  background.area <- as.GdkRectangle(background.area)
  cell.area <- as.GdkRectangle(cell.area)
  expose.area <- as.GdkRectangle(expose.area)
  

  w <- .RGtkCall("S_gtk_cell_renderer_class_render", object.class, object, window, widget, background.area, cell.area, expose.area, flags, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCellRendererClassActivate <-
function(object.class, object, event, widget, path, background.area, cell.area, flags)
{
  checkPtrType(object.class, "GtkCellRendererClass")
  checkPtrType(object, "GtkCellRenderer")
  checkPtrType(event, "GdkEvent")
  checkPtrType(widget, "GtkWidget")
  path <- as.character(path)
  background.area <- as.GdkRectangle(background.area)
  cell.area <- as.GdkRectangle(cell.area)
  

  w <- .RGtkCall("S_gtk_cell_renderer_class_activate", object.class, object, event, widget, path, background.area, cell.area, flags, PACKAGE = "RGtk2")

  return(w)
}

gtkCellRendererClassEditingCanceled <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkCellRendererClass")
  checkPtrType(object, "GtkCellRenderer")

  w <- .RGtkCall("S_gtk_cell_renderer_class_editing_canceled", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCellRendererClassEditingStarted <-
function(object.class, object, editable, path)
{
  checkPtrType(object.class, "GtkCellRendererClass")
  checkPtrType(object, "GtkCellRenderer")
  checkPtrType(editable, "GtkCellEditable")
  path <- as.character(path)

  w <- .RGtkCall("S_gtk_cell_renderer_class_editing_started", object.class, object, editable, path, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCellRendererClassStartEditing <-
function(object.class, object, event, widget, path, background.area, cell.area, flags)
{
  checkPtrType(object.class, "GtkCellRendererClass")
  checkPtrType(object, "GtkCellRenderer")
  checkPtrType(event, "GdkEvent")
  checkPtrType(widget, "GtkWidget")
  path <- as.character(path)
  background.area <- as.GdkRectangle(background.area)
  cell.area <- as.GdkRectangle(cell.area)
  

  w <- .RGtkCall("S_gtk_cell_renderer_class_start_editing", object.class, object, event, widget, path, background.area, cell.area, flags, PACKAGE = "RGtk2")

  return(w)
}

gtkCellRendererTextClassEdited <-
function(object.class, object, path, new.text)
{
  checkPtrType(object.class, "GtkCellRendererTextClass")
  checkPtrType(object, "GtkCellRendererText")
  path <- as.character(path)
  new.text <- as.character(new.text)

  w <- .RGtkCall("S_gtk_cell_renderer_text_class_edited", object.class, object, path, new.text, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCellRendererToggleClassToggled <-
function(object.class, object, path)
{
  checkPtrType(object.class, "GtkCellRendererToggleClass")
  checkPtrType(object, "GtkCellRendererToggle")
  path <- as.character(path)

  w <- .RGtkCall("S_gtk_cell_renderer_toggle_class_toggled", object.class, object, path, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCheckButtonClassDrawIndicator <-
function(object.class, object, area)
{
  checkPtrType(object.class, "GtkCheckButtonClass")
  checkPtrType(object, "GtkCheckButton")
  area <- as.GdkRectangle(area)

  w <- .RGtkCall("S_gtk_check_button_class_draw_indicator", object.class, object, area, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCheckMenuItemClassToggled <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkCheckMenuItemClass")
  checkPtrType(object, "GtkCheckMenuItem")

  w <- .RGtkCall("S_gtk_check_menu_item_class_toggled", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCheckMenuItemClassDrawIndicator <-
function(object.class, object, area)
{
  checkPtrType(object.class, "GtkCheckMenuItemClass")
  checkPtrType(object, "GtkCheckMenuItem")
  area <- as.GdkRectangle(area)

  w <- .RGtkCall("S_gtk_check_menu_item_class_draw_indicator", object.class, object, area, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCListClassSetScrollAdjustments <-
function(object.class, object, hadjustment, vadjustment)
{
  checkPtrType(object.class, "GtkCListClass")
  checkPtrType(object, "GtkCList")
  checkPtrType(hadjustment, "GtkAdjustment")
  checkPtrType(vadjustment, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_clist_class_set_scroll_adjustments", object.class, object, hadjustment, vadjustment, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCListClassRefresh <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkCListClass")
  checkPtrType(object, "GtkCList")

  w <- .RGtkCall("S_gtk_clist_class_refresh", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCListClassSelectRow <-
function(object.class, object, row, column, event)
{
  checkPtrType(object.class, "GtkCListClass")
  checkPtrType(object, "GtkCList")
  row <- as.integer(row)
  column <- as.integer(column)
  checkPtrType(event, "GdkEvent")

  w <- .RGtkCall("S_gtk_clist_class_select_row", object.class, object, row, column, event, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCListClassUnselectRow <-
function(object.class, object, row, column, event)
{
  checkPtrType(object.class, "GtkCListClass")
  checkPtrType(object, "GtkCList")
  row <- as.integer(row)
  column <- as.integer(column)
  checkPtrType(event, "GdkEvent")

  w <- .RGtkCall("S_gtk_clist_class_unselect_row", object.class, object, row, column, event, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCListClassRowMove <-
function(object.class, object, source.row, dest.row)
{
  checkPtrType(object.class, "GtkCListClass")
  checkPtrType(object, "GtkCList")
  source.row <- as.integer(source.row)
  dest.row <- as.integer(dest.row)

  w <- .RGtkCall("S_gtk_clist_class_row_move", object.class, object, source.row, dest.row, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCListClassClickColumn <-
function(object.class, object, column)
{
  checkPtrType(object.class, "GtkCListClass")
  checkPtrType(object, "GtkCList")
  column <- as.integer(column)

  w <- .RGtkCall("S_gtk_clist_class_click_column", object.class, object, column, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCListClassResizeColumn <-
function(object.class, object, column, width)
{
  checkPtrType(object.class, "GtkCListClass")
  checkPtrType(object, "GtkCList")
  column <- as.integer(column)
  width <- as.integer(width)

  w <- .RGtkCall("S_gtk_clist_class_resize_column", object.class, object, column, width, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCListClassToggleFocusRow <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkCListClass")
  checkPtrType(object, "GtkCList")

  w <- .RGtkCall("S_gtk_clist_class_toggle_focus_row", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCListClassSelectAll <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkCListClass")
  checkPtrType(object, "GtkCList")

  w <- .RGtkCall("S_gtk_clist_class_select_all", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCListClassUnselectAll <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkCListClass")
  checkPtrType(object, "GtkCList")

  w <- .RGtkCall("S_gtk_clist_class_unselect_all", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCListClassUndoSelection <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkCListClass")
  checkPtrType(object, "GtkCList")

  w <- .RGtkCall("S_gtk_clist_class_undo_selection", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCListClassStartSelection <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkCListClass")
  checkPtrType(object, "GtkCList")

  w <- .RGtkCall("S_gtk_clist_class_start_selection", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCListClassEndSelection <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkCListClass")
  checkPtrType(object, "GtkCList")

  w <- .RGtkCall("S_gtk_clist_class_end_selection", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCListClassExtendSelection <-
function(object.class, object, scroll.type, position, auto.start.selection)
{
  checkPtrType(object.class, "GtkCListClass")
  checkPtrType(object, "GtkCList")
  
  position <- as.numeric(position)
  auto.start.selection <- as.logical(auto.start.selection)

  w <- .RGtkCall("S_gtk_clist_class_extend_selection", object.class, object, scroll.type, position, auto.start.selection, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCListClassScrollHorizontal <-
function(object.class, object, scroll.type, position)
{
  checkPtrType(object.class, "GtkCListClass")
  checkPtrType(object, "GtkCList")
  
  position <- as.numeric(position)

  w <- .RGtkCall("S_gtk_clist_class_scroll_horizontal", object.class, object, scroll.type, position, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCListClassScrollVertical <-
function(object.class, object, scroll.type, position)
{
  checkPtrType(object.class, "GtkCListClass")
  checkPtrType(object, "GtkCList")
  
  position <- as.numeric(position)

  w <- .RGtkCall("S_gtk_clist_class_scroll_vertical", object.class, object, scroll.type, position, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCListClassToggleAddMode <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkCListClass")
  checkPtrType(object, "GtkCList")

  w <- .RGtkCall("S_gtk_clist_class_toggle_add_mode", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCListClassAbortColumnResize <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkCListClass")
  checkPtrType(object, "GtkCList")

  w <- .RGtkCall("S_gtk_clist_class_abort_column_resize", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCListClassResyncSelection <-
function(object.class, object, event)
{
  checkPtrType(object.class, "GtkCListClass")
  checkPtrType(object, "GtkCList")
  checkPtrType(event, "GdkEvent")

  w <- .RGtkCall("S_gtk_clist_class_resync_selection", object.class, object, event, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCListClassSelectionFind <-
function(object.class, object, row.number, row.list.element)
{
  checkPtrType(object.class, "GtkCListClass")
  checkPtrType(object, "GtkCList")
  row.number <- as.integer(row.number)
  row.list.element <- lapply(row.list.element, function(x) { x <- as.GList(x); x })

  w <- .RGtkCall("S_gtk_clist_class_selection_find", object.class, object, row.number, row.list.element, PACKAGE = "RGtk2")

  return(w)
}

gtkCListClassDrawRow <-
function(object.class, object, area, row, clist.row)
{
  checkPtrType(object.class, "GtkCListClass")
  checkPtrType(object, "GtkCList")
  area <- as.GdkRectangle(area)
  row <- as.integer(row)
  checkPtrType(clist.row, "GtkCListRow")

  w <- .RGtkCall("S_gtk_clist_class_draw_row", object.class, object, area, row, clist.row, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCListClassDrawDragHighlight <-
function(object.class, object, target.row, target.row.number, drag.pos)
{
  checkPtrType(object.class, "GtkCListClass")
  checkPtrType(object, "GtkCList")
  checkPtrType(target.row, "GtkCListRow")
  target.row.number <- as.integer(target.row.number)
  

  w <- .RGtkCall("S_gtk_clist_class_draw_drag_highlight", object.class, object, target.row, target.row.number, drag.pos, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCListClassClear <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkCListClass")
  checkPtrType(object, "GtkCList")

  w <- .RGtkCall("S_gtk_clist_class_clear", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCListClassFakeUnselectAll <-
function(object.class, object, row)
{
  checkPtrType(object.class, "GtkCListClass")
  checkPtrType(object, "GtkCList")
  row <- as.integer(row)

  w <- .RGtkCall("S_gtk_clist_class_fake_unselect_all", object.class, object, row, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCListClassSortList <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkCListClass")
  checkPtrType(object, "GtkCList")

  w <- .RGtkCall("S_gtk_clist_class_sort_list", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCListClassInsertRow <-
function(object.class, object, row, text)
{
  checkPtrType(object.class, "GtkCListClass")
  checkPtrType(object, "GtkCList")
  row <- as.integer(row)
  text <- as.list(as.character(text))

  w <- .RGtkCall("S_gtk_clist_class_insert_row", object.class, object, row, text, PACKAGE = "RGtk2")

  return(w)
}

gtkCListClassRemoveRow <-
function(object.class, object, row)
{
  checkPtrType(object.class, "GtkCListClass")
  checkPtrType(object, "GtkCList")
  row <- as.integer(row)

  w <- .RGtkCall("S_gtk_clist_class_remove_row", object.class, object, row, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCListClassSetCellContents <-
function(object.class, object, clist.row, column, type, text, spacing, pixmap, mask)
{
  checkPtrType(object.class, "GtkCListClass")
  checkPtrType(object, "GtkCList")
  checkPtrType(clist.row, "GtkCListRow")
  column <- as.integer(column)
  
  text <- as.character(text)
  spacing <- as.raw(spacing)
  checkPtrType(pixmap, "GdkPixmap")
  checkPtrType(mask, "GdkBitmap")

  w <- .RGtkCall("S_gtk_clist_class_set_cell_contents", object.class, object, clist.row, column, type, text, spacing, pixmap, mask, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCListClassCellSizeRequest <-
function(object.class, object, clist.row, column, requisition)
{
  checkPtrType(object.class, "GtkCListClass")
  checkPtrType(object, "GtkCList")
  checkPtrType(clist.row, "GtkCListRow")
  column <- as.integer(column)
  checkPtrType(requisition, "GtkRequisition")

  w <- .RGtkCall("S_gtk_clist_class_cell_size_request", object.class, object, clist.row, column, requisition, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkColorButtonClassColorSet <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkColorButtonClass")
  checkPtrType(object, "GtkColorButton")

  w <- .RGtkCall("S_gtk_color_button_class_color_set", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkColorSelectionClassColorChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkColorSelectionClass")
  checkPtrType(object, "GtkColorSelection")

  w <- .RGtkCall("S_gtk_color_selection_class_color_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkComboBoxClassChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkComboBoxClass")
  checkPtrType(object, "GtkComboBox")

  w <- .RGtkCall("S_gtk_combo_box_class_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkComboBoxClassGetActiveText <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkComboBoxClass")
  checkPtrType(object, "GtkComboBox")

  w <- .RGtkCall("S_gtk_combo_box_class_get_active_text", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkContainerClassAdd <-
function(object.class, object, widget)
{
  checkPtrType(object.class, "GtkContainerClass")
  checkPtrType(object, "GtkContainer")
  checkPtrType(widget, "GtkWidget")

  w <- .RGtkCall("S_gtk_container_class_add", object.class, object, widget, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkContainerClassRemove <-
function(object.class, object, widget)
{
  checkPtrType(object.class, "GtkContainerClass")
  checkPtrType(object, "GtkContainer")
  checkPtrType(widget, "GtkWidget")

  w <- .RGtkCall("S_gtk_container_class_remove", object.class, object, widget, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkContainerClassCheckResize <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkContainerClass")
  checkPtrType(object, "GtkContainer")

  w <- .RGtkCall("S_gtk_container_class_check_resize", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkContainerClassForall <-
function(object.class, object, include.internals, callback, callback.data)
{
  checkPtrType(object.class, "GtkContainerClass")
  checkPtrType(object, "GtkContainer")
  include.internals <- as.logical(include.internals)
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gtk_container_class_forall", object.class, object, include.internals, callback, callback.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkContainerClassSetFocusChild <-
function(object.class, object, widget)
{
  checkPtrType(object.class, "GtkContainerClass")
  checkPtrType(object, "GtkContainer")
  checkPtrType(widget, "GtkWidget")

  w <- .RGtkCall("S_gtk_container_class_set_focus_child", object.class, object, widget, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkContainerClassChildType <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkContainerClass")
  checkPtrType(object, "GtkContainer")

  w <- .RGtkCall("S_gtk_container_class_child_type", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkContainerClassCompositeName <-
function(object.class, object, child)
{
  checkPtrType(object.class, "GtkContainerClass")
  checkPtrType(object, "GtkContainer")
  checkPtrType(child, "GtkWidget")

  w <- .RGtkCall("S_gtk_container_class_composite_name", object.class, object, child, PACKAGE = "RGtk2")

  return(w)
}

gtkContainerClassSetChildProperty <-
function(object.class, object, child, property.id, value, pspec)
{
  checkPtrType(object.class, "GtkContainerClass")
  checkPtrType(object, "GtkContainer")
  checkPtrType(child, "GtkWidget")
  property.id <- as.numeric(property.id)
  
  pspec <- as.GParamSpec(pspec)

  w <- .RGtkCall("S_gtk_container_class_set_child_property", object.class, object, child, property.id, value, pspec, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkContainerClassGetChildProperty <-
function(object.class, object, child, property.id, value, pspec)
{
  checkPtrType(object.class, "GtkContainerClass")
  checkPtrType(object, "GtkContainer")
  checkPtrType(child, "GtkWidget")
  property.id <- as.numeric(property.id)
  
  pspec <- as.GParamSpec(pspec)

  w <- .RGtkCall("S_gtk_container_class_get_child_property", object.class, object, child, property.id, value, pspec, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCTreeClassTreeSelectRow <-
function(object.class, object, row, column)
{
  checkPtrType(object.class, "GtkCTreeClass")
  checkPtrType(object, "GtkCTree")
  checkPtrType(row, "GtkCTreeNode")
  column <- as.integer(column)

  w <- .RGtkCall("S_gtk_ctree_class_tree_select_row", object.class, object, row, column, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCTreeClassTreeUnselectRow <-
function(object.class, object, row, column)
{
  checkPtrType(object.class, "GtkCTreeClass")
  checkPtrType(object, "GtkCTree")
  checkPtrType(row, "GtkCTreeNode")
  column <- as.integer(column)

  w <- .RGtkCall("S_gtk_ctree_class_tree_unselect_row", object.class, object, row, column, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCTreeClassTreeExpand <-
function(object.class, object, node)
{
  checkPtrType(object.class, "GtkCTreeClass")
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")

  w <- .RGtkCall("S_gtk_ctree_class_tree_expand", object.class, object, node, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCTreeClassTreeCollapse <-
function(object.class, object, node)
{
  checkPtrType(object.class, "GtkCTreeClass")
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")

  w <- .RGtkCall("S_gtk_ctree_class_tree_collapse", object.class, object, node, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCTreeClassTreeMove <-
function(object.class, object, node, new.parent, new.sibling)
{
  checkPtrType(object.class, "GtkCTreeClass")
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")
  checkPtrType(new.parent, "GtkCTreeNode")
  checkPtrType(new.sibling, "GtkCTreeNode")

  w <- .RGtkCall("S_gtk_ctree_class_tree_move", object.class, object, node, new.parent, new.sibling, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCTreeClassChangeFocusRowExpansion <-
function(object.class, object, action)
{
  checkPtrType(object.class, "GtkCTreeClass")
  checkPtrType(object, "GtkCTree")
  

  w <- .RGtkCall("S_gtk_ctree_class_change_focus_row_expansion", object.class, object, action, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCurveClassCurveTypeChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkCurveClass")
  checkPtrType(object, "GtkCurve")

  w <- .RGtkCall("S_gtk_curve_class_curve_type_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkDialogClassResponse <-
function(object.class, object, response.id)
{
  checkPtrType(object.class, "GtkDialogClass")
  checkPtrType(object, "GtkDialog")
  response.id <- as.integer(response.id)

  w <- .RGtkCall("S_gtk_dialog_class_response", object.class, object, response.id, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkDialogClassClose <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkDialogClass")
  checkPtrType(object, "GtkDialog")

  w <- .RGtkCall("S_gtk_dialog_class_close", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkEditableIfaceInsertText <-
function(object.class, object, text, position)
{
  checkPtrType(object.class, "GtkEditableIface")
  checkPtrType(object, "GtkEditable")
  text <- as.character(text)
  position <- as.list(as.integer(position))

  w <- .RGtkCall("S_gtk_editable_iface_insert_text", object.class, object, text, position, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkEditableIfaceDeleteText <-
function(object.class, object, start.pos, end.pos)
{
  checkPtrType(object.class, "GtkEditableIface")
  checkPtrType(object, "GtkEditable")
  start.pos <- as.integer(start.pos)
  end.pos <- as.integer(end.pos)

  w <- .RGtkCall("S_gtk_editable_iface_delete_text", object.class, object, start.pos, end.pos, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkEditableIfaceChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkEditableIface")
  checkPtrType(object, "GtkEditable")

  w <- .RGtkCall("S_gtk_editable_iface_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkEditableIfaceDoInsertText <-
function(object.class, object, text, position)
{
  checkPtrType(object.class, "GtkEditableIface")
  checkPtrType(object, "GtkEditable")
  text <- as.character(text)
  position <- as.list(as.integer(position))

  w <- .RGtkCall("S_gtk_editable_iface_do_insert_text", object.class, object, text, position, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkEditableIfaceDoDeleteText <-
function(object.class, object, start.pos, end.pos)
{
  checkPtrType(object.class, "GtkEditableIface")
  checkPtrType(object, "GtkEditable")
  start.pos <- as.integer(start.pos)
  end.pos <- as.integer(end.pos)

  w <- .RGtkCall("S_gtk_editable_iface_do_delete_text", object.class, object, start.pos, end.pos, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkEditableIfaceGetChars <-
function(object.class, object, start.pos, end.pos)
{
  checkPtrType(object.class, "GtkEditableIface")
  checkPtrType(object, "GtkEditable")
  start.pos <- as.integer(start.pos)
  end.pos <- as.integer(end.pos)

  w <- .RGtkCall("S_gtk_editable_iface_get_chars", object.class, object, start.pos, end.pos, PACKAGE = "RGtk2")

  return(w)
}

gtkEditableIfaceSetSelectionBounds <-
function(object.class, object, start.pos, end.pos)
{
  checkPtrType(object.class, "GtkEditableIface")
  checkPtrType(object, "GtkEditable")
  start.pos <- as.integer(start.pos)
  end.pos <- as.integer(end.pos)

  w <- .RGtkCall("S_gtk_editable_iface_set_selection_bounds", object.class, object, start.pos, end.pos, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkEditableIfaceGetSelectionBounds <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkEditableIface")
  checkPtrType(object, "GtkEditable")

  w <- .RGtkCall("S_gtk_editable_iface_get_selection_bounds", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkEditableIfaceSetPosition <-
function(object.class, object, position)
{
  checkPtrType(object.class, "GtkEditableIface")
  checkPtrType(object, "GtkEditable")
  position <- as.integer(position)

  w <- .RGtkCall("S_gtk_editable_iface_set_position", object.class, object, position, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkEditableIfaceGetPosition <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkEditableIface")
  checkPtrType(object, "GtkEditable")

  w <- .RGtkCall("S_gtk_editable_iface_get_position", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkEntryClassPopulatePopup <-
function(object.class, object, menu)
{
  checkPtrType(object.class, "GtkEntryClass")
  checkPtrType(object, "GtkEntry")
  checkPtrType(menu, "GtkMenu")

  w <- .RGtkCall("S_gtk_entry_class_populate_popup", object.class, object, menu, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkEntryClassActivate <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkEntryClass")
  checkPtrType(object, "GtkEntry")

  w <- .RGtkCall("S_gtk_entry_class_activate", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkEntryClassMoveCursor <-
function(object.class, object, step, count, extend.selection)
{
  checkPtrType(object.class, "GtkEntryClass")
  checkPtrType(object, "GtkEntry")
  
  count <- as.integer(count)
  extend.selection <- as.logical(extend.selection)

  w <- .RGtkCall("S_gtk_entry_class_move_cursor", object.class, object, step, count, extend.selection, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkEntryClassInsertAtCursor <-
function(object.class, object, str)
{
  checkPtrType(object.class, "GtkEntryClass")
  checkPtrType(object, "GtkEntry")
  str <- as.character(str)

  w <- .RGtkCall("S_gtk_entry_class_insert_at_cursor", object.class, object, str, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkEntryClassDeleteFromCursor <-
function(object.class, object, type, count)
{
  checkPtrType(object.class, "GtkEntryClass")
  checkPtrType(object, "GtkEntry")
  
  count <- as.integer(count)

  w <- .RGtkCall("S_gtk_entry_class_delete_from_cursor", object.class, object, type, count, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkEntryClassBackspace <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkEntryClass")
  checkPtrType(object, "GtkEntry")

  w <- .RGtkCall("S_gtk_entry_class_backspace", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkEntryClassCutClipboard <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkEntryClass")
  checkPtrType(object, "GtkEntry")

  w <- .RGtkCall("S_gtk_entry_class_cut_clipboard", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkEntryClassCopyClipboard <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkEntryClass")
  checkPtrType(object, "GtkEntry")

  w <- .RGtkCall("S_gtk_entry_class_copy_clipboard", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkEntryClassPasteClipboard <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkEntryClass")
  checkPtrType(object, "GtkEntry")

  w <- .RGtkCall("S_gtk_entry_class_paste_clipboard", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkEntryClassToggleOverwrite <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkEntryClass")
  checkPtrType(object, "GtkEntry")

  w <- .RGtkCall("S_gtk_entry_class_toggle_overwrite", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkEntryCompletionClassMatchSelected <-
function(object.class, object, model, iter)
{
  checkPtrType(object.class, "GtkEntryCompletionClass")
  checkPtrType(object, "GtkEntryCompletion")
  checkPtrType(model, "GtkTreeModel")
  checkPtrType(iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_entry_completion_class_match_selected", object.class, object, model, iter, PACKAGE = "RGtk2")

  return(w)
}

gtkEntryCompletionClassActionActivated <-
function(object.class, object, index.)
{
  checkPtrType(object.class, "GtkEntryCompletionClass")
  checkPtrType(object, "GtkEntryCompletion")
  index. <- as.integer(index.)

  w <- .RGtkCall("S_gtk_entry_completion_class_action_activated", object.class, object, index., PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkEntryCompletionClassInsertPrefix <-
function(object.class, object, prefix)
{
  checkPtrType(object.class, "GtkEntryCompletionClass")
  checkPtrType(object, "GtkEntryCompletion")
  prefix <- as.character(prefix)

  w <- .RGtkCall("S_gtk_entry_completion_class_insert_prefix", object.class, object, prefix, PACKAGE = "RGtk2")

  return(w)
}

gtkExpanderClassActivate <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkExpanderClass")
  checkPtrType(object, "GtkExpander")

  w <- .RGtkCall("S_gtk_expander_class_activate", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkFontButtonClassFontSet <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkFontButtonClass")
  checkPtrType(object, "GtkFontButton")

  w <- .RGtkCall("S_gtk_font_button_class_font_set", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkFrameClassComputeChildAllocation <-
function(object.class, object, allocation)
{
  checkPtrType(object.class, "GtkFrameClass")
  checkPtrType(object, "GtkFrame")
  allocation <- as.GtkAllocation(allocation)

  w <- .RGtkCall("S_gtk_frame_class_compute_child_allocation", object.class, object, allocation, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkHandleBoxClassChildAttached <-
function(object.class, object, child)
{
  checkPtrType(object.class, "GtkHandleBoxClass")
  checkPtrType(object, "GtkHandleBox")
  checkPtrType(child, "GtkWidget")

  w <- .RGtkCall("S_gtk_handle_box_class_child_attached", object.class, object, child, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkHandleBoxClassChildDetached <-
function(object.class, object, child)
{
  checkPtrType(object.class, "GtkHandleBoxClass")
  checkPtrType(object, "GtkHandleBox")
  checkPtrType(child, "GtkWidget")

  w <- .RGtkCall("S_gtk_handle_box_class_child_detached", object.class, object, child, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkIconThemeClassChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkIconThemeClass")
  checkPtrType(object, "GtkIconTheme")

  w <- .RGtkCall("S_gtk_icon_theme_class_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkIconViewClassSetScrollAdjustments <-
function(object.class, object, hadjustment, vadjustment)
{
  checkPtrType(object.class, "GtkIconViewClass")
  checkPtrType(object, "GtkIconView")
  checkPtrType(hadjustment, "GtkAdjustment")
  checkPtrType(vadjustment, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_icon_view_class_set_scroll_adjustments", object.class, object, hadjustment, vadjustment, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkIconViewClassItemActivated <-
function(object.class, object, path)
{
  checkPtrType(object.class, "GtkIconViewClass")
  checkPtrType(object, "GtkIconView")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_icon_view_class_item_activated", object.class, object, path, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkIconViewClassSelectionChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkIconViewClass")
  checkPtrType(object, "GtkIconView")

  w <- .RGtkCall("S_gtk_icon_view_class_selection_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkIconViewClassSelectAll <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkIconViewClass")
  checkPtrType(object, "GtkIconView")

  w <- .RGtkCall("S_gtk_icon_view_class_select_all", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkIconViewClassUnselectAll <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkIconViewClass")
  checkPtrType(object, "GtkIconView")

  w <- .RGtkCall("S_gtk_icon_view_class_unselect_all", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkIconViewClassSelectCursorItem <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkIconViewClass")
  checkPtrType(object, "GtkIconView")

  w <- .RGtkCall("S_gtk_icon_view_class_select_cursor_item", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkIconViewClassToggleCursorItem <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkIconViewClass")
  checkPtrType(object, "GtkIconView")

  w <- .RGtkCall("S_gtk_icon_view_class_toggle_cursor_item", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkIconViewClassMoveCursor <-
function(object.class, object, step, count)
{
  checkPtrType(object.class, "GtkIconViewClass")
  checkPtrType(object, "GtkIconView")
  
  count <- as.integer(count)

  w <- .RGtkCall("S_gtk_icon_view_class_move_cursor", object.class, object, step, count, PACKAGE = "RGtk2")

  return(w)
}

gtkIconViewClassActivateCursorItem <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkIconViewClass")
  checkPtrType(object, "GtkIconView")

  w <- .RGtkCall("S_gtk_icon_view_class_activate_cursor_item", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkIMContextClassPreeditStart <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkIMContextClass")
  checkPtrType(object, "GtkIMContext")

  w <- .RGtkCall("S_gtk_imcontext_class_preedit_start", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkIMContextClassPreeditEnd <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkIMContextClass")
  checkPtrType(object, "GtkIMContext")

  w <- .RGtkCall("S_gtk_imcontext_class_preedit_end", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkIMContextClassPreeditChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkIMContextClass")
  checkPtrType(object, "GtkIMContext")

  w <- .RGtkCall("S_gtk_imcontext_class_preedit_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkIMContextClassCommit <-
function(object.class, object, str)
{
  checkPtrType(object.class, "GtkIMContextClass")
  checkPtrType(object, "GtkIMContext")
  str <- as.character(str)

  w <- .RGtkCall("S_gtk_imcontext_class_commit", object.class, object, str, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkIMContextClassRetrieveSurrounding <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkIMContextClass")
  checkPtrType(object, "GtkIMContext")

  w <- .RGtkCall("S_gtk_imcontext_class_retrieve_surrounding", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkIMContextClassDeleteSurrounding <-
function(object.class, object, offset, n.chars)
{
  checkPtrType(object.class, "GtkIMContextClass")
  checkPtrType(object, "GtkIMContext")
  offset <- as.integer(offset)
  n.chars <- as.integer(n.chars)

  w <- .RGtkCall("S_gtk_imcontext_class_delete_surrounding", object.class, object, offset, n.chars, PACKAGE = "RGtk2")

  return(w)
}

gtkIMContextClassSetClientWindow <-
function(object.class, object, window)
{
  checkPtrType(object.class, "GtkIMContextClass")
  checkPtrType(object, "GtkIMContext")
  checkPtrType(window, "GdkWindow")

  w <- .RGtkCall("S_gtk_imcontext_class_set_client_window", object.class, object, window, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkIMContextClassGetPreeditString <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkIMContextClass")
  checkPtrType(object, "GtkIMContext")

  w <- .RGtkCall("S_gtk_imcontext_class_get_preedit_string", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkIMContextClassFilterKeypress <-
function(object.class, object, event)
{
  checkPtrType(object.class, "GtkIMContextClass")
  checkPtrType(object, "GtkIMContext")
  checkPtrType(event, "GdkEventKey")

  w <- .RGtkCall("S_gtk_imcontext_class_filter_keypress", object.class, object, event, PACKAGE = "RGtk2")

  return(w)
}

gtkIMContextClassFocusIn <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkIMContextClass")
  checkPtrType(object, "GtkIMContext")

  w <- .RGtkCall("S_gtk_imcontext_class_focus_in", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkIMContextClassFocusOut <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkIMContextClass")
  checkPtrType(object, "GtkIMContext")

  w <- .RGtkCall("S_gtk_imcontext_class_focus_out", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkIMContextClassReset <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkIMContextClass")
  checkPtrType(object, "GtkIMContext")

  w <- .RGtkCall("S_gtk_imcontext_class_reset", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkIMContextClassSetCursorLocation <-
function(object.class, object, area)
{
  checkPtrType(object.class, "GtkIMContextClass")
  checkPtrType(object, "GtkIMContext")
  area <- as.GdkRectangle(area)

  w <- .RGtkCall("S_gtk_imcontext_class_set_cursor_location", object.class, object, area, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkIMContextClassSetUsePreedit <-
function(object.class, object, use.preedit)
{
  checkPtrType(object.class, "GtkIMContextClass")
  checkPtrType(object, "GtkIMContext")
  use.preedit <- as.logical(use.preedit)

  w <- .RGtkCall("S_gtk_imcontext_class_set_use_preedit", object.class, object, use.preedit, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkIMContextClassSetSurrounding <-
function(object.class, object, text, len, cursor.index)
{
  checkPtrType(object.class, "GtkIMContextClass")
  checkPtrType(object, "GtkIMContext")
  text <- as.character(text)
  len <- as.integer(len)
  cursor.index <- as.integer(cursor.index)

  w <- .RGtkCall("S_gtk_imcontext_class_set_surrounding", object.class, object, text, len, cursor.index, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkIMContextClassGetSurrounding <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkIMContextClass")
  checkPtrType(object, "GtkIMContext")

  w <- .RGtkCall("S_gtk_imcontext_class_get_surrounding", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkInputDialogClassEnableDevice <-
function(object.class, object, device)
{
  checkPtrType(object.class, "GtkInputDialogClass")
  checkPtrType(object, "GtkInputDialog")
  checkPtrType(device, "GdkDevice")

  w <- .RGtkCall("S_gtk_input_dialog_class_enable_device", object.class, object, device, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkInputDialogClassDisableDevice <-
function(object.class, object, device)
{
  checkPtrType(object.class, "GtkInputDialogClass")
  checkPtrType(object, "GtkInputDialog")
  checkPtrType(device, "GdkDevice")

  w <- .RGtkCall("S_gtk_input_dialog_class_disable_device", object.class, object, device, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkItemClassSelect <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkItemClass")
  checkPtrType(object, "GtkItem")

  w <- .RGtkCall("S_gtk_item_class_select", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkItemClassDeselect <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkItemClass")
  checkPtrType(object, "GtkItem")

  w <- .RGtkCall("S_gtk_item_class_deselect", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkItemClassToggle <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkItemClass")
  checkPtrType(object, "GtkItem")

  w <- .RGtkCall("S_gtk_item_class_toggle", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkLabelClassMoveCursor <-
function(object.class, object, step, count, extend.selection)
{
  checkPtrType(object.class, "GtkLabelClass")
  checkPtrType(object, "GtkLabel")
  
  count <- as.integer(count)
  extend.selection <- as.logical(extend.selection)

  w <- .RGtkCall("S_gtk_label_class_move_cursor", object.class, object, step, count, extend.selection, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkLabelClassCopyClipboard <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkLabelClass")
  checkPtrType(object, "GtkLabel")

  w <- .RGtkCall("S_gtk_label_class_copy_clipboard", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkLabelClassPopulatePopup <-
function(object.class, object, menu)
{
  checkPtrType(object.class, "GtkLabelClass")
  checkPtrType(object, "GtkLabel")
  checkPtrType(menu, "GtkMenu")

  w <- .RGtkCall("S_gtk_label_class_populate_popup", object.class, object, menu, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkLayoutClassSetScrollAdjustments <-
function(object.class, object, hadjustment, vadjustment)
{
  checkPtrType(object.class, "GtkLayoutClass")
  checkPtrType(object, "GtkLayout")
  checkPtrType(hadjustment, "GtkAdjustment")
  checkPtrType(vadjustment, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_layout_class_set_scroll_adjustments", object.class, object, hadjustment, vadjustment, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkListClassSelectionChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkListClass")
  checkPtrType(object, "GtkList")

  w <- .RGtkCall("S_gtk_list_class_selection_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkListClassSelectChild <-
function(object.class, object, child)
{
  checkPtrType(object.class, "GtkListClass")
  checkPtrType(object, "GtkList")
  checkPtrType(child, "GtkWidget")

  w <- .RGtkCall("S_gtk_list_class_select_child", object.class, object, child, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkListClassUnselectChild <-
function(object.class, object, child)
{
  checkPtrType(object.class, "GtkListClass")
  checkPtrType(object, "GtkList")
  checkPtrType(child, "GtkWidget")

  w <- .RGtkCall("S_gtk_list_class_unselect_child", object.class, object, child, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkListItemClassToggleFocusRow <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkListItemClass")
  checkPtrType(object, "GtkListItem")

  w <- .RGtkCall("S_gtk_list_item_class_toggle_focus_row", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkListItemClassSelectAll <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkListItemClass")
  checkPtrType(object, "GtkListItem")

  w <- .RGtkCall("S_gtk_list_item_class_select_all", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkListItemClassUnselectAll <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkListItemClass")
  checkPtrType(object, "GtkListItem")

  w <- .RGtkCall("S_gtk_list_item_class_unselect_all", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkListItemClassUndoSelection <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkListItemClass")
  checkPtrType(object, "GtkListItem")

  w <- .RGtkCall("S_gtk_list_item_class_undo_selection", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkListItemClassStartSelection <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkListItemClass")
  checkPtrType(object, "GtkListItem")

  w <- .RGtkCall("S_gtk_list_item_class_start_selection", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkListItemClassEndSelection <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkListItemClass")
  checkPtrType(object, "GtkListItem")

  w <- .RGtkCall("S_gtk_list_item_class_end_selection", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkListItemClassExtendSelection <-
function(object.class, object, scroll.type, position, auto.start.selection)
{
  checkPtrType(object.class, "GtkListItemClass")
  checkPtrType(object, "GtkListItem")
  
  position <- as.numeric(position)
  auto.start.selection <- as.logical(auto.start.selection)

  w <- .RGtkCall("S_gtk_list_item_class_extend_selection", object.class, object, scroll.type, position, auto.start.selection, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkListItemClassScrollHorizontal <-
function(object.class, object, scroll.type, position)
{
  checkPtrType(object.class, "GtkListItemClass")
  checkPtrType(object, "GtkListItem")
  
  position <- as.numeric(position)

  w <- .RGtkCall("S_gtk_list_item_class_scroll_horizontal", object.class, object, scroll.type, position, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkListItemClassScrollVertical <-
function(object.class, object, scroll.type, position)
{
  checkPtrType(object.class, "GtkListItemClass")
  checkPtrType(object, "GtkListItem")
  
  position <- as.numeric(position)

  w <- .RGtkCall("S_gtk_list_item_class_scroll_vertical", object.class, object, scroll.type, position, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkListItemClassToggleAddMode <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkListItemClass")
  checkPtrType(object, "GtkListItem")

  w <- .RGtkCall("S_gtk_list_item_class_toggle_add_mode", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkMenuItemClassActivate <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkMenuItemClass")
  checkPtrType(object, "GtkMenuItem")

  w <- .RGtkCall("S_gtk_menu_item_class_activate", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkMenuItemClassActivateItem <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkMenuItemClass")
  checkPtrType(object, "GtkMenuItem")

  w <- .RGtkCall("S_gtk_menu_item_class_activate_item", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkMenuItemClassToggleSizeRequest <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkMenuItemClass")
  checkPtrType(object, "GtkMenuItem")

  w <- .RGtkCall("S_gtk_menu_item_class_toggle_size_request", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkMenuItemClassToggleSizeAllocate <-
function(object.class, object, allocation)
{
  checkPtrType(object.class, "GtkMenuItemClass")
  checkPtrType(object, "GtkMenuItem")
  allocation <- as.integer(allocation)

  w <- .RGtkCall("S_gtk_menu_item_class_toggle_size_allocate", object.class, object, allocation, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkMenuShellClassDeactivate <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkMenuShellClass")
  checkPtrType(object, "GtkMenuShell")

  w <- .RGtkCall("S_gtk_menu_shell_class_deactivate", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkMenuShellClassSelectionDone <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkMenuShellClass")
  checkPtrType(object, "GtkMenuShell")

  w <- .RGtkCall("S_gtk_menu_shell_class_selection_done", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkMenuShellClassMoveCurrent <-
function(object.class, object, direction)
{
  checkPtrType(object.class, "GtkMenuShellClass")
  checkPtrType(object, "GtkMenuShell")
  

  w <- .RGtkCall("S_gtk_menu_shell_class_move_current", object.class, object, direction, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkMenuShellClassActivateCurrent <-
function(object.class, object, force.hide)
{
  checkPtrType(object.class, "GtkMenuShellClass")
  checkPtrType(object, "GtkMenuShell")
  force.hide <- as.logical(force.hide)

  w <- .RGtkCall("S_gtk_menu_shell_class_activate_current", object.class, object, force.hide, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkMenuShellClassCancel <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkMenuShellClass")
  checkPtrType(object, "GtkMenuShell")

  w <- .RGtkCall("S_gtk_menu_shell_class_cancel", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkMenuShellClassSelectItem <-
function(object.class, object, menu.item)
{
  checkPtrType(object.class, "GtkMenuShellClass")
  checkPtrType(object, "GtkMenuShell")
  checkPtrType(menu.item, "GtkWidget")

  w <- .RGtkCall("S_gtk_menu_shell_class_select_item", object.class, object, menu.item, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkMenuShellClassInsert <-
function(object.class, object, child, position)
{
  checkPtrType(object.class, "GtkMenuShellClass")
  checkPtrType(object, "GtkMenuShell")
  checkPtrType(child, "GtkWidget")
  position <- as.integer(position)

  w <- .RGtkCall("S_gtk_menu_shell_class_insert", object.class, object, child, position, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkMenuShellClassGetPopupDelay <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkMenuShellClass")
  checkPtrType(object, "GtkMenuShell")

  w <- .RGtkCall("S_gtk_menu_shell_class_get_popup_delay", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkMenuToolButtonClassShowMenu <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkMenuToolButtonClass")
  checkPtrType(object, "GtkMenuToolButton")

  w <- .RGtkCall("S_gtk_menu_tool_button_class_show_menu", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkNotebookClassSwitchPage <-
function(object.class, object, page, page.num)
{
  checkPtrType(object.class, "GtkNotebookClass")
  checkPtrType(object, "GtkNotebook")
  checkPtrType(page, "GtkNotebookPage")
  page.num <- as.numeric(page.num)

  w <- .RGtkCall("S_gtk_notebook_class_switch_page", object.class, object, page, page.num, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkNotebookClassSelectPage <-
function(object.class, object, move.focus)
{
  checkPtrType(object.class, "GtkNotebookClass")
  checkPtrType(object, "GtkNotebook")
  move.focus <- as.logical(move.focus)

  w <- .RGtkCall("S_gtk_notebook_class_select_page", object.class, object, move.focus, PACKAGE = "RGtk2")

  return(w)
}

gtkNotebookClassFocusTab <-
function(object.class, object, type)
{
  checkPtrType(object.class, "GtkNotebookClass")
  checkPtrType(object, "GtkNotebook")
  

  w <- .RGtkCall("S_gtk_notebook_class_focus_tab", object.class, object, type, PACKAGE = "RGtk2")

  return(w)
}

gtkNotebookClassChangeCurrentPage <-
function(object.class, object, offset)
{
  checkPtrType(object.class, "GtkNotebookClass")
  checkPtrType(object, "GtkNotebook")
  offset <- as.integer(offset)

  w <- .RGtkCall("S_gtk_notebook_class_change_current_page", object.class, object, offset, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkNotebookClassMoveFocusOut <-
function(object.class, object, direction)
{
  checkPtrType(object.class, "GtkNotebookClass")
  checkPtrType(object, "GtkNotebook")
  

  w <- .RGtkCall("S_gtk_notebook_class_move_focus_out", object.class, object, direction, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkOldEditableClassActivate <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkOldEditableClass")
  checkPtrType(object, "GtkOldEditable")

  w <- .RGtkCall("S_gtk_old_editable_class_activate", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkOldEditableClassSetEditable <-
function(object.class, object, is.editable)
{
  checkPtrType(object.class, "GtkOldEditableClass")
  checkPtrType(object, "GtkOldEditable")
  is.editable <- as.logical(is.editable)

  w <- .RGtkCall("S_gtk_old_editable_class_set_editable", object.class, object, is.editable, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkOldEditableClassMoveCursor <-
function(object.class, object, x, y)
{
  checkPtrType(object.class, "GtkOldEditableClass")
  checkPtrType(object, "GtkOldEditable")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_gtk_old_editable_class_move_cursor", object.class, object, x, y, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkOldEditableClassMoveWord <-
function(object.class, object, n)
{
  checkPtrType(object.class, "GtkOldEditableClass")
  checkPtrType(object, "GtkOldEditable")
  n <- as.integer(n)

  w <- .RGtkCall("S_gtk_old_editable_class_move_word", object.class, object, n, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkOldEditableClassMovePage <-
function(object.class, object, x, y)
{
  checkPtrType(object.class, "GtkOldEditableClass")
  checkPtrType(object, "GtkOldEditable")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_gtk_old_editable_class_move_page", object.class, object, x, y, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkOldEditableClassMoveToRow <-
function(object.class, object, row)
{
  checkPtrType(object.class, "GtkOldEditableClass")
  checkPtrType(object, "GtkOldEditable")
  row <- as.integer(row)

  w <- .RGtkCall("S_gtk_old_editable_class_move_to_row", object.class, object, row, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkOldEditableClassMoveToColumn <-
function(object.class, object, row)
{
  checkPtrType(object.class, "GtkOldEditableClass")
  checkPtrType(object, "GtkOldEditable")
  row <- as.integer(row)

  w <- .RGtkCall("S_gtk_old_editable_class_move_to_column", object.class, object, row, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkOldEditableClassKillChar <-
function(object.class, object, direction)
{
  checkPtrType(object.class, "GtkOldEditableClass")
  checkPtrType(object, "GtkOldEditable")
  direction <- as.integer(direction)

  w <- .RGtkCall("S_gtk_old_editable_class_kill_char", object.class, object, direction, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkOldEditableClassKillWord <-
function(object.class, object, direction)
{
  checkPtrType(object.class, "GtkOldEditableClass")
  checkPtrType(object, "GtkOldEditable")
  direction <- as.integer(direction)

  w <- .RGtkCall("S_gtk_old_editable_class_kill_word", object.class, object, direction, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkOldEditableClassKillLine <-
function(object.class, object, direction)
{
  checkPtrType(object.class, "GtkOldEditableClass")
  checkPtrType(object, "GtkOldEditable")
  direction <- as.integer(direction)

  w <- .RGtkCall("S_gtk_old_editable_class_kill_line", object.class, object, direction, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkOldEditableClassCutClipboard <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkOldEditableClass")
  checkPtrType(object, "GtkOldEditable")

  w <- .RGtkCall("S_gtk_old_editable_class_cut_clipboard", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkOldEditableClassCopyClipboard <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkOldEditableClass")
  checkPtrType(object, "GtkOldEditable")

  w <- .RGtkCall("S_gtk_old_editable_class_copy_clipboard", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkOldEditableClassPasteClipboard <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkOldEditableClass")
  checkPtrType(object, "GtkOldEditable")

  w <- .RGtkCall("S_gtk_old_editable_class_paste_clipboard", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkOldEditableClassUpdateText <-
function(object.class, object, start.pos, end.pos)
{
  checkPtrType(object.class, "GtkOldEditableClass")
  checkPtrType(object, "GtkOldEditable")
  start.pos <- as.integer(start.pos)
  end.pos <- as.integer(end.pos)

  w <- .RGtkCall("S_gtk_old_editable_class_update_text", object.class, object, start.pos, end.pos, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkOldEditableClassGetChars <-
function(object.class, object, start.pos, end.pos)
{
  checkPtrType(object.class, "GtkOldEditableClass")
  checkPtrType(object, "GtkOldEditable")
  start.pos <- as.integer(start.pos)
  end.pos <- as.integer(end.pos)

  w <- .RGtkCall("S_gtk_old_editable_class_get_chars", object.class, object, start.pos, end.pos, PACKAGE = "RGtk2")

  return(w)
}

gtkOldEditableClassSetSelection <-
function(object.class, object, start.pos, end.pos)
{
  checkPtrType(object.class, "GtkOldEditableClass")
  checkPtrType(object, "GtkOldEditable")
  start.pos <- as.integer(start.pos)
  end.pos <- as.integer(end.pos)

  w <- .RGtkCall("S_gtk_old_editable_class_set_selection", object.class, object, start.pos, end.pos, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkOldEditableClassSetPosition <-
function(object.class, object, position)
{
  checkPtrType(object.class, "GtkOldEditableClass")
  checkPtrType(object, "GtkOldEditable")
  position <- as.integer(position)

  w <- .RGtkCall("S_gtk_old_editable_class_set_position", object.class, object, position, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkOptionMenuClassChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkOptionMenuClass")
  checkPtrType(object, "GtkOptionMenu")

  w <- .RGtkCall("S_gtk_option_menu_class_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkPanedClassCycleChildFocus <-
function(object.class, object, reverse)
{
  checkPtrType(object.class, "GtkPanedClass")
  checkPtrType(object, "GtkPaned")
  reverse <- as.logical(reverse)

  w <- .RGtkCall("S_gtk_paned_class_cycle_child_focus", object.class, object, reverse, PACKAGE = "RGtk2")

  return(w)
}

gtkPanedClassToggleHandleFocus <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkPanedClass")
  checkPtrType(object, "GtkPaned")

  w <- .RGtkCall("S_gtk_paned_class_toggle_handle_focus", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkPanedClassMoveHandle <-
function(object.class, object, scroll)
{
  checkPtrType(object.class, "GtkPanedClass")
  checkPtrType(object, "GtkPaned")
  

  w <- .RGtkCall("S_gtk_paned_class_move_handle", object.class, object, scroll, PACKAGE = "RGtk2")

  return(w)
}

gtkPanedClassCycleHandleFocus <-
function(object.class, object, reverse)
{
  checkPtrType(object.class, "GtkPanedClass")
  checkPtrType(object, "GtkPaned")
  reverse <- as.logical(reverse)

  w <- .RGtkCall("S_gtk_paned_class_cycle_handle_focus", object.class, object, reverse, PACKAGE = "RGtk2")

  return(w)
}

gtkPanedClassAcceptPosition <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkPanedClass")
  checkPtrType(object, "GtkPaned")

  w <- .RGtkCall("S_gtk_paned_class_accept_position", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkPanedClassCancelPosition <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkPanedClass")
  checkPtrType(object, "GtkPaned")

  w <- .RGtkCall("S_gtk_paned_class_cancel_position", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkPlugClassEmbedded <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkPlugClass")
  checkPtrType(object, "GtkPlug")

  w <- .RGtkCall("S_gtk_plug_class_embedded", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkProgressClassPaint <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkProgressClass")
  checkPtrType(object, "GtkProgress")

  w <- .RGtkCall("S_gtk_progress_class_paint", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkProgressClassUpdate <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkProgressClass")
  checkPtrType(object, "GtkProgress")

  w <- .RGtkCall("S_gtk_progress_class_update", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkProgressClassActModeEnter <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkProgressClass")
  checkPtrType(object, "GtkProgress")

  w <- .RGtkCall("S_gtk_progress_class_act_mode_enter", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkRadioActionClassChanged <-
function(object.class, object, current)
{
  checkPtrType(object.class, "GtkRadioActionClass")
  checkPtrType(object, "GtkRadioAction")
  checkPtrType(current, "GtkRadioAction")

  w <- .RGtkCall("S_gtk_radio_action_class_changed", object.class, object, current, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkRadioButtonClassGroupChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkRadioButtonClass")
  checkPtrType(object, "GtkRadioButton")

  w <- .RGtkCall("S_gtk_radio_button_class_group_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkRadioMenuItemClassGroupChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkRadioMenuItemClass")
  checkPtrType(object, "GtkRadioMenuItem")

  w <- .RGtkCall("S_gtk_radio_menu_item_class_group_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkRangeClassValueChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkRangeClass")
  checkPtrType(object, "GtkRange")

  w <- .RGtkCall("S_gtk_range_class_value_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkRangeClassAdjustBounds <-
function(object.class, object, new.value)
{
  checkPtrType(object.class, "GtkRangeClass")
  checkPtrType(object, "GtkRange")
  new.value <- as.numeric(new.value)

  w <- .RGtkCall("S_gtk_range_class_adjust_bounds", object.class, object, new.value, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkRangeClassMoveSlider <-
function(object.class, object, scroll)
{
  checkPtrType(object.class, "GtkRangeClass")
  checkPtrType(object, "GtkRange")
  

  w <- .RGtkCall("S_gtk_range_class_move_slider", object.class, object, scroll, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkRangeClassGetRangeBorder <-
function(object.class, object, border.)
{
  checkPtrType(object.class, "GtkRangeClass")
  checkPtrType(object, "GtkRange")
  checkPtrType(border., "GtkBorder")

  w <- .RGtkCall("S_gtk_range_class_get_range_border", object.class, object, border., PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkRangeClassChangeValue <-
function(object.class, object, scroll, new.value)
{
  checkPtrType(object.class, "GtkRangeClass")
  checkPtrType(object, "GtkRange")
  
  new.value <- as.numeric(new.value)

  w <- .RGtkCall("S_gtk_range_class_change_value", object.class, object, scroll, new.value, PACKAGE = "RGtk2")

  return(w)
}

gtkRcStyleClassCreateRcStyle <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkRcStyleClass")
  checkPtrType(object, "GtkRcStyle")

  w <- .RGtkCall("S_gtk_rc_style_class_create_rc_style", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkRcStyleClassParse <-
function(object.class, object, settings, scanner)
{
  checkPtrType(object.class, "GtkRcStyleClass")
  checkPtrType(object, "GtkRcStyle")
  checkPtrType(settings, "GtkSettings")
  checkPtrType(scanner, "GScanner")

  w <- .RGtkCall("S_gtk_rc_style_class_parse", object.class, object, settings, scanner, PACKAGE = "RGtk2")

  return(w)
}

gtkRcStyleClassMerge <-
function(object.class, object, src)
{
  checkPtrType(object.class, "GtkRcStyleClass")
  checkPtrType(object, "GtkRcStyle")
  checkPtrType(src, "GtkRcStyle")

  w <- .RGtkCall("S_gtk_rc_style_class_merge", object.class, object, src, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkRcStyleClassCreateStyle <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkRcStyleClass")
  checkPtrType(object, "GtkRcStyle")

  w <- .RGtkCall("S_gtk_rc_style_class_create_style", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkRulerClassDrawTicks <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkRulerClass")
  checkPtrType(object, "GtkRuler")

  w <- .RGtkCall("S_gtk_ruler_class_draw_ticks", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkRulerClassDrawPos <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkRulerClass")
  checkPtrType(object, "GtkRuler")

  w <- .RGtkCall("S_gtk_ruler_class_draw_pos", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkScaleClassFormatValue <-
function(object.class, object, value)
{
  checkPtrType(object.class, "GtkScaleClass")
  checkPtrType(object, "GtkScale")
  value <- as.numeric(value)

  w <- .RGtkCall("S_gtk_scale_class_format_value", object.class, object, value, PACKAGE = "RGtk2")

  return(w)
}

gtkScaleClassDrawValue <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkScaleClass")
  checkPtrType(object, "GtkScale")

  w <- .RGtkCall("S_gtk_scale_class_draw_value", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkScaleClassGetLayoutOffsets <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkScaleClass")
  checkPtrType(object, "GtkScale")

  w <- .RGtkCall("S_gtk_scale_class_get_layout_offsets", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkScrolledWindowClassScrollChild <-
function(object.class, object, scroll, horizontal)
{
  checkPtrType(object.class, "GtkScrolledWindowClass")
  checkPtrType(object, "GtkScrolledWindow")
  
  horizontal <- as.logical(horizontal)

  w <- .RGtkCall("S_gtk_scrolled_window_class_scroll_child", object.class, object, scroll, horizontal, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkScrolledWindowClassMoveFocusOut <-
function(object.class, object, direction)
{
  checkPtrType(object.class, "GtkScrolledWindowClass")
  checkPtrType(object, "GtkScrolledWindow")
  

  w <- .RGtkCall("S_gtk_scrolled_window_class_move_focus_out", object.class, object, direction, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkSocketClassPlugAdded <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkSocketClass")
  checkPtrType(object, "GtkSocket")

  w <- .RGtkCall("S_gtk_socket_class_plug_added", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkSocketClassPlugRemoved <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkSocketClass")
  checkPtrType(object, "GtkSocket")

  w <- .RGtkCall("S_gtk_socket_class_plug_removed", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkSpinButtonClassInput <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkSpinButtonClass")
  checkPtrType(object, "GtkSpinButton")

  w <- .RGtkCall("S_gtk_spin_button_class_input", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkSpinButtonClassOutput <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkSpinButtonClass")
  checkPtrType(object, "GtkSpinButton")

  w <- .RGtkCall("S_gtk_spin_button_class_output", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkSpinButtonClassValueChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkSpinButtonClass")
  checkPtrType(object, "GtkSpinButton")

  w <- .RGtkCall("S_gtk_spin_button_class_value_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkSpinButtonClassChangeValue <-
function(object.class, object, scroll)
{
  checkPtrType(object.class, "GtkSpinButtonClass")
  checkPtrType(object, "GtkSpinButton")
  

  w <- .RGtkCall("S_gtk_spin_button_class_change_value", object.class, object, scroll, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkStatusbarClassTextPushed <-
function(object.class, object, context.id, text)
{
  checkPtrType(object.class, "GtkStatusbarClass")
  checkPtrType(object, "GtkStatusbar")
  context.id <- as.numeric(context.id)
  text <- as.character(text)

  w <- .RGtkCall("S_gtk_statusbar_class_text_pushed", object.class, object, context.id, text, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkStatusbarClassTextPopped <-
function(object.class, object, context.id, text)
{
  checkPtrType(object.class, "GtkStatusbarClass")
  checkPtrType(object, "GtkStatusbar")
  context.id <- as.numeric(context.id)
  text <- as.character(text)

  w <- .RGtkCall("S_gtk_statusbar_class_text_popped", object.class, object, context.id, text, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkStyleClassRealize <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkStyleClass")
  checkPtrType(object, "GtkStyle")

  w <- .RGtkCall("S_gtk_style_class_realize", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkStyleClassUnrealize <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkStyleClass")
  checkPtrType(object, "GtkStyle")

  w <- .RGtkCall("S_gtk_style_class_unrealize", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkStyleClassCopy <-
function(object.class, object, src)
{
  checkPtrType(object.class, "GtkStyleClass")
  checkPtrType(object, "GtkStyle")
  checkPtrType(src, "GtkStyle")

  w <- .RGtkCall("S_gtk_style_class_copy", object.class, object, src, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkStyleClassClone <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkStyleClass")
  checkPtrType(object, "GtkStyle")

  w <- .RGtkCall("S_gtk_style_class_clone", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkStyleClassInitFromRc <-
function(object.class, object, rc.style)
{
  checkPtrType(object.class, "GtkStyleClass")
  checkPtrType(object, "GtkStyle")
  checkPtrType(rc.style, "GtkRcStyle")

  w <- .RGtkCall("S_gtk_style_class_init_from_rc", object.class, object, rc.style, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkStyleClassSetBackground <-
function(object.class, object, window, state.type)
{
  checkPtrType(object.class, "GtkStyleClass")
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  

  w <- .RGtkCall("S_gtk_style_class_set_background", object.class, object, window, state.type, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkStyleClassRenderIcon <-
function(object.class, object, source, direction, state, size, widget, detail)
{
  checkPtrType(object.class, "GtkStyleClass")
  checkPtrType(object, "GtkStyle")
  checkPtrType(source, "GtkIconSource")
  
  
  
  checkPtrType(widget, "GtkWidget")
  detail <- as.character(detail)

  w <- .RGtkCall("S_gtk_style_class_render_icon", object.class, object, source, direction, state, size, widget, detail, PACKAGE = "RGtk2")

  return(w)
}

gtkStyleClassDrawHline <-
function(object.class, object, window, state.type, area, widget, detail, x1, x2, y)
{
  checkPtrType(object.class, "GtkStyleClass")
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  area <- as.GdkRectangle(area)
  checkPtrType(widget, "GtkWidget")
  detail <- as.character(detail)
  x1 <- as.integer(x1)
  x2 <- as.integer(x2)
  y <- as.integer(y)

  w <- .RGtkCall("S_gtk_style_class_draw_hline", object.class, object, window, state.type, area, widget, detail, x1, x2, y, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkStyleClassDrawVline <-
function(object.class, object, window, state.type, area, widget, detail, y1., y2., x)
{
  checkPtrType(object.class, "GtkStyleClass")
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  area <- as.GdkRectangle(area)
  checkPtrType(widget, "GtkWidget")
  detail <- as.character(detail)
  y1. <- as.integer(y1.)
  y2. <- as.integer(y2.)
  x <- as.integer(x)

  w <- .RGtkCall("S_gtk_style_class_draw_vline", object.class, object, window, state.type, area, widget, detail, y1., y2., x, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkStyleClassDrawShadow <-
function(object.class, object, window, state.type, shadow.type, area, widget, detail, x, y, width, height)
{
  checkPtrType(object.class, "GtkStyleClass")
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  area <- as.GdkRectangle(area)
  checkPtrType(widget, "GtkWidget")
  detail <- as.character(detail)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_style_class_draw_shadow", object.class, object, window, state.type, shadow.type, area, widget, detail, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkStyleClassDrawPolygon <-
function(object.class, object, window, state.type, shadow.type, area, widget, detail, point, npoints, fill)
{
  checkPtrType(object.class, "GtkStyleClass")
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  area <- as.GdkRectangle(area)
  checkPtrType(widget, "GtkWidget")
  detail <- as.character(detail)
  point <- as.GdkPoint(point)
  npoints <- as.integer(npoints)
  fill <- as.logical(fill)

  w <- .RGtkCall("S_gtk_style_class_draw_polygon", object.class, object, window, state.type, shadow.type, area, widget, detail, point, npoints, fill, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkStyleClassDrawArrow <-
function(object.class, object, window, state.type, shadow.type, area, widget, detail, arrow.type, fill, x, y, width, height)
{
  checkPtrType(object.class, "GtkStyleClass")
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  area <- as.GdkRectangle(area)
  checkPtrType(widget, "GtkWidget")
  detail <- as.character(detail)
  
  fill <- as.logical(fill)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_style_class_draw_arrow", object.class, object, window, state.type, shadow.type, area, widget, detail, arrow.type, fill, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkStyleClassDrawDiamond <-
function(object.class, object, window, state.type, shadow.type, area, widget, detail, x, y, width, height)
{
  checkPtrType(object.class, "GtkStyleClass")
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  area <- as.GdkRectangle(area)
  checkPtrType(widget, "GtkWidget")
  detail <- as.character(detail)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_style_class_draw_diamond", object.class, object, window, state.type, shadow.type, area, widget, detail, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkStyleClassDrawString <-
function(object.class, object, window, state.type, area, widget, detail, x, y, string)
{
  checkPtrType(object.class, "GtkStyleClass")
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  area <- as.GdkRectangle(area)
  checkPtrType(widget, "GtkWidget")
  detail <- as.character(detail)
  x <- as.integer(x)
  y <- as.integer(y)
  string <- as.character(string)

  w <- .RGtkCall("S_gtk_style_class_draw_string", object.class, object, window, state.type, area, widget, detail, x, y, string, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkStyleClassDrawBox <-
function(object.class, object, window, state.type, shadow.type, area, widget, detail, x, y, width, height)
{
  checkPtrType(object.class, "GtkStyleClass")
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  area <- as.GdkRectangle(area)
  checkPtrType(widget, "GtkWidget")
  detail <- as.character(detail)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_style_class_draw_box", object.class, object, window, state.type, shadow.type, area, widget, detail, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkStyleClassDrawFlatBox <-
function(object.class, object, window, state.type, shadow.type, area, widget, detail, x, y, width, height)
{
  checkPtrType(object.class, "GtkStyleClass")
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  area <- as.GdkRectangle(area)
  checkPtrType(widget, "GtkWidget")
  detail <- as.character(detail)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_style_class_draw_flat_box", object.class, object, window, state.type, shadow.type, area, widget, detail, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkStyleClassDrawCheck <-
function(object.class, object, window, state.type, shadow.type, area, widget, detail, x, y, width, height)
{
  checkPtrType(object.class, "GtkStyleClass")
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  area <- as.GdkRectangle(area)
  checkPtrType(widget, "GtkWidget")
  detail <- as.character(detail)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_style_class_draw_check", object.class, object, window, state.type, shadow.type, area, widget, detail, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkStyleClassDrawOption <-
function(object.class, object, window, state.type, shadow.type, area, widget, detail, x, y, width, height)
{
  checkPtrType(object.class, "GtkStyleClass")
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  area <- as.GdkRectangle(area)
  checkPtrType(widget, "GtkWidget")
  detail <- as.character(detail)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_style_class_draw_option", object.class, object, window, state.type, shadow.type, area, widget, detail, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkStyleClassDrawTab <-
function(object.class, object, window, state.type, shadow.type, area, widget, detail, x, y, width, height)
{
  checkPtrType(object.class, "GtkStyleClass")
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  area <- as.GdkRectangle(area)
  checkPtrType(widget, "GtkWidget")
  detail <- as.character(detail)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_style_class_draw_tab", object.class, object, window, state.type, shadow.type, area, widget, detail, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkStyleClassDrawShadowGap <-
function(object.class, object, window, state.type, shadow.type, area, widget, detail, x, y, width, height, gap.side, gap.x, gap.width)
{
  checkPtrType(object.class, "GtkStyleClass")
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  area <- as.GdkRectangle(area)
  checkPtrType(widget, "GtkWidget")
  detail <- as.character(detail)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)
  
  gap.x <- as.integer(gap.x)
  gap.width <- as.integer(gap.width)

  w <- .RGtkCall("S_gtk_style_class_draw_shadow_gap", object.class, object, window, state.type, shadow.type, area, widget, detail, x, y, width, height, gap.side, gap.x, gap.width, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkStyleClassDrawBoxGap <-
function(object.class, object, window, state.type, shadow.type, area, widget, detail, x, y, width, height, gap.side, gap.x, gap.width)
{
  checkPtrType(object.class, "GtkStyleClass")
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  area <- as.GdkRectangle(area)
  checkPtrType(widget, "GtkWidget")
  detail <- as.character(detail)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)
  
  gap.x <- as.integer(gap.x)
  gap.width <- as.integer(gap.width)

  w <- .RGtkCall("S_gtk_style_class_draw_box_gap", object.class, object, window, state.type, shadow.type, area, widget, detail, x, y, width, height, gap.side, gap.x, gap.width, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkStyleClassDrawExtension <-
function(object.class, object, window, state.type, shadow.type, area, widget, detail, x, y, width, height, gap.side)
{
  checkPtrType(object.class, "GtkStyleClass")
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  area <- as.GdkRectangle(area)
  checkPtrType(widget, "GtkWidget")
  detail <- as.character(detail)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)
  

  w <- .RGtkCall("S_gtk_style_class_draw_extension", object.class, object, window, state.type, shadow.type, area, widget, detail, x, y, width, height, gap.side, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkStyleClassDrawFocus <-
function(object.class, object, window, state.type, area, widget, detail, x, y, width, height)
{
  checkPtrType(object.class, "GtkStyleClass")
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  area <- as.GdkRectangle(area)
  checkPtrType(widget, "GtkWidget")
  detail <- as.character(detail)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_style_class_draw_focus", object.class, object, window, state.type, area, widget, detail, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkStyleClassDrawSlider <-
function(object.class, object, window, state.type, shadow.type, area, widget, detail, x, y, width, height, orientation)
{
  checkPtrType(object.class, "GtkStyleClass")
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  area <- as.GdkRectangle(area)
  checkPtrType(widget, "GtkWidget")
  detail <- as.character(detail)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)
  

  w <- .RGtkCall("S_gtk_style_class_draw_slider", object.class, object, window, state.type, shadow.type, area, widget, detail, x, y, width, height, orientation, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkStyleClassDrawHandle <-
function(object.class, object, window, state.type, shadow.type, area, widget, detail, x, y, width, height, orientation)
{
  checkPtrType(object.class, "GtkStyleClass")
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  area <- as.GdkRectangle(area)
  checkPtrType(widget, "GtkWidget")
  detail <- as.character(detail)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)
  

  w <- .RGtkCall("S_gtk_style_class_draw_handle", object.class, object, window, state.type, shadow.type, area, widget, detail, x, y, width, height, orientation, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkStyleClassDrawExpander <-
function(object.class, object, window, state.type, area, widget, detail, x, y, expander.style)
{
  checkPtrType(object.class, "GtkStyleClass")
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  area <- as.GdkRectangle(area)
  checkPtrType(widget, "GtkWidget")
  detail <- as.character(detail)
  x <- as.integer(x)
  y <- as.integer(y)
  

  w <- .RGtkCall("S_gtk_style_class_draw_expander", object.class, object, window, state.type, area, widget, detail, x, y, expander.style, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkStyleClassDrawLayout <-
function(object.class, object, window, state.type, use.text, area, widget, detail, x, y, layout)
{
  checkPtrType(object.class, "GtkStyleClass")
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  use.text <- as.logical(use.text)
  area <- as.GdkRectangle(area)
  checkPtrType(widget, "GtkWidget")
  detail <- as.character(detail)
  x <- as.integer(x)
  y <- as.integer(y)
  checkPtrType(layout, "PangoLayout")

  w <- .RGtkCall("S_gtk_style_class_draw_layout", object.class, object, window, state.type, use.text, area, widget, detail, x, y, layout, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkStyleClassDrawResizeGrip <-
function(object.class, object, window, state.type, area, widget, detail, edge, x, y, width, height)
{
  checkPtrType(object.class, "GtkStyleClass")
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  area <- as.GdkRectangle(area)
  checkPtrType(widget, "GtkWidget")
  detail <- as.character(detail)
  
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_style_class_draw_resize_grip", object.class, object, window, state.type, area, widget, detail, edge, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTextClassSetScrollAdjustments <-
function(object.class, object, hadjustment, vadjustment)
{
  checkPtrType(object.class, "GtkTextClass")
  checkPtrType(object, "GtkText")
  checkPtrType(hadjustment, "GtkAdjustment")
  checkPtrType(vadjustment, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_text_class_set_scroll_adjustments", object.class, object, hadjustment, vadjustment, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTextBufferClassInsertText <-
function(object.class, object, pos, text, length)
{
  checkPtrType(object.class, "GtkTextBufferClass")
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(pos, "GtkTextIter")
  text <- as.character(text)
  length <- as.integer(length)

  w <- .RGtkCall("S_gtk_text_buffer_class_insert_text", object.class, object, pos, text, length, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTextBufferClassInsertPixbuf <-
function(object.class, object, pos, pixbuf)
{
  checkPtrType(object.class, "GtkTextBufferClass")
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(pos, "GtkTextIter")
  checkPtrType(pixbuf, "GdkPixbuf")

  w <- .RGtkCall("S_gtk_text_buffer_class_insert_pixbuf", object.class, object, pos, pixbuf, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTextBufferClassInsertChildAnchor <-
function(object.class, object, pos, anchor)
{
  checkPtrType(object.class, "GtkTextBufferClass")
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(pos, "GtkTextIter")
  checkPtrType(anchor, "GtkTextChildAnchor")

  w <- .RGtkCall("S_gtk_text_buffer_class_insert_child_anchor", object.class, object, pos, anchor, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTextBufferClassDeleteRange <-
function(object.class, object, start, end)
{
  checkPtrType(object.class, "GtkTextBufferClass")
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(start, "GtkTextIter")
  checkPtrType(end, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_buffer_class_delete_range", object.class, object, start, end, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTextBufferClassChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkTextBufferClass")
  checkPtrType(object, "GtkTextBuffer")

  w <- .RGtkCall("S_gtk_text_buffer_class_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTextBufferClassModifiedChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkTextBufferClass")
  checkPtrType(object, "GtkTextBuffer")

  w <- .RGtkCall("S_gtk_text_buffer_class_modified_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTextBufferClassMarkSet <-
function(object.class, object, location, mark)
{
  checkPtrType(object.class, "GtkTextBufferClass")
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(location, "GtkTextIter")
  checkPtrType(mark, "GtkTextMark")

  w <- .RGtkCall("S_gtk_text_buffer_class_mark_set", object.class, object, location, mark, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTextBufferClassMarkDeleted <-
function(object.class, object, mark)
{
  checkPtrType(object.class, "GtkTextBufferClass")
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(mark, "GtkTextMark")

  w <- .RGtkCall("S_gtk_text_buffer_class_mark_deleted", object.class, object, mark, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTextBufferClassApplyTag <-
function(object.class, object, tag, start.char, end.char)
{
  checkPtrType(object.class, "GtkTextBufferClass")
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(tag, "GtkTextTag")
  checkPtrType(start.char, "GtkTextIter")
  checkPtrType(end.char, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_buffer_class_apply_tag", object.class, object, tag, start.char, end.char, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTextBufferClassRemoveTag <-
function(object.class, object, tag, start.char, end.char)
{
  checkPtrType(object.class, "GtkTextBufferClass")
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(tag, "GtkTextTag")
  checkPtrType(start.char, "GtkTextIter")
  checkPtrType(end.char, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_buffer_class_remove_tag", object.class, object, tag, start.char, end.char, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTextBufferClassBeginUserAction <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkTextBufferClass")
  checkPtrType(object, "GtkTextBuffer")

  w <- .RGtkCall("S_gtk_text_buffer_class_begin_user_action", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTextBufferClassEndUserAction <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkTextBufferClass")
  checkPtrType(object, "GtkTextBuffer")

  w <- .RGtkCall("S_gtk_text_buffer_class_end_user_action", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTextTagClassEvent <-
function(object.class, object, event.object, event, iter)
{
  checkPtrType(object.class, "GtkTextTagClass")
  checkPtrType(object, "GtkTextTag")
  checkPtrType(event.object, "GObject")
  checkPtrType(event, "GdkEvent")
  checkPtrType(iter, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_tag_class_event", object.class, object, event.object, event, iter, PACKAGE = "RGtk2")

  return(w)
}

gtkTextTagTableClassTagChanged <-
function(object.class, object, tag, size.changed)
{
  checkPtrType(object.class, "GtkTextTagTableClass")
  checkPtrType(object, "GtkTextTagTable")
  checkPtrType(tag, "GtkTextTag")
  size.changed <- as.logical(size.changed)

  w <- .RGtkCall("S_gtk_text_tag_table_class_tag_changed", object.class, object, tag, size.changed, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTextTagTableClassTagAdded <-
function(object.class, object, tag)
{
  checkPtrType(object.class, "GtkTextTagTableClass")
  checkPtrType(object, "GtkTextTagTable")
  checkPtrType(tag, "GtkTextTag")

  w <- .RGtkCall("S_gtk_text_tag_table_class_tag_added", object.class, object, tag, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTextTagTableClassTagRemoved <-
function(object.class, object, tag)
{
  checkPtrType(object.class, "GtkTextTagTableClass")
  checkPtrType(object, "GtkTextTagTable")
  checkPtrType(tag, "GtkTextTag")

  w <- .RGtkCall("S_gtk_text_tag_table_class_tag_removed", object.class, object, tag, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTextViewClassSetScrollAdjustments <-
function(object.class, object, hadjustment, vadjustment)
{
  checkPtrType(object.class, "GtkTextViewClass")
  checkPtrType(object, "GtkTextView")
  checkPtrType(hadjustment, "GtkAdjustment")
  checkPtrType(vadjustment, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_text_view_class_set_scroll_adjustments", object.class, object, hadjustment, vadjustment, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTextViewClassPopulatePopup <-
function(object.class, object, menu)
{
  checkPtrType(object.class, "GtkTextViewClass")
  checkPtrType(object, "GtkTextView")
  checkPtrType(menu, "GtkMenu")

  w <- .RGtkCall("S_gtk_text_view_class_populate_popup", object.class, object, menu, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTextViewClassMoveCursor <-
function(object.class, object, step, count, extend.selection)
{
  checkPtrType(object.class, "GtkTextViewClass")
  checkPtrType(object, "GtkTextView")
  
  count <- as.integer(count)
  extend.selection <- as.logical(extend.selection)

  w <- .RGtkCall("S_gtk_text_view_class_move_cursor", object.class, object, step, count, extend.selection, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTextViewClassPageHorizontally <-
function(object.class, object, count, extend.selection)
{
  checkPtrType(object.class, "GtkTextViewClass")
  checkPtrType(object, "GtkTextView")
  count <- as.integer(count)
  extend.selection <- as.logical(extend.selection)

  w <- .RGtkCall("S_gtk_text_view_class_page_horizontally", object.class, object, count, extend.selection, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTextViewClassSetAnchor <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkTextViewClass")
  checkPtrType(object, "GtkTextView")

  w <- .RGtkCall("S_gtk_text_view_class_set_anchor", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTextViewClassInsertAtCursor <-
function(object.class, object, str)
{
  checkPtrType(object.class, "GtkTextViewClass")
  checkPtrType(object, "GtkTextView")
  str <- as.character(str)

  w <- .RGtkCall("S_gtk_text_view_class_insert_at_cursor", object.class, object, str, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTextViewClassDeleteFromCursor <-
function(object.class, object, type, count)
{
  checkPtrType(object.class, "GtkTextViewClass")
  checkPtrType(object, "GtkTextView")
  
  count <- as.integer(count)

  w <- .RGtkCall("S_gtk_text_view_class_delete_from_cursor", object.class, object, type, count, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTextViewClassBackspace <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkTextViewClass")
  checkPtrType(object, "GtkTextView")

  w <- .RGtkCall("S_gtk_text_view_class_backspace", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTextViewClassCutClipboard <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkTextViewClass")
  checkPtrType(object, "GtkTextView")

  w <- .RGtkCall("S_gtk_text_view_class_cut_clipboard", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTextViewClassCopyClipboard <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkTextViewClass")
  checkPtrType(object, "GtkTextView")

  w <- .RGtkCall("S_gtk_text_view_class_copy_clipboard", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTextViewClassPasteClipboard <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkTextViewClass")
  checkPtrType(object, "GtkTextView")

  w <- .RGtkCall("S_gtk_text_view_class_paste_clipboard", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTextViewClassToggleOverwrite <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkTextViewClass")
  checkPtrType(object, "GtkTextView")

  w <- .RGtkCall("S_gtk_text_view_class_toggle_overwrite", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTextViewClassMoveFocus <-
function(object.class, object, direction)
{
  checkPtrType(object.class, "GtkTextViewClass")
  checkPtrType(object, "GtkTextView")
  

  w <- .RGtkCall("S_gtk_text_view_class_move_focus", object.class, object, direction, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTipsQueryClassStartQuery <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkTipsQueryClass")
  checkPtrType(object, "GtkTipsQuery")

  w <- .RGtkCall("S_gtk_tips_query_class_start_query", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTipsQueryClassStopQuery <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkTipsQueryClass")
  checkPtrType(object, "GtkTipsQuery")

  w <- .RGtkCall("S_gtk_tips_query_class_stop_query", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTipsQueryClassWidgetEntered <-
function(object.class, object, widget, tip.text, tip.private)
{
  checkPtrType(object.class, "GtkTipsQueryClass")
  checkPtrType(object, "GtkTipsQuery")
  checkPtrType(widget, "GtkWidget")
  tip.text <- as.character(tip.text)
  tip.private <- as.character(tip.private)

  w <- .RGtkCall("S_gtk_tips_query_class_widget_entered", object.class, object, widget, tip.text, tip.private, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTipsQueryClassWidgetSelected <-
function(object.class, object, widget, tip.text, tip.private, event)
{
  checkPtrType(object.class, "GtkTipsQueryClass")
  checkPtrType(object, "GtkTipsQuery")
  checkPtrType(widget, "GtkWidget")
  tip.text <- as.character(tip.text)
  tip.private <- as.character(tip.private)
  checkPtrType(event, "GdkEventButton")

  w <- .RGtkCall("S_gtk_tips_query_class_widget_selected", object.class, object, widget, tip.text, tip.private, event, PACKAGE = "RGtk2")

  return(w)
}

gtkToggleActionClassToggled <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkToggleActionClass")
  checkPtrType(object, "GtkToggleAction")

  w <- .RGtkCall("S_gtk_toggle_action_class_toggled", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkToggleButtonClassToggled <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkToggleButtonClass")
  checkPtrType(object, "GtkToggleButton")

  w <- .RGtkCall("S_gtk_toggle_button_class_toggled", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkToggleToolButtonClassToggled <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkToggleToolButtonClass")
  checkPtrType(object, "GtkToggleToolButton")

  w <- .RGtkCall("S_gtk_toggle_tool_button_class_toggled", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkToolbarClassOrientationChanged <-
function(object.class, object, orientation)
{
  checkPtrType(object.class, "GtkToolbarClass")
  checkPtrType(object, "GtkToolbar")
  

  w <- .RGtkCall("S_gtk_toolbar_class_orientation_changed", object.class, object, orientation, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkToolbarClassStyleChanged <-
function(object.class, object, style)
{
  checkPtrType(object.class, "GtkToolbarClass")
  checkPtrType(object, "GtkToolbar")
  

  w <- .RGtkCall("S_gtk_toolbar_class_style_changed", object.class, object, style, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkToolbarClassPopupContextMenu <-
function(object.class, object, x, y, button.number)
{
  checkPtrType(object.class, "GtkToolbarClass")
  checkPtrType(object, "GtkToolbar")
  x <- as.integer(x)
  y <- as.integer(y)
  button.number <- as.integer(button.number)

  w <- .RGtkCall("S_gtk_toolbar_class_popup_context_menu", object.class, object, x, y, button.number, PACKAGE = "RGtk2")

  return(w)
}

gtkToolButtonClassClicked <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkToolButtonClass")
  checkPtrType(object, "GtkToolButton")

  w <- .RGtkCall("S_gtk_tool_button_class_clicked", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkToolItemClassCreateMenuProxy <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkToolItemClass")
  checkPtrType(object, "GtkToolItem")

  w <- .RGtkCall("S_gtk_tool_item_class_create_menu_proxy", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkToolItemClassToolbarReconfigured <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkToolItemClass")
  checkPtrType(object, "GtkToolItem")

  w <- .RGtkCall("S_gtk_tool_item_class_toolbar_reconfigured", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkToolItemClassSetTooltip <-
function(object.class, object, tooltips, tip.text, tip.private)
{
  checkPtrType(object.class, "GtkToolItemClass")
  checkPtrType(object, "GtkToolItem")
  checkPtrType(tooltips, "GtkTooltips")
  tip.text <- as.character(tip.text)
  tip.private <- as.character(tip.private)

  w <- .RGtkCall("S_gtk_tool_item_class_set_tooltip", object.class, object, tooltips, tip.text, tip.private, PACKAGE = "RGtk2")

  return(w)
}

gtkTreeDragSourceIfaceRowDraggable <-
function(object.class, object, path)
{
  checkPtrType(object.class, "GtkTreeDragSourceIface")
  checkPtrType(object, "GtkTreeDragSource")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_drag_source_iface_row_draggable", object.class, object, path, PACKAGE = "RGtk2")

  return(w)
}

gtkTreeDragSourceIfaceDragDataGet <-
function(object.class, object, path, selection.data)
{
  checkPtrType(object.class, "GtkTreeDragSourceIface")
  checkPtrType(object, "GtkTreeDragSource")
  checkPtrType(path, "GtkTreePath")
  checkPtrType(selection.data, "GtkSelectionData")

  w <- .RGtkCall("S_gtk_tree_drag_source_iface_drag_data_get", object.class, object, path, selection.data, PACKAGE = "RGtk2")

  return(w)
}

gtkTreeDragSourceIfaceDragDataDelete <-
function(object.class, object, path)
{
  checkPtrType(object.class, "GtkTreeDragSourceIface")
  checkPtrType(object, "GtkTreeDragSource")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_drag_source_iface_drag_data_delete", object.class, object, path, PACKAGE = "RGtk2")

  return(w)
}

gtkTreeDragDestIfaceDragDataReceived <-
function(object.class, object, dest, selection.data)
{
  checkPtrType(object.class, "GtkTreeDragDestIface")
  checkPtrType(object, "GtkTreeDragDest")
  checkPtrType(dest, "GtkTreePath")
  checkPtrType(selection.data, "GtkSelectionData")

  w <- .RGtkCall("S_gtk_tree_drag_dest_iface_drag_data_received", object.class, object, dest, selection.data, PACKAGE = "RGtk2")

  return(w)
}

gtkTreeDragDestIfaceRowDropPossible <-
function(object.class, object, dest.path, selection.data)
{
  checkPtrType(object.class, "GtkTreeDragDestIface")
  checkPtrType(object, "GtkTreeDragDest")
  checkPtrType(dest.path, "GtkTreePath")
  checkPtrType(selection.data, "GtkSelectionData")

  w <- .RGtkCall("S_gtk_tree_drag_dest_iface_row_drop_possible", object.class, object, dest.path, selection.data, PACKAGE = "RGtk2")

  return(w)
}

gtkTreeSelectionClassChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkTreeSelectionClass")
  checkPtrType(object, "GtkTreeSelection")

  w <- .RGtkCall("S_gtk_tree_selection_class_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTreeClassSelectChild <-
function(object.class, object, child)
{
  checkPtrType(object.class, "GtkTreeClass")
  checkPtrType(object, "GtkTree")
  checkPtrType(child, "GtkWidget")

  w <- .RGtkCall("S_gtk_tree_class_select_child", object.class, object, child, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTreeClassUnselectChild <-
function(object.class, object, child)
{
  checkPtrType(object.class, "GtkTreeClass")
  checkPtrType(object, "GtkTree")
  checkPtrType(child, "GtkWidget")

  w <- .RGtkCall("S_gtk_tree_class_unselect_child", object.class, object, child, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTreeItemClassExpand <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkTreeItemClass")
  checkPtrType(object, "GtkTreeItem")

  w <- .RGtkCall("S_gtk_tree_item_class_expand", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTreeItemClassCollapse <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkTreeItemClass")
  checkPtrType(object, "GtkTreeItem")

  w <- .RGtkCall("S_gtk_tree_item_class_collapse", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTreeModelIfaceRowChanged <-
function(object.class, object, path, iter)
{
  checkPtrType(object.class, "GtkTreeModelIface")
  checkPtrType(object, "GtkTreeModel")
  checkPtrType(path, "GtkTreePath")
  checkPtrType(iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_model_iface_row_changed", object.class, object, path, iter, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTreeModelIfaceRowInserted <-
function(object.class, object, path, iter)
{
  checkPtrType(object.class, "GtkTreeModelIface")
  checkPtrType(object, "GtkTreeModel")
  checkPtrType(path, "GtkTreePath")
  checkPtrType(iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_model_iface_row_inserted", object.class, object, path, iter, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTreeModelIfaceRowHasChildToggled <-
function(object.class, object, path, iter)
{
  checkPtrType(object.class, "GtkTreeModelIface")
  checkPtrType(object, "GtkTreeModel")
  checkPtrType(path, "GtkTreePath")
  checkPtrType(iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_model_iface_row_has_child_toggled", object.class, object, path, iter, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTreeModelIfaceRowDeleted <-
function(object.class, object, path)
{
  checkPtrType(object.class, "GtkTreeModelIface")
  checkPtrType(object, "GtkTreeModel")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_model_iface_row_deleted", object.class, object, path, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTreeModelIfaceRowsReordered <-
function(object.class, object, path, iter, new.order)
{
  checkPtrType(object.class, "GtkTreeModelIface")
  checkPtrType(object, "GtkTreeModel")
  checkPtrType(path, "GtkTreePath")
  checkPtrType(iter, "GtkTreeIter")
  new.order <- as.list(as.integer(new.order))

  w <- .RGtkCall("S_gtk_tree_model_iface_rows_reordered", object.class, object, path, iter, new.order, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTreeModelIfaceGetFlags <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkTreeModelIface")
  checkPtrType(object, "GtkTreeModel")

  w <- .RGtkCall("S_gtk_tree_model_iface_get_flags", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkTreeModelIfaceGetNColumns <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkTreeModelIface")
  checkPtrType(object, "GtkTreeModel")

  w <- .RGtkCall("S_gtk_tree_model_iface_get_n_columns", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkTreeModelIfaceGetColumnType <-
function(object.class, object, index.)
{
  checkPtrType(object.class, "GtkTreeModelIface")
  checkPtrType(object, "GtkTreeModel")
  index. <- as.integer(index.)

  w <- .RGtkCall("S_gtk_tree_model_iface_get_column_type", object.class, object, index., PACKAGE = "RGtk2")

  return(w)
}

gtkTreeModelIfaceGetIter <-
function(object.class, object, path)
{
  checkPtrType(object.class, "GtkTreeModelIface")
  checkPtrType(object, "GtkTreeModel")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_model_iface_get_iter", object.class, object, path, PACKAGE = "RGtk2")

  return(w)
}

gtkTreeModelIfaceGetPath <-
function(object.class, object, iter)
{
  checkPtrType(object.class, "GtkTreeModelIface")
  checkPtrType(object, "GtkTreeModel")
  checkPtrType(iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_model_iface_get_path", object.class, object, iter, PACKAGE = "RGtk2")

  return(w)
}

gtkTreeModelIfaceGetValue <-
function(object.class, object, iter, column)
{
  checkPtrType(object.class, "GtkTreeModelIface")
  checkPtrType(object, "GtkTreeModel")
  checkPtrType(iter, "GtkTreeIter")
  column <- as.integer(column)

  w <- .RGtkCall("S_gtk_tree_model_iface_get_value", object.class, object, iter, column, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTreeModelIfaceIterNext <-
function(object.class, object, iter)
{
  checkPtrType(object.class, "GtkTreeModelIface")
  checkPtrType(object, "GtkTreeModel")
  checkPtrType(iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_model_iface_iter_next", object.class, object, iter, PACKAGE = "RGtk2")

  return(w)
}

gtkTreeModelIfaceIterChildren <-
function(object.class, object, parent)
{
  checkPtrType(object.class, "GtkTreeModelIface")
  checkPtrType(object, "GtkTreeModel")
  checkPtrType(parent, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_model_iface_iter_children", object.class, object, parent, PACKAGE = "RGtk2")

  return(w)
}

gtkTreeModelIfaceIterHasChild <-
function(object.class, object, iter)
{
  checkPtrType(object.class, "GtkTreeModelIface")
  checkPtrType(object, "GtkTreeModel")
  checkPtrType(iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_model_iface_iter_has_child", object.class, object, iter, PACKAGE = "RGtk2")

  return(w)
}

gtkTreeModelIfaceIterNChildren <-
function(object.class, object, iter)
{
  checkPtrType(object.class, "GtkTreeModelIface")
  checkPtrType(object, "GtkTreeModel")
  checkPtrType(iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_model_iface_iter_n_children", object.class, object, iter, PACKAGE = "RGtk2")

  return(w)
}

gtkTreeModelIfaceIterNthChild <-
function(object.class, object, parent, n)
{
  checkPtrType(object.class, "GtkTreeModelIface")
  checkPtrType(object, "GtkTreeModel")
  checkPtrType(parent, "GtkTreeIter")
  n <- as.integer(n)

  w <- .RGtkCall("S_gtk_tree_model_iface_iter_nth_child", object.class, object, parent, n, PACKAGE = "RGtk2")

  return(w)
}

gtkTreeModelIfaceIterParent <-
function(object.class, object, child)
{
  checkPtrType(object.class, "GtkTreeModelIface")
  checkPtrType(object, "GtkTreeModel")
  checkPtrType(child, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_model_iface_iter_parent", object.class, object, child, PACKAGE = "RGtk2")

  return(w)
}

gtkTreeModelIfaceRefNode <-
function(object.class, object, iter)
{
  checkPtrType(object.class, "GtkTreeModelIface")
  checkPtrType(object, "GtkTreeModel")
  checkPtrType(iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_model_iface_ref_node", object.class, object, iter, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTreeModelIfaceUnrefNode <-
function(object.class, object, iter)
{
  checkPtrType(object.class, "GtkTreeModelIface")
  checkPtrType(object, "GtkTreeModel")
  checkPtrType(iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_model_iface_unref_node", object.class, object, iter, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTreeSelectionClassChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkTreeSelectionClass")
  checkPtrType(object, "GtkTreeSelection")

  w <- .RGtkCall("S_gtk_tree_selection_class_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTreeSortableIfaceSortColumnChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkTreeSortableIface")
  checkPtrType(object, "GtkTreeSortable")

  w <- .RGtkCall("S_gtk_tree_sortable_iface_sort_column_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTreeSortableIfaceGetSortColumnId <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkTreeSortableIface")
  checkPtrType(object, "GtkTreeSortable")

  w <- .RGtkCall("S_gtk_tree_sortable_iface_get_sort_column_id", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkTreeSortableIfaceSetSortColumnId <-
function(object.class, object, sort.column.id, order)
{
  checkPtrType(object.class, "GtkTreeSortableIface")
  checkPtrType(object, "GtkTreeSortable")
  sort.column.id <- as.integer(sort.column.id)
  

  w <- .RGtkCall("S_gtk_tree_sortable_iface_set_sort_column_id", object.class, object, sort.column.id, order, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTreeSortableIfaceSetSortFunc <-
function(object.class, object, sort.column.id, func, data)
{
  checkPtrType(object.class, "GtkTreeSortableIface")
  checkPtrType(object, "GtkTreeSortable")
  sort.column.id <- as.integer(sort.column.id)
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_tree_sortable_iface_set_sort_func", object.class, object, sort.column.id, func, data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTreeSortableIfaceSetDefaultSortFunc <-
function(object.class, object, func, data)
{
  checkPtrType(object.class, "GtkTreeSortableIface")
  checkPtrType(object, "GtkTreeSortable")
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_tree_sortable_iface_set_default_sort_func", object.class, object, func, data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTreeSortableIfaceHasDefaultSortFunc <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkTreeSortableIface")
  checkPtrType(object, "GtkTreeSortable")

  w <- .RGtkCall("S_gtk_tree_sortable_iface_has_default_sort_func", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkTreeViewClassSetScrollAdjustments <-
function(object.class, object, hadjustment, vadjustment)
{
  checkPtrType(object.class, "GtkTreeViewClass")
  checkPtrType(object, "GtkTreeView")
  checkPtrType(hadjustment, "GtkAdjustment")
  checkPtrType(vadjustment, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_tree_view_class_set_scroll_adjustments", object.class, object, hadjustment, vadjustment, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTreeViewClassRowActivated <-
function(object.class, object, path, column)
{
  checkPtrType(object.class, "GtkTreeViewClass")
  checkPtrType(object, "GtkTreeView")
  checkPtrType(path, "GtkTreePath")
  checkPtrType(column, "GtkTreeViewColumn")

  w <- .RGtkCall("S_gtk_tree_view_class_row_activated", object.class, object, path, column, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTreeViewClassTestExpandRow <-
function(object.class, object, iter, path)
{
  checkPtrType(object.class, "GtkTreeViewClass")
  checkPtrType(object, "GtkTreeView")
  checkPtrType(iter, "GtkTreeIter")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_view_class_test_expand_row", object.class, object, iter, path, PACKAGE = "RGtk2")

  return(w)
}

gtkTreeViewClassTestCollapseRow <-
function(object.class, object, iter, path)
{
  checkPtrType(object.class, "GtkTreeViewClass")
  checkPtrType(object, "GtkTreeView")
  checkPtrType(iter, "GtkTreeIter")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_view_class_test_collapse_row", object.class, object, iter, path, PACKAGE = "RGtk2")

  return(w)
}

gtkTreeViewClassRowExpanded <-
function(object.class, object, iter, path)
{
  checkPtrType(object.class, "GtkTreeViewClass")
  checkPtrType(object, "GtkTreeView")
  checkPtrType(iter, "GtkTreeIter")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_view_class_row_expanded", object.class, object, iter, path, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTreeViewClassRowCollapsed <-
function(object.class, object, iter, path)
{
  checkPtrType(object.class, "GtkTreeViewClass")
  checkPtrType(object, "GtkTreeView")
  checkPtrType(iter, "GtkTreeIter")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_view_class_row_collapsed", object.class, object, iter, path, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTreeViewClassColumnsChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkTreeViewClass")
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_class_columns_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTreeViewClassCursorChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkTreeViewClass")
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_class_cursor_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkTreeViewClassMoveCursor <-
function(object.class, object, step, count)
{
  checkPtrType(object.class, "GtkTreeViewClass")
  checkPtrType(object, "GtkTreeView")
  
  count <- as.integer(count)

  w <- .RGtkCall("S_gtk_tree_view_class_move_cursor", object.class, object, step, count, PACKAGE = "RGtk2")

  return(w)
}

gtkTreeViewClassSelectAll <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkTreeViewClass")
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_class_select_all", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkTreeViewClassUnselectAll <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkTreeViewClass")
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_class_unselect_all", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkTreeViewClassSelectCursorRow <-
function(object.class, object, start.editing)
{
  checkPtrType(object.class, "GtkTreeViewClass")
  checkPtrType(object, "GtkTreeView")
  start.editing <- as.logical(start.editing)

  w <- .RGtkCall("S_gtk_tree_view_class_select_cursor_row", object.class, object, start.editing, PACKAGE = "RGtk2")

  return(w)
}

gtkTreeViewClassToggleCursorRow <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkTreeViewClass")
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_class_toggle_cursor_row", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkTreeViewClassExpandCollapseCursorRow <-
function(object.class, object, logical, expand, open.all)
{
  checkPtrType(object.class, "GtkTreeViewClass")
  checkPtrType(object, "GtkTreeView")
  logical <- as.logical(logical)
  expand <- as.logical(expand)
  open.all <- as.logical(open.all)

  w <- .RGtkCall("S_gtk_tree_view_class_expand_collapse_cursor_row", object.class, object, logical, expand, open.all, PACKAGE = "RGtk2")

  return(w)
}

gtkTreeViewClassSelectCursorParent <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkTreeViewClass")
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_class_select_cursor_parent", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkTreeViewClassStartInteractiveSearch <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkTreeViewClass")
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_class_start_interactive_search", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkTreeViewColumnClassClicked <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkTreeViewColumnClass")
  checkPtrType(object, "GtkTreeViewColumn")

  w <- .RGtkCall("S_gtk_tree_view_column_class_clicked", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkUIManagerClassAddWidget <-
function(object.class, object, widget)
{
  checkPtrType(object.class, "GtkUIManagerClass")
  checkPtrType(object, "GtkUIManager")
  checkPtrType(widget, "GtkWidget")

  w <- .RGtkCall("S_gtk_uimanager_class_add_widget", object.class, object, widget, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkUIManagerClassActionsChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkUIManagerClass")
  checkPtrType(object, "GtkUIManager")

  w <- .RGtkCall("S_gtk_uimanager_class_actions_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkUIManagerClassConnectProxy <-
function(object.class, object, action, proxy)
{
  checkPtrType(object.class, "GtkUIManagerClass")
  checkPtrType(object, "GtkUIManager")
  checkPtrType(action, "GtkAction")
  checkPtrType(proxy, "GtkWidget")

  w <- .RGtkCall("S_gtk_uimanager_class_connect_proxy", object.class, object, action, proxy, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkUIManagerClassDisconnectProxy <-
function(object.class, object, action, proxy)
{
  checkPtrType(object.class, "GtkUIManagerClass")
  checkPtrType(object, "GtkUIManager")
  checkPtrType(action, "GtkAction")
  checkPtrType(proxy, "GtkWidget")

  w <- .RGtkCall("S_gtk_uimanager_class_disconnect_proxy", object.class, object, action, proxy, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkUIManagerClassPreActivate <-
function(object.class, object, action)
{
  checkPtrType(object.class, "GtkUIManagerClass")
  checkPtrType(object, "GtkUIManager")
  checkPtrType(action, "GtkAction")

  w <- .RGtkCall("S_gtk_uimanager_class_pre_activate", object.class, object, action, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkUIManagerClassPostActivate <-
function(object.class, object, action)
{
  checkPtrType(object.class, "GtkUIManagerClass")
  checkPtrType(object, "GtkUIManager")
  checkPtrType(action, "GtkAction")

  w <- .RGtkCall("S_gtk_uimanager_class_post_activate", object.class, object, action, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkUIManagerClassGetWidget <-
function(object.class, object, path)
{
  checkPtrType(object.class, "GtkUIManagerClass")
  checkPtrType(object, "GtkUIManager")
  path <- as.character(path)

  w <- .RGtkCall("S_gtk_uimanager_class_get_widget", object.class, object, path, PACKAGE = "RGtk2")

  return(w)
}

gtkUIManagerClassGetAction <-
function(object.class, object, path)
{
  checkPtrType(object.class, "GtkUIManagerClass")
  checkPtrType(object, "GtkUIManager")
  path <- as.character(path)

  w <- .RGtkCall("S_gtk_uimanager_class_get_action", object.class, object, path, PACKAGE = "RGtk2")

  return(w)
}

gtkViewportClassSetScrollAdjustments <-
function(object.class, object, hadjustment, vadjustment)
{
  checkPtrType(object.class, "GtkViewportClass")
  checkPtrType(object, "GtkViewport")
  checkPtrType(hadjustment, "GtkAdjustment")
  checkPtrType(vadjustment, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_viewport_class_set_scroll_adjustments", object.class, object, hadjustment, vadjustment, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkWidgetClassDispatchChildPropertiesChanged <-
function(object.class, object, pspecs)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  pspecs <- lapply(pspecs, function(x) { x <- as.GParamSpec(x); x })

  w <- .RGtkCall("S_gtk_widget_class_dispatch_child_properties_changed", object.class, object, pspecs, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkWidgetClassShow <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_class_show", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkWidgetClassShowAll <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_class_show_all", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkWidgetClassHide <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_class_hide", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkWidgetClassHideAll <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_class_hide_all", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkWidgetClassMap <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_class_map", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkWidgetClassUnmap <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_class_unmap", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkWidgetClassRealize <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_class_realize", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkWidgetClassUnrealize <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_class_unrealize", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkWidgetClassSizeRequest <-
function(object.class, object, requisition)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(requisition, "GtkRequisition")

  w <- .RGtkCall("S_gtk_widget_class_size_request", object.class, object, requisition, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkWidgetClassSizeAllocate <-
function(object.class, object, allocation)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  allocation <- as.GtkAllocation(allocation)

  w <- .RGtkCall("S_gtk_widget_class_size_allocate", object.class, object, allocation, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkWidgetClassStateChanged <-
function(object.class, object, previous.state)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  

  w <- .RGtkCall("S_gtk_widget_class_state_changed", object.class, object, previous.state, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkWidgetClassParentSet <-
function(object.class, object, previous.parent)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(previous.parent, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_class_parent_set", object.class, object, previous.parent, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkWidgetClassHierarchyChanged <-
function(object.class, object, previous.toplevel)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(previous.toplevel, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_class_hierarchy_changed", object.class, object, previous.toplevel, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkWidgetClassStyleSet <-
function(object.class, object, previous.style)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(previous.style, "GtkStyle")

  w <- .RGtkCall("S_gtk_widget_class_style_set", object.class, object, previous.style, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkWidgetClassDirectionChanged <-
function(object.class, object, previous.direction)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  

  w <- .RGtkCall("S_gtk_widget_class_direction_changed", object.class, object, previous.direction, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkWidgetClassGrabNotify <-
function(object.class, object, was.grabbed)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  was.grabbed <- as.logical(was.grabbed)

  w <- .RGtkCall("S_gtk_widget_class_grab_notify", object.class, object, was.grabbed, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkWidgetClassChildNotify <-
function(object.class, object, pspec)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  pspec <- as.GParamSpec(pspec)

  w <- .RGtkCall("S_gtk_widget_class_child_notify", object.class, object, pspec, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkWidgetClassMnemonicActivate <-
function(object.class, object, group.cycling)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  group.cycling <- as.logical(group.cycling)

  w <- .RGtkCall("S_gtk_widget_class_mnemonic_activate", object.class, object, group.cycling, PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassGrabFocus <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_class_grab_focus", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkWidgetClassFocus <-
function(object.class, object, direction)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  

  w <- .RGtkCall("S_gtk_widget_class_focus", object.class, object, direction, PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassEvent <-
function(object.class, object, event)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(event, "GdkEvent")

  w <- .RGtkCall("S_gtk_widget_class_event", object.class, object, event, PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassButtonPressEvent <-
function(object.class, object, event)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(event, "GdkEventButton")

  w <- .RGtkCall("S_gtk_widget_class_button_press_event", object.class, object, event, PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassButtonReleaseEvent <-
function(object.class, object, event)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(event, "GdkEventButton")

  w <- .RGtkCall("S_gtk_widget_class_button_release_event", object.class, object, event, PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassScrollEvent <-
function(object.class, object, event)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(event, "GdkEventScroll")

  w <- .RGtkCall("S_gtk_widget_class_scroll_event", object.class, object, event, PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassMotionNotifyEvent <-
function(object.class, object, event)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(event, "GdkEventMotion")

  w <- .RGtkCall("S_gtk_widget_class_motion_notify_event", object.class, object, event, PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassDeleteEvent <-
function(object.class, object, event)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(event, "GdkEventAny")

  w <- .RGtkCall("S_gtk_widget_class_delete_event", object.class, object, event, PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassDestroyEvent <-
function(object.class, object, event)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(event, "GdkEventAny")

  w <- .RGtkCall("S_gtk_widget_class_destroy_event", object.class, object, event, PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassExposeEvent <-
function(object.class, object, event)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(event, "GdkEventExpose")

  w <- .RGtkCall("S_gtk_widget_class_expose_event", object.class, object, event, PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassKeyPressEvent <-
function(object.class, object, event)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(event, "GdkEventKey")

  w <- .RGtkCall("S_gtk_widget_class_key_press_event", object.class, object, event, PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassKeyReleaseEvent <-
function(object.class, object, event)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(event, "GdkEventKey")

  w <- .RGtkCall("S_gtk_widget_class_key_release_event", object.class, object, event, PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassEnterNotifyEvent <-
function(object.class, object, event)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(event, "GdkEventCrossing")

  w <- .RGtkCall("S_gtk_widget_class_enter_notify_event", object.class, object, event, PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassLeaveNotifyEvent <-
function(object.class, object, event)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(event, "GdkEventCrossing")

  w <- .RGtkCall("S_gtk_widget_class_leave_notify_event", object.class, object, event, PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassConfigureEvent <-
function(object.class, object, event)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(event, "GdkEventConfigure")

  w <- .RGtkCall("S_gtk_widget_class_configure_event", object.class, object, event, PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassFocusInEvent <-
function(object.class, object, event)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(event, "GdkEventFocus")

  w <- .RGtkCall("S_gtk_widget_class_focus_in_event", object.class, object, event, PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassFocusOutEvent <-
function(object.class, object, event)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(event, "GdkEventFocus")

  w <- .RGtkCall("S_gtk_widget_class_focus_out_event", object.class, object, event, PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassMapEvent <-
function(object.class, object, event)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(event, "GdkEventAny")

  w <- .RGtkCall("S_gtk_widget_class_map_event", object.class, object, event, PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassUnmapEvent <-
function(object.class, object, event)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(event, "GdkEventAny")

  w <- .RGtkCall("S_gtk_widget_class_unmap_event", object.class, object, event, PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassPropertyNotifyEvent <-
function(object.class, object, event)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(event, "GdkEventProperty")

  w <- .RGtkCall("S_gtk_widget_class_property_notify_event", object.class, object, event, PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassSelectionClearEvent <-
function(object.class, object, event)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(event, "GdkEventSelection")

  w <- .RGtkCall("S_gtk_widget_class_selection_clear_event", object.class, object, event, PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassSelectionRequestEvent <-
function(object.class, object, event)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(event, "GdkEventSelection")

  w <- .RGtkCall("S_gtk_widget_class_selection_request_event", object.class, object, event, PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassSelectionNotifyEvent <-
function(object.class, object, event)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(event, "GdkEventSelection")

  w <- .RGtkCall("S_gtk_widget_class_selection_notify_event", object.class, object, event, PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassProximityInEvent <-
function(object.class, object, event)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(event, "GdkEventProximity")

  w <- .RGtkCall("S_gtk_widget_class_proximity_in_event", object.class, object, event, PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassProximityOutEvent <-
function(object.class, object, event)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(event, "GdkEventProximity")

  w <- .RGtkCall("S_gtk_widget_class_proximity_out_event", object.class, object, event, PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassVisibilityNotifyEvent <-
function(object.class, object, event)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(event, "GdkEventVisibility")

  w <- .RGtkCall("S_gtk_widget_class_visibility_notify_event", object.class, object, event, PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassClientEvent <-
function(object.class, object, event)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(event, "GdkEventClient")

  w <- .RGtkCall("S_gtk_widget_class_client_event", object.class, object, event, PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassNoExposeEvent <-
function(object.class, object, event)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(event, "GdkEventAny")

  w <- .RGtkCall("S_gtk_widget_class_no_expose_event", object.class, object, event, PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassWindowStateEvent <-
function(object.class, object, event)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(event, "GdkEventWindowState")

  w <- .RGtkCall("S_gtk_widget_class_window_state_event", object.class, object, event, PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassSelectionGet <-
function(object.class, object, selection.data, info, time.)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(selection.data, "GtkSelectionData")
  info <- as.numeric(info)
  time. <- as.numeric(time.)

  w <- .RGtkCall("S_gtk_widget_class_selection_get", object.class, object, selection.data, info, time., PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkWidgetClassSelectionReceived <-
function(object.class, object, selection.data, time.)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(selection.data, "GtkSelectionData")
  time. <- as.numeric(time.)

  w <- .RGtkCall("S_gtk_widget_class_selection_received", object.class, object, selection.data, time., PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkWidgetClassDragBegin <-
function(object.class, object, context)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(context, "GdkDragContext")

  w <- .RGtkCall("S_gtk_widget_class_drag_begin", object.class, object, context, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkWidgetClassDragEnd <-
function(object.class, object, context)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(context, "GdkDragContext")

  w <- .RGtkCall("S_gtk_widget_class_drag_end", object.class, object, context, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkWidgetClassDragDataGet <-
function(object.class, object, context, selection.data, info, time.)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(context, "GdkDragContext")
  checkPtrType(selection.data, "GtkSelectionData")
  info <- as.numeric(info)
  time. <- as.numeric(time.)

  w <- .RGtkCall("S_gtk_widget_class_drag_data_get", object.class, object, context, selection.data, info, time., PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkWidgetClassDragDataDelete <-
function(object.class, object, context)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(context, "GdkDragContext")

  w <- .RGtkCall("S_gtk_widget_class_drag_data_delete", object.class, object, context, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkWidgetClassDragLeave <-
function(object.class, object, context, time.)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(context, "GdkDragContext")
  time. <- as.numeric(time.)

  w <- .RGtkCall("S_gtk_widget_class_drag_leave", object.class, object, context, time., PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkWidgetClassDragMotion <-
function(object.class, object, context, x, y, time.)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(context, "GdkDragContext")
  x <- as.integer(x)
  y <- as.integer(y)
  time. <- as.numeric(time.)

  w <- .RGtkCall("S_gtk_widget_class_drag_motion", object.class, object, context, x, y, time., PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassDragDrop <-
function(object.class, object, context, x, y, time.)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(context, "GdkDragContext")
  x <- as.integer(x)
  y <- as.integer(y)
  time. <- as.numeric(time.)

  w <- .RGtkCall("S_gtk_widget_class_drag_drop", object.class, object, context, x, y, time., PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassDragDataReceived <-
function(object.class, object, context, x, y, selection.data, info, time.)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(context, "GdkDragContext")
  x <- as.integer(x)
  y <- as.integer(y)
  checkPtrType(selection.data, "GtkSelectionData")
  info <- as.numeric(info)
  time. <- as.numeric(time.)

  w <- .RGtkCall("S_gtk_widget_class_drag_data_received", object.class, object, context, x, y, selection.data, info, time., PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkWidgetClassPopupMenu <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_class_popup_menu", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassShowHelp <-
function(object.class, object, help.type)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  

  w <- .RGtkCall("S_gtk_widget_class_show_help", object.class, object, help.type, PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassGetAccessible <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_class_get_accessible", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassScreenChanged <-
function(object.class, object, previous.screen)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(previous.screen, "GdkScreen")

  w <- .RGtkCall("S_gtk_widget_class_screen_changed", object.class, object, previous.screen, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkWidgetClassCanActivateAccel <-
function(object.class, object, signal.id)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  signal.id <- as.numeric(signal.id)

  w <- .RGtkCall("S_gtk_widget_class_can_activate_accel", object.class, object, signal.id, PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassGrabBrokenEvent <-
function(object.class, object, event)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")
  checkPtrType(event, "GdkEventGrabBroken")

  w <- .RGtkCall("S_gtk_widget_class_grab_broken_event", object.class, object, event, PACKAGE = "RGtk2")

  return(w)
}

gtkWindowClassSetFocus <-
function(object.class, object, focus)
{
  checkPtrType(object.class, "GtkWindowClass")
  checkPtrType(object, "GtkWindow")
  checkPtrType(focus, "GtkWidget")

  w <- .RGtkCall("S_gtk_window_class_set_focus", object.class, object, focus, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkWindowClassFrameEvent <-
function(object.class, object, event)
{
  checkPtrType(object.class, "GtkWindowClass")
  checkPtrType(object, "GtkWindow")
  checkPtrType(event, "GdkEvent")

  w <- .RGtkCall("S_gtk_window_class_frame_event", object.class, object, event, PACKAGE = "RGtk2")

  return(w)
}

gtkWindowClassActivateFocus <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkWindowClass")
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_class_activate_focus", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkWindowClassActivateDefault <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkWindowClass")
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_class_activate_default", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkWindowClassMoveFocus <-
function(object.class, object, direction)
{
  checkPtrType(object.class, "GtkWindowClass")
  checkPtrType(object, "GtkWindow")
  

  w <- .RGtkCall("S_gtk_window_class_move_focus", object.class, object, direction, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkWindowClassKeysChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkWindowClass")
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_class_keys_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkAssistantClassPrepare <-
function(object.class, object, page)
{
  checkPtrType(object.class, "GtkAssistantClass")
  checkPtrType(object, "GtkAssistant")
  checkPtrType(page, "GtkWidget")

  w <- .RGtkCall("S_gtk_assistant_class_prepare", object.class, object, page, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkAssistantClassApply <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkAssistantClass")
  checkPtrType(object, "GtkAssistant")

  w <- .RGtkCall("S_gtk_assistant_class_apply", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkAssistantClassClose <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkAssistantClass")
  checkPtrType(object, "GtkAssistant")

  w <- .RGtkCall("S_gtk_assistant_class_close", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkAssistantClassCancel <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkAssistantClass")
  checkPtrType(object, "GtkAssistant")

  w <- .RGtkCall("S_gtk_assistant_class_cancel", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCellRendererAccelClassAccelEdited <-
function(object.class, object, path.string, accel.key, accel.mods, hardware.keycode)
{
  checkPtrType(object.class, "GtkCellRendererAccelClass")
  checkPtrType(object, "GtkCellRendererAccel")
  path.string <- as.character(path.string)
  accel.key <- as.numeric(accel.key)
  
  hardware.keycode <- as.numeric(hardware.keycode)

  w <- .RGtkCall("S_gtk_cell_renderer_accel_class_accel_edited", object.class, object, path.string, accel.key, accel.mods, hardware.keycode, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkCellRendererAccelClassAccelCleared <-
function(object.class, object, path.string)
{
  checkPtrType(object.class, "GtkCellRendererAccelClass")
  checkPtrType(object, "GtkCellRendererAccel")
  path.string <- as.character(path.string)

  w <- .RGtkCall("S_gtk_cell_renderer_accel_class_accel_cleared", object.class, object, path.string, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkNotebookClassReorderTab <-
function(object.class, object, direction, move.to.last)
{
  checkPtrType(object.class, "GtkNotebookClass")
  checkPtrType(object, "GtkNotebook")
  
  move.to.last <- as.logical(move.to.last)

  w <- .RGtkCall("S_gtk_notebook_class_reorder_tab", object.class, object, direction, move.to.last, PACKAGE = "RGtk2")

  return(w)
}

gtkNotebookClassInsertPage <-
function(object.class, object, child, tab.label, menu.label, position)
{
  checkPtrType(object.class, "GtkNotebookClass")
  checkPtrType(object, "GtkNotebook")
  checkPtrType(child, "GtkWidget")
  checkPtrType(tab.label, "GtkWidget")
  checkPtrType(menu.label, "GtkWidget")
  position <- as.integer(position)

  w <- .RGtkCall("S_gtk_notebook_class_insert_page", object.class, object, child, tab.label, menu.label, position, PACKAGE = "RGtk2")

  return(w)
}

gtkPrintOperationClassDone <-
function(object.class, object, result)
{
  checkPtrType(object.class, "GtkPrintOperationClass")
  checkPtrType(object, "GtkPrintOperation")
  

  w <- .RGtkCall("S_gtk_print_operation_class_done", object.class, object, result, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkPrintOperationClassBeginPrint <-
function(object.class, object, context)
{
  checkPtrType(object.class, "GtkPrintOperationClass")
  checkPtrType(object, "GtkPrintOperation")
  checkPtrType(context, "GtkPrintContext")

  w <- .RGtkCall("S_gtk_print_operation_class_begin_print", object.class, object, context, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkPrintOperationClassPaginate <-
function(object.class, object, context)
{
  checkPtrType(object.class, "GtkPrintOperationClass")
  checkPtrType(object, "GtkPrintOperation")
  checkPtrType(context, "GtkPrintContext")

  w <- .RGtkCall("S_gtk_print_operation_class_paginate", object.class, object, context, PACKAGE = "RGtk2")

  return(w)
}

gtkPrintOperationClassRequestPageSetup <-
function(object.class, object, context, page.nr, setup)
{
  checkPtrType(object.class, "GtkPrintOperationClass")
  checkPtrType(object, "GtkPrintOperation")
  checkPtrType(context, "GtkPrintContext")
  page.nr <- as.integer(page.nr)
  checkPtrType(setup, "GtkPageSetup")

  w <- .RGtkCall("S_gtk_print_operation_class_request_page_setup", object.class, object, context, page.nr, setup, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkPrintOperationClassDrawPage <-
function(object.class, object, context, page.nr)
{
  checkPtrType(object.class, "GtkPrintOperationClass")
  checkPtrType(object, "GtkPrintOperation")
  checkPtrType(context, "GtkPrintContext")
  page.nr <- as.integer(page.nr)

  w <- .RGtkCall("S_gtk_print_operation_class_draw_page", object.class, object, context, page.nr, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkPrintOperationClassEndPrint <-
function(object.class, object, context)
{
  checkPtrType(object.class, "GtkPrintOperationClass")
  checkPtrType(object, "GtkPrintOperation")
  checkPtrType(context, "GtkPrintContext")

  w <- .RGtkCall("S_gtk_print_operation_class_end_print", object.class, object, context, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkPrintOperationClassStatusChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkPrintOperationClass")
  checkPtrType(object, "GtkPrintOperation")

  w <- .RGtkCall("S_gtk_print_operation_class_status_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkPrintOperationClassCreateCustomWidget <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkPrintOperationClass")
  checkPtrType(object, "GtkPrintOperation")

  w <- .RGtkCall("S_gtk_print_operation_class_create_custom_widget", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkPrintOperationClassCustomWidgetApply <-
function(object.class, object, widget)
{
  checkPtrType(object.class, "GtkPrintOperationClass")
  checkPtrType(object, "GtkPrintOperation")
  checkPtrType(widget, "GtkWidget")

  w <- .RGtkCall("S_gtk_print_operation_class_custom_widget_apply", object.class, object, widget, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkPrintOperationClassPreview <-
function(object.class, object, preview, context, parent)
{
  checkPtrType(object.class, "GtkPrintOperationClass")
  checkPtrType(object, "GtkPrintOperation")
  checkPtrType(preview, "GtkPrintOperationPreview")
  checkPtrType(context, "GtkPrintContext")
  checkPtrType(parent, "GtkWindow")

  w <- .RGtkCall("S_gtk_print_operation_class_preview", object.class, object, preview, context, parent, PACKAGE = "RGtk2")

  return(w)
}

gtkPrintOperationPreviewClassReady <-
function(object.class, object, context)
{
  checkPtrType(object.class, "GtkPrintOperationPreviewClass")
  checkPtrType(object, "GtkPrintOperationPreview")
  checkPtrType(context, "GtkPrintContext")

  w <- .RGtkCall("S_gtk_print_operation_preview_class_ready", object.class, object, context, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkPrintOperationPreviewClassGotPageSize <-
function(object.class, object, context, page.setup)
{
  checkPtrType(object.class, "GtkPrintOperationPreviewClass")
  checkPtrType(object, "GtkPrintOperationPreview")
  checkPtrType(context, "GtkPrintContext")
  checkPtrType(page.setup, "GtkPageSetup")

  w <- .RGtkCall("S_gtk_print_operation_preview_class_got_page_size", object.class, object, context, page.setup, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkPrintOperationPreviewClassRenderPage <-
function(object.class, object, page.nr)
{
  checkPtrType(object.class, "GtkPrintOperationPreviewClass")
  checkPtrType(object, "GtkPrintOperationPreview")
  page.nr <- as.integer(page.nr)

  w <- .RGtkCall("S_gtk_print_operation_preview_class_render_page", object.class, object, page.nr, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkPrintOperationPreviewClassIsSelected <-
function(object.class, object, page.nr)
{
  checkPtrType(object.class, "GtkPrintOperationPreviewClass")
  checkPtrType(object, "GtkPrintOperationPreview")
  page.nr <- as.integer(page.nr)

  w <- .RGtkCall("S_gtk_print_operation_preview_class_is_selected", object.class, object, page.nr, PACKAGE = "RGtk2")

  return(w)
}

gtkPrintOperationPreviewClassEndPreview <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkPrintOperationPreviewClass")
  checkPtrType(object, "GtkPrintOperationPreview")

  w <- .RGtkCall("S_gtk_print_operation_preview_class_end_preview", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkRecentChooserClassSetCurrentUri <-
function(object.class, object, uri, .errwarn = TRUE)
{
  checkPtrType(object.class, "GtkRecentChooserClass")
  checkPtrType(object, "GtkRecentChooser")
  uri <- as.character(uri)

  w <- .RGtkCall("S_gtk_recent_chooser_class_set_current_uri", object.class, object, uri, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gtkRecentChooserClassGetCurrentUri <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkRecentChooserClass")
  checkPtrType(object, "GtkRecentChooser")

  w <- .RGtkCall("S_gtk_recent_chooser_class_get_current_uri", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkRecentChooserClassSelectUri <-
function(object.class, object, uri, .errwarn = TRUE)
{
  checkPtrType(object.class, "GtkRecentChooserClass")
  checkPtrType(object, "GtkRecentChooser")
  uri <- as.character(uri)

  w <- .RGtkCall("S_gtk_recent_chooser_class_select_uri", object.class, object, uri, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gtkRecentChooserClassUnselectUri <-
function(object.class, object, uri)
{
  checkPtrType(object.class, "GtkRecentChooserClass")
  checkPtrType(object, "GtkRecentChooser")
  uri <- as.character(uri)

  w <- .RGtkCall("S_gtk_recent_chooser_class_unselect_uri", object.class, object, uri, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkRecentChooserClassSelectAll <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkRecentChooserClass")
  checkPtrType(object, "GtkRecentChooser")

  w <- .RGtkCall("S_gtk_recent_chooser_class_select_all", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkRecentChooserClassUnselectAll <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkRecentChooserClass")
  checkPtrType(object, "GtkRecentChooser")

  w <- .RGtkCall("S_gtk_recent_chooser_class_unselect_all", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkRecentChooserClassGetItems <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkRecentChooserClass")
  checkPtrType(object, "GtkRecentChooser")

  w <- .RGtkCall("S_gtk_recent_chooser_class_get_items", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkRecentChooserClassGetRecentManager <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkRecentChooserClass")
  checkPtrType(object, "GtkRecentChooser")

  w <- .RGtkCall("S_gtk_recent_chooser_class_get_recent_manager", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkRecentChooserClassAddFilter <-
function(object.class, object, filter)
{
  checkPtrType(object.class, "GtkRecentChooserClass")
  checkPtrType(object, "GtkRecentChooser")
  checkPtrType(filter, "GtkRecentFilter")

  w <- .RGtkCall("S_gtk_recent_chooser_class_add_filter", object.class, object, filter, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkRecentChooserClassRemoveFilter <-
function(object.class, object, filter)
{
  checkPtrType(object.class, "GtkRecentChooserClass")
  checkPtrType(object, "GtkRecentChooser")
  checkPtrType(filter, "GtkRecentFilter")

  w <- .RGtkCall("S_gtk_recent_chooser_class_remove_filter", object.class, object, filter, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkRecentChooserClassListFilters <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkRecentChooserClass")
  checkPtrType(object, "GtkRecentChooser")

  w <- .RGtkCall("S_gtk_recent_chooser_class_list_filters", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkRecentChooserClassSetSortFunc <-
function(object.class, object, sort.func, data)
{
  checkPtrType(object.class, "GtkRecentChooserClass")
  checkPtrType(object, "GtkRecentChooser")
  sort.func <- as.function(sort.func)
  

  w <- .RGtkCall("S_gtk_recent_chooser_class_set_sort_func", object.class, object, sort.func, data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkRecentChooserClassItemActivated <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkRecentChooserClass")
  checkPtrType(object, "GtkRecentChooser")

  w <- .RGtkCall("S_gtk_recent_chooser_class_item_activated", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkRecentChooserClassSelectionChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkRecentChooserClass")
  checkPtrType(object, "GtkRecentChooser")

  w <- .RGtkCall("S_gtk_recent_chooser_class_selection_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkRecentManagerClassChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkRecentManagerClass")
  checkPtrType(object, "GtkRecentManager")

  w <- .RGtkCall("S_gtk_recent_manager_class_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkSpinButtonClassWrapped <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkSpinButtonClass")
  checkPtrType(object, "GtkSpinButton")

  w <- .RGtkCall("S_gtk_spin_button_class_wrapped", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkStatusIconClassActivate <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkStatusIconClass")
  checkPtrType(object, "GtkStatusIcon")

  w <- .RGtkCall("S_gtk_status_icon_class_activate", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkStatusIconClassPopupMenu <-
function(object.class, object, button, activate.time)
{
  checkPtrType(object.class, "GtkStatusIconClass")
  checkPtrType(object, "GtkStatusIcon")
  button <- as.numeric(button)
  activate.time <- as.numeric(activate.time)

  w <- .RGtkCall("S_gtk_status_icon_class_popup_menu", object.class, object, button, activate.time, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkStatusIconClassSizeChanged <-
function(object.class, object, size)
{
  checkPtrType(object.class, "GtkStatusIconClass")
  checkPtrType(object, "GtkStatusIcon")
  size <- as.integer(size)

  w <- .RGtkCall("S_gtk_status_icon_class_size_changed", object.class, object, size, PACKAGE = "RGtk2")

  return(w)
}

gtkWidgetClassCompositedChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkWidgetClass")
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_class_composited_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkBuildableIfaceSetName <-
function(object.class, object, name)
{
  checkPtrType(object.class, "GtkBuildableIface")
  checkPtrType(object, "GtkBuildable")
  name <- as.character(name)

  w <- .RGtkCall("S_gtk_buildable_iface_set_name", object.class, object, name, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkBuildableIfaceGetName <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkBuildableIface")
  checkPtrType(object, "GtkBuildable")

  w <- .RGtkCall("S_gtk_buildable_iface_get_name", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkBuildableIfaceAddChild <-
function(object.class, object, builder, child, type)
{
  checkPtrType(object.class, "GtkBuildableIface")
  checkPtrType(object, "GtkBuildable")
  checkPtrType(builder, "GtkBuilder")
  checkPtrType(child, "GObject")
  type <- as.character(type)

  w <- .RGtkCall("S_gtk_buildable_iface_add_child", object.class, object, builder, child, type, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkBuildableIfaceSetBuildableProperty <-
function(object.class, object, builder, name, value)
{
  checkPtrType(object.class, "GtkBuildableIface")
  checkPtrType(object, "GtkBuildable")
  checkPtrType(builder, "GtkBuilder")
  name <- as.character(name)
  

  w <- .RGtkCall("S_gtk_buildable_iface_set_buildable_property", object.class, object, builder, name, value, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkBuildableIfaceConstructChild <-
function(object.class, object, builder, name)
{
  checkPtrType(object.class, "GtkBuildableIface")
  checkPtrType(object, "GtkBuildable")
  checkPtrType(builder, "GtkBuilder")
  name <- as.character(name)

  w <- .RGtkCall("S_gtk_buildable_iface_construct_child", object.class, object, builder, name, PACKAGE = "RGtk2")

  return(w)
}

gtkBuildableIfaceCustomTagStart <-
function(object.class, object, builder, child, tagname, parser)
{
  checkPtrType(object.class, "GtkBuildableIface")
  checkPtrType(object, "GtkBuildable")
  checkPtrType(builder, "GtkBuilder")
  checkPtrType(child, "GObject")
  tagname <- as.character(tagname)
  checkPtrType(parser, "GMarkupParser")

  w <- .RGtkCall("S_gtk_buildable_iface_custom_tag_start", object.class, object, builder, child, tagname, parser, PACKAGE = "RGtk2")

  return(w)
}

gtkBuildableIfaceCustomTagEnd <-
function(object.class, object, builder, child, tagname)
{
  checkPtrType(object.class, "GtkBuildableIface")
  checkPtrType(object, "GtkBuildable")
  checkPtrType(builder, "GtkBuilder")
  checkPtrType(child, "GObject")
  tagname <- as.character(tagname)

  w <- .RGtkCall("S_gtk_buildable_iface_custom_tag_end", object.class, object, builder, child, tagname, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkBuildableIfaceCustomFinished <-
function(object.class, object, builder, child, tagname, data)
{
  checkPtrType(object.class, "GtkBuildableIface")
  checkPtrType(object, "GtkBuildable")
  checkPtrType(builder, "GtkBuilder")
  checkPtrType(child, "GObject")
  tagname <- as.character(tagname)
  

  w <- .RGtkCall("S_gtk_buildable_iface_custom_finished", object.class, object, builder, child, tagname, data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkBuildableIfaceParserFinished <-
function(object.class, object, builder)
{
  checkPtrType(object.class, "GtkBuildableIface")
  checkPtrType(object, "GtkBuildable")
  checkPtrType(builder, "GtkBuilder")

  w <- .RGtkCall("S_gtk_buildable_iface_parser_finished", object.class, object, builder, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkBuildableIfaceGetInternalChild <-
function(object.class, object, builder, childname)
{
  checkPtrType(object.class, "GtkBuildableIface")
  checkPtrType(object, "GtkBuildable")
  checkPtrType(builder, "GtkBuilder")
  childname <- as.character(childname)

  w <- .RGtkCall("S_gtk_buildable_iface_get_internal_child", object.class, object, builder, childname, PACKAGE = "RGtk2")

  return(w)
}

gtkBuilderClassGetTypeFromName <-
function(object.class, object, type.name)
{
  checkPtrType(object.class, "GtkBuilderClass")
  checkPtrType(object, "GtkBuilder")
  type.name <- as.character(type.name)

  w <- .RGtkCall("S_gtk_builder_class_get_type_from_name", object.class, object, type.name, PACKAGE = "RGtk2")

  return(w)
}

gtkToolShellIfaceGetIconSize <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkToolShellIface")
  checkPtrType(object, "GtkToolShell")

  w <- .RGtkCall("S_gtk_tool_shell_iface_get_icon_size", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkToolShellIfaceGetOrientation <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkToolShellIface")
  checkPtrType(object, "GtkToolShell")

  w <- .RGtkCall("S_gtk_tool_shell_iface_get_orientation", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkToolShellIfaceGetStyle <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkToolShellIface")
  checkPtrType(object, "GtkToolShell")

  w <- .RGtkCall("S_gtk_tool_shell_iface_get_style", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkToolShellIfaceGetReliefStyle <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkToolShellIface")
  checkPtrType(object, "GtkToolShell")

  w <- .RGtkCall("S_gtk_tool_shell_iface_get_relief_style", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gtkToolShellIfaceRebuildMenu <-
function(object.class, object)
{
  checkPtrType(object.class, "GtkToolShellIface")
  checkPtrType(object, "GtkToolShell")

  w <- .RGtkCall("S_gtk_tool_shell_iface_rebuild_menu", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkActivatableIfaceUpdate <-
function(object.class, object, action, property.name)
{
  checkPtrType(object.class, "GtkActivatableIface")
  checkPtrType(object, "GtkActivatable")
  checkPtrType(action, "GtkAction")
  property.name <- as.character(property.name)

  w <- .RGtkCall("S_gtk_activatable_iface_update", object.class, object, action, property.name, PACKAGE = "RGtk2")

  return(invisible(w))
}

gtkActivatableIfaceSyncActionProperties <-
function(object.class, object, action)
{
  checkPtrType(object.class, "GtkActivatableIface")
  checkPtrType(object, "GtkActivatable")
  checkPtrType(action, "GtkAction")

  w <- .RGtkCall("S_gtk_activatable_iface_sync_action_properties", object.class, object, action, PACKAGE = "RGtk2")

  return(invisible(w))
}