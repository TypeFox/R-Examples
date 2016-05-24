
gdkNotifyStartupComplete <-
function()
{
  

  w <- .RGtkCall("S_gdk_notify_startup_complete", PACKAGE = "RGtk2")

  return(w)
} 


gdkGetProgramClass <-
function()
{
  

  w <- .RGtkCall("S_gdk_get_program_class", PACKAGE = "RGtk2")

  return(w)
} 


gdkSetProgramClass <-
function(program.class)
{
  program.class <- as.character(program.class)

  w <- .RGtkCall("S_gdk_set_program_class", program.class, PACKAGE = "RGtk2")

  return(w)
} 


gdkGetDisplay <-
function()
{
  

  w <- .RGtkCall("S_gdk_get_display", PACKAGE = "RGtk2")

  return(w)
} 


gdkPointerGrab <-
function(window, owner.events = FALSE, event.mask = 0, confine.to = NULL, cursor = NULL, time = "GDK_CURRENT_TIME")
{
  checkPtrType(window, "GdkWindow")
  owner.events <- as.logical(owner.events)
  
  if (!is.null( confine.to )) checkPtrType(confine.to, "GdkWindow")
  if (!is.null( cursor )) checkPtrType(cursor, "GdkCursor")
  time <- as.numeric(time)

  w <- .RGtkCall("S_gdk_pointer_grab", window, owner.events, event.mask, confine.to, cursor, time, PACKAGE = "RGtk2")

  return(w)
} 


gdkPointerUngrab <-
function(time = "GDK_CURRENT_TIME")
{
  time <- as.numeric(time)

  w <- .RGtkCall("S_gdk_pointer_ungrab", time, PACKAGE = "RGtk2")

  return(w)
} 


gdkKeyboardGrab <-
function(window, owner.events = FALSE, time = "GDK_CURRENT_TIME")
{
  checkPtrType(window, "GdkWindow")
  owner.events <- as.logical(owner.events)
  time <- as.numeric(time)

  w <- .RGtkCall("S_gdk_keyboard_grab", window, owner.events, time, PACKAGE = "RGtk2")

  return(w)
} 


gdkKeyboardUngrab <-
function(time = "GDK_CURRENT_TIME")
{
  time <- as.numeric(time)

  w <- .RGtkCall("S_gdk_keyboard_ungrab", time, PACKAGE = "RGtk2")

  return(w)
} 


gdkPointerIsGrabbed <-
function()
{
  

  w <- .RGtkCall("S_gdk_pointer_is_grabbed", PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenWidth <-
function()
{
  

  w <- .RGtkCall("S_gdk_screen_width", PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenHeight <-
function()
{
  

  w <- .RGtkCall("S_gdk_screen_height", PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenWidthMm <-
function()
{
  

  w <- .RGtkCall("S_gdk_screen_width_mm", PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenHeightMm <-
function()
{
  

  w <- .RGtkCall("S_gdk_screen_height_mm", PACKAGE = "RGtk2")

  return(w)
} 


gdkFlush <-
function()
{
  

  w <- .RGtkCall("S_gdk_flush", PACKAGE = "RGtk2")

  return(w)
} 


gdkBeep <-
function()
{
  

  w <- .RGtkCall("S_gdk_beep", PACKAGE = "RGtk2")

  return(w)
} 


gdkSetDoubleClickTime <-
function(msec)
{
  msec <- as.numeric(msec)

  w <- .RGtkCall("S_gdk_set_double_click_time", msec, PACKAGE = "RGtk2")

  return(w)
} 


gdkCairoCreate <-
function(drawable)
{
  checkPtrType(drawable, "GdkDrawable")

  w <- .RGtkCall("S_gdk_cairo_create", drawable, PACKAGE = "RGtk2")

  return(w)
} 


gdkCairoSetSourceColor <-
function(cr, color)
{
  checkPtrType(cr, "Cairo")
  color <- as.GdkColor(color)

  w <- .RGtkCall("S_gdk_cairo_set_source_color", cr, color, PACKAGE = "RGtk2")

  return(w)
} 


gdkCairoSetSourcePixbuf <-
function(cr, pixbuf, pixbuf.x, pixbuf.y)
{
  checkPtrType(cr, "Cairo")
  checkPtrType(pixbuf, "GdkPixbuf")
  pixbuf.x <- as.numeric(pixbuf.x)
  pixbuf.y <- as.numeric(pixbuf.y)

  w <- .RGtkCall("S_gdk_cairo_set_source_pixbuf", cr, pixbuf, pixbuf.x, pixbuf.y, PACKAGE = "RGtk2")

  return(w)
} 


gdkCairoRectangle <-
function(cr, rectangle)
{
  checkPtrType(cr, "Cairo")
  rectangle <- as.GdkRectangle(rectangle)

  w <- .RGtkCall("S_gdk_cairo_rectangle", cr, rectangle, PACKAGE = "RGtk2")

  return(w)
} 


gdkCairoRegion <-
function(cr, region)
{
  checkPtrType(cr, "Cairo")
  checkPtrType(region, "GdkRegion")

  w <- .RGtkCall("S_gdk_cairo_region", cr, region, PACKAGE = "RGtk2")

  return(w)
} 


gdkColormapGetType <-
function()
{
  

  w <- .RGtkCall("S_gdk_colormap_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gdkColormapNew <-
function(visual, allocate)
{
  checkPtrType(visual, "GdkVisual")
  allocate <- as.logical(allocate)

  w <- .RGtkCall("S_gdk_colormap_new", visual, allocate, PACKAGE = "RGtk2")

  return(w)
} 


gdkColormapGetSystem <-
function()
{
  

  w <- .RGtkCall("S_gdk_colormap_get_system", PACKAGE = "RGtk2")

  return(w)
} 


gdkColormapGetSystemSize <-
function()
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  

  w <- .RGtkCall("S_gdk_colormap_get_system_size", PACKAGE = "RGtk2")

  return(w)
} 


gdkColormapAllocColor <-
function(object, color, writeable, best.match)
{
  checkPtrType(object, "GdkColormap")
  color <- as.GdkColor(color)
  writeable <- as.logical(writeable)
  best.match <- as.logical(best.match)

  w <- .RGtkCall("S_gdk_colormap_alloc_color", object, color, writeable, best.match, PACKAGE = "RGtk2")

  return(w)
} 


gdkColormapFreeColors <-
function(object, colors)
{
  checkPtrType(object, "GdkColormap")
  colors <- lapply(colors, function(x) { x <- as.GdkColor(x); x })

  w <- .RGtkCall("S_gdk_colormap_free_colors", object, colors, PACKAGE = "RGtk2")

  return(w)
} 


gdkColormapQueryColor <-
function(object, pixel)
{
  checkPtrType(object, "GdkColormap")
  pixel <- as.numeric(pixel)

  w <- .RGtkCall("S_gdk_colormap_query_color", object, pixel, PACKAGE = "RGtk2")

  return(w)
} 


gdkColormapGetVisual <-
function(object)
{
  checkPtrType(object, "GdkColormap")

  w <- .RGtkCall("S_gdk_colormap_get_visual", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkColormapGetScreen <-
function(object)
{
  checkPtrType(object, "GdkColormap")

  w <- .RGtkCall("S_gdk_colormap_get_screen", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkColorParse <-
function(spec)
{
  spec <- as.character(spec)

  w <- .RGtkCall("S_gdk_color_parse", spec, PACKAGE = "RGtk2")

  return(w)
} 


gdkColorWhite <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  checkPtrType(object, "GdkColormap")

  w <- .RGtkCall("S_gdk_color_white", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkColorBlack <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  checkPtrType(object, "GdkColormap")

  w <- .RGtkCall("S_gdk_color_black", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkColorAlloc <-
function(object, color)
{
  if(getOption("depwarn"))
    .Deprecated("gdkColormapAllocColor", "RGtk2")

  checkPtrType(object, "GdkColormap")
  color <- as.GdkColor(color)

  w <- .RGtkCall("S_gdk_color_alloc", object, color, PACKAGE = "RGtk2")

  return(w)
} 


gdkColorChange <-
function(object, color)
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  checkPtrType(object, "GdkColormap")
  color <- as.GdkColor(color)

  w <- .RGtkCall("S_gdk_color_change", object, color, PACKAGE = "RGtk2")

  return(w)
} 


gdkCursorNew <-
function(cursor.type)
{
  

  w <- .RGtkCall("S_gdk_cursor_new", cursor.type, PACKAGE = "RGtk2")

  return(w)
} 


gdkCursorNewFromName <-
function(display, name)
{
  checkPtrType(display, "GdkDisplay")
  name <- as.character(name)

  w <- .RGtkCall("S_gdk_cursor_new_from_name", display, name, PACKAGE = "RGtk2")

  return(w)
} 


gdkCursorNewForDisplay <-
function(display, cursor.type)
{
  checkPtrType(display, "GdkDisplay")
  

  w <- .RGtkCall("S_gdk_cursor_new_for_display", display, cursor.type, PACKAGE = "RGtk2")

  return(w)
} 


gdkCursorNewFromPixmap <-
function(source, mask, fg, bg, x, y)
{
  checkPtrType(source, "GdkPixmap")
  checkPtrType(mask, "GdkPixmap")
  fg <- as.GdkColor(fg)
  bg <- as.GdkColor(bg)
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_gdk_cursor_new_from_pixmap", source, mask, fg, bg, x, y, PACKAGE = "RGtk2")

  return(w)
} 


gdkCursorNewFromPixbuf <-
function(display, source, x, y)
{
  checkPtrType(display, "GdkDisplay")
  checkPtrType(source, "GdkPixbuf")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_gdk_cursor_new_from_pixbuf", display, source, x, y, PACKAGE = "RGtk2")

  return(w)
} 


gdkCursorGetDisplay <-
function(object)
{
  checkPtrType(object, "GdkCursor")

  w <- .RGtkCall("S_gdk_cursor_get_display", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkCursorGetImage <-
function(object)
{
  checkPtrType(object, "GdkCursor")

  w <- .RGtkCall("S_gdk_cursor_get_image", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkDisplayGetType <-
function()
{
  

  w <- .RGtkCall("S_gdk_display_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gdkDisplayOpen <-
function(display.name)
{
  display.name <- as.character(display.name)

  w <- .RGtkCall("S_gdk_display_open", display.name, PACKAGE = "RGtk2")

  return(w)
} 


gdkDisplayGetName <-
function(object)
{
  checkPtrType(object, "GdkDisplay")

  w <- .RGtkCall("S_gdk_display_get_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkDisplayGetNScreens <-
function(object)
{
  checkPtrType(object, "GdkDisplay")

  w <- .RGtkCall("S_gdk_display_get_n_screens", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkDisplayGetScreen <-
function(object, screen.num)
{
  checkPtrType(object, "GdkDisplay")
  screen.num <- as.integer(screen.num)

  w <- .RGtkCall("S_gdk_display_get_screen", object, screen.num, PACKAGE = "RGtk2")

  return(w)
} 


gdkDisplayGetDefaultScreen <-
function(object)
{
  checkPtrType(object, "GdkDisplay")

  w <- .RGtkCall("S_gdk_display_get_default_screen", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkDisplayPointerUngrab <-
function(object, time. = "GDK_CURRENT_TIME")
{
  checkPtrType(object, "GdkDisplay")
  time. <- as.numeric(time.)

  w <- .RGtkCall("S_gdk_display_pointer_ungrab", object, time., PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDisplayKeyboardUngrab <-
function(object, time. = "GDK_CURRENT_TIME")
{
  checkPtrType(object, "GdkDisplay")
  time. <- as.numeric(time.)

  w <- .RGtkCall("S_gdk_display_keyboard_ungrab", object, time., PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDisplayPointerIsGrabbed <-
function(object)
{
  checkPtrType(object, "GdkDisplay")

  w <- .RGtkCall("S_gdk_display_pointer_is_grabbed", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkDisplayBeep <-
function(object)
{
  checkPtrType(object, "GdkDisplay")

  w <- .RGtkCall("S_gdk_display_beep", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDisplaySync <-
function(object)
{
  checkPtrType(object, "GdkDisplay")

  w <- .RGtkCall("S_gdk_display_sync", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDisplayClose <-
function(object)
{
  checkPtrType(object, "GdkDisplay")

  w <- .RGtkCall("S_gdk_display_close", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDisplayListDevices <-
function(object)
{
  checkPtrType(object, "GdkDisplay")

  w <- .RGtkCall("S_gdk_display_list_devices", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkDisplayGetEvent <-
function(object)
{
  checkPtrType(object, "GdkDisplay")

  w <- .RGtkCall("S_gdk_display_get_event", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkDisplayPeekEvent <-
function(object)
{
  checkPtrType(object, "GdkDisplay")

  w <- .RGtkCall("S_gdk_display_peek_event", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkDisplayPutEvent <-
function(object, event)
{
  checkPtrType(object, "GdkDisplay")
  checkPtrType(event, "GdkEvent")

  w <- .RGtkCall("S_gdk_display_put_event", object, event, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDisplayAddClientMessageFilter <-
function(object, message.type, func, data)
{
  checkPtrType(object, "GdkDisplay")
  message.type <- as.GdkAtom(message.type)
  func <- as.function(func)
  

  w <- .RGtkCall("S_gdk_display_add_client_message_filter", object, message.type, func, data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDisplaySetDoubleClickTime <-
function(object, msec)
{
  checkPtrType(object, "GdkDisplay")
  msec <- as.numeric(msec)

  w <- .RGtkCall("S_gdk_display_set_double_click_time", object, msec, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDisplayGetDefault <-
function()
{
  

  w <- .RGtkCall("S_gdk_display_get_default", PACKAGE = "RGtk2")

  return(w)
} 


gdkDisplayGetCorePointer <-
function(object)
{
  checkPtrType(object, "GdkDisplay")

  w <- .RGtkCall("S_gdk_display_get_core_pointer", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkDisplayGetPointer <-
function(object)
{
  checkPtrType(object, "GdkDisplay")

  w <- .RGtkCall("S_gdk_display_get_pointer", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDisplayGetWindowAtPointer <-
function(object)
{
  checkPtrType(object, "GdkDisplay")

  w <- .RGtkCall("S_gdk_display_get_window_at_pointer", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkDisplayWarpPointer <-
function(object, screen, x, y)
{
  checkPtrType(object, "GdkDisplay")
  checkPtrType(screen, "GdkScreen")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_gdk_display_warp_pointer", object, screen, x, y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDisplayStoreClipboard <-
function(object, clipboard.window, targets)
{
  checkPtrType(object, "GdkDisplay")
  checkPtrType(clipboard.window, "GdkWindow")
  targets <- lapply(targets, function(x) { x <- as.GdkAtom(x); x })

  w <- .RGtkCall("S_gdk_display_store_clipboard", object, clipboard.window, targets, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDisplaySupportsSelectionNotification <-
function(object)
{
  checkPtrType(object, "GdkDisplay")

  w <- .RGtkCall("S_gdk_display_supports_selection_notification", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkDisplayRequestSelectionNotification <-
function(object, selection)
{
  checkPtrType(object, "GdkDisplay")
  selection <- as.GdkAtom(selection)

  w <- .RGtkCall("S_gdk_display_request_selection_notification", object, selection, PACKAGE = "RGtk2")

  return(w)
} 


gdkDisplaySupportsClipboardPersistence <-
function(object)
{
  checkPtrType(object, "GdkDisplay")

  w <- .RGtkCall("S_gdk_display_supports_clipboard_persistence", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkDisplayManagerGetType <-
function()
{
  

  w <- .RGtkCall("S_gdk_display_manager_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gdkDisplayManagerGet <-
function()
{
  

  w <- .RGtkCall("S_gdk_display_manager_get", PACKAGE = "RGtk2")

  return(w)
} 


gdkDisplayManagerGetDefaultDisplay <-
function(object)
{
  checkPtrType(object, "GdkDisplayManager")

  w <- .RGtkCall("S_gdk_display_manager_get_default_display", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkDisplayManagerSetDefaultDisplay <-
function(object, display)
{
  checkPtrType(object, "GdkDisplayManager")
  checkPtrType(display, "GdkDisplay")

  w <- .RGtkCall("S_gdk_display_manager_set_default_display", object, display, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDisplayManagerListDisplays <-
function(object)
{
  checkPtrType(object, "GdkDisplayManager")

  w <- .RGtkCall("S_gdk_display_manager_list_displays", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkDisplayFlush <-
function(object)
{
  checkPtrType(object, "GdkDisplay")

  w <- .RGtkCall("S_gdk_display_flush", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDisplaySetDoubleClickDistance <-
function(object, distance)
{
  checkPtrType(object, "GdkDisplay")
  distance <- as.numeric(distance)

  w <- .RGtkCall("S_gdk_display_set_double_click_distance", object, distance, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDisplaySupportsCursorAlpha <-
function(object)
{
  checkPtrType(object, "GdkDisplay")

  w <- .RGtkCall("S_gdk_display_supports_cursor_alpha", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkDisplaySupportsCursorColor <-
function(object)
{
  checkPtrType(object, "GdkDisplay")

  w <- .RGtkCall("S_gdk_display_supports_cursor_color", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkDisplayGetDefaultCursorSize <-
function(object)
{
  checkPtrType(object, "GdkDisplay")

  w <- .RGtkCall("S_gdk_display_get_default_cursor_size", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkDisplayGetMaximalCursorSize <-
function(object)
{
  checkPtrType(object, "GdkDisplay")

  w <- .RGtkCall("S_gdk_display_get_maximal_cursor_size", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDisplayGetDefaultGroup <-
function(object)
{
  checkPtrType(object, "GdkDisplay")

  w <- .RGtkCall("S_gdk_display_get_default_group", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkDragContextGetType <-
function()
{
  

  w <- .RGtkCall("S_gdk_drag_context_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gdkDragContextNew <-
function()
{
  

  w <- .RGtkCall("S_gdk_drag_context_new", PACKAGE = "RGtk2")

  return(w)
} 


gdkDragStatus <-
function(object, action, time = "GDK_CURRENT_TIME")
{
  checkPtrType(object, "GdkDragContext")
  
  time <- as.numeric(time)

  w <- .RGtkCall("S_gdk_drag_status", object, action, time, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDropReply <-
function(object, ok, time = "GDK_CURRENT_TIME")
{
  checkPtrType(object, "GdkDragContext")
  ok <- as.logical(ok)
  time <- as.numeric(time)

  w <- .RGtkCall("S_gdk_drop_reply", object, ok, time, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDropFinish <-
function(object, success, time = "GDK_CURRENT_TIME")
{
  checkPtrType(object, "GdkDragContext")
  success <- as.logical(success)
  time <- as.numeric(time)

  w <- .RGtkCall("S_gdk_drop_finish", object, success, time, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDragGetSelection <-
function(object)
{
  checkPtrType(object, "GdkDragContext")

  w <- .RGtkCall("S_gdk_drag_get_selection", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkDragBegin <-
function(object, targets)
{
  checkPtrType(object, "GdkWindow")
  targets <- as.GList(targets)

  w <- .RGtkCall("S_gdk_drag_begin", object, targets, PACKAGE = "RGtk2")

  return(w)
} 


gdkDragGetProtocol <-
function(xid)
{
  xid <- as.numeric(xid)

  w <- .RGtkCall("S_gdk_drag_get_protocol", xid, PACKAGE = "RGtk2")

  return(w)
} 


gdkDragFindWindow <-
function(object, drag.window, x.root, y.root)
{
  checkPtrType(object, "GdkDragContext")
  checkPtrType(drag.window, "GdkWindow")
  x.root <- as.integer(x.root)
  y.root <- as.integer(y.root)

  w <- .RGtkCall("S_gdk_drag_find_window", object, drag.window, x.root, y.root, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDragGetProtocolForDisplay <-
function(display, xid)
{
  checkPtrType(display, "GdkDisplay")
  xid <- as.numeric(xid)

  w <- .RGtkCall("S_gdk_drag_get_protocol_for_display", display, xid, PACKAGE = "RGtk2")

  return(w)
} 


gdkDragFindWindowForScreen <-
function(object, drag.window, screen, x.root, y.root)
{
  checkPtrType(object, "GdkDragContext")
  checkPtrType(drag.window, "GdkWindow")
  checkPtrType(screen, "GdkScreen")
  x.root <- as.integer(x.root)
  y.root <- as.integer(y.root)

  w <- .RGtkCall("S_gdk_drag_find_window_for_screen", object, drag.window, screen, x.root, y.root, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDragMotion <-
function(object, dest.window, protocol, x.root, y.root, suggested.action, possible.actions, time)
{
  checkPtrType(object, "GdkDragContext")
  checkPtrType(dest.window, "GdkWindow")
  
  x.root <- as.integer(x.root)
  y.root <- as.integer(y.root)
  
  
  time <- as.numeric(time)

  w <- .RGtkCall("S_gdk_drag_motion", object, dest.window, protocol, x.root, y.root, suggested.action, possible.actions, time, PACKAGE = "RGtk2")

  return(w)
} 


gdkDragDrop <-
function(object, time)
{
  checkPtrType(object, "GdkDragContext")
  time <- as.numeric(time)

  w <- .RGtkCall("S_gdk_drag_drop", object, time, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDragAbort <-
function(object, time)
{
  checkPtrType(object, "GdkDragContext")
  time <- as.numeric(time)

  w <- .RGtkCall("S_gdk_drag_abort", object, time, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDragDropSucceeded <-
function(object)
{
  checkPtrType(object, "GdkDragContext")

  w <- .RGtkCall("S_gdk_drag_drop_succeeded", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkDrawableGetType <-
function()
{
  

  w <- .RGtkCall("S_gdk_drawable_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gdkDrawableSetData <-
function(object, key, data)
{
  checkPtrType(object, "GdkDrawable")
  key <- as.character(key)
  

  w <- .RGtkCall("S_gdk_drawable_set_data", object, key, data, PACKAGE = "RGtk2")

  return(w)
} 


gdkDrawableGetData <-
function(object, key)
{
  checkPtrType(object, "GdkDrawable")
  key <- as.character(key)

  w <- .RGtkCall("S_gdk_drawable_get_data", object, key, PACKAGE = "RGtk2")

  return(w)
} 


gdkDrawableGetSize <-
function(object)
{
  checkPtrType(object, "GdkDrawable")

  w <- .RGtkCall("S_gdk_drawable_get_size", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDrawableSetColormap <-
function(object, colormap)
{
  checkPtrType(object, "GdkDrawable")
  checkPtrType(colormap, "GdkColormap")

  w <- .RGtkCall("S_gdk_drawable_set_colormap", object, colormap, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDrawableGetColormap <-
function(object)
{
  checkPtrType(object, "GdkDrawable")

  w <- .RGtkCall("S_gdk_drawable_get_colormap", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkDrawableGetVisual <-
function(object)
{
  checkPtrType(object, "GdkDrawable")

  w <- .RGtkCall("S_gdk_drawable_get_visual", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkDrawableGetDepth <-
function(object)
{
  checkPtrType(object, "GdkDrawable")

  w <- .RGtkCall("S_gdk_drawable_get_depth", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkDrawableGetScreen <-
function(object)
{
  checkPtrType(object, "GdkDrawable")

  w <- .RGtkCall("S_gdk_drawable_get_screen", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkDrawableGetDisplay <-
function(object)
{
  checkPtrType(object, "GdkDrawable")

  w <- .RGtkCall("S_gdk_drawable_get_display", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkDrawPoint <-
function(object, gc, x, y)
{
  checkPtrType(object, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_gdk_draw_point", object, gc, x, y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDrawLine <-
function(object, gc, x1, y1, x2, y2)
{
  checkPtrType(object, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  x1 <- as.integer(x1)
  y1 <- as.integer(y1)
  x2 <- as.integer(x2)
  y2 <- as.integer(y2)

  w <- .RGtkCall("S_gdk_draw_line", object, gc, x1, y1, x2, y2, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDrawRectangle <-
function(object, gc, filled, x, y, width, height)
{
  checkPtrType(object, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  filled <- as.logical(filled)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gdk_draw_rectangle", object, gc, filled, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDrawArc <-
function(object, gc, filled, x, y, width, height, angle1, angle2)
{
  checkPtrType(object, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  filled <- as.logical(filled)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)
  angle1 <- as.integer(angle1)
  angle2 <- as.integer(angle2)

  w <- .RGtkCall("S_gdk_draw_arc", object, gc, filled, x, y, width, height, angle1, angle2, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDrawPolygon <-
function(object, gc, filled, points)
{
  checkPtrType(object, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  filled <- as.logical(filled)
  points <- lapply(points, function(x) { x <- as.GdkPoint(x); x })

  w <- .RGtkCall("S_gdk_draw_polygon", object, gc, filled, points, PACKAGE = "RGtk2")

  return(w)
} 


gdkDrawString <-
function(object, font, gc, x, y, string)
{
  if(getOption("depwarn"))
    .Deprecated("gdkDrawableDrawLayout", "RGtk2")

  checkPtrType(object, "GdkDrawable")
  checkPtrType(font, "GdkFont")
  checkPtrType(gc, "GdkGC")
  x <- as.integer(x)
  y <- as.integer(y)
  string <- as.character(string)

  w <- .RGtkCall("S_gdk_draw_string", object, font, gc, x, y, string, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDrawText <-
function(object, font, gc, x, y, text, text.length)
{
  if(getOption("depwarn"))
    .Deprecated("gdkDrawableDrawLayout", "RGtk2")

  checkPtrType(object, "GdkDrawable")
  checkPtrType(font, "GdkFont")
  checkPtrType(gc, "GdkGC")
  x <- as.integer(x)
  y <- as.integer(y)
  text <- as.character(text)
  text.length <- as.integer(text.length)

  w <- .RGtkCall("S_gdk_draw_text", object, font, gc, x, y, text, text.length, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDrawTextWc <-
function(object, font, gc, x, text)
{
  if(getOption("depwarn"))
    .Deprecated("gdkDrawableDrawLayout", "RGtk2")

  checkPtrType(object, "GdkDrawable")
  checkPtrType(font, "GdkFont")
  checkPtrType(gc, "GdkGC")
  x <- as.integer(x)
  text <- as.list(as.numeric(text))

  w <- .RGtkCall("S_gdk_draw_text_wc", object, font, gc, x, text, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDrawDrawable <-
function(object, gc, src, xsrc, ysrc, xdest, ydest, width, height)
{
  checkPtrType(object, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  checkPtrType(src, "GdkDrawable")
  xsrc <- as.integer(xsrc)
  ysrc <- as.integer(ysrc)
  xdest <- as.integer(xdest)
  ydest <- as.integer(ydest)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gdk_draw_drawable", object, gc, src, xsrc, ysrc, xdest, ydest, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDrawImage <-
function(object, gc, image, xsrc, ysrc, xdest, ydest, width, height)
{
  checkPtrType(object, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  checkPtrType(image, "GdkImage")
  xsrc <- as.integer(xsrc)
  ysrc <- as.integer(ysrc)
  xdest <- as.integer(xdest)
  ydest <- as.integer(ydest)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gdk_draw_image", object, gc, image, xsrc, ysrc, xdest, ydest, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDrawPoints <-
function(object, gc, points)
{
  checkPtrType(object, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  points <- lapply(points, function(x) { x <- as.GdkPoint(x); x })

  w <- .RGtkCall("S_gdk_draw_points", object, gc, points, PACKAGE = "RGtk2")

  return(w)
} 


gdkDrawSegments <-
function(object, gc, segs)
{
  checkPtrType(object, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  segs <- lapply(segs, function(x) { x <- as.GdkSegment(x); x })

  w <- .RGtkCall("S_gdk_draw_segments", object, gc, segs, PACKAGE = "RGtk2")

  return(w)
} 


gdkDrawLines <-
function(object, gc, points)
{
  checkPtrType(object, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  points <- lapply(points, function(x) { x <- as.GdkPoint(x); x })

  w <- .RGtkCall("S_gdk_draw_lines", object, gc, points, PACKAGE = "RGtk2")

  return(w)
} 


gdkDrawPixbuf <-
function(object, gc = NULL, pixbuf, src.x, src.y, dest.x, dest.y, width = -1, height = -1, dither = "GDK_RGB_DITHER_NORMAL", x.dither = 0, y.dither = 0)
{
  checkPtrType(object, "GdkDrawable")
  if (!is.null( gc )) checkPtrType(gc, "GdkGC")
  checkPtrType(pixbuf, "GdkPixbuf")
  src.x <- as.integer(src.x)
  src.y <- as.integer(src.y)
  dest.x <- as.integer(dest.x)
  dest.y <- as.integer(dest.y)
  width <- as.integer(width)
  height <- as.integer(height)
  
  x.dither <- as.integer(x.dither)
  y.dither <- as.integer(y.dither)

  w <- .RGtkCall("S_gdk_draw_pixbuf", object, gc, pixbuf, src.x, src.y, dest.x, dest.y, width, height, dither, x.dither, y.dither, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDrawGlyphs <-
function(object, gc, font, x, y, glyphs)
{
  checkPtrType(object, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  checkPtrType(font, "PangoFont")
  x <- as.integer(x)
  y <- as.integer(y)
  checkPtrType(glyphs, "PangoGlyphString")

  w <- .RGtkCall("S_gdk_draw_glyphs", object, gc, font, x, y, glyphs, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDrawLayoutLine <-
function(object, gc, x, y, line)
{
  checkPtrType(object, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  x <- as.integer(x)
  y <- as.integer(y)
  checkPtrType(line, "PangoLayoutLine")

  w <- .RGtkCall("S_gdk_draw_layout_line", object, gc, x, y, line, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDrawLayout <-
function(object, gc, x, y, layout)
{
  checkPtrType(object, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  x <- as.integer(x)
  y <- as.integer(y)
  checkPtrType(layout, "PangoLayout")

  w <- .RGtkCall("S_gdk_draw_layout", object, gc, x, y, layout, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDrawLayoutLineWithColors <-
function(drawable, gc, x, y, line, foreground, background)
{
  checkPtrType(drawable, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  x <- as.integer(x)
  y <- as.integer(y)
  checkPtrType(line, "PangoLayoutLine")
  foreground <- as.GdkColor(foreground)
  background <- as.GdkColor(background)

  w <- .RGtkCall("S_gdk_draw_layout_line_with_colors", drawable, gc, x, y, line, foreground, background, PACKAGE = "RGtk2")

  return(w)
} 


gdkDrawLayoutWithColors <-
function(drawable, gc, x, y, layout, foreground, background)
{
  if(getOption("depwarn"))
    .Deprecated("gdkDrawableDrawLayout", "RGtk2")

  checkPtrType(drawable, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  x <- as.integer(x)
  y <- as.integer(y)
  checkPtrType(layout, "PangoLayout")
  foreground <- as.GdkColor(foreground)
  background <- as.GdkColor(background)

  w <- .RGtkCall("S_gdk_draw_layout_with_colors", drawable, gc, x, y, layout, foreground, background, PACKAGE = "RGtk2")

  return(w)
} 


gdkDrawGlyphsTransformed <-
function(drawable, gc, matrix, font, x, y, glyphs)
{
  checkPtrType(drawable, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  checkPtrType(matrix, "PangoMatrix")
  checkPtrType(font, "PangoFont")
  x <- as.integer(x)
  y <- as.integer(y)
  checkPtrType(glyphs, "PangoGlyphString")

  w <- .RGtkCall("S_gdk_draw_glyphs_transformed", drawable, gc, matrix, font, x, y, glyphs, PACKAGE = "RGtk2")

  return(w)
} 


gdkDrawTrapezoids <-
function(drawable, gc, trapezoids)
{
  checkPtrType(drawable, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  trapezoids <- lapply(trapezoids, function(x) { x <- as.GdkTrapezoid(x); x })

  w <- .RGtkCall("S_gdk_draw_trapezoids", drawable, gc, trapezoids, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDrawableGetImage <-
function(object, x, y, width, height)
{
  checkPtrType(object, "GdkDrawable")
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gdk_drawable_get_image", object, x, y, width, height, PACKAGE = "RGtk2")

  return(w)
} 


gdkDrawableCopyToImage <-
function(object, image = NULL, src.x, src.y, dest.x, dest.y, width, height)
{
  checkPtrType(object, "GdkDrawable")
  if (!is.null( image )) checkPtrType(image, "GdkImage")
  src.x <- as.integer(src.x)
  src.y <- as.integer(src.y)
  dest.x <- as.integer(dest.x)
  dest.y <- as.integer(dest.y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gdk_drawable_copy_to_image", object, image, src.x, src.y, dest.x, dest.y, width, height, PACKAGE = "RGtk2")

  return(w)
} 


gdkDrawableGetClipRegion <-
function(object)
{
  checkPtrType(object, "GdkDrawable")

  w <- .RGtkCall("S_gdk_drawable_get_clip_region", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkDrawableGetVisibleRegion <-
function(object)
{
  checkPtrType(object, "GdkDrawable")

  w <- .RGtkCall("S_gdk_drawable_get_visible_region", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkEventGetType <-
function()
{
  

  w <- .RGtkCall("S_gdk_event_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gdkEventsPending <-
function()
{
  

  w <- .RGtkCall("S_gdk_events_pending", PACKAGE = "RGtk2")

  return(w)
} 


gdkEventGet <-
function()
{
  

  w <- .RGtkCall("S_gdk_event_get", PACKAGE = "RGtk2")

  return(w)
} 


gdkEventPeek <-
function()
{
  

  w <- .RGtkCall("S_gdk_event_peek", PACKAGE = "RGtk2")

  return(w)
} 


gdkEventGetGraphicsExpose <-
function(window)
{
  checkPtrType(window, "GdkWindow")

  w <- .RGtkCall("S_gdk_event_get_graphics_expose", window, PACKAGE = "RGtk2")

  return(w)
} 


gdkEventPut <-
function(object)
{
  checkPtrType(object, "GdkEvent")

  w <- .RGtkCall("S_gdk_event_put", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkEventNew <-
function(type)
{
  

  w <- .RGtkCall("S_gdk_event_new", type, PACKAGE = "RGtk2")

  return(w)
} 


gdkEventCopy <-
function(object)
{
  checkPtrType(object, "GdkEvent")

  w <- .RGtkCall("S_gdk_event_copy", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkEventGetTime <-
function(object)
{
  checkPtrType(object, "GdkEvent")

  w <- .RGtkCall("S_gdk_event_get_time", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkEventGetState <-
function(object)
{
  checkPtrType(object, "GdkEvent")

  w <- .RGtkCall("S_gdk_event_get_state", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkEventGetCoords <-
function(object)
{
  checkPtrType(object, "GdkEvent")

  w <- .RGtkCall("S_gdk_event_get_coords", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkEventGetRootCoords <-
function(object)
{
  checkPtrType(object, "GdkEvent")

  w <- .RGtkCall("S_gdk_event_get_root_coords", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkEventGetAxis <-
function(object, axis.use)
{
  checkPtrType(object, "GdkEvent")
  

  w <- .RGtkCall("S_gdk_event_get_axis", object, axis.use, PACKAGE = "RGtk2")

  return(w)
} 


gdkEventHandlerSet <-
function(func, data)
{
  func <- as.function(func)
  

  w <- .RGtkCall("S_gdk_event_handler_set", func, data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkEventSetScreen <-
function(object, screen)
{
  checkPtrType(object, "GdkEvent")
  checkPtrType(screen, "GdkScreen")

  w <- .RGtkCall("S_gdk_event_set_screen", object, screen, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkEventGetScreen <-
function(object)
{
  checkPtrType(object, "GdkEvent")

  w <- .RGtkCall("S_gdk_event_get_screen", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkSetShowEvents <-
function(show.events)
{
  show.events <- as.logical(show.events)

  w <- .RGtkCall("S_gdk_set_show_events", show.events, PACKAGE = "RGtk2")

  return(w)
} 


gdkGetShowEvents <-
function()
{
  

  w <- .RGtkCall("S_gdk_get_show_events", PACKAGE = "RGtk2")

  return(w)
} 


gdkAddClientMessageFilter <-
function(message.type, func, data)
{
  message.type <- as.GdkAtom(message.type)
  func <- as.function(func)
  

  w <- .RGtkCall("S_gdk_add_client_message_filter", message.type, func, data, PACKAGE = "RGtk2")

  return(w)
} 


gdkSettingGet <-
function(name)
{
  name <- as.character(name)

  w <- .RGtkCall("S_gdk_setting_get", name, PACKAGE = "RGtk2")

  return(w)
} 


gdkFontId <-
function(object)
{
  checkPtrType(object, "GdkFont")

  w <- .RGtkCall("S_gdk_font_id", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkFontLoadForDisplay <-
function(display, font.name)
{
  checkPtrType(display, "GdkDisplay")
  font.name <- as.character(font.name)

  w <- .RGtkCall("S_gdk_font_load_for_display", display, font.name, PACKAGE = "RGtk2")

  return(w)
} 


gdkFontsetLoadForDisplay <-
function(display, fontset.name)
{
  checkPtrType(display, "GdkDisplay")
  fontset.name <- as.character(fontset.name)

  w <- .RGtkCall("S_gdk_fontset_load_for_display", display, fontset.name, PACKAGE = "RGtk2")

  return(w)
} 


gdkFontFromDescriptionForDisplay <-
function(display, font.desc)
{
  checkPtrType(display, "GdkDisplay")
  checkPtrType(font.desc, "PangoFontDescription")

  w <- .RGtkCall("S_gdk_font_from_description_for_display", display, font.desc, PACKAGE = "RGtk2")

  return(w)
} 


gdkFontLoad <-
function(font.name)
{
  font.name <- as.character(font.name)

  w <- .RGtkCall("S_gdk_font_load", font.name, PACKAGE = "RGtk2")

  return(w)
} 


gdkFontsetLoad <-
function(fontset.name)
{
  fontset.name <- as.character(fontset.name)

  w <- .RGtkCall("S_gdk_fontset_load", fontset.name, PACKAGE = "RGtk2")

  return(w)
} 


gdkFontFromDescription <-
function(font.desc)
{
  checkPtrType(font.desc, "PangoFontDescription")

  w <- .RGtkCall("S_gdk_font_from_description", font.desc, PACKAGE = "RGtk2")

  return(w)
} 


gdkStringWidth <-
function(object, string)
{
  checkPtrType(object, "GdkFont")
  string <- as.character(string)

  w <- .RGtkCall("S_gdk_string_width", object, string, PACKAGE = "RGtk2")

  return(w)
} 


gdkTextWidth <-
function(object, text, text.length = -1)
{
  checkPtrType(object, "GdkFont")
  text <- as.character(text)
  text.length <- as.integer(text.length)

  w <- .RGtkCall("S_gdk_text_width", object, text, text.length, PACKAGE = "RGtk2")

  return(w)
} 


gdkTextWidthWc <-
function(object, text)
{
  checkPtrType(object, "GdkFont")
  text <- as.list(as.numeric(text))

  w <- .RGtkCall("S_gdk_text_width_wc", object, text, PACKAGE = "RGtk2")

  return(w)
} 


gdkCharWidth <-
function(object, character)
{
  checkPtrType(object, "GdkFont")
  character <- as.character(character)

  w <- .RGtkCall("S_gdk_char_width", object, character, PACKAGE = "RGtk2")

  return(w)
} 


gdkCharWidthWc <-
function(object, character)
{
  checkPtrType(object, "GdkFont")
  character <- as.numeric(character)

  w <- .RGtkCall("S_gdk_char_width_wc", object, character, PACKAGE = "RGtk2")

  return(w)
} 


gdkStringMeasure <-
function(object, string)
{
  checkPtrType(object, "GdkFont")
  string <- as.character(string)

  w <- .RGtkCall("S_gdk_string_measure", object, string, PACKAGE = "RGtk2")

  return(w)
} 


gdkTextMeasure <-
function(object, text, text.length = -1)
{
  checkPtrType(object, "GdkFont")
  text <- as.character(text)
  text.length <- as.integer(text.length)

  w <- .RGtkCall("S_gdk_text_measure", object, text, text.length, PACKAGE = "RGtk2")

  return(w)
} 


gdkCharMeasure <-
function(object, character)
{
  checkPtrType(object, "GdkFont")
  character <- as.character(character)

  w <- .RGtkCall("S_gdk_char_measure", object, character, PACKAGE = "RGtk2")

  return(w)
} 


gdkStringHeight <-
function(object, string)
{
  checkPtrType(object, "GdkFont")
  string <- as.character(string)

  w <- .RGtkCall("S_gdk_string_height", object, string, PACKAGE = "RGtk2")

  return(w)
} 


gdkTextHeight <-
function(object, text, text.length = -1)
{
  checkPtrType(object, "GdkFont")
  text <- as.character(text)
  text.length <- as.integer(text.length)

  w <- .RGtkCall("S_gdk_text_height", object, text, text.length, PACKAGE = "RGtk2")

  return(w)
} 


gdkCharHeight <-
function(object, character)
{
  checkPtrType(object, "GdkFont")
  character <- as.character(character)

  w <- .RGtkCall("S_gdk_char_height", object, character, PACKAGE = "RGtk2")

  return(w)
} 


gdkTextExtentsWc <-
function(object, text)
{
  checkPtrType(object, "GdkFont")
  text <- as.list(as.numeric(text))

  w <- .RGtkCall("S_gdk_text_extents_wc", object, text, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkStringExtents <-
function(object, string)
{
  checkPtrType(object, "GdkFont")
  string <- as.character(string)

  w <- .RGtkCall("S_gdk_string_extents", object, string, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkFontGetDisplay <-
function(object)
{
  checkPtrType(object, "GdkFont")

  w <- .RGtkCall("S_gdk_font_get_display", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkGCGetType <-
function()
{
  

  w <- .RGtkCall("S_gdk_gc_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gdkGCNew <-
function(drawable)
{
  checkPtrType(drawable, "GdkDrawable")

  w <- .RGtkCall("S_gdk_gc_new", drawable, PACKAGE = "RGtk2")

  return(w)
} 


gdkGCNewWithValues <-
function(object, values)
{
  checkPtrType(object, "GdkDrawable")
  values <- as.GdkGCValues(values)

  w <- .RGtkCall("S_gdk_gc_new_with_values", object, values, PACKAGE = "RGtk2")

  return(w)
} 


gdkGCGetValues <-
function(object)
{
  checkPtrType(object, "GdkGC")

  w <- .RGtkCall("S_gdk_gc_get_values", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkGCSetValues <-
function(object, values)
{
  checkPtrType(object, "GdkGC")
  values <- as.GdkGCValues(values)

  w <- .RGtkCall("S_gdk_gc_set_values", object, values, PACKAGE = "RGtk2")

  return(w)
} 


gdkGCSetForeground <-
function(object, color)
{
  checkPtrType(object, "GdkGC")
  color <- as.GdkColor(color)

  w <- .RGtkCall("S_gdk_gc_set_foreground", object, color, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkGCSetBackground <-
function(object, color)
{
  checkPtrType(object, "GdkGC")
  color <- as.GdkColor(color)

  w <- .RGtkCall("S_gdk_gc_set_background", object, color, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkGCSetFont <-
function(object, font)
{
  checkPtrType(object, "GdkGC")
  checkPtrType(font, "GdkFont")

  w <- .RGtkCall("S_gdk_gc_set_font", object, font, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkGCSetFunction <-
function(object, fun)
{
  checkPtrType(object, "GdkGC")
  

  w <- .RGtkCall("S_gdk_gc_set_function", object, fun, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkGCSetFill <-
function(object, fill)
{
  checkPtrType(object, "GdkGC")
  

  w <- .RGtkCall("S_gdk_gc_set_fill", object, fill, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkGCSetTile <-
function(object, tile)
{
  checkPtrType(object, "GdkGC")
  checkPtrType(tile, "GdkPixmap")

  w <- .RGtkCall("S_gdk_gc_set_tile", object, tile, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkGCSetStipple <-
function(object, stipple)
{
  checkPtrType(object, "GdkGC")
  checkPtrType(stipple, "GdkPixmap")

  w <- .RGtkCall("S_gdk_gc_set_stipple", object, stipple, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkGCSetTsOrigin <-
function(object, x, y)
{
  checkPtrType(object, "GdkGC")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_gdk_gc_set_ts_origin", object, x, y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkGCSetClipOrigin <-
function(object, x, y)
{
  checkPtrType(object, "GdkGC")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_gdk_gc_set_clip_origin", object, x, y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkGCSetClipMask <-
function(object, mask)
{
  checkPtrType(object, "GdkGC")
  checkPtrType(mask, "GdkBitmap")

  w <- .RGtkCall("S_gdk_gc_set_clip_mask", object, mask, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkGCSetClipRectangle <-
function(object, rectangle)
{
  checkPtrType(object, "GdkGC")
  rectangle <- as.GdkRectangle(rectangle)

  w <- .RGtkCall("S_gdk_gc_set_clip_rectangle", object, rectangle, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkGCSetClipRegion <-
function(object, region)
{
  checkPtrType(object, "GdkGC")
  checkPtrType(region, "GdkRegion")

  w <- .RGtkCall("S_gdk_gc_set_clip_region", object, region, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkGCSetSubwindow <-
function(object, mode)
{
  checkPtrType(object, "GdkGC")
  

  w <- .RGtkCall("S_gdk_gc_set_subwindow", object, mode, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkGCSetExposures <-
function(object, exposures)
{
  checkPtrType(object, "GdkGC")
  exposures <- as.logical(exposures)

  w <- .RGtkCall("S_gdk_gc_set_exposures", object, exposures, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkGCSetLineAttributes <-
function(object, line.width, line.style, cap.style, join.style)
{
  checkPtrType(object, "GdkGC")
  line.width <- as.integer(line.width)
  
  
  

  w <- .RGtkCall("S_gdk_gc_set_line_attributes", object, line.width, line.style, cap.style, join.style, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkGCSetDashes <-
function(object, dash.list)
{
  checkPtrType(object, "GdkGC")
  dash.list <- as.list(as.raw(dash.list))

  w <- .RGtkCall("S_gdk_gc_set_dashes", object, dash.list, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkGCOffset <-
function(object, x.offset, y.offset)
{
  checkPtrType(object, "GdkGC")
  x.offset <- as.integer(x.offset)
  y.offset <- as.integer(y.offset)

  w <- .RGtkCall("S_gdk_gc_offset", object, x.offset, y.offset, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkGCCopy <-
function(object, src.gc)
{
  checkPtrType(object, "GdkGC")
  checkPtrType(src.gc, "GdkGC")

  w <- .RGtkCall("S_gdk_gc_copy", object, src.gc, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkGCSetColormap <-
function(object, colormap)
{
  checkPtrType(object, "GdkGC")
  checkPtrType(colormap, "GdkColormap")

  w <- .RGtkCall("S_gdk_gc_set_colormap", object, colormap, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkGCGetColormap <-
function(object)
{
  checkPtrType(object, "GdkGC")

  w <- .RGtkCall("S_gdk_gc_get_colormap", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkGCSetRgbFgColor <-
function(object, color)
{
  checkPtrType(object, "GdkGC")
  color <- as.GdkColor(color)

  w <- .RGtkCall("S_gdk_gc_set_rgb_fg_color", object, color, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkGCSetRgbBgColor <-
function(object, color)
{
  checkPtrType(object, "GdkGC")
  color <- as.GdkColor(color)

  w <- .RGtkCall("S_gdk_gc_set_rgb_bg_color", object, color, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkGCGetScreen <-
function(object)
{
  checkPtrType(object, "GdkGC")

  w <- .RGtkCall("S_gdk_gc_get_screen", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkImageGetType <-
function()
{
  

  w <- .RGtkCall("S_gdk_image_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gdkImageNew <-
function(type, visual, width, height)
{
  
  checkPtrType(visual, "GdkVisual")
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gdk_image_new", type, visual, width, height, PACKAGE = "RGtk2")

  return(w)
} 


gdkImageGet <-
function(object, x, y, width, height)
{
  if(getOption("depwarn"))
    .Deprecated("gdkDrawableGetImage", "RGtk2")

  checkPtrType(object, "GdkDrawable")
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gdk_image_get", object, x, y, width, height, PACKAGE = "RGtk2")

  return(w)
} 


gdkImagePutPixel <-
function(object, x, y, pixel)
{
  checkPtrType(object, "GdkImage")
  x <- as.integer(x)
  y <- as.integer(y)
  pixel <- as.numeric(pixel)

  w <- .RGtkCall("S_gdk_image_put_pixel", object, x, y, pixel, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkImageGetPixel <-
function(object, x, y)
{
  checkPtrType(object, "GdkImage")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_gdk_image_get_pixel", object, x, y, PACKAGE = "RGtk2")

  return(w)
} 


gdkImageSetColormap <-
function(object, colormap)
{
  checkPtrType(object, "GdkImage")
  checkPtrType(colormap, "GdkColormap")

  w <- .RGtkCall("S_gdk_image_set_colormap", object, colormap, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkImageGetColormap <-
function(object)
{
  checkPtrType(object, "GdkImage")

  w <- .RGtkCall("S_gdk_image_get_colormap", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkDeviceGetType <-
function()
{
  

  w <- .RGtkCall("S_gdk_device_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gdkDevicesList <-
function()
{
  

  w <- .RGtkCall("S_gdk_devices_list", PACKAGE = "RGtk2")

  return(w)
} 


gdkDeviceSetSource <-
function(object, source)
{
  checkPtrType(object, "GdkDevice")
  

  w <- .RGtkCall("S_gdk_device_set_source", object, source, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDeviceSetMode <-
function(object, mode)
{
  checkPtrType(object, "GdkDevice")
  

  w <- .RGtkCall("S_gdk_device_set_mode", object, mode, PACKAGE = "RGtk2")

  return(w)
} 


gdkDeviceSetKey <-
function(object, index, keyval, modifiers)
{
  checkPtrType(object, "GdkDevice")
  index <- as.numeric(index)
  keyval <- as.numeric(keyval)
  

  w <- .RGtkCall("S_gdk_device_set_key", object, index, keyval, modifiers, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDeviceSetAxisUse <-
function(object, index, use)
{
  checkPtrType(object, "GdkDevice")
  index <- as.numeric(index)
  

  w <- .RGtkCall("S_gdk_device_set_axis_use", object, index, use, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDeviceGetState <-
function(object, window)
{
  checkPtrType(object, "GdkDevice")
  checkPtrType(window, "GdkWindow")

  w <- .RGtkCall("S_gdk_device_get_state", object, window, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDeviceGetHistory <-
function(object, window, start, stop)
{
  checkPtrType(object, "GdkDevice")
  checkPtrType(window, "GdkWindow")
  start <- as.numeric(start)
  stop <- as.numeric(stop)

  w <- .RGtkCall("S_gdk_device_get_history", object, window, start, stop, PACKAGE = "RGtk2")

  return(w)
} 


gdkDeviceGetAxis <-
function(object, axes, use)
{
  checkPtrType(object, "GdkDevice")
  axes <- as.list(as.numeric(axes))
  

  w <- .RGtkCall("S_gdk_device_get_axis", object, axes, use, PACKAGE = "RGtk2")

  return(w)
} 


gdkInputSetExtensionEvents <-
function(object, mask, mode)
{
  checkPtrType(object, "GdkWindow")
  mask <- as.integer(mask)
  

  w <- .RGtkCall("S_gdk_input_set_extension_events", object, mask, mode, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkDeviceGetCorePointer <-
function()
{
  

  w <- .RGtkCall("S_gdk_device_get_core_pointer", PACKAGE = "RGtk2")

  return(w)
} 


gdkKeymapGetType <-
function()
{
  

  w <- .RGtkCall("S_gdk_keymap_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gdkKeymapGetDefault <-
function()
{
  

  w <- .RGtkCall("S_gdk_keymap_get_default", PACKAGE = "RGtk2")

  return(w)
} 


gdkKeymapGetForDisplay <-
function(display)
{
  checkPtrType(display, "GdkDisplay")

  w <- .RGtkCall("S_gdk_keymap_get_for_display", display, PACKAGE = "RGtk2")

  return(w)
} 


gdkKeymapLookupKey <-
function(object, key)
{
  checkPtrType(object, "GdkKeymap")
  key <- as.GdkKeymapKey(key)

  w <- .RGtkCall("S_gdk_keymap_lookup_key", object, key, PACKAGE = "RGtk2")

  return(w)
} 


gdkKeymapTranslateKeyboardState <-
function(object, hardware.keycode, state, group)
{
  checkPtrType(object, "GdkKeymap")
  hardware.keycode <- as.numeric(hardware.keycode)
  
  group <- as.integer(group)

  w <- .RGtkCall("S_gdk_keymap_translate_keyboard_state", object, hardware.keycode, state, group, PACKAGE = "RGtk2")

  return(w)
} 


gdkKeymapGetEntriesForKeyval <-
function(object, keyval)
{
  checkPtrType(object, "GdkKeymap")
  keyval <- as.numeric(keyval)

  w <- .RGtkCall("S_gdk_keymap_get_entries_for_keyval", object, keyval, PACKAGE = "RGtk2")

  return(w)
} 


gdkKeymapGetEntriesForKeycode <-
function(object, hardware.keycode)
{
  checkPtrType(object, "GdkKeymap")
  hardware.keycode <- as.numeric(hardware.keycode)

  w <- .RGtkCall("S_gdk_keymap_get_entries_for_keycode", object, hardware.keycode, PACKAGE = "RGtk2")

  return(w)
} 


gdkKeymapGetDirection <-
function(object)
{
  checkPtrType(object, "GdkKeymap")

  w <- .RGtkCall("S_gdk_keymap_get_direction", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkKeyvalName <-
function(keyval)
{
  keyval <- as.numeric(keyval)

  w <- .RGtkCall("S_gdk_keyval_name", keyval, PACKAGE = "RGtk2")

  return(w)
} 


gdkKeyvalFromName <-
function(keyval.name)
{
  keyval.name <- as.character(keyval.name)

  w <- .RGtkCall("S_gdk_keyval_from_name", keyval.name, PACKAGE = "RGtk2")

  return(w)
} 


gdkKeyvalConvertCase <-
function(symbol)
{
  symbol <- as.numeric(symbol)

  w <- .RGtkCall("S_gdk_keyval_convert_case", symbol, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkKeyvalToUpper <-
function(keyval)
{
  keyval <- as.numeric(keyval)

  w <- .RGtkCall("S_gdk_keyval_to_upper", keyval, PACKAGE = "RGtk2")

  return(w)
} 


gdkKeyvalToLower <-
function(keyval)
{
  keyval <- as.numeric(keyval)

  w <- .RGtkCall("S_gdk_keyval_to_lower", keyval, PACKAGE = "RGtk2")

  return(w)
} 


gdkKeyvalIsUpper <-
function(keyval)
{
  keyval <- as.numeric(keyval)

  w <- .RGtkCall("S_gdk_keyval_is_upper", keyval, PACKAGE = "RGtk2")

  return(w)
} 


gdkKeyvalIsLower <-
function(keyval)
{
  keyval <- as.numeric(keyval)

  w <- .RGtkCall("S_gdk_keyval_is_lower", keyval, PACKAGE = "RGtk2")

  return(w)
} 


gdkKeyvalToUnicode <-
function(keyval)
{
  keyval <- as.numeric(keyval)

  w <- .RGtkCall("S_gdk_keyval_to_unicode", keyval, PACKAGE = "RGtk2")

  return(w)
} 


gdkUnicodeToKeyval <-
function(wc)
{
  wc <- as.numeric(wc)

  w <- .RGtkCall("S_gdk_unicode_to_keyval", wc, PACKAGE = "RGtk2")

  return(w)
} 


gdkPangoRendererGetType <-
function()
{
  

  w <- .RGtkCall("S_gdk_pango_renderer_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gdkPangoRendererNew <-
function(screen)
{
  checkPtrType(screen, "GdkScreen")

  w <- .RGtkCall("S_gdk_pango_renderer_new", screen, PACKAGE = "RGtk2")

  return(w)
} 


gdkPangoRendererGetDefault <-
function(screen)
{
  checkPtrType(screen, "GdkScreen")

  w <- .RGtkCall("S_gdk_pango_renderer_get_default", screen, PACKAGE = "RGtk2")

  return(w)
} 


gdkPangoRendererSetDrawable <-
function(object, drawable = NULL)
{
  checkPtrType(object, "GdkPangoRenderer")
  if (!is.null( drawable )) checkPtrType(drawable, "GdkDrawable")

  w <- .RGtkCall("S_gdk_pango_renderer_set_drawable", object, drawable, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkPangoRendererSetGc <-
function(object, gc = NULL)
{
  checkPtrType(object, "GdkPangoRenderer")
  if (!is.null( gc )) checkPtrType(gc, "GdkGC")

  w <- .RGtkCall("S_gdk_pango_renderer_set_gc", object, gc, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkPangoRendererSetStipple <-
function(object, part, stipple)
{
  checkPtrType(object, "GdkPangoRenderer")
  
  checkPtrType(stipple, "GdkBitmap")

  w <- .RGtkCall("S_gdk_pango_renderer_set_stipple", object, part, stipple, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkPangoRendererSetOverrideColor <-
function(object, part, color = NULL)
{
  checkPtrType(object, "GdkPangoRenderer")
  
  if (!is.null( color )) color <- as.GdkColor(color)

  w <- .RGtkCall("S_gdk_pango_renderer_set_override_color", object, part, color, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkPangoContextGetForScreen <-
function(screen)
{
  checkPtrType(screen, "GdkScreen")

  w <- .RGtkCall("S_gdk_pango_context_get_for_screen", screen, PACKAGE = "RGtk2")

  return(w)
} 


gdkPangoContextGet <-
function()
{
  

  w <- .RGtkCall("S_gdk_pango_context_get", PACKAGE = "RGtk2")

  return(w)
} 


gdkPangoContextSetColormap <-
function(context, colormap)
{
  checkPtrType(context, "PangoContext")
  checkPtrType(colormap, "GdkColormap")

  w <- .RGtkCall("S_gdk_pango_context_set_colormap", context, colormap, PACKAGE = "RGtk2")

  return(w)
} 


gdkPangoLayoutLineGetClipRegion <-
function(line, x.origin, index.ranges)
{
  checkPtrType(line, "PangoLayoutLine")
  x.origin <- as.integer(x.origin)
  index.ranges <- as.list(as.integer(index.ranges))

  w <- .RGtkCall("S_gdk_pango_layout_line_get_clip_region", line, x.origin, index.ranges, PACKAGE = "RGtk2")

  return(w)
} 


gdkPangoLayoutGetClipRegion <-
function(layout, x.origin, index.ranges)
{
  checkPtrType(layout, "PangoLayout")
  x.origin <- as.integer(x.origin)
  index.ranges <- as.list(as.integer(index.ranges))

  w <- .RGtkCall("S_gdk_pango_layout_get_clip_region", layout, x.origin, index.ranges, PACKAGE = "RGtk2")

  return(w)
} 


gdkPangoAttrStippleNew <-
function(stipple)
{
  checkPtrType(stipple, "GdkBitmap")

  w <- .RGtkCall("S_gdk_pango_attr_stipple_new", stipple, PACKAGE = "RGtk2")

  return(w)
} 


gdkPangoAttrEmbossedNew <-
function(embossed)
{
  embossed <- as.logical(embossed)

  w <- .RGtkCall("S_gdk_pango_attr_embossed_new", embossed, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufRenderThresholdAlpha <-
function(object, bitmap, src.x, src.y, dest.x, dest.y, width = -1, height = -1, alpha.threshold)
{
  checkPtrType(object, "GdkPixbuf")
  checkPtrType(bitmap, "GdkBitmap")
  src.x <- as.integer(src.x)
  src.y <- as.integer(src.y)
  dest.x <- as.integer(dest.x)
  dest.y <- as.integer(dest.y)
  width <- as.integer(width)
  height <- as.integer(height)
  alpha.threshold <- as.integer(alpha.threshold)

  w <- .RGtkCall("S_gdk_pixbuf_render_threshold_alpha", object, bitmap, src.x, src.y, dest.x, dest.y, width, height, alpha.threshold, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkPixbufRenderToDrawable <-
function(object, drawable, gc, src.x = 0, src.y = 0, dest.x, dest.y, width = -1, height = -1, dither = "GDK_RGB_DITHER_NORMAL", x.dither = 0, y.dither = 0)
{
  if(getOption("depwarn"))
    .Deprecated("gdkDrawableDrawPixbuf", "RGtk2")

  checkPtrType(object, "GdkPixbuf")
  checkPtrType(drawable, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  src.x <- as.integer(src.x)
  src.y <- as.integer(src.y)
  dest.x <- as.integer(dest.x)
  dest.y <- as.integer(dest.y)
  width <- as.integer(width)
  height <- as.integer(height)
  
  x.dither <- as.integer(x.dither)
  y.dither <- as.integer(y.dither)

  w <- .RGtkCall("S_gdk_pixbuf_render_to_drawable", object, drawable, gc, src.x, src.y, dest.x, dest.y, width, height, dither, x.dither, y.dither, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkPixbufRenderToDrawableAlpha <-
function(object, drawable, src.x = 0, src.y = 0, dest.x, dest.y, width = -1, height = -1, alpha.mode = NULL, alpha.threshold = NULL, dither = "GDK_RGB_DITHER_NORMAL", x.dither = 0, y.dither = 0)
{
  if(getOption("depwarn"))
    .Deprecated("gdkDrawableDrawPixbuf", "RGtk2")

  checkPtrType(object, "GdkPixbuf")
  checkPtrType(drawable, "GdkDrawable")
  src.x <- as.integer(src.x)
  src.y <- as.integer(src.y)
  dest.x <- as.integer(dest.x)
  dest.y <- as.integer(dest.y)
  width <- as.integer(width)
  height <- as.integer(height)
  
  alpha.threshold <- as.integer(alpha.threshold)
  
  x.dither <- as.integer(x.dither)
  y.dither <- as.integer(y.dither)

  w <- .RGtkCall("S_gdk_pixbuf_render_to_drawable_alpha", object, drawable, src.x, src.y, dest.x, dest.y, width, height, alpha.mode, alpha.threshold, dither, x.dither, y.dither, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkPixbufRenderPixmapAndMask <-
function(object, alpha.threshold = 127)
{
  checkPtrType(object, "GdkPixbuf")
  alpha.threshold <- as.integer(alpha.threshold)

  w <- .RGtkCall("S_gdk_pixbuf_render_pixmap_and_mask", object, alpha.threshold, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkPixbufRenderPixmapAndMaskForColormap <-
function(object, colormap, alpha.threshold = 127)
{
  checkPtrType(object, "GdkPixbuf")
  checkPtrType(colormap, "GdkColormap")
  alpha.threshold <- as.integer(alpha.threshold)

  w <- .RGtkCall("S_gdk_pixbuf_render_pixmap_and_mask_for_colormap", object, colormap, alpha.threshold, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkPixbufGetFromDrawable <-
function(dest = NULL, src, cmap = NULL, src.x, src.y, dest.x, dest.y, width, height)
{
  if (!is.null( dest )) checkPtrType(dest, "GdkPixbuf")
  checkPtrType(src, "GdkDrawable")
  if (!is.null( cmap )) checkPtrType(cmap, "GdkColormap")
  src.x <- as.integer(src.x)
  src.y <- as.integer(src.y)
  dest.x <- as.integer(dest.x)
  dest.y <- as.integer(dest.y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gdk_pixbuf_get_from_drawable", dest, src, cmap, src.x, src.y, dest.x, dest.y, width, height, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufGetFromImage <-
function(src, cmap, src.x, src.y, dest.x, dest.y, width, height)
{
  checkPtrType(src, "GdkImage")
  checkPtrType(cmap, "GdkColormap")
  src.x <- as.integer(src.x)
  src.y <- as.integer(src.y)
  dest.x <- as.integer(dest.x)
  dest.y <- as.integer(dest.y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gdk_pixbuf_get_from_image", src, cmap, src.x, src.y, dest.x, dest.y, width, height, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixmapGetType <-
function()
{
  

  w <- .RGtkCall("S_gdk_pixmap_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gdkPixmapNew <-
function(drawable = NULL, width, height, depth = -1)
{
  if (!is.null( drawable )) checkPtrType(drawable, "GdkDrawable")
  width <- as.integer(width)
  height <- as.integer(height)
  depth <- as.integer(depth)

  w <- .RGtkCall("S_gdk_pixmap_new", drawable, width, height, depth, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixmapCreateFromData <-
function(drawable = NULL, data, height, depth, fg, bg)
{
  if (!is.null( drawable )) checkPtrType(drawable, "GdkDrawable")
  data <- as.list(as.raw(data))
  height <- as.integer(height)
  depth <- as.integer(depth)
  fg <- as.GdkColor(fg)
  bg <- as.GdkColor(bg)

  w <- .RGtkCall("S_gdk_pixmap_create_from_data", drawable, data, height, depth, fg, bg, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixmapCreateFromXpm <-
function(drawable, transparent.color, filename)
{
  checkPtrType(drawable, "GdkDrawable")
  transparent.color <- as.GdkColor(transparent.color)
  filename <- as.character(filename)

  w <- .RGtkCall("S_gdk_pixmap_create_from_xpm", drawable, transparent.color, filename, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixmapColormapCreateFromXpm <-
function(drawable, colormap, transparent.color, filename)
{
  checkPtrType(drawable, "GdkDrawable")
  checkPtrType(colormap, "GdkColormap")
  transparent.color <- as.GdkColor(transparent.color)
  filename <- as.character(filename)

  w <- .RGtkCall("S_gdk_pixmap_colormap_create_from_xpm", drawable, colormap, transparent.color, filename, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixmapCreateFromXpmD <-
function(drawable, transparent.color, data)
{
  checkPtrType(drawable, "GdkDrawable")
  transparent.color <- as.GdkColor(transparent.color)
  data <- as.list(as.character(data))

  w <- .RGtkCall("S_gdk_pixmap_create_from_xpm_d", drawable, transparent.color, data, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixmapColormapCreateFromXpmD <-
function(drawable, colormap, transparent.color, data)
{
  checkPtrType(drawable, "GdkDrawable")
  checkPtrType(colormap, "GdkColormap")
  transparent.color <- as.GdkColor(transparent.color)
  data <- as.list(as.character(data))

  w <- .RGtkCall("S_gdk_pixmap_colormap_create_from_xpm_d", drawable, colormap, transparent.color, data, PACKAGE = "RGtk2")

  return(w)
} 


gdkAtomName <-
function(atom)
{
  atom <- as.GdkAtom(atom)

  w <- .RGtkCall("S_gdk_atom_name", atom, PACKAGE = "RGtk2")

  return(w)
} 


gdkAtomIntern <-
function(atom.name, only.if.exists = FALSE)
{
  atom.name <- as.character(atom.name)
  only.if.exists <- as.logical(only.if.exists)

  w <- .RGtkCall("S_gdk_atom_intern", atom.name, only.if.exists, PACKAGE = "RGtk2")

  return(w)
} 


gdkPropertyGet <-
function(object, property, type, offset, length, pdelete)
{
  checkPtrType(object, "GdkWindow")
  property <- as.GdkAtom(property)
  type <- as.GdkAtom(type)
  offset <- as.numeric(offset)
  length <- as.numeric(length)
  pdelete <- as.integer(pdelete)

  w <- .RGtkCall("S_gdk_property_get", object, property, type, offset, length, pdelete, PACKAGE = "RGtk2")

  return(w)
} 


gdkPropertyChange <-
function(object, property, type, format, mode, data)
{
  checkPtrType(object, "GdkWindow")
  property <- as.GdkAtom(property)
  type <- as.GdkAtom(type)
  format <- as.integer(format)
  
  data <- as.list(as.raw(data))

  w <- .RGtkCall("S_gdk_property_change", object, property, type, format, mode, data, PACKAGE = "RGtk2")

  return(w)
} 


gdkPropertyDelete <-
function(object, property)
{
  checkPtrType(object, "GdkWindow")
  property <- as.GdkAtom(property)

  w <- .RGtkCall("S_gdk_property_delete", object, property, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkRgbXpixelFromRgb <-
function(rgb)
{
  rgb <- as.numeric(rgb)

  w <- .RGtkCall("S_gdk_rgb_xpixel_from_rgb", rgb, PACKAGE = "RGtk2")

  return(w)
} 


gdkRgbGcSetForeground <-
function(gc, rgb)
{
  checkPtrType(gc, "GdkGC")
  rgb <- as.numeric(rgb)

  w <- .RGtkCall("S_gdk_rgb_gc_set_foreground", gc, rgb, PACKAGE = "RGtk2")

  return(w)
} 


gdkRgbGcSetBackground <-
function(gc, rgb)
{
  checkPtrType(gc, "GdkGC")
  rgb <- as.numeric(rgb)

  w <- .RGtkCall("S_gdk_rgb_gc_set_background", gc, rgb, PACKAGE = "RGtk2")

  return(w)
} 


gdkDrawRgbImageDithalign <-
function(object, gc, x, y, width, height, dith, rgb.buf, xdith, ydith)
{
  checkPtrType(object, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)
  
  rgb.buf <- as.list(as.raw(rgb.buf))
  xdith <- as.integer(xdith)
  ydith <- as.integer(ydith)

  w <- .RGtkCall("S_gdk_draw_rgb_image_dithalign", object, gc, x, y, width, height, dith, rgb.buf, xdith, ydith, PACKAGE = "RGtk2")

  return(w)
} 


gdkDrawRgb32Image <-
function(object, gc, x, y, width, height, dith, buf)
{
  checkPtrType(object, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)
  
  buf <- as.list(as.raw(buf))

  w <- .RGtkCall("S_gdk_draw_rgb_32_image", object, gc, x, y, width, height, dith, buf, PACKAGE = "RGtk2")

  return(w)
} 


gdkDrawRgb32ImageDithalign <-
function(object, gc, x, y, width, height, dith, buf, xdith, ydith)
{
  checkPtrType(object, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)
  
  buf <- as.list(as.raw(buf))
  xdith <- as.integer(xdith)
  ydith <- as.integer(ydith)

  w <- .RGtkCall("S_gdk_draw_rgb_32_image_dithalign", object, gc, x, y, width, height, dith, buf, xdith, ydith, PACKAGE = "RGtk2")

  return(w)
} 


gdkRgbFindColor <-
function(colormap, color)
{
  checkPtrType(colormap, "GdkColormap")
  color <- as.GdkColor(color)

  w <- .RGtkCall("S_gdk_rgb_find_color", colormap, color, PACKAGE = "RGtk2")

  return(w)
} 


gdkRgbColormapDitherable <-
function(colormap)
{
  checkPtrType(colormap, "GdkColormap")

  w <- .RGtkCall("S_gdk_rgb_colormap_ditherable", colormap, PACKAGE = "RGtk2")

  return(w)
} 


gdkDrawGrayImage <-
function(object, gc, x, y, width, height, dith, buf)
{
  checkPtrType(object, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)
  
  buf <- as.list(as.raw(buf))

  w <- .RGtkCall("S_gdk_draw_gray_image", object, gc, x, y, width, height, dith, buf, PACKAGE = "RGtk2")

  return(w)
} 


gdkRgbCmapNew <-
function(colors)
{
  colors <- as.list(as.numeric(colors))

  w <- .RGtkCall("S_gdk_rgb_cmap_new", colors, PACKAGE = "RGtk2")

  return(w)
} 


gdkDrawIndexedImage <-
function(object, gc, x, y, width, height, dith, buf, cmap)
{
  checkPtrType(object, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)
  
  buf <- as.list(as.raw(buf))
  cmap <- as.GdkRgbCmap(cmap)

  w <- .RGtkCall("S_gdk_draw_indexed_image", object, gc, x, y, width, height, dith, buf, cmap, PACKAGE = "RGtk2")

  return(w)
} 


gdkRgbDitherable <-
function()
{
  

  w <- .RGtkCall("S_gdk_rgb_ditherable", PACKAGE = "RGtk2")

  return(w)
} 


gdkRgbSetVerbose <-
function(verbose)
{
  verbose <- as.logical(verbose)

  w <- .RGtkCall("S_gdk_rgb_set_verbose", verbose, PACKAGE = "RGtk2")

  return(w)
} 


gdkRgbSetInstall <-
function(install)
{
  install <- as.logical(install)

  w <- .RGtkCall("S_gdk_rgb_set_install", install, PACKAGE = "RGtk2")

  return(w)
} 


gdkRgbSetMinColors <-
function(min.colors)
{
  min.colors <- as.integer(min.colors)

  w <- .RGtkCall("S_gdk_rgb_set_min_colors", min.colors, PACKAGE = "RGtk2")

  return(w)
} 


gdkRgbGetColormap <-
function()
{
  

  w <- .RGtkCall("S_gdk_rgb_get_colormap", PACKAGE = "RGtk2")

  return(w)
} 


gdkRgbGetCmap <-
function()
{
  if(getOption("depwarn"))
    .Deprecated("gdkRgbGetColormap", "RGtk2")

  

  w <- .RGtkCall("S_gdk_rgb_get_cmap", PACKAGE = "RGtk2")

  return(w)
} 


gdkRgbGetVisual <-
function()
{
  

  w <- .RGtkCall("S_gdk_rgb_get_visual", PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenGetType <-
function()
{
  

  w <- .RGtkCall("S_gdk_screen_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenGetDefaultColormap <-
function(object)
{
  checkPtrType(object, "GdkScreen")

  w <- .RGtkCall("S_gdk_screen_get_default_colormap", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenSetDefaultColormap <-
function(object, colormap)
{
  checkPtrType(object, "GdkScreen")
  checkPtrType(colormap, "GdkColormap")

  w <- .RGtkCall("S_gdk_screen_set_default_colormap", object, colormap, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkScreenGetSystemColormap <-
function(object)
{
  checkPtrType(object, "GdkScreen")

  w <- .RGtkCall("S_gdk_screen_get_system_colormap", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenGetSystemVisual <-
function(object)
{
  checkPtrType(object, "GdkScreen")

  w <- .RGtkCall("S_gdk_screen_get_system_visual", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenGetRgbColormap <-
function(object)
{
  checkPtrType(object, "GdkScreen")

  w <- .RGtkCall("S_gdk_screen_get_rgb_colormap", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenGetRgbaColormap <-
function(object)
{
  checkPtrType(object, "GdkScreen")

  w <- .RGtkCall("S_gdk_screen_get_rgba_colormap", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenGetRgbVisual <-
function(object)
{
  checkPtrType(object, "GdkScreen")

  w <- .RGtkCall("S_gdk_screen_get_rgb_visual", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenGetRgbaVisual <-
function(object)
{
  checkPtrType(object, "GdkScreen")

  w <- .RGtkCall("S_gdk_screen_get_rgba_visual", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenGetRootWindow <-
function(object)
{
  checkPtrType(object, "GdkScreen")

  w <- .RGtkCall("S_gdk_screen_get_root_window", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenGetDisplay <-
function(object)
{
  checkPtrType(object, "GdkScreen")

  w <- .RGtkCall("S_gdk_screen_get_display", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenGetNumber <-
function(object)
{
  checkPtrType(object, "GdkScreen")

  w <- .RGtkCall("S_gdk_screen_get_number", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenGetWidth <-
function(object)
{
  checkPtrType(object, "GdkScreen")

  w <- .RGtkCall("S_gdk_screen_get_width", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenGetHeight <-
function(object)
{
  checkPtrType(object, "GdkScreen")

  w <- .RGtkCall("S_gdk_screen_get_height", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenGetWidthMm <-
function(object)
{
  checkPtrType(object, "GdkScreen")

  w <- .RGtkCall("S_gdk_screen_get_width_mm", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenGetHeightMm <-
function(object)
{
  checkPtrType(object, "GdkScreen")

  w <- .RGtkCall("S_gdk_screen_get_height_mm", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenListVisuals <-
function(object)
{
  checkPtrType(object, "GdkScreen")

  w <- .RGtkCall("S_gdk_screen_list_visuals", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenGetToplevelWindows <-
function(object)
{
  checkPtrType(object, "GdkScreen")

  w <- .RGtkCall("S_gdk_screen_get_toplevel_windows", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenMakeDisplayName <-
function(object)
{
  checkPtrType(object, "GdkScreen")

  w <- .RGtkCall("S_gdk_screen_make_display_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenGetNMonitors <-
function(object)
{
  checkPtrType(object, "GdkScreen")

  w <- .RGtkCall("S_gdk_screen_get_n_monitors", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenGetMonitorGeometry <-
function(object, monitor.num)
{
  checkPtrType(object, "GdkScreen")
  monitor.num <- as.integer(monitor.num)

  w <- .RGtkCall("S_gdk_screen_get_monitor_geometry", object, monitor.num, PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenGetMonitorAtPoint <-
function(object, x, y)
{
  checkPtrType(object, "GdkScreen")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_gdk_screen_get_monitor_at_point", object, x, y, PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenGetMonitorAtWindow <-
function(object, window)
{
  checkPtrType(object, "GdkScreen")
  checkPtrType(window, "GdkWindow")

  w <- .RGtkCall("S_gdk_screen_get_monitor_at_window", object, window, PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenBroadcastClientMessage <-
function(object, event)
{
  checkPtrType(object, "GdkScreen")
  checkPtrType(event, "GdkEvent")

  w <- .RGtkCall("S_gdk_screen_broadcast_client_message", object, event, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkScreenGetDefault <-
function()
{
  

  w <- .RGtkCall("S_gdk_screen_get_default", PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenGetSetting <-
function(object, name)
{
  checkPtrType(object, "GdkScreen")
  name <- as.character(name)

  w <- .RGtkCall("S_gdk_screen_get_setting", object, name, PACKAGE = "RGtk2")

  return(w)
} 


gdkSpawnCommandLineOnScreen <-
function(screen, command.line, .errwarn = TRUE)
{
  checkPtrType(screen, "GdkScreen")
  command.line <- as.character(command.line)

  w <- .RGtkCall("S_gdk_spawn_command_line_on_screen", screen, command.line, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkAlternativeDialogButtonOrder <-
function(object)
{
  checkPtrType(object, "GdkScreen")

  w <- .RGtkCall("S_gtk_alternative_dialog_button_order", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkSelectionOwnerGet <-
function(selection)
{
  selection <- as.GdkAtom(selection)

  w <- .RGtkCall("S_gdk_selection_owner_get", selection, PACKAGE = "RGtk2")

  return(w)
} 


gdkSelectionOwnerGetForDisplay <-
function(display, selection)
{
  checkPtrType(display, "GdkDisplay")
  selection <- as.GdkAtom(selection)

  w <- .RGtkCall("S_gdk_selection_owner_get_for_display", display, selection, PACKAGE = "RGtk2")

  return(w)
} 


gdkSelectionPropertyGet <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_selection_property_get", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkVisualGetBestDepth <-
function()
{
  

  w <- .RGtkCall("S_gdk_visual_get_best_depth", PACKAGE = "RGtk2")

  return(w)
} 


gdkVisualGetBestType <-
function()
{
  

  w <- .RGtkCall("S_gdk_visual_get_best_type", PACKAGE = "RGtk2")

  return(w)
} 


gdkVisualGetSystem <-
function()
{
  

  w <- .RGtkCall("S_gdk_visual_get_system", PACKAGE = "RGtk2")

  return(w)
} 


gdkVisualGetBest <-
function()
{
  

  w <- .RGtkCall("S_gdk_visual_get_best", PACKAGE = "RGtk2")

  return(w)
} 


gdkVisualGetBestWithDepth <-
function(depth)
{
  depth <- as.integer(depth)

  w <- .RGtkCall("S_gdk_visual_get_best_with_depth", depth, PACKAGE = "RGtk2")

  return(w)
} 


gdkVisualGetBestWithType <-
function(visual.type)
{
  

  w <- .RGtkCall("S_gdk_visual_get_best_with_type", visual.type, PACKAGE = "RGtk2")

  return(w)
} 


gdkVisualGetBestWithBoth <-
function(depth, visual.type)
{
  depth <- as.integer(depth)
  

  w <- .RGtkCall("S_gdk_visual_get_best_with_both", depth, visual.type, PACKAGE = "RGtk2")

  return(w)
} 


gdkQueryDepths <-
function()
{
  

  w <- .RGtkCall("S_gdk_query_depths", PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkQueryVisualTypes <-
function()
{
  

  w <- .RGtkCall("S_gdk_query_visual_types", PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkListVisuals <-
function()
{
  

  w <- .RGtkCall("S_gdk_list_visuals", PACKAGE = "RGtk2")

  return(w)
} 


gdkVisualGetScreen <-
function(object)
{
  checkPtrType(object, "GdkVisual")

  w <- .RGtkCall("S_gdk_visual_get_screen", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkWindowObjectGetType <-
function()
{
  

  w <- .RGtkCall("S_gdk_window_object_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gdkWindowSetKeepAbove <-
function(object, setting)
{
  checkPtrType(object, "GdkWindow")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gdk_window_set_keep_above", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowSetKeepBelow <-
function(object, setting)
{
  checkPtrType(object, "GdkWindow")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gdk_window_set_keep_below", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowDestroy <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_destroy", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowGetWindowType <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_get_window_type", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkWindowAtPointer <-
function()
{
  

  w <- .RGtkCall("S_gdk_window_at_pointer", PACKAGE = "RGtk2")

  return(w)
} 


gdkWindowShow <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_show", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowShowUnraised <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_show_unraised", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowHide <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_hide", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowWithdraw <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_withdraw", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowMove <-
function(object, x, y)
{
  checkPtrType(object, "GdkWindow")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_gdk_window_move", object, x, y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowResize <-
function(object, width, height)
{
  checkPtrType(object, "GdkWindow")
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gdk_window_resize", object, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowMoveResize <-
function(object, x, y, width, height)
{
  checkPtrType(object, "GdkWindow")
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gdk_window_move_resize", object, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowMoveRegion <-
function(object, region, x, y)
{
  checkPtrType(object, "GdkWindow")
  checkPtrType(region, "GdkRegion")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_gdk_window_move_region", object, region, x, y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowReparent <-
function(object, new.parent, x, y)
{
  checkPtrType(object, "GdkWindow")
  checkPtrType(new.parent, "GdkWindow")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_gdk_window_reparent", object, new.parent, x, y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowClear <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_clear", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowClearArea <-
function(object, x, y, width, height)
{
  checkPtrType(object, "GdkWindow")
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gdk_window_clear_area", object, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowClearAreaE <-
function(object, x, y, width, height)
{
  checkPtrType(object, "GdkWindow")
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gdk_window_clear_area_e", object, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowRaise <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_raise", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowLower <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_lower", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowFocus <-
function(object, timestamp = "GDK_CURRENT_TIME")
{
  checkPtrType(object, "GdkWindow")
  timestamp <- as.numeric(timestamp)

  w <- .RGtkCall("S_gdk_window_focus", object, timestamp, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowSetUserData <-
function(object, user.data = NULL)
{
  checkPtrType(object, "GdkWindow")
  if (!is.null( user.data )) checkPtrType(user.data, "GtkWidget")

  w <- .RGtkCall("S_gdk_window_set_user_data", object, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowGetUserData <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_get_user_data", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkWindowSetOverrideRedirect <-
function(object, override.redirect)
{
  checkPtrType(object, "GdkWindow")
  override.redirect <- as.logical(override.redirect)

  w <- .RGtkCall("S_gdk_window_set_override_redirect", object, override.redirect, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowAddFilter <-
function(object, fun, data)
{
  checkPtrType(object, "GdkWindow")
  fun <- as.function(fun)
  

  w <- .RGtkCall("S_gdk_window_add_filter", object, fun, data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowRemoveFilter <-
function(object, fun, data)
{
  checkPtrType(object, "GdkWindow")
  fun <- as.function(fun)
  

  w <- .RGtkCall("S_gdk_window_remove_filter", object, fun, data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowScroll <-
function(object, dx, dy)
{
  checkPtrType(object, "GdkWindow")
  dx <- as.integer(dx)
  dy <- as.integer(dy)

  w <- .RGtkCall("S_gdk_window_scroll", object, dx, dy, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowShapeCombineMask <-
function(object, shape.mask = NULL, offset.x, offset.y)
{
  checkPtrType(object, "GdkWindow")
  if (!is.null( shape.mask )) checkPtrType(shape.mask, "GdkBitmap")
  offset.x <- as.integer(offset.x)
  offset.y <- as.integer(offset.y)

  w <- .RGtkCall("S_gdk_window_shape_combine_mask", object, shape.mask, offset.x, offset.y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowShapeCombineRegion <-
function(object, shape.region = NULL, offset.x, offset.y)
{
  checkPtrType(object, "GdkWindow")
  if (!is.null( shape.region )) checkPtrType(shape.region, "GdkRegion")
  offset.x <- as.integer(offset.x)
  offset.y <- as.integer(offset.y)

  w <- .RGtkCall("S_gdk_window_shape_combine_region", object, shape.region, offset.x, offset.y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowSetChildShapes <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_set_child_shapes", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowMergeChildShapes <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_merge_child_shapes", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowIsVisible <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_is_visible", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkWindowIsViewable <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_is_viewable", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkWindowGetState <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_get_state", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkWindowSetStaticGravities <-
function(object, use.static)
{
  checkPtrType(object, "GdkWindow")
  use.static <- as.logical(use.static)

  w <- .RGtkCall("S_gdk_window_set_static_gravities", object, use.static, PACKAGE = "RGtk2")

  return(w)
} 


gdkWindowSetHints <-
function(object, x, y, min.width, min.height, max.width, max.height, flags)
{
  checkPtrType(object, "GdkWindow")
  x <- as.integer(x)
  y <- as.integer(y)
  min.width <- as.integer(min.width)
  min.height <- as.integer(min.height)
  max.width <- as.integer(max.width)
  max.height <- as.integer(max.height)
  flags <- as.integer(flags)

  w <- .RGtkCall("S_gdk_window_set_hints", object, x, y, min.width, min.height, max.width, max.height, flags, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowSetTypeHint <-
function(object, hint)
{
  checkPtrType(object, "GdkWindow")
  

  w <- .RGtkCall("S_gdk_window_set_type_hint", object, hint, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowSetModalHint <-
function(object, modal)
{
  checkPtrType(object, "GdkWindow")
  modal <- as.logical(modal)

  w <- .RGtkCall("S_gdk_window_set_modal_hint", object, modal, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowSetSkipTaskbarHint <-
function(object, modal)
{
  checkPtrType(object, "GdkWindow")
  modal <- as.logical(modal)

  w <- .RGtkCall("S_gdk_window_set_skip_taskbar_hint", object, modal, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowSetSkipPagerHint <-
function(object, modal)
{
  checkPtrType(object, "GdkWindow")
  modal <- as.logical(modal)

  w <- .RGtkCall("S_gdk_window_set_skip_pager_hint", object, modal, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowSetUrgencyHint <-
function(object, urgent)
{
  checkPtrType(object, "GdkWindow")
  urgent <- as.logical(urgent)

  w <- .RGtkCall("S_gdk_window_set_urgency_hint", object, urgent, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowSetGeometryHints <-
function(object, geometry)
{
  checkPtrType(object, "GdkWindow")
  geometry <- as.GdkGeometry(geometry)

  w <- .RGtkCall("S_gdk_window_set_geometry_hints", object, geometry, PACKAGE = "RGtk2")

  return(w)
} 


gdkWindowBeginPaintRect <-
function(object, rectangle)
{
  checkPtrType(object, "GdkWindow")
  rectangle <- as.GdkRectangle(rectangle)

  w <- .RGtkCall("S_gdk_window_begin_paint_rect", object, rectangle, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowBeginPaintRegion <-
function(object, region)
{
  checkPtrType(object, "GdkWindow")
  checkPtrType(region, "GdkRegion")

  w <- .RGtkCall("S_gdk_window_begin_paint_region", object, region, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowEndPaint <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_end_paint", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowSetTitle <-
function(object, title)
{
  checkPtrType(object, "GdkWindow")
  title <- as.character(title)

  w <- .RGtkCall("S_gdk_window_set_title", object, title, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowSetRole <-
function(object, role)
{
  checkPtrType(object, "GdkWindow")
  role <- as.character(role)

  w <- .RGtkCall("S_gdk_window_set_role", object, role, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowSetTransientFor <-
function(object, leader)
{
  checkPtrType(object, "GdkWindow")
  checkPtrType(leader, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_set_transient_for", object, leader, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowSetBackground <-
function(object, color)
{
  checkPtrType(object, "GdkWindow")
  color <- as.GdkColor(color)

  w <- .RGtkCall("S_gdk_window_set_background", object, color, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowSetBackPixmap <-
function(object, pixmap = NULL, parent.relative)
{
  checkPtrType(object, "GdkWindow")
  if (!is.null( pixmap )) checkPtrType(pixmap, "GdkPixmap")
  parent.relative <- as.logical(parent.relative)

  w <- .RGtkCall("S_gdk_window_set_back_pixmap", object, pixmap, parent.relative, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowSetCursor <-
function(object, cursor = NULL)
{
  checkPtrType(object, "GdkWindow")
  if (!is.null( cursor )) checkPtrType(cursor, "GdkCursor")

  w <- .RGtkCall("S_gdk_window_set_cursor", object, cursor, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowGetGeometry <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_get_geometry", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowGetPosition <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_get_position", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowGetOrigin <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_get_origin", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkWindowGetDeskrelativeOrigin <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_get_deskrelative_origin", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkWindowGetRootOrigin <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_get_root_origin", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowGetFrameExtents <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_get_frame_extents", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkWindowGetPointer <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_get_pointer", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkWindowGetParent <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_get_parent", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkWindowGetToplevel <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_get_toplevel", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkWindowGetChildren <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_get_children", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkWindowPeekChildren <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_peek_children", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkWindowGetEvents <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_get_events", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkWindowSetEvents <-
function(object, event.mask)
{
  checkPtrType(object, "GdkWindow")
  

  w <- .RGtkCall("S_gdk_window_set_events", object, event.mask, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowSetIconList <-
function(object, pixbufs)
{
  checkPtrType(object, "GdkWindow")
  pixbufs <- lapply(pixbufs, function(x) { x <- as.GList(x); x })

  w <- .RGtkCall("S_gdk_window_set_icon_list", object, pixbufs, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowSetIcon <-
function(object, icon.window, pixmap, mask)
{
  checkPtrType(object, "GdkWindow")
  checkPtrType(icon.window, "GdkWindow")
  checkPtrType(pixmap, "GdkPixmap")
  checkPtrType(mask, "GdkBitmap")

  w <- .RGtkCall("S_gdk_window_set_icon", object, icon.window, pixmap, mask, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowSetIconName <-
function(object, name)
{
  checkPtrType(object, "GdkWindow")
  name <- as.character(name)

  w <- .RGtkCall("S_gdk_window_set_icon_name", object, name, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowSetGroup <-
function(object, leader)
{
  checkPtrType(object, "GdkWindow")
  checkPtrType(leader, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_set_group", object, leader, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowGetGroup <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_get_group", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkWindowSetDecorations <-
function(object, decorations)
{
  checkPtrType(object, "GdkWindow")
  

  w <- .RGtkCall("S_gdk_window_set_decorations", object, decorations, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowGetDecorations <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_get_decorations", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkWindowSetFunctions <-
function(object, functions)
{
  checkPtrType(object, "GdkWindow")
  

  w <- .RGtkCall("S_gdk_window_set_functions", object, functions, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowGetToplevels <-
function()
{
  

  w <- .RGtkCall("S_gdk_window_get_toplevels", PACKAGE = "RGtk2")

  return(w)
} 


gdkWindowIconify <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_iconify", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowDeiconify <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_deiconify", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowStick <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_stick", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowUnstick <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_unstick", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowMaximize <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_maximize", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowUnmaximize <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_unmaximize", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowFullscreen <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_fullscreen", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowUnfullscreen <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_unfullscreen", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowRegisterDnd <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_register_dnd", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowBeginResizeDrag <-
function(object, edge, button, root.x, root.y, timestamp)
{
  checkPtrType(object, "GdkWindow")
  
  button <- as.integer(button)
  root.x <- as.integer(root.x)
  root.y <- as.integer(root.y)
  timestamp <- as.numeric(timestamp)

  w <- .RGtkCall("S_gdk_window_begin_resize_drag", object, edge, button, root.x, root.y, timestamp, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowBeginMoveDrag <-
function(object, button, root.x, root.y, timestamp)
{
  checkPtrType(object, "GdkWindow")
  button <- as.integer(button)
  root.x <- as.integer(root.x)
  root.y <- as.integer(root.y)
  timestamp <- as.numeric(timestamp)

  w <- .RGtkCall("S_gdk_window_begin_move_drag", object, button, root.x, root.y, timestamp, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowInvalidateRect <-
function(object, rect = NULL, invalidate.children)
{
  checkPtrType(object, "GdkWindow")
  if (!is.null( rect )) rect <- as.GdkRectangle(rect)
  invalidate.children <- as.logical(invalidate.children)

  w <- .RGtkCall("S_gdk_window_invalidate_rect", object, rect, invalidate.children, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowInvalidateRegion <-
function(object, region, invalidate.children)
{
  checkPtrType(object, "GdkWindow")
  checkPtrType(region, "GdkRegion")
  invalidate.children <- as.logical(invalidate.children)

  w <- .RGtkCall("S_gdk_window_invalidate_region", object, region, invalidate.children, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowGetUpdateArea <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_get_update_area", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkWindowFreezeUpdates <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_freeze_updates", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowThawUpdates <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_thaw_updates", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowProcessAllUpdates <-
function()
{
  

  w <- .RGtkCall("S_gdk_window_process_all_updates", PACKAGE = "RGtk2")

  return(w)
} 


gdkWindowProcessUpdates <-
function(object, update.children)
{
  checkPtrType(object, "GdkWindow")
  update.children <- as.logical(update.children)

  w <- .RGtkCall("S_gdk_window_process_updates", object, update.children, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowSetDebugUpdates <-
function(setting)
{
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gdk_window_set_debug_updates", setting, PACKAGE = "RGtk2")

  return(w)
} 


gdkWindowGetInternalPaintInfo <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_get_internal_paint_info", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkGetDefaultRootWindow <-
function()
{
  

  w <- .RGtkCall("S_gdk_get_default_root_window", PACKAGE = "RGtk2")

  return(w)
} 


gdkWindowSetAcceptFocus <-
function(object, accept.focus)
{
  checkPtrType(object, "GdkWindow")
  accept.focus <- as.logical(accept.focus)

  w <- .RGtkCall("S_gdk_window_set_accept_focus", object, accept.focus, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowSetFocusOnMap <-
function(object, focus.on.map)
{
  checkPtrType(object, "GdkWindow")
  focus.on.map <- as.logical(focus.on.map)

  w <- .RGtkCall("S_gdk_window_set_focus_on_map", object, focus.on.map, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowEnableSynchronizedConfigure <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_enable_synchronized_configure", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowConfigureFinished <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_configure_finished", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDragFinish <-
function(object, success, del, time = "GDK_CURRENT_TIME")
{
  checkPtrType(object, "GdkDragContext")
  success <- as.logical(success)
  del <- as.logical(del)
  time <- as.numeric(time)

  w <- .RGtkCall("S_gtk_drag_finish", object, success, del, time, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDragSetIconName <-
function(object, icon.name, hot.x, hot.y)
{
  checkPtrType(object, "GdkDragContext")
  icon.name <- as.character(icon.name)
  hot.x <- as.integer(hot.x)
  hot.y <- as.integer(hot.y)

  w <- .RGtkCall("S_gtk_drag_set_icon_name", object, icon.name, hot.x, hot.y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDragSetIconWidget <-
function(object, widget, hot.x, hot.y)
{
  checkPtrType(object, "GdkDragContext")
  checkPtrType(widget, "GtkWidget")
  hot.x <- as.integer(hot.x)
  hot.y <- as.integer(hot.y)

  w <- .RGtkCall("S_gtk_drag_set_icon_widget", object, widget, hot.x, hot.y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDragSetIconPixmap <-
function(object, colormap, pixmap, mask, hot.x, hot.y)
{
  checkPtrType(object, "GdkDragContext")
  checkPtrType(colormap, "GdkColormap")
  checkPtrType(pixmap, "GdkPixmap")
  checkPtrType(mask, "GdkBitmap")
  hot.x <- as.integer(hot.x)
  hot.y <- as.integer(hot.y)

  w <- .RGtkCall("S_gtk_drag_set_icon_pixmap", object, colormap, pixmap, mask, hot.x, hot.y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDragSetIconPixbuf <-
function(object, pixbuf, hot.x, hot.y)
{
  checkPtrType(object, "GdkDragContext")
  checkPtrType(pixbuf, "GdkPixbuf")
  hot.x <- as.integer(hot.x)
  hot.y <- as.integer(hot.y)

  w <- .RGtkCall("S_gtk_drag_set_icon_pixbuf", object, pixbuf, hot.x, hot.y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDragSetIconStock <-
function(object, stock.id, hot.x, hot.y)
{
  checkPtrType(object, "GdkDragContext")
  stock.id <- as.character(stock.id)
  hot.x <- as.integer(hot.x)
  hot.y <- as.integer(hot.y)

  w <- .RGtkCall("S_gtk_drag_set_icon_stock", object, stock.id, hot.x, hot.y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDragSetIconDefault <-
function(object)
{
  checkPtrType(object, "GdkDragContext")

  w <- .RGtkCall("S_gtk_drag_set_icon_default", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkPixbufGetColorspace <-
function(object)
{
  checkPtrType(object, "GdkPixbuf")

  w <- .RGtkCall("S_gdk_pixbuf_get_colorspace", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufGetNChannels <-
function(object)
{
  checkPtrType(object, "GdkPixbuf")

  w <- .RGtkCall("S_gdk_pixbuf_get_n_channels", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufGetHasAlpha <-
function(object)
{
  checkPtrType(object, "GdkPixbuf")

  w <- .RGtkCall("S_gdk_pixbuf_get_has_alpha", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufGetBitsPerSample <-
function(object)
{
  checkPtrType(object, "GdkPixbuf")

  w <- .RGtkCall("S_gdk_pixbuf_get_bits_per_sample", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufGetPixels <-
function(object)
{
  checkPtrType(object, "GdkPixbuf")

  w <- .RGtkCall("S_gdk_pixbuf_get_pixels", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufGetWidth <-
function(object)
{
  checkPtrType(object, "GdkPixbuf")

  w <- .RGtkCall("S_gdk_pixbuf_get_width", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufGetHeight <-
function(object)
{
  checkPtrType(object, "GdkPixbuf")

  w <- .RGtkCall("S_gdk_pixbuf_get_height", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufGetRowstride <-
function(object)
{
  checkPtrType(object, "GdkPixbuf")

  w <- .RGtkCall("S_gdk_pixbuf_get_rowstride", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufNew <-
function(colorspace, has.alpha, bits.per.sample, width, height)
{
  
  has.alpha <- as.logical(has.alpha)
  bits.per.sample <- as.integer(bits.per.sample)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gdk_pixbuf_new", colorspace, has.alpha, bits.per.sample, width, height, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufCopy <-
function(object)
{
  checkPtrType(object, "GdkPixbuf")

  w <- .RGtkCall("S_gdk_pixbuf_copy", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufNewFromFile <-
function(filename, .errwarn = TRUE)
{
  filename <- as.character(filename)

  w <- .RGtkCall("S_gdk_pixbuf_new_from_file", filename, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gdkPixbufNewFromFileAtSize <-
function(filename, width, height, .errwarn = TRUE)
{
  filename <- as.character(filename)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gdk_pixbuf_new_from_file_at_size", filename, width, height, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gdkPixbufNewFromFileAtScale <-
function(filename, width, height, preserve.aspect.ratio, .errwarn = TRUE)
{
  filename <- as.character(filename)
  width <- as.integer(width)
  height <- as.integer(height)
  preserve.aspect.ratio <- as.logical(preserve.aspect.ratio)

  w <- .RGtkCall("S_gdk_pixbuf_new_from_file_at_scale", filename, width, height, preserve.aspect.ratio, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gdkPixbufNewFromXpmData <-
function(data)
{
  data <- as.list(as.character(data))

  w <- .RGtkCall("S_gdk_pixbuf_new_from_xpm_data", data, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufNewSubpixbuf <-
function(object, src.x, src.y, width, height)
{
  checkPtrType(object, "GdkPixbuf")
  src.x <- as.integer(src.x)
  src.y <- as.integer(src.y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gdk_pixbuf_new_subpixbuf", object, src.x, src.y, width, height, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufFill <-
function(object, pixel)
{
  checkPtrType(object, "GdkPixbuf")
  pixel <- as.numeric(pixel)

  w <- .RGtkCall("S_gdk_pixbuf_fill", object, pixel, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkPixbufSavev <-
function(object, filename, type, option.keys, option.values, .errwarn = TRUE)
{
  checkPtrType(object, "GdkPixbuf")
  filename <- as.character(filename)
  type <- as.character(type)
  option.keys <- as.list(as.character(option.keys))
  option.values <- as.list(as.character(option.values))

  w <- .RGtkCall("S_gdk_pixbuf_savev", object, filename, type, option.keys, option.values, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gdkPixbufSaveToCallbackv <-
function(object, save.func, user.data, type, option.keys, option.values, .errwarn = TRUE)
{
  checkPtrType(object, "GdkPixbuf")
  save.func <- as.function(save.func)
  
  type <- as.character(type)
  option.keys <- as.list(as.character(option.keys))
  option.values <- as.list(as.character(option.values))

  w <- .RGtkCall("S_gdk_pixbuf_save_to_callbackv", object, save.func, user.data, type, option.keys, option.values, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gdkPixbufSaveToBufferv <-
function(object, type, option.keys, option.values, .errwarn = TRUE)
{
  checkPtrType(object, "GdkPixbuf")
  type <- as.character(type)
  option.keys <- as.list(as.character(option.keys))
  option.values <- as.list(as.character(option.values))

  w <- .RGtkCall("S_gdk_pixbuf_save_to_bufferv", object, type, option.keys, option.values, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(invisible(w))
} 


gdkPixbufAddAlpha <-
function(object, substitute.color, r, g, b)
{
  checkPtrType(object, "GdkPixbuf")
  substitute.color <- as.logical(substitute.color)
  r <- as.raw(r)
  g <- as.raw(g)
  b <- as.raw(b)

  w <- .RGtkCall("S_gdk_pixbuf_add_alpha", object, substitute.color, r, g, b, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufCopyArea <-
function(object, src.x, src.y, width, height, dest.pixbuf, dest.x, dest.y)
{
  checkPtrType(object, "GdkPixbuf")
  src.x <- as.integer(src.x)
  src.y <- as.integer(src.y)
  width <- as.integer(width)
  height <- as.integer(height)
  checkPtrType(dest.pixbuf, "GdkPixbuf")
  dest.x <- as.integer(dest.x)
  dest.y <- as.integer(dest.y)

  w <- .RGtkCall("S_gdk_pixbuf_copy_area", object, src.x, src.y, width, height, dest.pixbuf, dest.x, dest.y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkPixbufSaturateAndPixelate <-
function(object, dest, saturation, pixelate)
{
  checkPtrType(object, "GdkPixbuf")
  checkPtrType(dest, "GdkPixbuf")
  saturation <- as.numeric(saturation)
  pixelate <- as.logical(pixelate)

  w <- .RGtkCall("S_gdk_pixbuf_saturate_and_pixelate", object, dest, saturation, pixelate, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkPixbufScale <-
function(object, dest, dest.x, dest.y, dest.width, dest.height, offset.x, offset.y, scale.x, scale.y, interp.type)
{
  checkPtrType(object, "GdkPixbuf")
  checkPtrType(dest, "GdkPixbuf")
  dest.x <- as.integer(dest.x)
  dest.y <- as.integer(dest.y)
  dest.width <- as.integer(dest.width)
  dest.height <- as.integer(dest.height)
  offset.x <- as.numeric(offset.x)
  offset.y <- as.numeric(offset.y)
  scale.x <- as.numeric(scale.x)
  scale.y <- as.numeric(scale.y)
  

  w <- .RGtkCall("S_gdk_pixbuf_scale", object, dest, dest.x, dest.y, dest.width, dest.height, offset.x, offset.y, scale.x, scale.y, interp.type, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkPixbufComposite <-
function(object, dest, dest.x, dest.y, dest.width, dest.height, offset.x, offset.y, scale.x, scale.y, interp.type, overall.alpha)
{
  checkPtrType(object, "GdkPixbuf")
  checkPtrType(dest, "GdkPixbuf")
  dest.x <- as.integer(dest.x)
  dest.y <- as.integer(dest.y)
  dest.width <- as.integer(dest.width)
  dest.height <- as.integer(dest.height)
  offset.x <- as.numeric(offset.x)
  offset.y <- as.numeric(offset.y)
  scale.x <- as.numeric(scale.x)
  scale.y <- as.numeric(scale.y)
  
  overall.alpha <- as.integer(overall.alpha)

  w <- .RGtkCall("S_gdk_pixbuf_composite", object, dest, dest.x, dest.y, dest.width, dest.height, offset.x, offset.y, scale.x, scale.y, interp.type, overall.alpha, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkPixbufCompositeColor <-
function(object, dest, dest.x, dest.y, dest.width, dest.height, offset.x, offset.y, scale.x, scale.y, interp.type, overall.alpha, check.x, check.y, check.size, color1, color2)
{
  checkPtrType(object, "GdkPixbuf")
  checkPtrType(dest, "GdkPixbuf")
  dest.x <- as.integer(dest.x)
  dest.y <- as.integer(dest.y)
  dest.width <- as.integer(dest.width)
  dest.height <- as.integer(dest.height)
  offset.x <- as.numeric(offset.x)
  offset.y <- as.numeric(offset.y)
  scale.x <- as.numeric(scale.x)
  scale.y <- as.numeric(scale.y)
  
  overall.alpha <- as.integer(overall.alpha)
  check.x <- as.integer(check.x)
  check.y <- as.integer(check.y)
  check.size <- as.integer(check.size)
  color1 <- as.numeric(color1)
  color2 <- as.numeric(color2)

  w <- .RGtkCall("S_gdk_pixbuf_composite_color", object, dest, dest.x, dest.y, dest.width, dest.height, offset.x, offset.y, scale.x, scale.y, interp.type, overall.alpha, check.x, check.y, check.size, color1, color2, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkPixbufRotateSimple <-
function(object, angle)
{
  checkPtrType(object, "GdkPixbuf")
  

  w <- .RGtkCall("S_gdk_pixbuf_rotate_simple", object, angle, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufFlip <-
function(object, horizontal)
{
  checkPtrType(object, "GdkPixbuf")
  horizontal <- as.logical(horizontal)

  w <- .RGtkCall("S_gdk_pixbuf_flip", object, horizontal, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufScaleSimple <-
function(object, dest.width, dest.height, interp.type)
{
  checkPtrType(object, "GdkPixbuf")
  dest.width <- as.integer(dest.width)
  dest.height <- as.integer(dest.height)
  

  w <- .RGtkCall("S_gdk_pixbuf_scale_simple", object, dest.width, dest.height, interp.type, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufCompositeColorSimple <-
function(object, dest.width, dest.height, interp.type, overall.alpha, check.size, color1, color2)
{
  checkPtrType(object, "GdkPixbuf")
  dest.width <- as.integer(dest.width)
  dest.height <- as.integer(dest.height)
  
  overall.alpha <- as.integer(overall.alpha)
  check.size <- as.integer(check.size)
  color1 <- as.numeric(color1)
  color2 <- as.numeric(color2)

  w <- .RGtkCall("S_gdk_pixbuf_composite_color_simple", object, dest.width, dest.height, interp.type, overall.alpha, check.size, color1, color2, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufAnimationGetType <-
function()
{
  

  w <- .RGtkCall("S_gdk_pixbuf_animation_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufAnimationNewFromFile <-
function(filename, .errwarn = TRUE)
{
  filename <- as.character(filename)

  w <- .RGtkCall("S_gdk_pixbuf_animation_new_from_file", filename, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gdkPixbufAnimationGetWidth <-
function(object)
{
  checkPtrType(object, "GdkPixbufAnimation")

  w <- .RGtkCall("S_gdk_pixbuf_animation_get_width", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufAnimationGetHeight <-
function(object)
{
  checkPtrType(object, "GdkPixbufAnimation")

  w <- .RGtkCall("S_gdk_pixbuf_animation_get_height", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufAnimationIsStaticImage <-
function(object)
{
  checkPtrType(object, "GdkPixbufAnimation")

  w <- .RGtkCall("S_gdk_pixbuf_animation_is_static_image", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufAnimationGetStaticImage <-
function(object)
{
  checkPtrType(object, "GdkPixbufAnimation")

  w <- .RGtkCall("S_gdk_pixbuf_animation_get_static_image", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufAnimationGetIter <-
function(object, start.time)
{
  checkPtrType(object, "GdkPixbufAnimation")
  start.time <- as.GTimeVal(start.time)

  w <- .RGtkCall("S_gdk_pixbuf_animation_get_iter", object, start.time, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufAnimationIterGetType <-
function()
{
  

  w <- .RGtkCall("S_gdk_pixbuf_animation_iter_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufAnimationIterGetDelayTime <-
function(object)
{
  checkPtrType(object, "GdkPixbufAnimationIter")

  w <- .RGtkCall("S_gdk_pixbuf_animation_iter_get_delay_time", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufAnimationIterGetPixbuf <-
function(object)
{
  checkPtrType(object, "GdkPixbufAnimationIter")

  w <- .RGtkCall("S_gdk_pixbuf_animation_iter_get_pixbuf", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufAnimationIterOnCurrentlyLoadingFrame <-
function(object)
{
  checkPtrType(object, "GdkPixbufAnimationIter")

  w <- .RGtkCall("S_gdk_pixbuf_animation_iter_on_currently_loading_frame", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufAnimationIterAdvance <-
function(object, current.time)
{
  checkPtrType(object, "GdkPixbufAnimationIter")
  current.time <- as.GTimeVal(current.time)

  w <- .RGtkCall("S_gdk_pixbuf_animation_iter_advance", object, current.time, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufSimpleAnimNew <-
function(width, height, rate)
{
  width <- as.integer(width)
  height <- as.integer(height)
  rate <- as.numeric(rate)

  w <- .RGtkCall("S_gdk_pixbuf_simple_anim_new", width, height, rate, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufSimpleAnimAddFrame <-
function(object, pixbuf)
{
  checkPtrType(object, "GdkPixbufSimpleAnim")
  checkPtrType(pixbuf, "GdkPixbuf")

  w <- .RGtkCall("S_gdk_pixbuf_simple_anim_add_frame", object, pixbuf, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkPixbufGetOption <-
function(object, key)
{
  checkPtrType(object, "GdkPixbuf")
  key <- as.character(key)

  w <- .RGtkCall("S_gdk_pixbuf_get_option", object, key, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufSetOption <-
function(object, key, value)
{
  checkPtrType(object, "GdkPixbuf")
  key <- as.character(key)
  value <- as.character(value)

  w <- .RGtkCall("S_gdk_pixbuf_set_option", object, key, value, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufGetFormats <-
function()
{
  

  w <- .RGtkCall("S_gdk_pixbuf_get_formats", PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufGetFileInfo <-
function(filename)
{
  filename <- as.character(filename)

  w <- .RGtkCall("S_gdk_pixbuf_get_file_info", filename, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufFormatGetName <-
function(object)
{
  checkPtrType(object, "GdkPixbufFormat")

  w <- .RGtkCall("S_gdk_pixbuf_format_get_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufFormatIsScalable <-
function(object)
{
  checkPtrType(object, "GdkPixbufFormat")

  w <- .RGtkCall("S_gdk_pixbuf_format_is_scalable", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufFormatIsDisabled <-
function(object)
{
  checkPtrType(object, "GdkPixbufFormat")

  w <- .RGtkCall("S_gdk_pixbuf_format_is_disabled", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufFormatSetDisabled <-
function(object, disabled)
{
  checkPtrType(object, "GdkPixbufFormat")
  disabled <- as.logical(disabled)

  w <- .RGtkCall("S_gdk_pixbuf_format_set_disabled", object, disabled, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkPixbufFormatGetLicense <-
function(object)
{
  checkPtrType(object, "GdkPixbufFormat")

  w <- .RGtkCall("S_gdk_pixbuf_format_get_license", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufFormatGetDescription <-
function(object)
{
  checkPtrType(object, "GdkPixbufFormat")

  w <- .RGtkCall("S_gdk_pixbuf_format_get_description", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufFormatGetMimeTypes <-
function(object)
{
  checkPtrType(object, "GdkPixbufFormat")

  w <- .RGtkCall("S_gdk_pixbuf_format_get_mime_types", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufFormatGetExtensions <-
function(object)
{
  checkPtrType(object, "GdkPixbufFormat")

  w <- .RGtkCall("S_gdk_pixbuf_format_get_extensions", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufFormatIsWritable <-
function(object)
{
  checkPtrType(object, "GdkPixbufFormat")

  w <- .RGtkCall("S_gdk_pixbuf_format_is_writable", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufLoaderGetType <-
function()
{
  

  w <- .RGtkCall("S_gdk_pixbuf_loader_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufLoaderNew <-
function()
{
  

  w <- .RGtkCall("S_gdk_pixbuf_loader_new", PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufLoaderNewWithType <-
function(image.type, .errwarn = TRUE)
{
  image.type <- as.character(image.type)

  w <- .RGtkCall("S_gdk_pixbuf_loader_new_with_type", image.type, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gdkPixbufLoaderNewWithMimeType <-
function(mime.type, .errwarn = TRUE)
{
  mime.type <- as.character(mime.type)

  w <- .RGtkCall("S_gdk_pixbuf_loader_new_with_mime_type", mime.type, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gdkPixbufLoaderWrite <-
function(object, buf, .errwarn = TRUE)
{
  checkPtrType(object, "GdkPixbufLoader")
  buf <- as.list(as.raw(buf))

  w <- .RGtkCall("S_gdk_pixbuf_loader_write", object, buf, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gdkPixbufLoaderGetPixbuf <-
function(object)
{
  checkPtrType(object, "GdkPixbufLoader")

  w <- .RGtkCall("S_gdk_pixbuf_loader_get_pixbuf", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufLoaderGetAnimation <-
function(object)
{
  checkPtrType(object, "GdkPixbufLoader")

  w <- .RGtkCall("S_gdk_pixbuf_loader_get_animation", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufLoaderClose <-
function(object, .errwarn = TRUE)
{
  checkPtrType(object, "GdkPixbufLoader")

  w <- .RGtkCall("S_gdk_pixbuf_loader_close", object, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gdkPixbufLoaderSetSize <-
function(object, width, height)
{
  checkPtrType(object, "GdkPixbufLoader")
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gdk_pixbuf_loader_set_size", object, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkPixbufLoaderGetFormat <-
function(object)
{
  checkPtrType(object, "GdkPixbufLoader")

  w <- .RGtkCall("S_gdk_pixbuf_loader_get_format", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkRectangleIntersect <-
function(src1, src2)
{
  src1 <- as.GdkRectangle(src1)
  src2 <- as.GdkRectangle(src2)

  w <- .RGtkCall("S_gdk_rectangle_intersect", src1, src2, PACKAGE = "RGtk2")

  return(w)
} 


gdkRectangleUnion <-
function(src1, src2)
{
  src1 <- as.GdkRectangle(src1)
  src2 <- as.GdkRectangle(src2)

  w <- .RGtkCall("S_gdk_rectangle_union", src1, src2, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkRegionNew <-
function()
{
  

  w <- .RGtkCall("S_gdk_region_new", PACKAGE = "RGtk2")

  return(w)
} 


gdkRegionPolygon <-
function(points, fill.rule)
{
  points <- lapply(points, function(x) { x <- as.GdkPoint(x); x })
  

  w <- .RGtkCall("S_gdk_region_polygon", points, fill.rule, PACKAGE = "RGtk2")

  return(w)
} 


gdkRegionCopy <-
function(object)
{
  checkPtrType(object, "GdkRegion")

  w <- .RGtkCall("S_gdk_region_copy", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkRegionRectangle <-
function(rectangle)
{
  rectangle <- as.GdkRectangle(rectangle)

  w <- .RGtkCall("S_gdk_region_rectangle", rectangle, PACKAGE = "RGtk2")

  return(w)
} 


gdkRegionGetClipbox <-
function(object)
{
  checkPtrType(object, "GdkRegion")

  w <- .RGtkCall("S_gdk_region_get_clipbox", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkRegionGetRectangles <-
function(object)
{
  checkPtrType(object, "GdkRegion")

  w <- .RGtkCall("S_gdk_region_get_rectangles", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkRegionEmpty <-
function(object)
{
  checkPtrType(object, "GdkRegion")

  w <- .RGtkCall("S_gdk_region_empty", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkRegionEqual <-
function(object, region2)
{
  checkPtrType(object, "GdkRegion")
  checkPtrType(region2, "GdkRegion")

  w <- .RGtkCall("S_gdk_region_equal", object, region2, PACKAGE = "RGtk2")

  return(w)
} 


gdkRegionPointIn <-
function(object, x, y)
{
  checkPtrType(object, "GdkRegion")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_gdk_region_point_in", object, x, y, PACKAGE = "RGtk2")

  return(w)
} 


gdkRegionRectIn <-
function(object, rect)
{
  checkPtrType(object, "GdkRegion")
  rect <- as.GdkRectangle(rect)

  w <- .RGtkCall("S_gdk_region_rect_in", object, rect, PACKAGE = "RGtk2")

  return(w)
} 


gdkRegionOffset <-
function(object, dx, dy)
{
  checkPtrType(object, "GdkRegion")
  dx <- as.integer(dx)
  dy <- as.integer(dy)

  w <- .RGtkCall("S_gdk_region_offset", object, dx, dy, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkRegionShrink <-
function(object, dx, dy)
{
  checkPtrType(object, "GdkRegion")
  dx <- as.integer(dx)
  dy <- as.integer(dy)

  w <- .RGtkCall("S_gdk_region_shrink", object, dx, dy, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkRegionUnionWithRect <-
function(object, rect)
{
  checkPtrType(object, "GdkRegion")
  rect <- as.GdkRectangle(rect)

  w <- .RGtkCall("S_gdk_region_union_with_rect", object, rect, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkRegionIntersect <-
function(object, source2)
{
  checkPtrType(object, "GdkRegion")
  checkPtrType(source2, "GdkRegion")

  w <- .RGtkCall("S_gdk_region_intersect", object, source2, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkRegionUnion <-
function(object, source2)
{
  checkPtrType(object, "GdkRegion")
  checkPtrType(source2, "GdkRegion")

  w <- .RGtkCall("S_gdk_region_union", object, source2, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkRegionSubtract <-
function(object, source2)
{
  checkPtrType(object, "GdkRegion")
  checkPtrType(source2, "GdkRegion")

  w <- .RGtkCall("S_gdk_region_subtract", object, source2, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkRegionXor <-
function(object, source2)
{
  checkPtrType(object, "GdkRegion")
  checkPtrType(source2, "GdkRegion")

  w <- .RGtkCall("S_gdk_region_xor", object, source2, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkRegionSpansIntersectForeach <-
function(object, spans, sorted, fun, data)
{
  checkPtrType(object, "GdkRegion")
  spans <- lapply(spans, function(x) { x <- as.GdkSpan(x); x })
  sorted <- as.logical(sorted)
  fun <- as.function(fun)
  

  w <- .RGtkCall("S_gdk_region_spans_intersect_foreach", object, spans, sorted, fun, data, PACKAGE = "RGtk2")

  return(w)
} 


gdkCairoSetSourcePixmap <-
function(cr, pixmap, pixmap.x, pixmap.y)
{
  checkPtrType(cr, "Cairo")
  checkPtrType(pixmap, "GdkPixmap")
  pixmap.x <- as.numeric(pixmap.x)
  pixmap.y <- as.numeric(pixmap.y)

  w <- .RGtkCall("S_gdk_cairo_set_source_pixmap", cr, pixmap, pixmap.x, pixmap.y, PACKAGE = "RGtk2")

  return(w)
} 


gdkDisplaySupportsShapes <-
function(object)
{
  checkPtrType(object, "GdkDisplay")

  w <- .RGtkCall("S_gdk_display_supports_shapes", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkDisplaySupportsInputShapes <-
function(object)
{
  checkPtrType(object, "GdkDisplay")

  w <- .RGtkCall("S_gdk_display_supports_input_shapes", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkAtomInternStaticString <-
function(atom.name)
{
  atom.name <- as.character(atom.name)

  w <- .RGtkCall("S_gdk_atom_intern_static_string", atom.name, PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenIsComposited <-
function(object)
{
  checkPtrType(object, "GdkScreen")

  w <- .RGtkCall("S_gdk_screen_is_composited", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenSetFontOptions <-
function(object, options)
{
  checkPtrType(object, "GdkScreen")
  checkPtrType(options, "CairoFontOptions")

  w <- .RGtkCall("S_gdk_screen_set_font_options", object, options, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkScreenGetFontOptions <-
function(object)
{
  checkPtrType(object, "GdkScreen")

  w <- .RGtkCall("S_gdk_screen_get_font_options", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenSetResolution <-
function(object, dpi)
{
  checkPtrType(object, "GdkScreen")
  dpi <- as.numeric(dpi)

  w <- .RGtkCall("S_gdk_screen_set_resolution", object, dpi, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkScreenGetResolution <-
function(object)
{
  checkPtrType(object, "GdkScreen")

  w <- .RGtkCall("S_gdk_screen_get_resolution", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenGetActiveWindow <-
function(object)
{
  checkPtrType(object, "GdkScreen")

  w <- .RGtkCall("S_gdk_screen_get_active_window", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenGetWindowStack <-
function(object)
{
  checkPtrType(object, "GdkScreen")

  w <- .RGtkCall("S_gdk_screen_get_window_stack", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkWindowInputShapeCombineMask <-
function(object, mask, x, y)
{
  checkPtrType(object, "GdkWindow")
  checkPtrType(mask, "GdkBitmap")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_gdk_window_input_shape_combine_mask", object, mask, x, y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowInputShapeCombineRegion <-
function(object, shape.region, offset.x, offset.y)
{
  checkPtrType(object, "GdkWindow")
  checkPtrType(shape.region, "GdkRegion")
  offset.x <- as.integer(offset.x)
  offset.y <- as.integer(offset.y)

  w <- .RGtkCall("S_gdk_window_input_shape_combine_region", object, shape.region, offset.x, offset.y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowSetChildInputShapes <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_set_child_input_shapes", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowMergeChildInputShapes <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_merge_child_input_shapes", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowGetTypeHint <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_get_type_hint", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkColorToString <-
function(object)
{
  object <- as.GdkColor(object)

  w <- .RGtkCall("S_gdk_color_to_string", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkDisplaySupportsComposite <-
function(object)
{
  checkPtrType(object, "GdkDisplay")

  w <- .RGtkCall("S_gdk_display_supports_composite", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkEventRequestMotions <-
function(event)
{
  checkPtrType(event, "GdkEventMotion")

  w <- .RGtkCall("S_gdk_event_request_motions", event, PACKAGE = "RGtk2")

  return(w)
} 


gdkKeymapHaveBidiLayouts <-
function(object)
{
  checkPtrType(object, "GdkKeymap")

  w <- .RGtkCall("S_gdk_keymap_have_bidi_layouts", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkPangoAttrEmbossColorNew <-
function(color)
{
  color <- as.GdkColor(color)

  w <- .RGtkCall("S_gdk_pango_attr_emboss_color_new", color, PACKAGE = "RGtk2")

  return(w)
} 


gdkWindowSetComposited <-
function(object, composited)
{
  checkPtrType(object, "GdkWindow")
  composited <- as.logical(composited)

  w <- .RGtkCall("S_gdk_window_set_composited", object, composited, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowSetStartupId <-
function(object, startup.id)
{
  checkPtrType(object, "GdkWindow")
  startup.id <- as.character(startup.id)

  w <- .RGtkCall("S_gdk_window_set_startup_id", object, startup.id, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowBeep <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_beep", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowSetOpacity <-
function(object, opacity)
{
  checkPtrType(object, "GdkWindow")
  opacity <- as.numeric(opacity)

  w <- .RGtkCall("S_gdk_window_set_opacity", object, opacity, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkNotifyStartupCompleteWithId <-
function(id)
{
  id <- as.character(id)

  w <- .RGtkCall("S_gdk_notify_startup_complete_with_id", id, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufApplyEmbeddedOrientation <-
function(object)
{
  checkPtrType(object, "GdkPixbuf")

  w <- .RGtkCall("S_gdk_pixbuf_apply_embedded_orientation", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkAppLaunchContextGetType <-
function()
{
  

  w <- .RGtkCall("S_gdk_app_launch_context_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gdkAppLaunchContextNew <-
function()
{
  

  w <- .RGtkCall("S_gdk_app_launch_context_new", PACKAGE = "RGtk2")

  return(w)
} 


gdkAppLaunchContextSetDisplay <-
function(object, display)
{
  checkPtrType(object, "GdkAppLaunchContext")
  checkPtrType(display, "GdkDisplay")

  w <- .RGtkCall("S_gdk_app_launch_context_set_display", object, display, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkAppLaunchContextSetScreen <-
function(object, screen)
{
  checkPtrType(object, "GdkAppLaunchContext")
  checkPtrType(screen, "GdkScreen")

  w <- .RGtkCall("S_gdk_app_launch_context_set_screen", object, screen, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkAppLaunchContextSetDesktop <-
function(object, desktop)
{
  checkPtrType(object, "GdkAppLaunchContext")
  desktop <- as.integer(desktop)

  w <- .RGtkCall("S_gdk_app_launch_context_set_desktop", object, desktop, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkAppLaunchContextSetTimestamp <-
function(object, timestamp)
{
  checkPtrType(object, "GdkAppLaunchContext")
  timestamp <- as.numeric(timestamp)

  w <- .RGtkCall("S_gdk_app_launch_context_set_timestamp", object, timestamp, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkAppLaunchContextSetIcon <-
function(object, icon = NULL)
{
  checkPtrType(object, "GdkAppLaunchContext")
  if (!is.null( icon )) checkPtrType(icon, "GIcon")

  w <- .RGtkCall("S_gdk_app_launch_context_set_icon", object, icon, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkAppLaunchContextSetIconName <-
function(object, icon.name = NULL)
{
  checkPtrType(object, "GdkAppLaunchContext")
  if (!is.null( icon.name )) icon.name <- as.character(icon.name)

  w <- .RGtkCall("S_gdk_app_launch_context_set_icon_name", object, icon.name, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkScreenGetMonitorWidthMm <-
function(object, monitor.num)
{
  checkPtrType(object, "GdkScreen")
  monitor.num <- as.integer(monitor.num)

  w <- .RGtkCall("S_gdk_screen_get_monitor_width_mm", object, monitor.num, PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenGetMonitorHeightMm <-
function(object, monitor.num)
{
  checkPtrType(object, "GdkScreen")
  monitor.num <- as.integer(monitor.num)

  w <- .RGtkCall("S_gdk_screen_get_monitor_height_mm", object, monitor.num, PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenGetMonitorPlugName <-
function(object, monitor.num)
{
  checkPtrType(object, "GdkScreen")
  monitor.num <- as.integer(monitor.num)

  w <- .RGtkCall("S_gdk_screen_get_monitor_plug_name", object, monitor.num, PACKAGE = "RGtk2")

  return(w)
} 


gdkWindowRedirectToDrawable <-
function(object, drawable, src.x, src.y, dest.x, dest.y, width, height)
{
  checkPtrType(object, "GdkWindow")
  checkPtrType(drawable, "GdkDrawable")
  src.x <- as.integer(src.x)
  src.y <- as.integer(src.y)
  dest.x <- as.integer(dest.x)
  dest.y <- as.integer(dest.y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gdk_window_redirect_to_drawable", object, drawable, src.x, src.y, dest.x, dest.y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowRemoveRedirection <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_remove_redirection", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkPixbufNewFromStream <-
function(stream, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(stream, "GInputStream")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gdk_pixbuf_new_from_stream", stream, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gdkPixbufNewFromStreamAtScale <-
function(stream, width = -1, height = -1, preserve.aspect.ratio = 1, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(stream, "GInputStream")
  width <- as.integer(width)
  height <- as.integer(height)
  preserve.aspect.ratio <- as.logical(preserve.aspect.ratio)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gdk_pixbuf_new_from_stream_at_scale", stream, width, height, preserve.aspect.ratio, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gdkTestRenderSync <-
function(window)
{
  checkPtrType(window, "GdkWindow")

  w <- .RGtkCall("S_gdk_test_render_sync", window, PACKAGE = "RGtk2")

  return(w)
} 


gdkTestSimulateKey <-
function(window, x, y, keyval, modifiers, key.pressrelease)
{
  checkPtrType(window, "GdkWindow")
  x <- as.integer(x)
  y <- as.integer(y)
  keyval <- as.numeric(keyval)
  
  

  w <- .RGtkCall("S_gdk_test_simulate_key", window, x, y, keyval, modifiers, key.pressrelease, PACKAGE = "RGtk2")

  return(w)
} 


gdkTestSimulateButton <-
function(window, x, y, button, modifiers, button.pressrelease)
{
  checkPtrType(window, "GdkWindow")
  x <- as.integer(x)
  y <- as.integer(y)
  button <- as.numeric(button)
  
  

  w <- .RGtkCall("S_gdk_test_simulate_button", window, x, y, button, modifiers, button.pressrelease, PACKAGE = "RGtk2")

  return(w)
} 


gdkPixbufSaveToStream <-
function(object, stream, type, cancellable, .errwarn = TRUE)
{
  checkPtrType(object, "GdkPixbuf")
  checkPtrType(stream, "GOutputStream")
  type <- as.character(type)
  checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gdk_pixbuf_save_to_stream", object, stream, type, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gdkKeymapGetCapsLockState <-
function(object)
{
  checkPtrType(object, "GdkKeymap")

  w <- .RGtkCall("S_gdk_keymap_get_caps_lock_state", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkCairoResetClip <-
function(cr, drawable)
{
  checkPtrType(cr, "Cairo")
  checkPtrType(drawable, "GdkDrawable")

  w <- .RGtkCall("S_gdk_cairo_reset_clip", cr, drawable, PACKAGE = "RGtk2")

  return(w)
} 


gdkOffscreenWindowGetPixmap <-
function(window)
{
  checkPtrType(window, "GdkWindow")

  w <- .RGtkCall("S_gdk_offscreen_window_get_pixmap", window, PACKAGE = "RGtk2")

  return(w)
} 


gdkOffscreenWindowSetEmbedder <-
function(window, embedder)
{
  checkPtrType(window, "GdkWindow")
  checkPtrType(embedder, "GdkWindow")

  w <- .RGtkCall("S_gdk_offscreen_window_set_embedder", window, embedder, PACKAGE = "RGtk2")

  return(w)
} 


gdkOffscreenWindowGetEmbedder <-
function(window)
{
  checkPtrType(window, "GdkWindow")

  w <- .RGtkCall("S_gdk_offscreen_window_get_embedder", window, PACKAGE = "RGtk2")

  return(w)
} 


gdkRegionRectEqual <-
function(object, rectangle)
{
  checkPtrType(object, "GdkRegion")
  rectangle <- as.GdkRectangle(rectangle)

  w <- .RGtkCall("S_gdk_region_rect_equal", object, rectangle, PACKAGE = "RGtk2")

  return(w)
} 


gdkWindowEnsureNative <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_ensure_native", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkWindowFlush <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_flush", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowGeometryChanged <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_geometry_changed", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowGetCursor <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_get_cursor", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkWindowRestack <-
function(object, sibling, above)
{
  checkPtrType(object, "GdkWindow")
  checkPtrType(sibling, "GdkWindow")
  above <- as.logical(above)

  w <- .RGtkCall("S_gdk_window_restack", object, sibling, above, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkWindowIsDestroyed <-
function(object)
{
  checkPtrType(object, "GdkWindow")

  w <- .RGtkCall("S_gdk_window_is_destroyed", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkWindowGetRootCoords <-
function(object, x, y)
{
  checkPtrType(object, "GdkWindow")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_gdk_window_get_root_coords", object, x, y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkPixbufSimpleAnimSetLoop <-
function(object, loop)
{
  checkPtrType(object, "GdkPixbufSimpleAnim")
  loop <- as.logical(loop)

  w <- .RGtkCall("S_gdk_pixbuf_simple_anim_set_loop", object, loop, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gdkPixbufSimpleAnimGetLoop <-
function(object)
{
  checkPtrType(object, "GdkPixbufSimpleAnim")

  w <- .RGtkCall("S_gdk_pixbuf_simple_anim_get_loop", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkKeymapAddVirtualModifiers <-
function(object)
{
  checkPtrType(object, "GdkKeymap")

  w <- .RGtkCall("S_gdk_keymap_add_virtual_modifiers", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkKeymapMapVirtualModifiers <-
function(object)
{
  checkPtrType(object, "GdkKeymap")

  w <- .RGtkCall("S_gdk_keymap_map_virtual_modifiers", object, PACKAGE = "RGtk2")

  return(w)
} 


gdkScreenGetPrimaryMonitor <-
function(object)
{
  checkPtrType(object, "GdkScreen")

  w <- .RGtkCall("S_gdk_screen_get_primary_monitor", object, PACKAGE = "RGtk2")

  return(w)
} 

