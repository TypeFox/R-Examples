
gtkObjectFlags <-
function(object)
{
  checkPtrType(object, "GtkObject")

  w <- .RGtkCall("S_GTK_OBJECT_FLAGS", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetSetFlags <-
function(wid, flags)
{
  checkPtrType(wid, "GtkWidget")
  

  w <- .RGtkCall("S_GTK_WIDGET_SET_FLAGS", wid, flags, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetUnsetFlags <-
function(wid, flags)
{
  checkPtrType(wid, "GtkWidget")
  

  w <- .RGtkCall("S_GTK_WIDGET_UNSET_FLAGS", wid, flags, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetIsSensitive <-
function(wid)
{
  checkPtrType(wid, "GtkWidget")

  w <- .RGtkCall("S_GTK_WIDGET_IS_SENSITIVE", wid, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetState <-
function(wid)
{
  checkPtrType(wid, "GtkWidget")

  w <- .RGtkCall("S_GTK_WIDGET_STATE", wid, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetSavedState <-
function(wid)
{
  checkPtrType(wid, "GtkWidget")

  w <- .RGtkCall("S_GTK_WIDGET_SAVED_STATE", wid, PACKAGE = "RGtk2")

  return(w)
} 


gtkCTreeRow <-
function(node)
{
  checkPtrType(node, "GtkCTreeNode")

  w <- .RGtkCall("S_GTK_CTREE_ROW", node, PACKAGE = "RGtk2")

  return(w)
} 


gtkAboutDialogGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_about_dialog_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkAboutDialogNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_about_dialog_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkAboutDialogGetName <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("getProgramName", "RGtk2")

  checkPtrType(object, "GtkAboutDialog")

  w <- .RGtkCall("S_gtk_about_dialog_get_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkAboutDialogSetName <-
function(object, name = NULL)
{
  if(getOption("depwarn"))
    .Deprecated("setProgramName", "RGtk2")

  checkPtrType(object, "GtkAboutDialog")
  if (!is.null( name )) name <- as.character(name)

  w <- .RGtkCall("S_gtk_about_dialog_set_name", object, name, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAboutDialogGetVersion <-
function(object)
{
  checkPtrType(object, "GtkAboutDialog")

  w <- .RGtkCall("S_gtk_about_dialog_get_version", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkAboutDialogSetVersion <-
function(object, version = NULL)
{
  checkPtrType(object, "GtkAboutDialog")
  if (!is.null( version )) version <- as.character(version)

  w <- .RGtkCall("S_gtk_about_dialog_set_version", object, version, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAboutDialogGetCopyright <-
function(object)
{
  checkPtrType(object, "GtkAboutDialog")

  w <- .RGtkCall("S_gtk_about_dialog_get_copyright", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkAboutDialogSetCopyright <-
function(object, copyright = NULL)
{
  checkPtrType(object, "GtkAboutDialog")
  if (!is.null( copyright )) copyright <- as.character(copyright)

  w <- .RGtkCall("S_gtk_about_dialog_set_copyright", object, copyright, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAboutDialogGetComments <-
function(object)
{
  checkPtrType(object, "GtkAboutDialog")

  w <- .RGtkCall("S_gtk_about_dialog_get_comments", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkAboutDialogSetComments <-
function(object, comments = NULL)
{
  checkPtrType(object, "GtkAboutDialog")
  if (!is.null( comments )) comments <- as.character(comments)

  w <- .RGtkCall("S_gtk_about_dialog_set_comments", object, comments, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAboutDialogGetLicense <-
function(object)
{
  checkPtrType(object, "GtkAboutDialog")

  w <- .RGtkCall("S_gtk_about_dialog_get_license", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkAboutDialogSetLicense <-
function(object, license = NULL)
{
  checkPtrType(object, "GtkAboutDialog")
  if (!is.null( license )) license <- as.character(license)

  w <- .RGtkCall("S_gtk_about_dialog_set_license", object, license, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAboutDialogGetWrapLicense <-
function(object)
{
  checkPtrType(object, "GtkAboutDialog")

  w <- .RGtkCall("S_gtk_about_dialog_get_wrap_license", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkAboutDialogSetWrapLicense <-
function(object, wrap.license)
{
  checkPtrType(object, "GtkAboutDialog")
  wrap.license <- as.logical(wrap.license)

  w <- .RGtkCall("S_gtk_about_dialog_set_wrap_license", object, wrap.license, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAboutDialogGetWebsite <-
function(object)
{
  checkPtrType(object, "GtkAboutDialog")

  w <- .RGtkCall("S_gtk_about_dialog_get_website", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkAboutDialogSetWebsite <-
function(object, website = NULL)
{
  checkPtrType(object, "GtkAboutDialog")
  if (!is.null( website )) website <- as.character(website)

  w <- .RGtkCall("S_gtk_about_dialog_set_website", object, website, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAboutDialogGetWebsiteLabel <-
function(object)
{
  checkPtrType(object, "GtkAboutDialog")

  w <- .RGtkCall("S_gtk_about_dialog_get_website_label", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkAboutDialogSetWebsiteLabel <-
function(object, website.label = NULL)
{
  checkPtrType(object, "GtkAboutDialog")
  if (!is.null( website.label )) website.label <- as.character(website.label)

  w <- .RGtkCall("S_gtk_about_dialog_set_website_label", object, website.label, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAboutDialogGetAuthors <-
function(object)
{
  checkPtrType(object, "GtkAboutDialog")

  w <- .RGtkCall("S_gtk_about_dialog_get_authors", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkAboutDialogSetAuthors <-
function(object, authors)
{
  checkPtrType(object, "GtkAboutDialog")
  authors <- as.list(as.character(authors))

  w <- .RGtkCall("S_gtk_about_dialog_set_authors", object, authors, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAboutDialogGetDocumenters <-
function(object)
{
  checkPtrType(object, "GtkAboutDialog")

  w <- .RGtkCall("S_gtk_about_dialog_get_documenters", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkAboutDialogSetDocumenters <-
function(object, documenters)
{
  checkPtrType(object, "GtkAboutDialog")
  documenters <- as.list(as.character(documenters))

  w <- .RGtkCall("S_gtk_about_dialog_set_documenters", object, documenters, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAboutDialogGetArtists <-
function(object)
{
  checkPtrType(object, "GtkAboutDialog")

  w <- .RGtkCall("S_gtk_about_dialog_get_artists", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkAboutDialogSetArtists <-
function(object, artists)
{
  checkPtrType(object, "GtkAboutDialog")
  artists <- as.list(as.character(artists))

  w <- .RGtkCall("S_gtk_about_dialog_set_artists", object, artists, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAboutDialogGetTranslatorCredits <-
function(object)
{
  checkPtrType(object, "GtkAboutDialog")

  w <- .RGtkCall("S_gtk_about_dialog_get_translator_credits", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkAboutDialogSetTranslatorCredits <-
function(object, translator.credits = NULL)
{
  checkPtrType(object, "GtkAboutDialog")
  if (!is.null( translator.credits )) translator.credits <- as.character(translator.credits)

  w <- .RGtkCall("S_gtk_about_dialog_set_translator_credits", object, translator.credits, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAboutDialogGetLogo <-
function(object)
{
  checkPtrType(object, "GtkAboutDialog")

  w <- .RGtkCall("S_gtk_about_dialog_get_logo", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkAboutDialogSetLogo <-
function(object, logo = NULL)
{
  checkPtrType(object, "GtkAboutDialog")
  if (!is.null( logo )) checkPtrType(logo, "GdkPixbuf")

  w <- .RGtkCall("S_gtk_about_dialog_set_logo", object, logo, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAboutDialogGetLogoIconName <-
function(object)
{
  checkPtrType(object, "GtkAboutDialog")

  w <- .RGtkCall("S_gtk_about_dialog_get_logo_icon_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkAboutDialogSetLogoIconName <-
function(object, icon.name = NULL)
{
  checkPtrType(object, "GtkAboutDialog")
  if (!is.null( icon.name )) icon.name <- as.character(icon.name)

  w <- .RGtkCall("S_gtk_about_dialog_set_logo_icon_name", object, icon.name, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAboutDialogSetEmailHook <-
function(func, data = NULL)
{
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_about_dialog_set_email_hook", func, data, PACKAGE = "RGtk2")

  return(w)
} 


gtkAboutDialogSetUrlHook <-
function(func, data = NULL)
{
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_about_dialog_set_url_hook", func, data, PACKAGE = "RGtk2")

  return(w)
} 


gtkAccelGroupGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_accel_group_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkAccelGroupNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_accel_group_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkAccelGroupLock <-
function(object)
{
  checkPtrType(object, "GtkAccelGroup")

  w <- .RGtkCall("S_gtk_accel_group_lock", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAccelGroupUnlock <-
function(object)
{
  checkPtrType(object, "GtkAccelGroup")

  w <- .RGtkCall("S_gtk_accel_group_unlock", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAccelGroupConnect <-
function(object, accel.key, accel.mods, accel.flags, closure)
{
  checkPtrType(object, "GtkAccelGroup")
  accel.key <- as.numeric(accel.key)
  
  
  closure <- as.GClosure(closure)

  w <- .RGtkCall("S_gtk_accel_group_connect", object, accel.key, accel.mods, accel.flags, closure, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAccelGroupConnectByPath <-
function(object, accel.path, closure)
{
  checkPtrType(object, "GtkAccelGroup")
  accel.path <- as.character(accel.path)
  closure <- as.GClosure(closure)

  w <- .RGtkCall("S_gtk_accel_group_connect_by_path", object, accel.path, closure, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAccelGroupDisconnect <-
function(object, closure)
{
  checkPtrType(object, "GtkAccelGroup")
  closure <- as.GClosure(closure)

  w <- .RGtkCall("S_gtk_accel_group_disconnect", object, closure, PACKAGE = "RGtk2")

  return(w)
} 


gtkAccelGroupDisconnectKey <-
function(object, accel.key, accel.mods)
{
  checkPtrType(object, "GtkAccelGroup")
  accel.key <- as.numeric(accel.key)
  

  w <- .RGtkCall("S_gtk_accel_group_disconnect_key", object, accel.key, accel.mods, PACKAGE = "RGtk2")

  return(w)
} 


gtkAccelGroupActivate <-
function(object, accel.quark, acceleratable, accel.key, accel.mods)
{
  checkPtrType(object, "GtkAccelGroup")
  accel.quark <- as.GQuark(accel.quark)
  checkPtrType(acceleratable, "GObject")
  accel.key <- as.numeric(accel.key)
  

  w <- .RGtkCall("S_gtk_accel_group_activate", object, accel.quark, acceleratable, accel.key, accel.mods, PACKAGE = "RGtk2")

  return(w)
} 


gtkAccelGroupsActivate <-
function(object, accel.key, accel.mods)
{
  checkPtrType(object, "GObject")
  accel.key <- as.numeric(accel.key)
  

  w <- .RGtkCall("S_gtk_accel_groups_activate", object, accel.key, accel.mods, PACKAGE = "RGtk2")

  return(w)
} 


gtkAccelGroupsFromObject <-
function(object)
{
  checkPtrType(object, "GObject")

  w <- .RGtkCall("S_gtk_accel_groups_from_object", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkAccelGroupFind <-
function(object, find.func, data = NULL)
{
  checkPtrType(object, "GtkAccelGroup")
  find.func <- as.function(find.func)
  

  w <- .RGtkCall("S_gtk_accel_group_find", object, find.func, data, PACKAGE = "RGtk2")

  return(w)
} 


gtkAccelGroupFromAccelClosure <-
function(closure)
{
  closure <- as.GClosure(closure)

  w <- .RGtkCall("S_gtk_accel_group_from_accel_closure", closure, PACKAGE = "RGtk2")

  return(w)
} 


gtkAcceleratorValid <-
function(keyval, modifiers)
{
  keyval <- as.numeric(keyval)
  

  w <- .RGtkCall("S_gtk_accelerator_valid", keyval, modifiers, PACKAGE = "RGtk2")

  return(w)
} 


gtkAcceleratorParse <-
function(accelerator)
{
  accelerator <- as.character(accelerator)

  w <- .RGtkCall("S_gtk_accelerator_parse", accelerator, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAcceleratorName <-
function(accelerator.key, accelerator.mods)
{
  accelerator.key <- as.numeric(accelerator.key)
  

  w <- .RGtkCall("S_gtk_accelerator_name", accelerator.key, accelerator.mods, PACKAGE = "RGtk2")

  return(w)
} 


gtkAcceleratorSetDefaultModMask <-
function(default.mod.mask)
{
  

  w <- .RGtkCall("S_gtk_accelerator_set_default_mod_mask", default.mod.mask, PACKAGE = "RGtk2")

  return(w)
} 


gtkAcceleratorGetDefaultModMask <-
function()
{
  

  w <- .RGtkCall("S_gtk_accelerator_get_default_mod_mask", PACKAGE = "RGtk2")

  return(w)
} 


gtkAcceleratorGetLabel <-
function(accelerator.key, accelerator.mods)
{
  accelerator.key <- as.numeric(accelerator.key)
  

  w <- .RGtkCall("S_gtk_accelerator_get_label", accelerator.key, accelerator.mods, PACKAGE = "RGtk2")

  return(w)
} 


gtkAccelLabelGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_accel_label_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkAccelLabelNew <-
function(string = NULL, show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_accel_label_new", string, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkAccelLabelGetAccelWidget <-
function(object)
{
  checkPtrType(object, "GtkAccelLabel")

  w <- .RGtkCall("S_gtk_accel_label_get_accel_widget", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkAccelLabelGetAccelWidth <-
function(object)
{
  checkPtrType(object, "GtkAccelLabel")

  w <- .RGtkCall("S_gtk_accel_label_get_accel_width", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkAccelLabelSetAccelWidget <-
function(object, accel.widget)
{
  checkPtrType(object, "GtkAccelLabel")
  checkPtrType(accel.widget, "GtkWidget")

  w <- .RGtkCall("S_gtk_accel_label_set_accel_widget", object, accel.widget, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAccelLabelSetAccelClosure <-
function(object, accel.closure)
{
  checkPtrType(object, "GtkAccelLabel")
  accel.closure <- as.GClosure(accel.closure)

  w <- .RGtkCall("S_gtk_accel_label_set_accel_closure", object, accel.closure, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAccelLabelRefetch <-
function(object)
{
  checkPtrType(object, "GtkAccelLabel")

  w <- .RGtkCall("S_gtk_accel_label_refetch", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkAccelMapAddEntry <-
function(accel.path, accel.key, accel.mods)
{
  accel.path <- as.character(accel.path)
  accel.key <- as.numeric(accel.key)
  

  w <- .RGtkCall("S_gtk_accel_map_add_entry", accel.path, accel.key, accel.mods, PACKAGE = "RGtk2")

  return(w)
} 


gtkAccelMapLookupEntry <-
function(accel.path)
{
  accel.path <- as.character(accel.path)

  w <- .RGtkCall("S_gtk_accel_map_lookup_entry", accel.path, PACKAGE = "RGtk2")

  return(w)
} 


gtkAccelMapChangeEntry <-
function(accel.path, accel.key, accel.mods, replace)
{
  accel.path <- as.character(accel.path)
  accel.key <- as.numeric(accel.key)
  
  replace <- as.logical(replace)

  w <- .RGtkCall("S_gtk_accel_map_change_entry", accel.path, accel.key, accel.mods, replace, PACKAGE = "RGtk2")

  return(w)
} 


gtkAccelMapLoad <-
function(file.name)
{
  file.name <- as.character(file.name)

  w <- .RGtkCall("S_gtk_accel_map_load", file.name, PACKAGE = "RGtk2")

  return(w)
} 


gtkAccelMapSave <-
function(file.name)
{
  file.name <- as.character(file.name)

  w <- .RGtkCall("S_gtk_accel_map_save", file.name, PACKAGE = "RGtk2")

  return(w)
} 


gtkAccelMapForeach <-
function(data = NULL, foreach.func)
{
  
  foreach.func <- as.function(foreach.func)

  w <- .RGtkCall("S_gtk_accel_map_foreach", data, foreach.func, PACKAGE = "RGtk2")

  return(w)
} 


gtkAccelMapLoadFd <-
function(fd)
{
  fd <- as.integer(fd)

  w <- .RGtkCall("S_gtk_accel_map_load_fd", fd, PACKAGE = "RGtk2")

  return(w)
} 


gtkAccelMapLoadScanner <-
function(scanner)
{
  checkPtrType(scanner, "GScanner")

  w <- .RGtkCall("S_gtk_accel_map_load_scanner", scanner, PACKAGE = "RGtk2")

  return(w)
} 


gtkAccelMapSaveFd <-
function(fd)
{
  fd <- as.integer(fd)

  w <- .RGtkCall("S_gtk_accel_map_save_fd", fd, PACKAGE = "RGtk2")

  return(w)
} 


gtkAccelMapLockPath <-
function(accel.path)
{
  accel.path <- as.character(accel.path)

  w <- .RGtkCall("S_gtk_accel_map_lock_path", accel.path, PACKAGE = "RGtk2")

  return(w)
} 


gtkAccelMapUnlockPath <-
function(accel.path)
{
  accel.path <- as.character(accel.path)

  w <- .RGtkCall("S_gtk_accel_map_unlock_path", accel.path, PACKAGE = "RGtk2")

  return(w)
} 


gtkAccelMapAddFilter <-
function(filter.pattern)
{
  filter.pattern <- as.character(filter.pattern)

  w <- .RGtkCall("S_gtk_accel_map_add_filter", filter.pattern, PACKAGE = "RGtk2")

  return(w)
} 


gtkAccelMapForeachUnfiltered <-
function(data = NULL, foreach.func)
{
  
  foreach.func <- as.function(foreach.func)

  w <- .RGtkCall("S_gtk_accel_map_foreach_unfiltered", data, foreach.func, PACKAGE = "RGtk2")

  return(w)
} 


gtkAccelMapGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_accel_map_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkAccelMapGet <-
function()
{
  

  w <- .RGtkCall("S_gtk_accel_map_get", PACKAGE = "RGtk2")

  return(w)
} 


gtkAccessibleGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_accessible_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkAccessibleConnectWidgetDestroyed <-
function(object)
{
  checkPtrType(object, "GtkAccessible")

  w <- .RGtkCall("S_gtk_accessible_connect_widget_destroyed", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActionGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_action_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkActionNew <-
function(name = NULL, label = NULL, tooltip = NULL, stock.id = NULL)
{
  

  w <- .RGtkCall("S_gtk_action_new", name, label, tooltip, stock.id, PACKAGE = "RGtk2")

  return(w)
} 


gtkActionGetName <-
function(object)
{
  checkPtrType(object, "GtkAction")

  w <- .RGtkCall("S_gtk_action_get_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkActionIsSensitive <-
function(object)
{
  checkPtrType(object, "GtkAction")

  w <- .RGtkCall("S_gtk_action_is_sensitive", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkActionGetSensitive <-
function(object)
{
  checkPtrType(object, "GtkAction")

  w <- .RGtkCall("S_gtk_action_get_sensitive", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkActionIsVisible <-
function(object)
{
  checkPtrType(object, "GtkAction")

  w <- .RGtkCall("S_gtk_action_is_visible", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkActionGetVisible <-
function(object)
{
  checkPtrType(object, "GtkAction")

  w <- .RGtkCall("S_gtk_action_get_visible", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkActionActivate <-
function(object)
{
  checkPtrType(object, "GtkAction")

  w <- .RGtkCall("S_gtk_action_activate", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActionCreateIcon <-
function(object, icon.size)
{
  checkPtrType(object, "GtkAction")
  

  w <- .RGtkCall("S_gtk_action_create_icon", object, icon.size, PACKAGE = "RGtk2")

  return(w)
} 


gtkActionCreateMenuItem <-
function(object)
{
  checkPtrType(object, "GtkAction")

  w <- .RGtkCall("S_gtk_action_create_menu_item", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkActionCreateToolItem <-
function(object)
{
  checkPtrType(object, "GtkAction")

  w <- .RGtkCall("S_gtk_action_create_tool_item", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkActionConnectProxy <-
function(object, proxy)
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  checkPtrType(object, "GtkAction")
  checkPtrType(proxy, "GtkWidget")

  w <- .RGtkCall("S_gtk_action_connect_proxy", object, proxy, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActionDisconnectProxy <-
function(object, proxy)
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  checkPtrType(object, "GtkAction")
  checkPtrType(proxy, "GtkWidget")

  w <- .RGtkCall("S_gtk_action_disconnect_proxy", object, proxy, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActionGetProxies <-
function(object)
{
  checkPtrType(object, "GtkAction")

  w <- .RGtkCall("S_gtk_action_get_proxies", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkActionConnectAccelerator <-
function(object)
{
  checkPtrType(object, "GtkAction")

  w <- .RGtkCall("S_gtk_action_connect_accelerator", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActionDisconnectAccelerator <-
function(object)
{
  checkPtrType(object, "GtkAction")

  w <- .RGtkCall("S_gtk_action_disconnect_accelerator", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActionGetAccelPath <-
function(object)
{
  checkPtrType(object, "GtkAction")

  w <- .RGtkCall("S_gtk_action_get_accel_path", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkActionGetAccelClosure <-
function(object)
{
  checkPtrType(object, "GtkAction")

  w <- .RGtkCall("S_gtk_action_get_accel_closure", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkActionBlockActivateFrom <-
function(object, proxy)
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  checkPtrType(object, "GtkAction")
  checkPtrType(proxy, "GtkWidget")

  w <- .RGtkCall("S_gtk_action_block_activate_from", object, proxy, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActionUnblockActivateFrom <-
function(object, proxy)
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  checkPtrType(object, "GtkAction")
  checkPtrType(proxy, "GtkWidget")

  w <- .RGtkCall("S_gtk_action_unblock_activate_from", object, proxy, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActionSetAccelPath <-
function(object, accel.path)
{
  checkPtrType(object, "GtkAction")
  accel.path <- as.character(accel.path)

  w <- .RGtkCall("S_gtk_action_set_accel_path", object, accel.path, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActionSetAccelGroup <-
function(object, accel.group)
{
  checkPtrType(object, "GtkAction")
  checkPtrType(accel.group, "GtkAccelGroup")

  w <- .RGtkCall("S_gtk_action_set_accel_group", object, accel.group, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActionSetSensitive <-
function(object, sensitive)
{
  checkPtrType(object, "GtkAction")
  sensitive <- as.logical(sensitive)

  w <- .RGtkCall("S_gtk_action_set_sensitive", object, sensitive, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActionSetVisible <-
function(object, visible)
{
  checkPtrType(object, "GtkAction")
  visible <- as.logical(visible)

  w <- .RGtkCall("S_gtk_action_set_visible", object, visible, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActionGroupGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_action_group_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkActionGroupNew <-
function(name = NULL)
{
  

  w <- .RGtkCall("S_gtk_action_group_new", name, PACKAGE = "RGtk2")

  return(w)
} 


gtkActionGroupGetName <-
function(object)
{
  checkPtrType(object, "GtkActionGroup")

  w <- .RGtkCall("S_gtk_action_group_get_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkActionGroupGetSensitive <-
function(object)
{
  checkPtrType(object, "GtkActionGroup")

  w <- .RGtkCall("S_gtk_action_group_get_sensitive", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkActionGroupSetSensitive <-
function(object, sensitive)
{
  checkPtrType(object, "GtkActionGroup")
  sensitive <- as.logical(sensitive)

  w <- .RGtkCall("S_gtk_action_group_set_sensitive", object, sensitive, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActionGroupGetVisible <-
function(object)
{
  checkPtrType(object, "GtkActionGroup")

  w <- .RGtkCall("S_gtk_action_group_get_visible", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkActionGroupSetVisible <-
function(object, visible)
{
  checkPtrType(object, "GtkActionGroup")
  visible <- as.logical(visible)

  w <- .RGtkCall("S_gtk_action_group_set_visible", object, visible, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActionGroupGetAction <-
function(object, action.name)
{
  checkPtrType(object, "GtkActionGroup")
  action.name <- as.character(action.name)

  w <- .RGtkCall("S_gtk_action_group_get_action", object, action.name, PACKAGE = "RGtk2")

  return(w)
} 


gtkActionGroupListActions <-
function(object)
{
  checkPtrType(object, "GtkActionGroup")

  w <- .RGtkCall("S_gtk_action_group_list_actions", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkActionGroupAddAction <-
function(object, action)
{
  checkPtrType(object, "GtkActionGroup")
  checkPtrType(action, "GtkAction")

  w <- .RGtkCall("S_gtk_action_group_add_action", object, action, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActionGroupAddActionWithAccel <-
function(object, action, accelerator = NULL)
{
  checkPtrType(object, "GtkActionGroup")
  checkPtrType(action, "GtkAction")
  if (!is.null( accelerator )) accelerator <- as.character(accelerator)

  w <- .RGtkCall("S_gtk_action_group_add_action_with_accel", object, action, accelerator, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActionGroupRemoveAction <-
function(object, action)
{
  checkPtrType(object, "GtkActionGroup")
  checkPtrType(action, "GtkAction")

  w <- .RGtkCall("S_gtk_action_group_remove_action", object, action, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActionGroupAddActions <-
function(object, entries, user.data = NULL)
{
  checkPtrType(object, "GtkActionGroup")
  entries <- lapply(entries, function(x) { x <- as.GtkActionEntry(x); x })
  

  w <- .RGtkCall("S_gtk_action_group_add_actions", object, entries, user.data, PACKAGE = "RGtk2")

  return(w)
} 


gtkActionGroupAddToggleActions <-
function(object, entries, user.data = NULL)
{
  checkPtrType(object, "GtkActionGroup")
  entries <- lapply(entries, function(x) { x <- as.GtkToggleActionEntry(x); x })
  

  w <- .RGtkCall("S_gtk_action_group_add_toggle_actions", object, entries, user.data, PACKAGE = "RGtk2")

  return(w)
} 


gtkActionGroupAddRadioActions <-
function(object, entries, value, on.change = NULL, user.data = NULL)
{
  checkPtrType(object, "GtkActionGroup")
  entries <- lapply(entries, function(x) { x <- as.GtkRadioActionEntry(x); x })
  value <- as.integer(value)
  if (!is.null( on.change )) on.change <- as.function(on.change)
  

  w <- .RGtkCall("S_gtk_action_group_add_radio_actions", object, entries, value, on.change, user.data, PACKAGE = "RGtk2")

  return(w)
} 


gtkActionGroupAddActionsFull <-
function(object, entries, user.data = NULL)
{
  checkPtrType(object, "GtkActionGroup")
  entries <- lapply(entries, function(x) { x <- as.GtkActionEntry(x); x })
  

  w <- .RGtkCall("S_gtk_action_group_add_actions_full", object, entries, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActionGroupAddToggleActionsFull <-
function(object, entries, user.data = NULL)
{
  checkPtrType(object, "GtkActionGroup")
  entries <- lapply(entries, function(x) { x <- as.GtkToggleActionEntry(x); x })
  

  w <- .RGtkCall("S_gtk_action_group_add_toggle_actions_full", object, entries, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActionGroupAddRadioActionsFull <-
function(object, entries, value, on.change = NULL, user.data = NULL)
{
  checkPtrType(object, "GtkActionGroup")
  entries <- lapply(entries, function(x) { x <- as.GtkRadioActionEntry(x); x })
  value <- as.integer(value)
  if (!is.null( on.change )) on.change <- as.function(on.change)
  

  w <- .RGtkCall("S_gtk_action_group_add_radio_actions_full", object, entries, value, on.change, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActionGroupSetTranslateFunc <-
function(object, func, data = NULL)
{
  checkPtrType(object, "GtkActionGroup")
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_action_group_set_translate_func", object, func, data, PACKAGE = "RGtk2")

  return(w)
} 


gtkActionGroupSetTranslationDomain <-
function(object, domain)
{
  checkPtrType(object, "GtkActionGroup")
  domain <- as.character(domain)

  w <- .RGtkCall("S_gtk_action_group_set_translation_domain", object, domain, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActionGroupTranslateString <-
function(object, string)
{
  checkPtrType(object, "GtkActionGroup")
  string <- as.character(string)

  w <- .RGtkCall("S_gtk_action_group_translate_string", object, string, PACKAGE = "RGtk2")

  return(w)
} 


gtkAdjustmentGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_adjustment_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkAdjustmentNew <-
function(value = NULL, lower = NULL, upper = NULL, step.incr = NULL, page.incr = NULL, page.size = NULL)
{
  

  w <- .RGtkCall("S_gtk_adjustment_new", value, lower, upper, step.incr, page.incr, page.size, PACKAGE = "RGtk2")

  return(w)
} 


gtkAdjustmentChanged <-
function(object)
{
  checkPtrType(object, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_adjustment_changed", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAdjustmentValueChanged <-
function(object)
{
  checkPtrType(object, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_adjustment_value_changed", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAdjustmentClampPage <-
function(object, lower, upper)
{
  checkPtrType(object, "GtkAdjustment")
  lower <- as.numeric(lower)
  upper <- as.numeric(upper)

  w <- .RGtkCall("S_gtk_adjustment_clamp_page", object, lower, upper, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAdjustmentGetValue <-
function(object)
{
  checkPtrType(object, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_adjustment_get_value", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkAdjustmentSetValue <-
function(object, value)
{
  checkPtrType(object, "GtkAdjustment")
  value <- as.numeric(value)

  w <- .RGtkCall("S_gtk_adjustment_set_value", object, value, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAlignmentGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_alignment_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkAlignmentNew <-
function(xalign = NULL, yalign = NULL, xscale = NULL, yscale = NULL, show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_alignment_new", xalign, yalign, xscale, yscale, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkAlignmentSet <-
function(object, xalign, yalign, xscale, yscale)
{
  checkPtrType(object, "GtkAlignment")
  xalign <- as.numeric(xalign)
  yalign <- as.numeric(yalign)
  xscale <- as.numeric(xscale)
  yscale <- as.numeric(yscale)

  w <- .RGtkCall("S_gtk_alignment_set", object, xalign, yalign, xscale, yscale, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAlignmentSetPadding <-
function(object, padding.top, padding.bottom, padding.left, padding.right)
{
  checkPtrType(object, "GtkAlignment")
  padding.top <- as.numeric(padding.top)
  padding.bottom <- as.numeric(padding.bottom)
  padding.left <- as.numeric(padding.left)
  padding.right <- as.numeric(padding.right)

  w <- .RGtkCall("S_gtk_alignment_set_padding", object, padding.top, padding.bottom, padding.left, padding.right, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAlignmentGetPadding <-
function(object)
{
  checkPtrType(object, "GtkAlignment")

  w <- .RGtkCall("S_gtk_alignment_get_padding", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkArrowGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_arrow_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkArrowNew <-
function(arrow.type = NULL, shadow.type = NULL, show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_arrow_new", arrow.type, shadow.type, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkArrowSet <-
function(object, arrow.type, shadow.type)
{
  checkPtrType(object, "GtkArrow")
  
  

  w <- .RGtkCall("S_gtk_arrow_set", object, arrow.type, shadow.type, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAspectFrameGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_aspect_frame_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkAspectFrameNew <-
function(label = NULL, xalign = NULL, yalign = NULL, ratio = NULL, obey.child = NULL, show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_aspect_frame_new", label, xalign, yalign, ratio, obey.child, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkAspectFrameSet <-
function(object, xalign = 0, yalign = 0, ratio = 1, obey.child = 1)
{
  checkPtrType(object, "GtkAspectFrame")
  xalign <- as.numeric(xalign)
  yalign <- as.numeric(yalign)
  ratio <- as.numeric(ratio)
  obey.child <- as.logical(obey.child)

  w <- .RGtkCall("S_gtk_aspect_frame_set", object, xalign, yalign, ratio, obey.child, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkButtonBoxGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_button_box_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkButtonBoxGetLayout <-
function(object)
{
  checkPtrType(object, "GtkButtonBox")

  w <- .RGtkCall("S_gtk_button_box_get_layout", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkButtonBoxSetLayout <-
function(object, layout.style)
{
  checkPtrType(object, "GtkButtonBox")
  

  w <- .RGtkCall("S_gtk_button_box_set_layout", object, layout.style, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkButtonBoxGetChildSecondary <-
function(object, child)
{
  checkPtrType(object, "GtkButtonBox")
  checkPtrType(child, "GtkWidget")

  w <- .RGtkCall("S_gtk_button_box_get_child_secondary", object, child, PACKAGE = "RGtk2")

  return(w)
} 


gtkButtonBoxSetChildSecondary <-
function(object, child, is.secondary)
{
  checkPtrType(object, "GtkButtonBox")
  checkPtrType(child, "GtkWidget")
  is.secondary <- as.logical(is.secondary)

  w <- .RGtkCall("S_gtk_button_box_set_child_secondary", object, child, is.secondary, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkButtonBoxSetChildSize <-
function(object, min.width, min.height)
{
  if(getOption("depwarn"))
    .Deprecated("style properties child-min-width/-height", "RGtk2")

  checkPtrType(object, "GtkButtonBox")
  min.width <- as.integer(min.width)
  min.height <- as.integer(min.height)

  w <- .RGtkCall("S_gtk_button_box_set_child_size", object, min.width, min.height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkButtonBoxSetChildIpadding <-
function(object, ipad.x, ipad.y)
{
  if(getOption("depwarn"))
    .Deprecated("style properties child-internal-pad-x/-y", "RGtk2")

  checkPtrType(object, "GtkButtonBox")
  ipad.x <- as.integer(ipad.x)
  ipad.y <- as.integer(ipad.y)

  w <- .RGtkCall("S_gtk_button_box_set_child_ipadding", object, ipad.x, ipad.y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkButtonBoxGetChildSize <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("style properties child-min-width/-height", "RGtk2")

  checkPtrType(object, "GtkButtonBox")

  w <- .RGtkCall("S_gtk_button_box_get_child_size", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkButtonBoxGetChildIpadding <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("style properties child-internal-pad-x/-y", "RGtk2")

  checkPtrType(object, "GtkButtonBox")

  w <- .RGtkCall("S_gtk_button_box_get_child_ipadding", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkBinGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_bin_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkBinGetChild <-
function(object)
{
  checkPtrType(object, "GtkBin")

  w <- .RGtkCall("S_gtk_bin_get_child", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkBoxGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_box_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkBoxPackStart <-
function(object, child, expand = TRUE, fill = TRUE, padding = 0)
{
  checkPtrType(object, "GtkBox")
  checkPtrType(child, "GtkWidget")
  expand <- as.logical(expand)
  fill <- as.logical(fill)
  padding <- as.numeric(padding)

  w <- .RGtkCall("S_gtk_box_pack_start", object, child, expand, fill, padding, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkBoxPackEnd <-
function(object, child, expand = TRUE, fill = TRUE, padding = 0)
{
  checkPtrType(object, "GtkBox")
  checkPtrType(child, "GtkWidget")
  expand <- as.logical(expand)
  fill <- as.logical(fill)
  padding <- as.numeric(padding)

  w <- .RGtkCall("S_gtk_box_pack_end", object, child, expand, fill, padding, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkBoxPackStartDefaults <-
function(object, widget)
{
  if(getOption("depwarn"))
    .Deprecated("packStart", "RGtk2")

  checkPtrType(object, "GtkBox")
  checkPtrType(widget, "GtkWidget")

  w <- .RGtkCall("S_gtk_box_pack_start_defaults", object, widget, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkBoxPackEndDefaults <-
function(object, widget)
{
  if(getOption("depwarn"))
    .Deprecated("packEnd", "RGtk2")

  checkPtrType(object, "GtkBox")
  checkPtrType(widget, "GtkWidget")

  w <- .RGtkCall("S_gtk_box_pack_end_defaults", object, widget, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkBoxSetHomogeneous <-
function(object, homogeneous)
{
  checkPtrType(object, "GtkBox")
  homogeneous <- as.logical(homogeneous)

  w <- .RGtkCall("S_gtk_box_set_homogeneous", object, homogeneous, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkBoxGetHomogeneous <-
function(object)
{
  checkPtrType(object, "GtkBox")

  w <- .RGtkCall("S_gtk_box_get_homogeneous", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkBoxSetSpacing <-
function(object, spacing)
{
  checkPtrType(object, "GtkBox")
  spacing <- as.integer(spacing)

  w <- .RGtkCall("S_gtk_box_set_spacing", object, spacing, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkBoxGetSpacing <-
function(object)
{
  checkPtrType(object, "GtkBox")

  w <- .RGtkCall("S_gtk_box_get_spacing", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkBoxReorderChild <-
function(object, child, position)
{
  checkPtrType(object, "GtkBox")
  checkPtrType(child, "GtkWidget")
  position <- as.integer(position)

  w <- .RGtkCall("S_gtk_box_reorder_child", object, child, position, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkBoxQueryChildPacking <-
function(object, child)
{
  checkPtrType(object, "GtkBox")
  checkPtrType(child, "GtkWidget")

  w <- .RGtkCall("S_gtk_box_query_child_packing", object, child, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkBoxSetChildPacking <-
function(object, child, expand, fill, padding, pack.type)
{
  checkPtrType(object, "GtkBox")
  checkPtrType(child, "GtkWidget")
  expand <- as.logical(expand)
  fill <- as.logical(fill)
  padding <- as.numeric(padding)
  

  w <- .RGtkCall("S_gtk_box_set_child_packing", object, child, expand, fill, padding, pack.type, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkButtonGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_button_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkButtonNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_button_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkButtonNewWithLabel <-
function(label, show = TRUE)
{
  label <- as.character(label)

  w <- .RGtkCall("S_gtk_button_new_with_label", label, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkButtonNewFromStock <-
function(stock.id, show = TRUE)
{
  stock.id <- as.character(stock.id)

  w <- .RGtkCall("S_gtk_button_new_from_stock", stock.id, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkButtonNewWithMnemonic <-
function(label, show = TRUE)
{
  label <- as.character(label)

  w <- .RGtkCall("S_gtk_button_new_with_mnemonic", label, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkButtonPressed <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("'button-press-event' signal", "RGtk2")

  checkPtrType(object, "GtkButton")

  w <- .RGtkCall("S_gtk_button_pressed", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkButtonReleased <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("'button-release-event' signal", "RGtk2")

  checkPtrType(object, "GtkButton")

  w <- .RGtkCall("S_gtk_button_released", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkButtonClicked <-
function(object)
{
  checkPtrType(object, "GtkButton")

  w <- .RGtkCall("S_gtk_button_clicked", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkButtonEnter <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("'enter-notify-event' signal", "RGtk2")

  checkPtrType(object, "GtkButton")

  w <- .RGtkCall("S_gtk_button_enter", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkButtonLeave <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("'leave-notify-event' signal", "RGtk2")

  checkPtrType(object, "GtkButton")

  w <- .RGtkCall("S_gtk_button_leave", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkButtonSetRelief <-
function(object, newstyle)
{
  checkPtrType(object, "GtkButton")
  

  w <- .RGtkCall("S_gtk_button_set_relief", object, newstyle, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkButtonGetRelief <-
function(object)
{
  checkPtrType(object, "GtkButton")

  w <- .RGtkCall("S_gtk_button_get_relief", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkButtonSetLabel <-
function(object, label)
{
  checkPtrType(object, "GtkButton")
  label <- as.character(label)

  w <- .RGtkCall("S_gtk_button_set_label", object, label, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkButtonGetLabel <-
function(object)
{
  checkPtrType(object, "GtkButton")

  w <- .RGtkCall("S_gtk_button_get_label", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkButtonSetUseUnderline <-
function(object, use.underline)
{
  checkPtrType(object, "GtkButton")
  use.underline <- as.logical(use.underline)

  w <- .RGtkCall("S_gtk_button_set_use_underline", object, use.underline, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkButtonGetUseUnderline <-
function(object)
{
  checkPtrType(object, "GtkButton")

  w <- .RGtkCall("S_gtk_button_get_use_underline", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkButtonSetUseStock <-
function(object, use.stock)
{
  checkPtrType(object, "GtkButton")
  use.stock <- as.logical(use.stock)

  w <- .RGtkCall("S_gtk_button_set_use_stock", object, use.stock, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkButtonGetUseStock <-
function(object)
{
  checkPtrType(object, "GtkButton")

  w <- .RGtkCall("S_gtk_button_get_use_stock", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkButtonSetFocusOnClick <-
function(object, focus.on.click)
{
  checkPtrType(object, "GtkButton")
  focus.on.click <- as.logical(focus.on.click)

  w <- .RGtkCall("S_gtk_button_set_focus_on_click", object, focus.on.click, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkButtonGetFocusOnClick <-
function(object)
{
  checkPtrType(object, "GtkButton")

  w <- .RGtkCall("S_gtk_button_get_focus_on_click", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkButtonSetAlignment <-
function(object, xalign, yalign)
{
  checkPtrType(object, "GtkButton")
  xalign <- as.numeric(xalign)
  yalign <- as.numeric(yalign)

  w <- .RGtkCall("S_gtk_button_set_alignment", object, xalign, yalign, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkButtonGetAlignment <-
function(object)
{
  checkPtrType(object, "GtkButton")

  w <- .RGtkCall("S_gtk_button_get_alignment", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkButtonSetImage <-
function(object, image)
{
  checkPtrType(object, "GtkButton")
  checkPtrType(image, "GtkWidget")

  w <- .RGtkCall("S_gtk_button_set_image", object, image, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkButtonGetImage <-
function(object)
{
  checkPtrType(object, "GtkButton")

  w <- .RGtkCall("S_gtk_button_get_image", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkCalendarGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_calendar_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkCalendarNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_calendar_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkCalendarSelectMonth <-
function(object, month, year)
{
  checkPtrType(object, "GtkCalendar")
  month <- as.numeric(month)
  year <- as.numeric(year)

  w <- .RGtkCall("S_gtk_calendar_select_month", object, month, year, PACKAGE = "RGtk2")

  return(w)
} 


gtkCalendarSelectDay <-
function(object, day)
{
  checkPtrType(object, "GtkCalendar")
  day <- as.numeric(day)

  w <- .RGtkCall("S_gtk_calendar_select_day", object, day, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCalendarMarkDay <-
function(object, day)
{
  checkPtrType(object, "GtkCalendar")
  day <- as.numeric(day)

  w <- .RGtkCall("S_gtk_calendar_mark_day", object, day, PACKAGE = "RGtk2")

  return(w)
} 


gtkCalendarUnmarkDay <-
function(object, day)
{
  checkPtrType(object, "GtkCalendar")
  day <- as.numeric(day)

  w <- .RGtkCall("S_gtk_calendar_unmark_day", object, day, PACKAGE = "RGtk2")

  return(w)
} 


gtkCalendarClearMarks <-
function(object)
{
  checkPtrType(object, "GtkCalendar")

  w <- .RGtkCall("S_gtk_calendar_clear_marks", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCalendarSetDisplayOptions <-
function(object, flags)
{
  checkPtrType(object, "GtkCalendar")
  

  w <- .RGtkCall("S_gtk_calendar_set_display_options", object, flags, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCalendarGetDisplayOptions <-
function(object)
{
  checkPtrType(object, "GtkCalendar")

  w <- .RGtkCall("S_gtk_calendar_get_display_options", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkCalendarDisplayOptions <-
function(object, flags)
{
  if(getOption("depwarn"))
    .Deprecated("setDisplayOptions", "RGtk2")

  checkPtrType(object, "GtkCalendar")
  

  w <- .RGtkCall("S_gtk_calendar_display_options", object, flags, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCalendarGetDate <-
function(object)
{
  checkPtrType(object, "GtkCalendar")

  w <- .RGtkCall("S_gtk_calendar_get_date", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCalendarFreeze <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  checkPtrType(object, "GtkCalendar")

  w <- .RGtkCall("S_gtk_calendar_freeze", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCalendarThaw <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  checkPtrType(object, "GtkCalendar")

  w <- .RGtkCall("S_gtk_calendar_thaw", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCellEditableGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_cell_editable_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkCellEditableStartEditing <-
function(object, event = NULL)
{
  checkPtrType(object, "GtkCellEditable")
  if (!is.null( event )) checkPtrType(event, "GdkEvent")

  w <- .RGtkCall("S_gtk_cell_editable_start_editing", object, event, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCellEditableEditingDone <-
function(object)
{
  checkPtrType(object, "GtkCellEditable")

  w <- .RGtkCall("S_gtk_cell_editable_editing_done", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCellEditableRemoveWidget <-
function(object)
{
  checkPtrType(object, "GtkCellEditable")

  w <- .RGtkCall("S_gtk_cell_editable_remove_widget", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCellLayoutGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_cell_layout_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkCellLayoutPackStart <-
function(object, cell, expand = TRUE)
{
  checkPtrType(object, "GtkCellLayout")
  checkPtrType(cell, "GtkCellRenderer")
  expand <- as.logical(expand)

  w <- .RGtkCall("S_gtk_cell_layout_pack_start", object, cell, expand, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCellLayoutPackEnd <-
function(object, cell, expand = TRUE)
{
  checkPtrType(object, "GtkCellLayout")
  checkPtrType(cell, "GtkCellRenderer")
  expand <- as.logical(expand)

  w <- .RGtkCall("S_gtk_cell_layout_pack_end", object, cell, expand, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCellLayoutClear <-
function(object)
{
  checkPtrType(object, "GtkCellLayout")

  w <- .RGtkCall("S_gtk_cell_layout_clear", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCellLayoutAddAttribute <-
function(object, cell, attribute, column)
{
  checkPtrType(object, "GtkCellLayout")
  checkPtrType(cell, "GtkCellRenderer")
  attribute <- as.character(attribute)
  column <- as.integer(column)

  w <- .RGtkCall("S_gtk_cell_layout_add_attribute", object, cell, attribute, column, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCellLayoutSetCellDataFunc <-
function(object, cell, func, func.data = NULL)
{
  checkPtrType(object, "GtkCellLayout")
  checkPtrType(cell, "GtkCellRenderer")
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_cell_layout_set_cell_data_func", object, cell, func, func.data, PACKAGE = "RGtk2")

  return(w)
} 


gtkCellLayoutClearAttributes <-
function(object, cell)
{
  checkPtrType(object, "GtkCellLayout")
  checkPtrType(cell, "GtkCellRenderer")

  w <- .RGtkCall("S_gtk_cell_layout_clear_attributes", object, cell, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCellLayoutReorder <-
function(object, cell, position)
{
  checkPtrType(object, "GtkCellLayout")
  checkPtrType(cell, "GtkCellRenderer")
  position <- as.integer(position)

  w <- .RGtkCall("S_gtk_cell_layout_reorder", object, cell, position, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCellRendererGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_cell_renderer_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkCellRendererGetSize <-
function(object, widget, cell.area = NULL)
{
  checkPtrType(object, "GtkCellRenderer")
  checkPtrType(widget, "GtkWidget")
  if (!is.null( cell.area )) cell.area <- as.GdkRectangle(cell.area)

  w <- .RGtkCall("S_gtk_cell_renderer_get_size", object, widget, cell.area, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCellRendererRender <-
function(object, window, widget, background.area, cell.area, expose.area, flags)
{
  checkPtrType(object, "GtkCellRenderer")
  checkPtrType(window, "GdkWindow")
  checkPtrType(widget, "GtkWidget")
  background.area <- as.GdkRectangle(background.area)
  cell.area <- as.GdkRectangle(cell.area)
  expose.area <- as.GdkRectangle(expose.area)
  

  w <- .RGtkCall("S_gtk_cell_renderer_render", object, window, widget, background.area, cell.area, expose.area, flags, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCellRendererActivate <-
function(object, event, widget, path, background.area, cell.area, flags)
{
  checkPtrType(object, "GtkCellRenderer")
  checkPtrType(event, "GdkEvent")
  checkPtrType(widget, "GtkWidget")
  path <- as.character(path)
  background.area <- as.GdkRectangle(background.area)
  cell.area <- as.GdkRectangle(cell.area)
  

  w <- .RGtkCall("S_gtk_cell_renderer_activate", object, event, widget, path, background.area, cell.area, flags, PACKAGE = "RGtk2")

  return(w)
} 


gtkCellRendererStartEditing <-
function(object, event, widget, path, background.area, cell.area, flags)
{
  checkPtrType(object, "GtkCellRenderer")
  checkPtrType(event, "GdkEvent")
  checkPtrType(widget, "GtkWidget")
  path <- as.character(path)
  background.area <- as.GdkRectangle(background.area)
  cell.area <- as.GdkRectangle(cell.area)
  

  w <- .RGtkCall("S_gtk_cell_renderer_start_editing", object, event, widget, path, background.area, cell.area, flags, PACKAGE = "RGtk2")

  return(w)
} 


gtkCellRendererSetFixedSize <-
function(object, width, height)
{
  checkPtrType(object, "GtkCellRenderer")
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_cell_renderer_set_fixed_size", object, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCellRendererGetFixedSize <-
function(object)
{
  checkPtrType(object, "GtkCellRenderer")

  w <- .RGtkCall("S_gtk_cell_renderer_get_fixed_size", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCellRendererEditingCanceled <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("stopEditing", "RGtk2")

  checkPtrType(object, "GtkCellRenderer")

  w <- .RGtkCall("S_gtk_cell_renderer_editing_canceled", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCellRendererStopEditing <-
function(object, canceled)
{
  checkPtrType(object, "GtkCellRenderer")
  canceled <- as.logical(canceled)

  w <- .RGtkCall("S_gtk_cell_renderer_stop_editing", object, canceled, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCellRendererComboGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_cell_renderer_combo_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkCellRendererComboNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_cell_renderer_combo_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkCellRendererPixbufGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_cell_renderer_pixbuf_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkCellRendererPixbufNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_cell_renderer_pixbuf_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkCellRendererProgressGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_cell_renderer_progress_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkCellRendererProgressNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_cell_renderer_progress_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkCellRendererTextGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_cell_renderer_text_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkCellRendererTextNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_cell_renderer_text_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkCellRendererTextSetFixedHeightFromFont <-
function(object, number.of.rows)
{
  checkPtrType(object, "GtkCellRendererText")
  number.of.rows <- as.integer(number.of.rows)

  w <- .RGtkCall("S_gtk_cell_renderer_text_set_fixed_height_from_font", object, number.of.rows, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCellRendererToggleGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_cell_renderer_toggle_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkCellRendererToggleNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_cell_renderer_toggle_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkCellRendererToggleGetRadio <-
function(object)
{
  checkPtrType(object, "GtkCellRendererToggle")

  w <- .RGtkCall("S_gtk_cell_renderer_toggle_get_radio", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkCellRendererToggleSetRadio <-
function(object, radio)
{
  checkPtrType(object, "GtkCellRendererToggle")
  radio <- as.logical(radio)

  w <- .RGtkCall("S_gtk_cell_renderer_toggle_set_radio", object, radio, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCellRendererToggleGetActive <-
function(object)
{
  checkPtrType(object, "GtkCellRendererToggle")

  w <- .RGtkCall("S_gtk_cell_renderer_toggle_get_active", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkCellRendererToggleSetActive <-
function(object, setting)
{
  checkPtrType(object, "GtkCellRendererToggle")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_cell_renderer_toggle_set_active", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCellViewGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_cell_view_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkCellViewNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_cell_view_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkCellViewNewWithText <-
function(text)
{
  text <- as.character(text)

  w <- .RGtkCall("S_gtk_cell_view_new_with_text", text, PACKAGE = "RGtk2")

  return(w)
} 


gtkCellViewNewWithMarkup <-
function(markup)
{
  markup <- as.character(markup)

  w <- .RGtkCall("S_gtk_cell_view_new_with_markup", markup, PACKAGE = "RGtk2")

  return(w)
} 


gtkCellViewNewWithPixbuf <-
function(pixbuf)
{
  checkPtrType(pixbuf, "GdkPixbuf")

  w <- .RGtkCall("S_gtk_cell_view_new_with_pixbuf", pixbuf, PACKAGE = "RGtk2")

  return(w)
} 


gtkCellViewSetModel <-
function(object, model = NULL)
{
  checkPtrType(object, "GtkCellView")
  if (!is.null( model )) checkPtrType(model, "GtkTreeModel")

  w <- .RGtkCall("S_gtk_cell_view_set_model", object, model, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCellViewSetDisplayedRow <-
function(object, path = NULL)
{
  checkPtrType(object, "GtkCellView")
  if (!is.null( path )) checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_cell_view_set_displayed_row", object, path, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCellViewGetDisplayedRow <-
function(object)
{
  checkPtrType(object, "GtkCellView")

  w <- .RGtkCall("S_gtk_cell_view_get_displayed_row", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkCellViewGetSizeOfRow <-
function(object, path)
{
  checkPtrType(object, "GtkCellView")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_cell_view_get_size_of_row", object, path, PACKAGE = "RGtk2")

  return(w)
} 


gtkCellViewSetBackgroundColor <-
function(object, color)
{
  checkPtrType(object, "GtkCellView")
  color <- as.GdkColor(color)

  w <- .RGtkCall("S_gtk_cell_view_set_background_color", object, color, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCellViewGetCellRenderers <-
function(object)
{
  checkPtrType(object, "GtkCellView")

  w <- .RGtkCall("S_gtk_cell_view_get_cell_renderers", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkClipboardSetCanStore <-
function(object, targets)
{
  checkPtrType(object, "GtkClipboard")
  targets <- lapply(targets, function(x) { x <- as.GtkTargetEntry(x); x })

  w <- .RGtkCall("S_gtk_clipboard_set_can_store", object, targets, PACKAGE = "RGtk2")

  return(w)
} 


gtkClipboardStore <-
function(object)
{
  checkPtrType(object, "GtkClipboard")

  w <- .RGtkCall("S_gtk_clipboard_store", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCheckButtonGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_check_button_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkCheckButtonNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_check_button_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkCheckButtonNewWithLabel <-
function(label, show = TRUE)
{
  label <- as.character(label)

  w <- .RGtkCall("S_gtk_check_button_new_with_label", label, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkCheckButtonNewWithMnemonic <-
function(label, show = TRUE)
{
  label <- as.character(label)

  w <- .RGtkCall("S_gtk_check_button_new_with_mnemonic", label, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkCheckMenuItemGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_check_menu_item_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkCheckMenuItemNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_check_menu_item_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkCheckMenuItemNewWithLabel <-
function(label, show = TRUE)
{
  label <- as.character(label)

  w <- .RGtkCall("S_gtk_check_menu_item_new_with_label", label, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkCheckMenuItemNewWithMnemonic <-
function(label, show = TRUE)
{
  label <- as.character(label)

  w <- .RGtkCall("S_gtk_check_menu_item_new_with_mnemonic", label, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkCheckMenuItemSetActive <-
function(object, is.active)
{
  checkPtrType(object, "GtkCheckMenuItem")
  is.active <- as.logical(is.active)

  w <- .RGtkCall("S_gtk_check_menu_item_set_active", object, is.active, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCheckMenuItemGetActive <-
function(object)
{
  checkPtrType(object, "GtkCheckMenuItem")

  w <- .RGtkCall("S_gtk_check_menu_item_get_active", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkCheckMenuItemToggled <-
function(object)
{
  checkPtrType(object, "GtkCheckMenuItem")

  w <- .RGtkCall("S_gtk_check_menu_item_toggled", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCheckMenuItemSetInconsistent <-
function(object, setting)
{
  checkPtrType(object, "GtkCheckMenuItem")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_check_menu_item_set_inconsistent", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCheckMenuItemGetInconsistent <-
function(object)
{
  checkPtrType(object, "GtkCheckMenuItem")

  w <- .RGtkCall("S_gtk_check_menu_item_get_inconsistent", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkCheckMenuItemSetDrawAsRadio <-
function(object, draw.as.radio)
{
  checkPtrType(object, "GtkCheckMenuItem")
  draw.as.radio <- as.logical(draw.as.radio)

  w <- .RGtkCall("S_gtk_check_menu_item_set_draw_as_radio", object, draw.as.radio, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCheckMenuItemGetDrawAsRadio <-
function(object)
{
  checkPtrType(object, "GtkCheckMenuItem")

  w <- .RGtkCall("S_gtk_check_menu_item_get_draw_as_radio", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkCheckMenuItemSetShowToggle <-
function(object, always)
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  checkPtrType(object, "GtkCheckMenuItem")
  always <- as.logical(always)

  w <- .RGtkCall("S_gtk_check_menu_item_set_show_toggle", object, always, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCheckMenuItemSetState <-
function(object, is.active)
{
  if(getOption("depwarn"))
    .Deprecated("gtkCheckMenuItemSetActive", "RGtk2")

  checkPtrType(object, "GtkCheckMenuItem")
  is.active <- as.logical(is.active)

  w <- .RGtkCall("S_gtk_check_menu_item_set_state", object, is.active, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkClipboardGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_clipboard_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkClipboardGetForDisplay <-
function(display, selection = "GDK_SELECTION_CLIPBOARD")
{
  checkPtrType(display, "GdkDisplay")
  selection <- as.GdkAtom(selection)

  w <- .RGtkCall("S_gtk_clipboard_get_for_display", display, selection, PACKAGE = "RGtk2")

  return(w)
} 


gtkClipboardGet <-
function(selection = "GDK_SELECTION_CLIPBOARD")
{
  selection <- as.GdkAtom(selection)

  w <- .RGtkCall("S_gtk_clipboard_get", selection, PACKAGE = "RGtk2")

  return(w)
} 


gtkClipboardGetDisplay <-
function(object)
{
  checkPtrType(object, "GtkClipboard")

  w <- .RGtkCall("S_gtk_clipboard_get_display", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkClipboardSetWithData <-
function(object, targets, get.func, user.data = NULL)
{
  checkPtrType(object, "GtkClipboard")
  targets <- lapply(targets, function(x) { x <- as.GtkTargetEntry(x); x })
  get.func <- as.function(get.func)
  

  w <- .RGtkCall("S_gtk_clipboard_set_with_data", object, targets, get.func, user.data, PACKAGE = "RGtk2")

  return(w)
} 


gtkClipboardGetOwner <-
function(object)
{
  checkPtrType(object, "GtkClipboard")

  w <- .RGtkCall("S_gtk_clipboard_get_owner", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkClipboardClear <-
function(object)
{
  checkPtrType(object, "GtkClipboard")

  w <- .RGtkCall("S_gtk_clipboard_clear", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkClipboardSetText <-
function(object, text, len = -1)
{
  checkPtrType(object, "GtkClipboard")
  text <- as.character(text)
  len <- as.integer(len)

  w <- .RGtkCall("S_gtk_clipboard_set_text", object, text, len, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkClipboardSetImage <-
function(object, pixbuf)
{
  checkPtrType(object, "GtkClipboard")
  checkPtrType(pixbuf, "GdkPixbuf")

  w <- .RGtkCall("S_gtk_clipboard_set_image", object, pixbuf, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkClipboardRequestContents <-
function(object, target, callback, user.data = NULL)
{
  checkPtrType(object, "GtkClipboard")
  target <- as.GdkAtom(target)
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gtk_clipboard_request_contents", object, target, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkClipboardRequestImage <-
function(object, callback, user.data = NULL)
{
  checkPtrType(object, "GtkClipboard")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gtk_clipboard_request_image", object, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkClipboardRequestText <-
function(object, callback, user.data = NULL)
{
  checkPtrType(object, "GtkClipboard")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gtk_clipboard_request_text", object, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkClipboardRequestTargets <-
function(object, callback, user.data = NULL)
{
  checkPtrType(object, "GtkClipboard")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gtk_clipboard_request_targets", object, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkClipboardWaitForContents <-
function(object, target)
{
  checkPtrType(object, "GtkClipboard")
  target <- as.GdkAtom(target)

  w <- .RGtkCall("S_gtk_clipboard_wait_for_contents", object, target, PACKAGE = "RGtk2")

  return(w)
} 


gtkClipboardWaitForImage <-
function(object)
{
  checkPtrType(object, "GtkClipboard")

  w <- .RGtkCall("S_gtk_clipboard_wait_for_image", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkClipboardWaitForText <-
function(object)
{
  checkPtrType(object, "GtkClipboard")

  w <- .RGtkCall("S_gtk_clipboard_wait_for_text", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkClipboardWaitIsImageAvailable <-
function(object)
{
  checkPtrType(object, "GtkClipboard")

  w <- .RGtkCall("S_gtk_clipboard_wait_is_image_available", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkClipboardWaitIsTextAvailable <-
function(object)
{
  checkPtrType(object, "GtkClipboard")

  w <- .RGtkCall("S_gtk_clipboard_wait_is_text_available", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkClipboardWaitIsTargetAvailable <-
function(object, target)
{
  checkPtrType(object, "GtkClipboard")
  target <- as.GdkAtom(target)

  w <- .RGtkCall("S_gtk_clipboard_wait_is_target_available", object, target, PACKAGE = "RGtk2")

  return(w)
} 


gtkClipboardWaitForTargets <-
function(object)
{
  checkPtrType(object, "GtkClipboard")

  w <- .RGtkCall("S_gtk_clipboard_wait_for_targets", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkCListGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_clist_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkCListNew <-
function(columns = 1, show = TRUE)
{
  if(getOption("depwarn"))
    .Deprecated("GtkListStore/GtkTreeView", "RGtk2")

  columns <- as.integer(columns)

  w <- .RGtkCall("S_gtk_clist_new", columns, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkCListNewWithTitles <-
function(columns = 1, titles, show = TRUE)
{
  if(getOption("depwarn"))
    .Deprecated("GtkListStore/GtkTreeView", "RGtk2")

  columns <- as.integer(columns)
  titles <- as.list(as.character(titles))

  w <- .RGtkCall("S_gtk_clist_new_with_titles", columns, titles, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkCListSetHadjustment <-
function(object, adjustment)
{
  checkPtrType(object, "GtkCList")
  checkPtrType(adjustment, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_clist_set_hadjustment", object, adjustment, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListSetVadjustment <-
function(object, adjustment)
{
  checkPtrType(object, "GtkCList")
  checkPtrType(adjustment, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_clist_set_vadjustment", object, adjustment, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListGetHadjustment <-
function(object)
{
  checkPtrType(object, "GtkCList")

  w <- .RGtkCall("S_gtk_clist_get_hadjustment", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkCListGetVadjustment <-
function(object)
{
  checkPtrType(object, "GtkCList")

  w <- .RGtkCall("S_gtk_clist_get_vadjustment", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkCListSetShadowType <-
function(object, type)
{
  checkPtrType(object, "GtkCList")
  

  w <- .RGtkCall("S_gtk_clist_set_shadow_type", object, type, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListSetSelectionMode <-
function(object, mode)
{
  checkPtrType(object, "GtkCList")
  

  w <- .RGtkCall("S_gtk_clist_set_selection_mode", object, mode, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListSetReorderable <-
function(object, reorderable)
{
  checkPtrType(object, "GtkCList")
  reorderable <- as.logical(reorderable)

  w <- .RGtkCall("S_gtk_clist_set_reorderable", object, reorderable, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListSetUseDragIcons <-
function(object, use.icons)
{
  checkPtrType(object, "GtkCList")
  use.icons <- as.logical(use.icons)

  w <- .RGtkCall("S_gtk_clist_set_use_drag_icons", object, use.icons, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListSetButtonActions <-
function(object, button, button.actions)
{
  checkPtrType(object, "GtkCList")
  button <- as.numeric(button)
  button.actions <- as.raw(button.actions)

  w <- .RGtkCall("S_gtk_clist_set_button_actions", object, button, button.actions, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListFreeze <-
function(object)
{
  checkPtrType(object, "GtkCList")

  w <- .RGtkCall("S_gtk_clist_freeze", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListThaw <-
function(object)
{
  checkPtrType(object, "GtkCList")

  w <- .RGtkCall("S_gtk_clist_thaw", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListColumnTitlesShow <-
function(object)
{
  checkPtrType(object, "GtkCList")

  w <- .RGtkCall("S_gtk_clist_column_titles_show", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListColumnTitlesHide <-
function(object)
{
  checkPtrType(object, "GtkCList")

  w <- .RGtkCall("S_gtk_clist_column_titles_hide", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListColumnTitleActive <-
function(object, column)
{
  checkPtrType(object, "GtkCList")
  column <- as.integer(column)

  w <- .RGtkCall("S_gtk_clist_column_title_active", object, column, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListColumnTitlePassive <-
function(object, column)
{
  checkPtrType(object, "GtkCList")
  column <- as.integer(column)

  w <- .RGtkCall("S_gtk_clist_column_title_passive", object, column, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListColumnTitlesActive <-
function(object)
{
  checkPtrType(object, "GtkCList")

  w <- .RGtkCall("S_gtk_clist_column_titles_active", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListColumnTitlesPassive <-
function(object)
{
  checkPtrType(object, "GtkCList")

  w <- .RGtkCall("S_gtk_clist_column_titles_passive", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListSetColumnTitle <-
function(object, column, title)
{
  checkPtrType(object, "GtkCList")
  column <- as.integer(column)
  title <- as.character(title)

  w <- .RGtkCall("S_gtk_clist_set_column_title", object, column, title, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListGetColumnTitle <-
function(object, column)
{
  checkPtrType(object, "GtkCList")
  column <- as.integer(column)

  w <- .RGtkCall("S_gtk_clist_get_column_title", object, column, PACKAGE = "RGtk2")

  return(w)
} 


gtkCListSetColumnWidget <-
function(object, column, widget)
{
  checkPtrType(object, "GtkCList")
  column <- as.integer(column)
  checkPtrType(widget, "GtkWidget")

  w <- .RGtkCall("S_gtk_clist_set_column_widget", object, column, widget, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListGetColumnWidget <-
function(object, column)
{
  checkPtrType(object, "GtkCList")
  column <- as.integer(column)

  w <- .RGtkCall("S_gtk_clist_get_column_widget", object, column, PACKAGE = "RGtk2")

  return(w)
} 


gtkCListSetColumnJustification <-
function(object, column, justification)
{
  checkPtrType(object, "GtkCList")
  column <- as.integer(column)
  

  w <- .RGtkCall("S_gtk_clist_set_column_justification", object, column, justification, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListSetColumnVisibility <-
function(object, column, visible)
{
  checkPtrType(object, "GtkCList")
  column <- as.integer(column)
  visible <- as.logical(visible)

  w <- .RGtkCall("S_gtk_clist_set_column_visibility", object, column, visible, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListSetColumnResizeable <-
function(object, column, resizeable)
{
  checkPtrType(object, "GtkCList")
  column <- as.integer(column)
  resizeable <- as.logical(resizeable)

  w <- .RGtkCall("S_gtk_clist_set_column_resizeable", object, column, resizeable, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListSetColumnAutoResize <-
function(object, column, auto.resize)
{
  checkPtrType(object, "GtkCList")
  column <- as.integer(column)
  auto.resize <- as.logical(auto.resize)

  w <- .RGtkCall("S_gtk_clist_set_column_auto_resize", object, column, auto.resize, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListColumnsAutosize <-
function(object)
{
  checkPtrType(object, "GtkCList")

  w <- .RGtkCall("S_gtk_clist_columns_autosize", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkCListOptimalColumnWidth <-
function(object, column)
{
  checkPtrType(object, "GtkCList")
  column <- as.integer(column)

  w <- .RGtkCall("S_gtk_clist_optimal_column_width", object, column, PACKAGE = "RGtk2")

  return(w)
} 


gtkCListSetColumnWidth <-
function(object, column, width)
{
  checkPtrType(object, "GtkCList")
  column <- as.integer(column)
  width <- as.integer(width)

  w <- .RGtkCall("S_gtk_clist_set_column_width", object, column, width, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListSetColumnMinWidth <-
function(object, column, min.width)
{
  checkPtrType(object, "GtkCList")
  column <- as.integer(column)
  min.width <- as.integer(min.width)

  w <- .RGtkCall("S_gtk_clist_set_column_min_width", object, column, min.width, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListSetColumnMaxWidth <-
function(object, column, max.width)
{
  checkPtrType(object, "GtkCList")
  column <- as.integer(column)
  max.width <- as.integer(max.width)

  w <- .RGtkCall("S_gtk_clist_set_column_max_width", object, column, max.width, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListSetRowHeight <-
function(object, height)
{
  checkPtrType(object, "GtkCList")
  height <- as.numeric(height)

  w <- .RGtkCall("S_gtk_clist_set_row_height", object, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListMoveto <-
function(object, row, column, row.align, col.align)
{
  checkPtrType(object, "GtkCList")
  row <- as.integer(row)
  column <- as.integer(column)
  row.align <- as.numeric(row.align)
  col.align <- as.numeric(col.align)

  w <- .RGtkCall("S_gtk_clist_moveto", object, row, column, row.align, col.align, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListRowIsVisible <-
function(object, row)
{
  checkPtrType(object, "GtkCList")
  row <- as.integer(row)

  w <- .RGtkCall("S_gtk_clist_row_is_visible", object, row, PACKAGE = "RGtk2")

  return(w)
} 


gtkCListGetCellType <-
function(object, row, column)
{
  checkPtrType(object, "GtkCList")
  row <- as.integer(row)
  column <- as.integer(column)

  w <- .RGtkCall("S_gtk_clist_get_cell_type", object, row, column, PACKAGE = "RGtk2")

  return(w)
} 


gtkCListSetText <-
function(object, row, column, text)
{
  checkPtrType(object, "GtkCList")
  row <- as.integer(row)
  column <- as.integer(column)
  text <- as.character(text)

  w <- .RGtkCall("S_gtk_clist_set_text", object, row, column, text, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListGetText <-
function(object, row, column)
{
  checkPtrType(object, "GtkCList")
  row <- as.integer(row)
  column <- as.integer(column)

  w <- .RGtkCall("S_gtk_clist_get_text", object, row, column, PACKAGE = "RGtk2")

  return(w)
} 


gtkCListSetPixmap <-
function(object, row, column, pixmap, mask = NULL)
{
  checkPtrType(object, "GtkCList")
  row <- as.integer(row)
  column <- as.integer(column)
  checkPtrType(pixmap, "GdkPixmap")
  if (!is.null( mask )) checkPtrType(mask, "GdkBitmap")

  w <- .RGtkCall("S_gtk_clist_set_pixmap", object, row, column, pixmap, mask, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListGetPixmap <-
function(object, row, column)
{
  checkPtrType(object, "GtkCList")
  row <- as.integer(row)
  column <- as.integer(column)

  w <- .RGtkCall("S_gtk_clist_get_pixmap", object, row, column, PACKAGE = "RGtk2")

  return(w)
} 


gtkCListSetPixtext <-
function(object, row, column, text, spacing, pixmap, mask)
{
  checkPtrType(object, "GtkCList")
  row <- as.integer(row)
  column <- as.integer(column)
  text <- as.character(text)
  spacing <- as.raw(spacing)
  checkPtrType(pixmap, "GdkPixmap")
  checkPtrType(mask, "GdkBitmap")

  w <- .RGtkCall("S_gtk_clist_set_pixtext", object, row, column, text, spacing, pixmap, mask, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListGetPixtext <-
function(object, row, column)
{
  checkPtrType(object, "GtkCList")
  row <- as.integer(row)
  column <- as.integer(column)

  w <- .RGtkCall("S_gtk_clist_get_pixtext", object, row, column, PACKAGE = "RGtk2")

  return(w)
} 


gtkCListSetForeground <-
function(object, row, color)
{
  checkPtrType(object, "GtkCList")
  row <- as.integer(row)
  color <- as.GdkColor(color)

  w <- .RGtkCall("S_gtk_clist_set_foreground", object, row, color, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListSetBackground <-
function(object, row, color)
{
  checkPtrType(object, "GtkCList")
  row <- as.integer(row)
  color <- as.GdkColor(color)

  w <- .RGtkCall("S_gtk_clist_set_background", object, row, color, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListSetCellStyle <-
function(object, row, column, style)
{
  checkPtrType(object, "GtkCList")
  row <- as.integer(row)
  column <- as.integer(column)
  checkPtrType(style, "GtkStyle")

  w <- .RGtkCall("S_gtk_clist_set_cell_style", object, row, column, style, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListGetCellStyle <-
function(object, row, column)
{
  checkPtrType(object, "GtkCList")
  row <- as.integer(row)
  column <- as.integer(column)

  w <- .RGtkCall("S_gtk_clist_get_cell_style", object, row, column, PACKAGE = "RGtk2")

  return(w)
} 


gtkCListSetRowStyle <-
function(object, row, style)
{
  checkPtrType(object, "GtkCList")
  row <- as.integer(row)
  checkPtrType(style, "GtkStyle")

  w <- .RGtkCall("S_gtk_clist_set_row_style", object, row, style, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListGetRowStyle <-
function(object, row)
{
  checkPtrType(object, "GtkCList")
  row <- as.integer(row)

  w <- .RGtkCall("S_gtk_clist_get_row_style", object, row, PACKAGE = "RGtk2")

  return(w)
} 


gtkCListSetShift <-
function(object, row, column, vertical, horizontal)
{
  checkPtrType(object, "GtkCList")
  row <- as.integer(row)
  column <- as.integer(column)
  vertical <- as.integer(vertical)
  horizontal <- as.integer(horizontal)

  w <- .RGtkCall("S_gtk_clist_set_shift", object, row, column, vertical, horizontal, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListSetSelectable <-
function(object, row, selectable)
{
  checkPtrType(object, "GtkCList")
  row <- as.integer(row)
  selectable <- as.logical(selectable)

  w <- .RGtkCall("S_gtk_clist_set_selectable", object, row, selectable, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListGetSelectable <-
function(object, row)
{
  checkPtrType(object, "GtkCList")
  row <- as.integer(row)

  w <- .RGtkCall("S_gtk_clist_get_selectable", object, row, PACKAGE = "RGtk2")

  return(w)
} 


gtkCListPrepend <-
function(object, text)
{
  checkPtrType(object, "GtkCList")
  text <- as.list(as.character(text))

  w <- .RGtkCall("S_gtk_clist_prepend", object, text, PACKAGE = "RGtk2")

  return(w)
} 


gtkCListAppend <-
function(object, text)
{
  checkPtrType(object, "GtkCList")
  text <- as.list(as.character(text))

  w <- .RGtkCall("S_gtk_clist_append", object, text, PACKAGE = "RGtk2")

  return(w)
} 


gtkCListInsert <-
function(object, row, text)
{
  checkPtrType(object, "GtkCList")
  row <- as.integer(row)
  text <- as.list(as.character(text))

  w <- .RGtkCall("S_gtk_clist_insert", object, row, text, PACKAGE = "RGtk2")

  return(w)
} 


gtkCListRemove <-
function(object, row)
{
  checkPtrType(object, "GtkCList")
  row <- as.integer(row)

  w <- .RGtkCall("S_gtk_clist_remove", object, row, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListSetRowData <-
function(object, row, data = NULL)
{
  checkPtrType(object, "GtkCList")
  row <- as.integer(row)
  

  w <- .RGtkCall("S_gtk_clist_set_row_data", object, row, data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListSetRowDataFull <-
function(object, row, data = NULL)
{
  checkPtrType(object, "GtkCList")
  row <- as.integer(row)
  

  w <- .RGtkCall("S_gtk_clist_set_row_data_full", object, row, data, PACKAGE = "RGtk2")

  return(w)
} 


gtkCListGetRowData <-
function(object, row)
{
  checkPtrType(object, "GtkCList")
  row <- as.integer(row)

  w <- .RGtkCall("S_gtk_clist_get_row_data", object, row, PACKAGE = "RGtk2")

  return(w)
} 


gtkCListFindRowFromData <-
function(object, data)
{
  checkPtrType(object, "GtkCList")
  

  w <- .RGtkCall("S_gtk_clist_find_row_from_data", object, data, PACKAGE = "RGtk2")

  return(w)
} 


gtkCListSelectRow <-
function(object, row, column)
{
  checkPtrType(object, "GtkCList")
  row <- as.integer(row)
  column <- as.integer(column)

  w <- .RGtkCall("S_gtk_clist_select_row", object, row, column, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListUnselectRow <-
function(object, row, column)
{
  checkPtrType(object, "GtkCList")
  row <- as.integer(row)
  column <- as.integer(column)

  w <- .RGtkCall("S_gtk_clist_unselect_row", object, row, column, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListUndoSelection <-
function(object)
{
  checkPtrType(object, "GtkCList")

  w <- .RGtkCall("S_gtk_clist_undo_selection", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListClear <-
function(object)
{
  checkPtrType(object, "GtkCList")

  w <- .RGtkCall("S_gtk_clist_clear", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListGetSelectionInfo <-
function(object, x, y)
{
  checkPtrType(object, "GtkCList")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_gtk_clist_get_selection_info", object, x, y, PACKAGE = "RGtk2")

  return(w)
} 


gtkCListSelectAll <-
function(object)
{
  checkPtrType(object, "GtkCList")

  w <- .RGtkCall("S_gtk_clist_select_all", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListUnselectAll <-
function(object)
{
  checkPtrType(object, "GtkCList")

  w <- .RGtkCall("S_gtk_clist_unselect_all", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListSwapRows <-
function(object, row1, row2)
{
  checkPtrType(object, "GtkCList")
  row1 <- as.integer(row1)
  row2 <- as.integer(row2)

  w <- .RGtkCall("S_gtk_clist_swap_rows", object, row1, row2, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListRowMove <-
function(object, source.row, dest.row)
{
  checkPtrType(object, "GtkCList")
  source.row <- as.integer(source.row)
  dest.row <- as.integer(dest.row)

  w <- .RGtkCall("S_gtk_clist_row_move", object, source.row, dest.row, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListSetSortColumn <-
function(object, column)
{
  checkPtrType(object, "GtkCList")
  column <- as.integer(column)

  w <- .RGtkCall("S_gtk_clist_set_sort_column", object, column, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListSetSortType <-
function(object, sort.type)
{
  checkPtrType(object, "GtkCList")
  

  w <- .RGtkCall("S_gtk_clist_set_sort_type", object, sort.type, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListSort <-
function(object)
{
  checkPtrType(object, "GtkCList")

  w <- .RGtkCall("S_gtk_clist_sort", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCListSetAutoSort <-
function(object, auto.sort)
{
  checkPtrType(object, "GtkCList")
  auto.sort <- as.logical(auto.sort)

  w <- .RGtkCall("S_gtk_clist_set_auto_sort", object, auto.sort, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkColorButtonGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_color_button_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkColorButtonNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_color_button_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkColorButtonNewWithColor <-
function(color, show = TRUE)
{
  color <- as.GdkColor(color)

  w <- .RGtkCall("S_gtk_color_button_new_with_color", color, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkColorButtonSetColor <-
function(object, color)
{
  checkPtrType(object, "GtkColorButton")
  color <- as.GdkColor(color)

  w <- .RGtkCall("S_gtk_color_button_set_color", object, color, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkColorButtonSetAlpha <-
function(object, alpha)
{
  checkPtrType(object, "GtkColorButton")
  alpha <- as.integer(alpha)

  w <- .RGtkCall("S_gtk_color_button_set_alpha", object, alpha, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkColorButtonGetColor <-
function(object)
{
  checkPtrType(object, "GtkColorButton")

  w <- .RGtkCall("S_gtk_color_button_get_color", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkColorButtonGetAlpha <-
function(object)
{
  checkPtrType(object, "GtkColorButton")

  w <- .RGtkCall("S_gtk_color_button_get_alpha", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkColorButtonSetUseAlpha <-
function(object, use.alpha)
{
  checkPtrType(object, "GtkColorButton")
  use.alpha <- as.logical(use.alpha)

  w <- .RGtkCall("S_gtk_color_button_set_use_alpha", object, use.alpha, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkColorButtonGetUseAlpha <-
function(object)
{
  checkPtrType(object, "GtkColorButton")

  w <- .RGtkCall("S_gtk_color_button_get_use_alpha", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkColorButtonSetTitle <-
function(object, title)
{
  checkPtrType(object, "GtkColorButton")
  title <- as.character(title)

  w <- .RGtkCall("S_gtk_color_button_set_title", object, title, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkColorButtonGetTitle <-
function(object)
{
  checkPtrType(object, "GtkColorButton")

  w <- .RGtkCall("S_gtk_color_button_get_title", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkColorSelectionGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_color_selection_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkColorSelectionNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_color_selection_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkColorSelectionGetHasOpacityControl <-
function(object)
{
  checkPtrType(object, "GtkColorSelection")

  w <- .RGtkCall("S_gtk_color_selection_get_has_opacity_control", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkColorSelectionSetHasOpacityControl <-
function(object, has.opacity)
{
  checkPtrType(object, "GtkColorSelection")
  has.opacity <- as.logical(has.opacity)

  w <- .RGtkCall("S_gtk_color_selection_set_has_opacity_control", object, has.opacity, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkColorSelectionGetHasPalette <-
function(object)
{
  checkPtrType(object, "GtkColorSelection")

  w <- .RGtkCall("S_gtk_color_selection_get_has_palette", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkColorSelectionSetHasPalette <-
function(object, has.palette)
{
  checkPtrType(object, "GtkColorSelection")
  has.palette <- as.logical(has.palette)

  w <- .RGtkCall("S_gtk_color_selection_set_has_palette", object, has.palette, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkColorSelectionSetCurrentColor <-
function(object, color)
{
  checkPtrType(object, "GtkColorSelection")
  color <- as.GdkColor(color)

  w <- .RGtkCall("S_gtk_color_selection_set_current_color", object, color, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkColorSelectionSetCurrentAlpha <-
function(object, alpha)
{
  checkPtrType(object, "GtkColorSelection")
  alpha <- as.integer(alpha)

  w <- .RGtkCall("S_gtk_color_selection_set_current_alpha", object, alpha, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkColorSelectionGetCurrentColor <-
function(object)
{
  checkPtrType(object, "GtkColorSelection")

  w <- .RGtkCall("S_gtk_color_selection_get_current_color", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkColorSelectionGetCurrentAlpha <-
function(object)
{
  checkPtrType(object, "GtkColorSelection")

  w <- .RGtkCall("S_gtk_color_selection_get_current_alpha", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkColorSelectionSetPreviousColor <-
function(object, color)
{
  checkPtrType(object, "GtkColorSelection")
  color <- as.GdkColor(color)

  w <- .RGtkCall("S_gtk_color_selection_set_previous_color", object, color, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkColorSelectionSetPreviousAlpha <-
function(object, alpha)
{
  checkPtrType(object, "GtkColorSelection")
  alpha <- as.integer(alpha)

  w <- .RGtkCall("S_gtk_color_selection_set_previous_alpha", object, alpha, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkColorSelectionGetPreviousColor <-
function(object, color)
{
  checkPtrType(object, "GtkColorSelection")
  color <- as.GdkColor(color)

  w <- .RGtkCall("S_gtk_color_selection_get_previous_color", object, color, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkColorSelectionGetPreviousAlpha <-
function(object)
{
  checkPtrType(object, "GtkColorSelection")

  w <- .RGtkCall("S_gtk_color_selection_get_previous_alpha", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkColorSelectionIsAdjusting <-
function(object)
{
  checkPtrType(object, "GtkColorSelection")

  w <- .RGtkCall("S_gtk_color_selection_is_adjusting", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkColorSelectionPaletteFromString <-
function(str)
{
  str <- as.character(str)

  w <- .RGtkCall("S_gtk_color_selection_palette_from_string", str, PACKAGE = "RGtk2")

  return(w)
} 


gtkColorSelectionPaletteToString <-
function(colors)
{
  colors <- lapply(colors, function(x) { x <- as.GdkColor(x); x })

  w <- .RGtkCall("S_gtk_color_selection_palette_to_string", colors, PACKAGE = "RGtk2")

  return(w)
} 


gtkColorSelectionSetChangePaletteHook <-
function(func)
{
  if(getOption("depwarn"))
    .Deprecated("setChangePaletteWithScreenHook", "RGtk2")

  func <- as.function(func)

  w <- .RGtkCall("S_gtk_color_selection_set_change_palette_hook", func, PACKAGE = "RGtk2")

  return(w)
} 


gtkColorSelectionSetChangePaletteWithScreenHook <-
function(func)
{
  func <- as.function(func)

  w <- .RGtkCall("S_gtk_color_selection_set_change_palette_with_screen_hook", func, PACKAGE = "RGtk2")

  return(w)
} 


gtkColorSelectionSetColor <-
function(object, color)
{
  if(getOption("depwarn"))
    .Deprecated("gtkColorSelectionSetCurrentColor", "RGtk2")

  checkPtrType(object, "GtkColorSelection")
  color <- as.list(as.numeric(color))

  w <- .RGtkCall("S_gtk_color_selection_set_color", object, color, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkColorSelectionGetColor <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("gtkColorSelectionGetCurrentColor", "RGtk2")

  checkPtrType(object, "GtkColorSelection")

  w <- .RGtkCall("S_gtk_color_selection_get_color", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkColorSelectionSetUpdatePolicy <-
function(object, policy)
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  checkPtrType(object, "GtkColorSelection")
  

  w <- .RGtkCall("S_gtk_color_selection_set_update_policy", object, policy, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkColorSelectionDialogGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_color_selection_dialog_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkColorSelectionDialogNew <-
function(title = NULL, show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_color_selection_dialog_new", title, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkComboGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_combo_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkComboNew <-
function(show = TRUE)
{
  if(getOption("depwarn"))
    .Deprecated("GtkComboBoxEntry", "RGtk2")

  

  w <- .RGtkCall("S_gtk_combo_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkComboSetValueInList <-
function(object, val, ok.if.empty)
{
  checkPtrType(object, "GtkCombo")
  val <- as.logical(val)
  ok.if.empty <- as.logical(ok.if.empty)

  w <- .RGtkCall("S_gtk_combo_set_value_in_list", object, val, ok.if.empty, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkComboSetUseArrows <-
function(object, val)
{
  checkPtrType(object, "GtkCombo")
  val <- as.logical(val)

  w <- .RGtkCall("S_gtk_combo_set_use_arrows", object, val, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkComboSetUseArrowsAlways <-
function(object, val)
{
  checkPtrType(object, "GtkCombo")
  val <- as.logical(val)

  w <- .RGtkCall("S_gtk_combo_set_use_arrows_always", object, val, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkComboSetCaseSensitive <-
function(object, val)
{
  checkPtrType(object, "GtkCombo")
  val <- as.logical(val)

  w <- .RGtkCall("S_gtk_combo_set_case_sensitive", object, val, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkComboSetItemString <-
function(object, item, item.value)
{
  checkPtrType(object, "GtkCombo")
  checkPtrType(item, "GtkItem")
  item.value <- as.character(item.value)

  w <- .RGtkCall("S_gtk_combo_set_item_string", object, item, item.value, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkComboSetPopdownStrings <-
function(object, strings)
{
  checkPtrType(object, "GtkCombo")
  strings <- as.GList(strings)

  w <- .RGtkCall("S_gtk_combo_set_popdown_strings", object, strings, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkComboDisableActivate <-
function(object)
{
  checkPtrType(object, "GtkCombo")

  w <- .RGtkCall("S_gtk_combo_disable_activate", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkComboBoxGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_combo_box_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkComboBoxNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_combo_box_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkComboBoxNewWithModel <-
function(model, show = TRUE)
{
  checkPtrType(model, "GtkTreeModel")

  w <- .RGtkCall("S_gtk_combo_box_new_with_model", model, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkComboBoxSetWrapWidth <-
function(object, width)
{
  checkPtrType(object, "GtkComboBox")
  width <- as.integer(width)

  w <- .RGtkCall("S_gtk_combo_box_set_wrap_width", object, width, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkComboBoxSetRowSpanColumn <-
function(object, row.span)
{
  checkPtrType(object, "GtkComboBox")
  row.span <- as.integer(row.span)

  w <- .RGtkCall("S_gtk_combo_box_set_row_span_column", object, row.span, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkComboBoxSetColumnSpanColumn <-
function(object, column.span)
{
  checkPtrType(object, "GtkComboBox")
  column.span <- as.integer(column.span)

  w <- .RGtkCall("S_gtk_combo_box_set_column_span_column", object, column.span, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkComboBoxGetActive <-
function(object)
{
  checkPtrType(object, "GtkComboBox")

  w <- .RGtkCall("S_gtk_combo_box_get_active", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkComboBoxSetActive <-
function(object, index)
{
  checkPtrType(object, "GtkComboBox")
  index <- as.integer(index)

  w <- .RGtkCall("S_gtk_combo_box_set_active", object, index, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkComboBoxGetActiveIter <-
function(object)
{
  checkPtrType(object, "GtkComboBox")

  w <- .RGtkCall("S_gtk_combo_box_get_active_iter", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkComboBoxSetActiveIter <-
function(object, iter)
{
  checkPtrType(object, "GtkComboBox")
  checkPtrType(iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_combo_box_set_active_iter", object, iter, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkComboBoxSetModel <-
function(object, model = NULL)
{
  checkPtrType(object, "GtkComboBox")
  if (!is.null( model )) checkPtrType(model, "GtkTreeModel")

  w <- .RGtkCall("S_gtk_combo_box_set_model", object, model, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkComboBoxGetModel <-
function(object)
{
  checkPtrType(object, "GtkComboBox")

  w <- .RGtkCall("S_gtk_combo_box_get_model", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkComboBoxNewText <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_combo_box_new_text", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkComboBoxAppendText <-
function(object, text)
{
  checkPtrType(object, "GtkComboBox")
  text <- as.character(text)

  w <- .RGtkCall("S_gtk_combo_box_append_text", object, text, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkComboBoxInsertText <-
function(object, position, text)
{
  checkPtrType(object, "GtkComboBox")
  position <- as.integer(position)
  text <- as.character(text)

  w <- .RGtkCall("S_gtk_combo_box_insert_text", object, position, text, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkComboBoxPrependText <-
function(object, text)
{
  checkPtrType(object, "GtkComboBox")
  text <- as.character(text)

  w <- .RGtkCall("S_gtk_combo_box_prepend_text", object, text, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkComboBoxRemoveText <-
function(object, position)
{
  checkPtrType(object, "GtkComboBox")
  position <- as.integer(position)

  w <- .RGtkCall("S_gtk_combo_box_remove_text", object, position, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkComboBoxPopup <-
function(object)
{
  checkPtrType(object, "GtkComboBox")

  w <- .RGtkCall("S_gtk_combo_box_popup", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkComboBoxPopdown <-
function(object)
{
  checkPtrType(object, "GtkComboBox")

  w <- .RGtkCall("S_gtk_combo_box_popdown", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkComboBoxGetWrapWidth <-
function(object)
{
  checkPtrType(object, "GtkComboBox")

  w <- .RGtkCall("S_gtk_combo_box_get_wrap_width", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkComboBoxGetRowSpanColumn <-
function(object)
{
  checkPtrType(object, "GtkComboBox")

  w <- .RGtkCall("S_gtk_combo_box_get_row_span_column", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkComboBoxGetColumnSpanColumn <-
function(object)
{
  checkPtrType(object, "GtkComboBox")

  w <- .RGtkCall("S_gtk_combo_box_get_column_span_column", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkComboBoxSetAddTearoffs <-
function(object, add.tearoffs)
{
  checkPtrType(object, "GtkComboBox")
  add.tearoffs <- as.logical(add.tearoffs)

  w <- .RGtkCall("S_gtk_combo_box_set_add_tearoffs", object, add.tearoffs, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkComboBoxGetAddTearoffs <-
function(object)
{
  checkPtrType(object, "GtkComboBox")

  w <- .RGtkCall("S_gtk_combo_box_get_add_tearoffs", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkComboBoxGetFocusOnClick <-
function(object)
{
  checkPtrType(object, "GtkComboBox")

  w <- .RGtkCall("S_gtk_combo_box_get_focus_on_click", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkComboBoxSetFocusOnClick <-
function(object, focus.on.click)
{
  checkPtrType(object, "GtkComboBox")
  focus.on.click <- as.logical(focus.on.click)

  w <- .RGtkCall("S_gtk_combo_box_set_focus_on_click", object, focus.on.click, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkComboBoxSetRowSeparatorFunc <-
function(object, func, data = NULL)
{
  checkPtrType(object, "GtkComboBox")
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_combo_box_set_row_separator_func", object, func, data, PACKAGE = "RGtk2")

  return(w)
} 


gtkComboBoxGetRowSeparatorFunc <-
function(object)
{
  checkPtrType(object, "GtkComboBox")

  w <- .RGtkCall("S_gtk_combo_box_get_row_separator_func", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkComboBoxGetActiveText <-
function(object)
{
  checkPtrType(object, "GtkComboBox")

  w <- .RGtkCall("S_gtk_combo_box_get_active_text", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkComboBoxGetPopupAccessible <-
function(object)
{
  checkPtrType(object, "GtkComboBox")

  w <- .RGtkCall("S_gtk_combo_box_get_popup_accessible", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkComboBoxEntryGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_combo_box_entry_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkComboBoxEntryNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_combo_box_entry_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkComboBoxEntryNewWithModel <-
function(model, text.column, show = TRUE)
{
  checkPtrType(model, "GtkTreeModel")
  text.column <- as.integer(text.column)

  w <- .RGtkCall("S_gtk_combo_box_entry_new_with_model", model, text.column, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkComboBoxEntrySetTextColumn <-
function(object, text.column)
{
  checkPtrType(object, "GtkComboBoxEntry")
  text.column <- as.integer(text.column)

  w <- .RGtkCall("S_gtk_combo_box_entry_set_text_column", object, text.column, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkComboBoxEntryGetTextColumn <-
function(object)
{
  checkPtrType(object, "GtkComboBoxEntry")

  w <- .RGtkCall("S_gtk_combo_box_entry_get_text_column", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkComboBoxEntryNewText <-
function()
{
  

  w <- .RGtkCall("S_gtk_combo_box_entry_new_text", PACKAGE = "RGtk2")

  return(w)
} 


gtkContainerGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_container_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkContainerSetBorderWidth <-
function(object, border.width)
{
  checkPtrType(object, "GtkContainer")
  border.width <- as.numeric(border.width)

  w <- .RGtkCall("S_gtk_container_set_border_width", object, border.width, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkContainerGetBorderWidth <-
function(object)
{
  checkPtrType(object, "GtkContainer")

  w <- .RGtkCall("S_gtk_container_get_border_width", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkContainerAdd <-
function(object, widget)
{
  checkPtrType(object, "GtkContainer")
  checkPtrType(widget, "GtkWidget")

  w <- .RGtkCall("S_gtk_container_add", object, widget, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkContainerRemove <-
function(object, widget)
{
  checkPtrType(object, "GtkContainer")
  checkPtrType(widget, "GtkWidget")

  w <- .RGtkCall("S_gtk_container_remove", object, widget, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkContainerSetResizeMode <-
function(object, resize.mode)
{
  checkPtrType(object, "GtkContainer")
  

  w <- .RGtkCall("S_gtk_container_set_resize_mode", object, resize.mode, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkContainerGetResizeMode <-
function(object)
{
  checkPtrType(object, "GtkContainer")

  w <- .RGtkCall("S_gtk_container_get_resize_mode", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkContainerCheckResize <-
function(object)
{
  checkPtrType(object, "GtkContainer")

  w <- .RGtkCall("S_gtk_container_check_resize", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkContainerForeach <-
function(object, callback, callback.data = NULL)
{
  checkPtrType(object, "GtkContainer")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gtk_container_foreach", object, callback, callback.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkContainerForeachFull <-
function(object, callback, callback.data = NULL)
{
  checkPtrType(object, "GtkContainer")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gtk_container_foreach_full", object, callback, callback.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkContainerGetChildren <-
function(object)
{
  checkPtrType(object, "GtkContainer")

  w <- .RGtkCall("S_gtk_container_get_children", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkContainerChildren <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("gtkContainerGetChildren", "RGtk2")

  checkPtrType(object, "GtkContainer")

  w <- .RGtkCall("S_gtk_container_children", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkContainerPropagateExpose <-
function(object, child, event)
{
  checkPtrType(object, "GtkContainer")
  checkPtrType(child, "GtkWidget")
  checkPtrType(event, "GdkEventExpose")

  w <- .RGtkCall("S_gtk_container_propagate_expose", object, child, event, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkContainerSetFocusChain <-
function(object, focusable.widgets)
{
  checkPtrType(object, "GtkContainer")
  focusable.widgets <- as.GList(focusable.widgets)

  w <- .RGtkCall("S_gtk_container_set_focus_chain", object, focusable.widgets, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkContainerGetFocusChain <-
function(object)
{
  checkPtrType(object, "GtkContainer")

  w <- .RGtkCall("S_gtk_container_get_focus_chain", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkContainerUnsetFocusChain <-
function(object)
{
  checkPtrType(object, "GtkContainer")

  w <- .RGtkCall("S_gtk_container_unset_focus_chain", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkContainerSetReallocateRedraws <-
function(object, needs.redraws)
{
  checkPtrType(object, "GtkContainer")
  needs.redraws <- as.logical(needs.redraws)

  w <- .RGtkCall("S_gtk_container_set_reallocate_redraws", object, needs.redraws, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkContainerSetFocusChild <-
function(object, child)
{
  checkPtrType(object, "GtkContainer")
  checkPtrType(child, "GtkWidget")

  w <- .RGtkCall("S_gtk_container_set_focus_child", object, child, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkContainerSetFocusVadjustment <-
function(object, adjustment)
{
  checkPtrType(object, "GtkContainer")
  checkPtrType(adjustment, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_container_set_focus_vadjustment", object, adjustment, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkContainerGetFocusVadjustment <-
function(object)
{
  checkPtrType(object, "GtkContainer")

  w <- .RGtkCall("S_gtk_container_get_focus_vadjustment", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkContainerSetFocusHadjustment <-
function(object, adjustment)
{
  checkPtrType(object, "GtkContainer")
  checkPtrType(adjustment, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_container_set_focus_hadjustment", object, adjustment, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkContainerGetFocusHadjustment <-
function(object)
{
  checkPtrType(object, "GtkContainer")

  w <- .RGtkCall("S_gtk_container_get_focus_hadjustment", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkContainerResizeChildren <-
function(object)
{
  checkPtrType(object, "GtkContainer")

  w <- .RGtkCall("S_gtk_container_resize_children", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkContainerChildType <-
function(object)
{
  checkPtrType(object, "GtkContainer")

  w <- .RGtkCall("S_gtk_container_child_type", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkContainerClassInstallChildProperty <-
function(cclass, property.id, pspec)
{
  checkPtrType(cclass, "GtkContainerClass")
  property.id <- as.numeric(property.id)
  pspec <- as.GParamSpec(pspec)

  w <- .RGtkCall("S_gtk_container_class_install_child_property", cclass, property.id, pspec, PACKAGE = "RGtk2")

  return(w)
} 


gtkContainerClassFindChildProperty <-
function(cclass, property.name)
{
  checkPtrType(cclass, "GObjectClass")
  property.name <- as.character(property.name)

  w <- .RGtkCall("S_gtk_container_class_find_child_property", cclass, property.name, PACKAGE = "RGtk2")

  return(w)
} 


gtkContainerClassListChildProperties <-
function(cclass)
{
  checkPtrType(cclass, "GObjectClass")

  w <- .RGtkCall("S_gtk_container_class_list_child_properties", cclass, PACKAGE = "RGtk2")

  return(w)
} 


gtkContainerChildSetProperty <-
function(object, child, property.name, value)
{
  checkPtrType(object, "GtkContainer")
  checkPtrType(child, "GtkWidget")
  property.name <- as.character(property.name)
  

  w <- .RGtkCall("S_gtk_container_child_set_property", object, child, property.name, value, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkContainerChildGetProperty <-
function(object, child, property.name)
{
  checkPtrType(object, "GtkContainer")
  checkPtrType(child, "GtkWidget")
  property.name <- as.character(property.name)

  w <- .RGtkCall("S_gtk_container_child_get_property", object, child, property.name, PACKAGE = "RGtk2")

  return(w)
} 


gtkContainerForall <-
function(object, callback, callback.data = NULL)
{
  checkPtrType(object, "GtkContainer")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gtk_container_forall", object, callback, callback.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_ctree_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkCTreeNewWithTitles <-
function(columns = 1, tree.column = 0, titles, show = TRUE)
{
  if(getOption("depwarn"))
    .Deprecated("GtkTreeStore/GtkTreeView", "RGtk2")

  columns <- as.integer(columns)
  tree.column <- as.integer(tree.column)
  titles <- as.list(as.character(titles))

  w <- .RGtkCall("S_gtk_ctree_new_with_titles", columns, tree.column, titles, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkCTreeNew <-
function(columns = 1, tree.column = 0, show = TRUE)
{
  if(getOption("depwarn"))
    .Deprecated("GtkTreeStore/GtkTreeView", "RGtk2")

  columns <- as.integer(columns)
  tree.column <- as.integer(tree.column)

  w <- .RGtkCall("S_gtk_ctree_new", columns, tree.column, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkCTreeInsertNode <-
function(object, parent, sibling, text, spacing = 5, pixmap.closed = NULL, mask.closed = NULL, pixmap.opened = NULL, mask.opened = NULL, is.leaf = 1, expanded = 0)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(parent, "GtkCTreeNode")
  checkPtrType(sibling, "GtkCTreeNode")
  text <- as.list(as.character(text))
  spacing <- as.raw(spacing)
  if (!is.null( pixmap.closed )) checkPtrType(pixmap.closed, "GdkPixmap")
  if (!is.null( mask.closed )) checkPtrType(mask.closed, "GdkBitmap")
  if (!is.null( pixmap.opened )) checkPtrType(pixmap.opened, "GdkPixmap")
  if (!is.null( mask.opened )) checkPtrType(mask.opened, "GdkBitmap")
  is.leaf <- as.logical(is.leaf)
  expanded <- as.logical(expanded)

  w <- .RGtkCall("S_gtk_ctree_insert_node", object, parent, sibling, text, spacing, pixmap.closed, mask.closed, pixmap.opened, mask.opened, is.leaf, expanded, PACKAGE = "RGtk2")

  return(w)
} 


gtkCTreeRemoveNode <-
function(object, node)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")

  w <- .RGtkCall("S_gtk_ctree_remove_node", object, node, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeInsertGnode <-
function(object, parent, sibling, gnode, func, data = NULL)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(parent, "GtkCTreeNode")
  checkPtrType(sibling, "GtkCTreeNode")
  checkPtrType(gnode, "GNode")
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_ctree_insert_gnode", object, parent, sibling, gnode, func, data, PACKAGE = "RGtk2")

  return(w)
} 


gtkCTreeExportToGnode <-
function(object, parent, sibling, node, func, data = NULL)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(parent, "GNode")
  checkPtrType(sibling, "GNode")
  checkPtrType(node, "GtkCTreeNode")
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_ctree_export_to_gnode", object, parent, sibling, node, func, data, PACKAGE = "RGtk2")

  return(w)
} 


gtkCTreePostRecursive <-
function(object, node, func, data = NULL)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_ctree_post_recursive", object, node, func, data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreePostRecursiveToDepth <-
function(object, node, depth, func, data = NULL)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")
  depth <- as.integer(depth)
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_ctree_post_recursive_to_depth", object, node, depth, func, data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreePreRecursive <-
function(object, node, func, data = NULL)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_ctree_pre_recursive", object, node, func, data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreePreRecursiveToDepth <-
function(object, node, depth, func, data = NULL)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")
  depth <- as.integer(depth)
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_ctree_pre_recursive_to_depth", object, node, depth, func, data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeIsViewable <-
function(object, node)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")

  w <- .RGtkCall("S_gtk_ctree_is_viewable", object, node, PACKAGE = "RGtk2")

  return(w)
} 


gtkCTreeLast <-
function(object, node)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")

  w <- .RGtkCall("S_gtk_ctree_last", object, node, PACKAGE = "RGtk2")

  return(w)
} 


gtkCTreeFindNodePtr <-
function(object, ctree.row)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(ctree.row, "GtkCTreeRow")

  w <- .RGtkCall("S_gtk_ctree_find_node_ptr", object, ctree.row, PACKAGE = "RGtk2")

  return(w)
} 


gtkCTreeNodeNth <-
function(object, row)
{
  checkPtrType(object, "GtkCTree")
  row <- as.numeric(row)

  w <- .RGtkCall("S_gtk_ctree_node_nth", object, row, PACKAGE = "RGtk2")

  return(w)
} 


gtkCTreeFind <-
function(object, node, child)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")
  checkPtrType(child, "GtkCTreeNode")

  w <- .RGtkCall("S_gtk_ctree_find", object, node, child, PACKAGE = "RGtk2")

  return(w)
} 


gtkCTreeIsAncestor <-
function(object, node, child)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")
  checkPtrType(child, "GtkCTreeNode")

  w <- .RGtkCall("S_gtk_ctree_is_ancestor", object, node, child, PACKAGE = "RGtk2")

  return(w)
} 


gtkCTreeFindByRowData <-
function(object, node, data = NULL)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")
  

  w <- .RGtkCall("S_gtk_ctree_find_by_row_data", object, node, data, PACKAGE = "RGtk2")

  return(w)
} 


gtkCTreeFindAllByRowData <-
function(object, node, data = NULL)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")
  

  w <- .RGtkCall("S_gtk_ctree_find_all_by_row_data", object, node, data, PACKAGE = "RGtk2")

  return(w)
} 


gtkCTreeFindByRowDataCustom <-
function(object, node, data = NULL, func)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")
  
  func <- as.function(func)

  w <- .RGtkCall("S_gtk_ctree_find_by_row_data_custom", object, node, data, func, PACKAGE = "RGtk2")

  return(w)
} 


gtkCTreeFindAllByRowDataCustom <-
function(object, node, data = NULL, func)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")
  
  func <- as.function(func)

  w <- .RGtkCall("S_gtk_ctree_find_all_by_row_data_custom", object, node, data, func, PACKAGE = "RGtk2")

  return(w)
} 


gtkCTreeIsHotSpot <-
function(object, x, y)
{
  checkPtrType(object, "GtkCTree")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_gtk_ctree_is_hot_spot", object, x, y, PACKAGE = "RGtk2")

  return(w)
} 


gtkCTreeMove <-
function(object, node, new.parent = NULL, new.sibling = NULL)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")
  if (!is.null( new.parent )) checkPtrType(new.parent, "GtkCTreeNode")
  if (!is.null( new.sibling )) checkPtrType(new.sibling, "GtkCTreeNode")

  w <- .RGtkCall("S_gtk_ctree_move", object, node, new.parent, new.sibling, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeExpand <-
function(object, node)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")

  w <- .RGtkCall("S_gtk_ctree_expand", object, node, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeExpandRecursive <-
function(object, node)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")

  w <- .RGtkCall("S_gtk_ctree_expand_recursive", object, node, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeExpandToDepth <-
function(object, node, depth)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")
  depth <- as.integer(depth)

  w <- .RGtkCall("S_gtk_ctree_expand_to_depth", object, node, depth, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeCollapse <-
function(object, node)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")

  w <- .RGtkCall("S_gtk_ctree_collapse", object, node, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeCollapseRecursive <-
function(object, node)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")

  w <- .RGtkCall("S_gtk_ctree_collapse_recursive", object, node, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeCollapseToDepth <-
function(object, node, depth)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")
  depth <- as.integer(depth)

  w <- .RGtkCall("S_gtk_ctree_collapse_to_depth", object, node, depth, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeToggleExpansion <-
function(object, node)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")

  w <- .RGtkCall("S_gtk_ctree_toggle_expansion", object, node, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeToggleExpansionRecursive <-
function(object, node)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")

  w <- .RGtkCall("S_gtk_ctree_toggle_expansion_recursive", object, node, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeSelect <-
function(object, node)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")

  w <- .RGtkCall("S_gtk_ctree_select", object, node, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeSelectRecursive <-
function(object, node)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")

  w <- .RGtkCall("S_gtk_ctree_select_recursive", object, node, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeUnselect <-
function(object, node)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")

  w <- .RGtkCall("S_gtk_ctree_unselect", object, node, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeUnselectRecursive <-
function(object, node)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")

  w <- .RGtkCall("S_gtk_ctree_unselect_recursive", object, node, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeRealSelectRecursive <-
function(object, node, state)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")
  state <- as.integer(state)

  w <- .RGtkCall("S_gtk_ctree_real_select_recursive", object, node, state, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeNodeSetText <-
function(object, node, column, text)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")
  column <- as.integer(column)
  text <- as.character(text)

  w <- .RGtkCall("S_gtk_ctree_node_set_text", object, node, column, text, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeNodeSetPixmap <-
function(object, node, column, pixmap, mask = NULL)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")
  column <- as.integer(column)
  checkPtrType(pixmap, "GdkPixmap")
  if (!is.null( mask )) checkPtrType(mask, "GdkBitmap")

  w <- .RGtkCall("S_gtk_ctree_node_set_pixmap", object, node, column, pixmap, mask, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeNodeSetPixtext <-
function(object, node, column, text, spacing, pixmap, mask = NULL)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")
  column <- as.integer(column)
  text <- as.character(text)
  spacing <- as.raw(spacing)
  checkPtrType(pixmap, "GdkPixmap")
  if (!is.null( mask )) checkPtrType(mask, "GdkBitmap")

  w <- .RGtkCall("S_gtk_ctree_node_set_pixtext", object, node, column, text, spacing, pixmap, mask, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeSetNodeInfo <-
function(object, node, text, spacing, pixmap.closed = NULL, mask.closed = NULL, pixmap.opened = NULL, mask.opened = NULL, is.leaf, expanded)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")
  text <- as.character(text)
  spacing <- as.raw(spacing)
  if (!is.null( pixmap.closed )) checkPtrType(pixmap.closed, "GdkPixmap")
  if (!is.null( mask.closed )) checkPtrType(mask.closed, "GdkBitmap")
  if (!is.null( pixmap.opened )) checkPtrType(pixmap.opened, "GdkPixmap")
  if (!is.null( mask.opened )) checkPtrType(mask.opened, "GdkBitmap")
  is.leaf <- as.logical(is.leaf)
  expanded <- as.logical(expanded)

  w <- .RGtkCall("S_gtk_ctree_set_node_info", object, node, text, spacing, pixmap.closed, mask.closed, pixmap.opened, mask.opened, is.leaf, expanded, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeNodeSetShift <-
function(object, node, column, vertical, horizontal)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")
  column <- as.integer(column)
  vertical <- as.integer(vertical)
  horizontal <- as.integer(horizontal)

  w <- .RGtkCall("S_gtk_ctree_node_set_shift", object, node, column, vertical, horizontal, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeNodeSetSelectable <-
function(object, node, selectable)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")
  selectable <- as.logical(selectable)

  w <- .RGtkCall("S_gtk_ctree_node_set_selectable", object, node, selectable, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeNodeGetSelectable <-
function(object, node)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")

  w <- .RGtkCall("S_gtk_ctree_node_get_selectable", object, node, PACKAGE = "RGtk2")

  return(w)
} 


gtkCTreeNodeGetCellType <-
function(object, node, column)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")
  column <- as.integer(column)

  w <- .RGtkCall("S_gtk_ctree_node_get_cell_type", object, node, column, PACKAGE = "RGtk2")

  return(w)
} 


gtkCTreeNodeGetText <-
function(object, node, column)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")
  column <- as.integer(column)

  w <- .RGtkCall("S_gtk_ctree_node_get_text", object, node, column, PACKAGE = "RGtk2")

  return(w)
} 


gtkCTreeNodeGetPixmap <-
function(object, node, column)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")
  column <- as.integer(column)

  w <- .RGtkCall("S_gtk_ctree_node_get_pixmap", object, node, column, PACKAGE = "RGtk2")

  return(w)
} 


gtkCTreeNodeGetPixtext <-
function(object, node, column)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")
  column <- as.integer(column)

  w <- .RGtkCall("S_gtk_ctree_node_get_pixtext", object, node, column, PACKAGE = "RGtk2")

  return(w)
} 


gtkCTreeGetNodeInfo <-
function(object, node)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")

  w <- .RGtkCall("S_gtk_ctree_get_node_info", object, node, PACKAGE = "RGtk2")

  return(w)
} 


gtkCTreeNodeSetRowStyle <-
function(object, node, style)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")
  checkPtrType(style, "GtkStyle")

  w <- .RGtkCall("S_gtk_ctree_node_set_row_style", object, node, style, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeNodeGetRowStyle <-
function(object, node)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")

  w <- .RGtkCall("S_gtk_ctree_node_get_row_style", object, node, PACKAGE = "RGtk2")

  return(w)
} 


gtkCTreeNodeSetCellStyle <-
function(object, node, column, style)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")
  column <- as.integer(column)
  checkPtrType(style, "GtkStyle")

  w <- .RGtkCall("S_gtk_ctree_node_set_cell_style", object, node, column, style, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeNodeGetCellStyle <-
function(object, node, column)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")
  column <- as.integer(column)

  w <- .RGtkCall("S_gtk_ctree_node_get_cell_style", object, node, column, PACKAGE = "RGtk2")

  return(w)
} 


gtkCTreeNodeSetForeground <-
function(object, node, color)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")
  color <- as.GdkColor(color)

  w <- .RGtkCall("S_gtk_ctree_node_set_foreground", object, node, color, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeNodeSetBackground <-
function(object, node, color)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")
  color <- as.GdkColor(color)

  w <- .RGtkCall("S_gtk_ctree_node_set_background", object, node, color, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeNodeSetRowData <-
function(object, node, data)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")
  

  w <- .RGtkCall("S_gtk_ctree_node_set_row_data", object, node, data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeNodeSetRowDataFull <-
function(object, node, data)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")
  

  w <- .RGtkCall("S_gtk_ctree_node_set_row_data_full", object, node, data, PACKAGE = "RGtk2")

  return(w)
} 


gtkCTreeNodeGetRowData <-
function(object, node)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")

  w <- .RGtkCall("S_gtk_ctree_node_get_row_data", object, node, PACKAGE = "RGtk2")

  return(w)
} 


gtkCTreeNodeMoveto <-
function(object, node, column, row.align, col.align)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")
  column <- as.integer(column)
  row.align <- as.numeric(row.align)
  col.align <- as.numeric(col.align)

  w <- .RGtkCall("S_gtk_ctree_node_moveto", object, node, column, row.align, col.align, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeNodeIsVisible <-
function(object, node)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")

  w <- .RGtkCall("S_gtk_ctree_node_is_visible", object, node, PACKAGE = "RGtk2")

  return(w)
} 


gtkCTreeSetIndent <-
function(object, indent)
{
  checkPtrType(object, "GtkCTree")
  indent <- as.integer(indent)

  w <- .RGtkCall("S_gtk_ctree_set_indent", object, indent, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeSetSpacing <-
function(object, spacing)
{
  checkPtrType(object, "GtkCTree")
  spacing <- as.integer(spacing)

  w <- .RGtkCall("S_gtk_ctree_set_spacing", object, spacing, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeSetShowStub <-
function(object, show.stub)
{
  checkPtrType(object, "GtkCTree")
  show.stub <- as.logical(show.stub)

  w <- .RGtkCall("S_gtk_ctree_set_show_stub", object, show.stub, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeSetLineStyle <-
function(object, line.style)
{
  checkPtrType(object, "GtkCTree")
  

  w <- .RGtkCall("S_gtk_ctree_set_line_style", object, line.style, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeSetExpanderStyle <-
function(object, expander.style)
{
  checkPtrType(object, "GtkCTree")
  

  w <- .RGtkCall("S_gtk_ctree_set_expander_style", object, expander.style, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeSortNode <-
function(object, node)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")

  w <- .RGtkCall("S_gtk_ctree_sort_node", object, node, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeSortRecursive <-
function(object, node)
{
  checkPtrType(object, "GtkCTree")
  checkPtrType(node, "GtkCTreeNode")

  w <- .RGtkCall("S_gtk_ctree_sort_recursive", object, node, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCTreeNodeGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_ctree_node_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkCurveGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_curve_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkCurveNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_curve_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkCurveReset <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  checkPtrType(object, "GtkCurve")

  w <- .RGtkCall("S_gtk_curve_reset", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCurveSetGamma <-
function(object, gamma)
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  checkPtrType(object, "GtkCurve")
  gamma <- as.numeric(gamma)

  w <- .RGtkCall("S_gtk_curve_set_gamma", object, gamma, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCurveSetRange <-
function(object, min.x, max.x, min.y, max.y)
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  checkPtrType(object, "GtkCurve")
  min.x <- as.numeric(min.x)
  max.x <- as.numeric(max.x)
  min.y <- as.numeric(min.y)
  max.y <- as.numeric(max.y)

  w <- .RGtkCall("S_gtk_curve_set_range", object, min.x, max.x, min.y, max.y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCurveGetVector <-
function(object, veclen)
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  checkPtrType(object, "GtkCurve")
  veclen <- as.integer(veclen)

  w <- .RGtkCall("S_gtk_curve_get_vector", object, veclen, PACKAGE = "RGtk2")

  return(w)
} 


gtkCurveSetVector <-
function(object, vector)
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  checkPtrType(object, "GtkCurve")
  vector <- as.list(as.numeric(vector))

  w <- .RGtkCall("S_gtk_curve_set_vector", object, vector, PACKAGE = "RGtk2")

  return(w)
} 


gtkCurveSetCurveType <-
function(object, type)
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  checkPtrType(object, "GtkCurve")
  

  w <- .RGtkCall("S_gtk_curve_set_curve_type", object, type, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDialogGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_dialog_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkDialogNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_dialog_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkDialogAddActionWidget <-
function(object, child, response.id)
{
  checkPtrType(object, "GtkDialog")
  checkPtrType(child, "GtkWidget")
  response.id <- as.integer(response.id)

  w <- .RGtkCall("S_gtk_dialog_add_action_widget", object, child, response.id, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDialogAddButton <-
function(object, button.text, response.id)
{
  checkPtrType(object, "GtkDialog")
  button.text <- as.character(button.text)
  response.id <- as.integer(response.id)

  w <- .RGtkCall("S_gtk_dialog_add_button", object, button.text, response.id, PACKAGE = "RGtk2")

  return(w)
} 


gtkDialogSetResponseSensitive <-
function(object, response.id, setting)
{
  checkPtrType(object, "GtkDialog")
  response.id <- as.integer(response.id)
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_dialog_set_response_sensitive", object, response.id, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDialogSetDefaultResponse <-
function(object, response.id)
{
  checkPtrType(object, "GtkDialog")
  response.id <- as.integer(response.id)

  w <- .RGtkCall("S_gtk_dialog_set_default_response", object, response.id, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDialogGetResponseForWidget <-
function(object, widget)
{
  checkPtrType(object, "GtkDialog")
  checkPtrType(widget, "GtkWidget")

  w <- .RGtkCall("S_gtk_dialog_get_response_for_widget", object, widget, PACKAGE = "RGtk2")

  return(w)
} 


gtkDialogSetHasSeparator <-
function(object, setting)
{
  checkPtrType(object, "GtkDialog")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_dialog_set_has_separator", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDialogGetHasSeparator <-
function(object)
{
  checkPtrType(object, "GtkDialog")

  w <- .RGtkCall("S_gtk_dialog_get_has_separator", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkDialogResponse <-
function(object, response.id)
{
  checkPtrType(object, "GtkDialog")
  response.id <- as.integer(response.id)

  w <- .RGtkCall("S_gtk_dialog_response", object, response.id, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDialogRun <-
function(object)
{
  checkPtrType(object, "GtkDialog")

  w <- .RGtkCall("S_gtk_dialog_run", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkDialogSetAlternativeButtonOrderFromArray <-
function(object, new.order)
{
  checkPtrType(object, "GtkDialog")
  new.order <- as.list(as.integer(new.order))

  w <- .RGtkCall("S_gtk_dialog_set_alternative_button_order_from_array", object, new.order, PACKAGE = "RGtk2")

  return(w)
} 


gtkDragCheckThreshold <-
function(object, start.x, start.y, current.x, current.y)
{
  checkPtrType(object, "GtkWidget")
  start.x <- as.integer(start.x)
  start.y <- as.integer(start.y)
  current.x <- as.integer(current.x)
  current.y <- as.integer(current.y)

  w <- .RGtkCall("S_gtk_drag_check_threshold", object, start.x, start.y, current.x, current.y, PACKAGE = "RGtk2")

  return(w)
} 


gtkDragGetData <-
function(object, context, target, time = "GDK_CURRENT_TIME")
{
  checkPtrType(object, "GtkWidget")
  checkPtrType(context, "GdkDragContext")
  target <- as.GdkAtom(target)
  time <- as.numeric(time)

  w <- .RGtkCall("S_gtk_drag_get_data", object, context, target, time, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDragHighlight <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_drag_highlight", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDragUnhighlight <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_drag_unhighlight", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDragDestSet <-
function(object, flags, targets, actions)
{
  checkPtrType(object, "GtkWidget")
  
  targets <- lapply(targets, function(x) { x <- as.GtkTargetEntry(x); x })
  

  w <- .RGtkCall("S_gtk_drag_dest_set", object, flags, targets, actions, PACKAGE = "RGtk2")

  return(w)
} 


gtkDragDestSetProxy <-
function(object, proxy.window, protocol, use.coordinates)
{
  checkPtrType(object, "GtkWidget")
  checkPtrType(proxy.window, "GdkWindow")
  
  use.coordinates <- as.logical(use.coordinates)

  w <- .RGtkCall("S_gtk_drag_dest_set_proxy", object, proxy.window, protocol, use.coordinates, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDragDestUnset <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_drag_dest_unset", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDragDestFindTarget <-
function(object, context, target.list)
{
  checkPtrType(object, "GtkWidget")
  checkPtrType(context, "GdkDragContext")
  checkPtrType(target.list, "GtkTargetList")

  w <- .RGtkCall("S_gtk_drag_dest_find_target", object, context, target.list, PACKAGE = "RGtk2")

  return(w)
} 


gtkDragDestGetTargetList <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_drag_dest_get_target_list", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkDragDestSetTargetList <-
function(object, target.list)
{
  checkPtrType(object, "GtkWidget")
  checkPtrType(target.list, "GtkTargetList")

  w <- .RGtkCall("S_gtk_drag_dest_set_target_list", object, target.list, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDragSourceSet <-
function(object, start.button.mask, targets, actions)
{
  checkPtrType(object, "GtkWidget")
  
  targets <- lapply(targets, function(x) { x <- as.GtkTargetEntry(x); x })
  

  w <- .RGtkCall("S_gtk_drag_source_set", object, start.button.mask, targets, actions, PACKAGE = "RGtk2")

  return(w)
} 


gtkDragSourceUnset <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_drag_source_unset", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDragSourceSetIcon <-
function(object, colormap, pixmap, mask = NULL)
{
  checkPtrType(object, "GtkWidget")
  checkPtrType(colormap, "GdkColormap")
  checkPtrType(pixmap, "GdkPixmap")
  if (!is.null( mask )) checkPtrType(mask, "GdkBitmap")

  w <- .RGtkCall("S_gtk_drag_source_set_icon", object, colormap, pixmap, mask, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDragSourceSetIconPixbuf <-
function(object, pixbuf)
{
  checkPtrType(object, "GtkWidget")
  checkPtrType(pixbuf, "GdkPixbuf")

  w <- .RGtkCall("S_gtk_drag_source_set_icon_pixbuf", object, pixbuf, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDragSourceSetIconStock <-
function(object, stock.id)
{
  checkPtrType(object, "GtkWidget")
  stock.id <- as.character(stock.id)

  w <- .RGtkCall("S_gtk_drag_source_set_icon_stock", object, stock.id, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDragSourceGetTargetList <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_drag_source_get_target_list", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkDragSourceSetTargetList <-
function(object, target.list)
{
  checkPtrType(object, "GtkWidget")
  checkPtrType(target.list, "GtkTargetList")

  w <- .RGtkCall("S_gtk_drag_source_set_target_list", object, target.list, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDragBegin <-
function(object, targets, actions, button, event)
{
  checkPtrType(object, "GtkWidget")
  checkPtrType(targets, "GtkTargetList")
  
  button <- as.integer(button)
  checkPtrType(event, "GdkEvent")

  w <- .RGtkCall("S_gtk_drag_begin", object, targets, actions, button, event, PACKAGE = "RGtk2")

  return(w)
} 


gtkDragSetDefaultIcon <-
function(colormap, pixmap, mask, hot.x, hot.y)
{
  if(getOption("depwarn"))
    .Deprecated("a different stock pixbuf for GTK_STOCK_DND", "RGtk2")

  checkPtrType(colormap, "GdkColormap")
  checkPtrType(pixmap, "GdkPixmap")
  checkPtrType(mask, "GdkBitmap")
  hot.x <- as.integer(hot.x)
  hot.y <- as.integer(hot.y)

  w <- .RGtkCall("S_gtk_drag_set_default_icon", colormap, pixmap, mask, hot.x, hot.y, PACKAGE = "RGtk2")

  return(w)
} 


gtkDragDestAddTextTargets <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_drag_dest_add_text_targets", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDragDestAddImageTargets <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_drag_dest_add_image_targets", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDragDestAddUriTargets <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_drag_dest_add_uri_targets", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDragSourceAddTextTargets <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_drag_source_add_text_targets", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDragSourceAddImageTargets <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_drag_source_add_image_targets", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDragSourceAddUriTargets <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_drag_source_add_uri_targets", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTargetListAddTextTargets <-
function(list, info)
{
  checkPtrType(list, "GtkTargetList")
  info <- as.numeric(info)

  w <- .RGtkCall("S_gtk_target_list_add_text_targets", list, info, PACKAGE = "RGtk2")

  return(w)
} 


gtkTargetListAddImageTargets <-
function(list, info, writable)
{
  checkPtrType(list, "GtkTargetList")
  info <- as.numeric(info)
  writable <- as.logical(writable)

  w <- .RGtkCall("S_gtk_target_list_add_image_targets", list, info, writable, PACKAGE = "RGtk2")

  return(w)
} 


gtkTargetListAddUriTargets <-
function(list, info)
{
  checkPtrType(list, "GtkTargetList")
  info <- as.numeric(info)

  w <- .RGtkCall("S_gtk_target_list_add_uri_targets", list, info, PACKAGE = "RGtk2")

  return(w)
} 


gtkDragGetSourceWidget <-
function(context)
{
  checkPtrType(context, "GdkDragContext")

  w <- .RGtkCall("S_gtk_drag_get_source_widget", context, PACKAGE = "RGtk2")

  return(w)
} 


gtkDragSourceSetIconName <-
function(widget, icon.name)
{
  checkPtrType(widget, "GtkWidget")
  icon.name <- as.character(icon.name)

  w <- .RGtkCall("S_gtk_drag_source_set_icon_name", widget, icon.name, PACKAGE = "RGtk2")

  return(w)
} 


gtkDrawingAreaGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_drawing_area_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkDrawingAreaNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_drawing_area_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkDrawingAreaSize <-
function(object, width, height)
{
  if(getOption("depwarn"))
    .Deprecated("gtkWidgetSetSizeRequest", "RGtk2")

  checkPtrType(object, "GtkDrawingArea")
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_drawing_area_size", object, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEditableGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_editable_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkEditableSelectRegion <-
function(object, start, end)
{
  checkPtrType(object, "GtkEditable")
  start <- as.integer(start)
  end <- as.integer(end)

  w <- .RGtkCall("S_gtk_editable_select_region", object, start, end, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEditableGetSelectionBounds <-
function(object)
{
  checkPtrType(object, "GtkEditable")

  w <- .RGtkCall("S_gtk_editable_get_selection_bounds", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEditableInsertText <-
function(object, new.text, position = 0)
{
  checkPtrType(object, "GtkEditable")
  new.text <- as.character(new.text)
  position <- as.list(as.integer(position))

  w <- .RGtkCall("S_gtk_editable_insert_text", object, new.text, position, PACKAGE = "RGtk2")

  return(w)
} 


gtkEditableDeleteText <-
function(object, start.pos, end.pos)
{
  checkPtrType(object, "GtkEditable")
  start.pos <- as.integer(start.pos)
  end.pos <- as.integer(end.pos)

  w <- .RGtkCall("S_gtk_editable_delete_text", object, start.pos, end.pos, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEditableGetChars <-
function(object, start.pos, end.pos)
{
  checkPtrType(object, "GtkEditable")
  start.pos <- as.integer(start.pos)
  end.pos <- as.integer(end.pos)

  w <- .RGtkCall("S_gtk_editable_get_chars", object, start.pos, end.pos, PACKAGE = "RGtk2")

  return(w)
} 


gtkEditableCutClipboard <-
function(object)
{
  checkPtrType(object, "GtkEditable")

  w <- .RGtkCall("S_gtk_editable_cut_clipboard", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEditableCopyClipboard <-
function(object)
{
  checkPtrType(object, "GtkEditable")

  w <- .RGtkCall("S_gtk_editable_copy_clipboard", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEditablePasteClipboard <-
function(object)
{
  checkPtrType(object, "GtkEditable")

  w <- .RGtkCall("S_gtk_editable_paste_clipboard", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEditableDeleteSelection <-
function(object)
{
  checkPtrType(object, "GtkEditable")

  w <- .RGtkCall("S_gtk_editable_delete_selection", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEditableSetPosition <-
function(object, position)
{
  checkPtrType(object, "GtkEditable")
  position <- as.integer(position)

  w <- .RGtkCall("S_gtk_editable_set_position", object, position, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEditableGetPosition <-
function(object)
{
  checkPtrType(object, "GtkEditable")

  w <- .RGtkCall("S_gtk_editable_get_position", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEditableSetEditable <-
function(object, is.editable)
{
  checkPtrType(object, "GtkEditable")
  is.editable <- as.logical(is.editable)

  w <- .RGtkCall("S_gtk_editable_set_editable", object, is.editable, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEditableGetEditable <-
function(object)
{
  checkPtrType(object, "GtkEditable")

  w <- .RGtkCall("S_gtk_editable_get_editable", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_entry_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_entry_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkEntryNewWithMaxLength <-
function(max = 0, show = TRUE)
{
  if(getOption("depwarn"))
    .Deprecated("gtkEntryNew", "RGtk2")

  max <- as.integer(max)

  w <- .RGtkCall("S_gtk_entry_new_with_max_length", max, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkEntrySetVisibility <-
function(object, visible)
{
  checkPtrType(object, "GtkEntry")
  visible <- as.logical(visible)

  w <- .RGtkCall("S_gtk_entry_set_visibility", object, visible, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryGetVisibility <-
function(object)
{
  checkPtrType(object, "GtkEntry")

  w <- .RGtkCall("S_gtk_entry_get_visibility", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntrySetInvisibleChar <-
function(object, ch)
{
  checkPtrType(object, "GtkEntry")
  ch <- as.numeric(ch)

  w <- .RGtkCall("S_gtk_entry_set_invisible_char", object, ch, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryGetInvisibleChar <-
function(object)
{
  checkPtrType(object, "GtkEntry")

  w <- .RGtkCall("S_gtk_entry_get_invisible_char", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntrySetHasFrame <-
function(object, setting)
{
  checkPtrType(object, "GtkEntry")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_entry_set_has_frame", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryGetHasFrame <-
function(object)
{
  checkPtrType(object, "GtkEntry")

  w <- .RGtkCall("S_gtk_entry_get_has_frame", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntrySetMaxLength <-
function(object, max)
{
  checkPtrType(object, "GtkEntry")
  max <- as.integer(max)

  w <- .RGtkCall("S_gtk_entry_set_max_length", object, max, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryGetMaxLength <-
function(object)
{
  checkPtrType(object, "GtkEntry")

  w <- .RGtkCall("S_gtk_entry_get_max_length", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntrySetActivatesDefault <-
function(object, setting)
{
  checkPtrType(object, "GtkEntry")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_entry_set_activates_default", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryGetActivatesDefault <-
function(object)
{
  checkPtrType(object, "GtkEntry")

  w <- .RGtkCall("S_gtk_entry_get_activates_default", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntrySetWidthChars <-
function(object, n.chars)
{
  checkPtrType(object, "GtkEntry")
  n.chars <- as.integer(n.chars)

  w <- .RGtkCall("S_gtk_entry_set_width_chars", object, n.chars, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryGetWidthChars <-
function(object)
{
  checkPtrType(object, "GtkEntry")

  w <- .RGtkCall("S_gtk_entry_get_width_chars", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntrySetText <-
function(object, text)
{
  checkPtrType(object, "GtkEntry")
  text <- as.character(text)

  w <- .RGtkCall("S_gtk_entry_set_text", object, text, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryGetText <-
function(object)
{
  checkPtrType(object, "GtkEntry")

  w <- .RGtkCall("S_gtk_entry_get_text", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryGetLayout <-
function(object)
{
  checkPtrType(object, "GtkEntry")

  w <- .RGtkCall("S_gtk_entry_get_layout", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryGetLayoutOffsets <-
function(object)
{
  checkPtrType(object, "GtkEntry")

  w <- .RGtkCall("S_gtk_entry_get_layout_offsets", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryLayoutIndexToTextIndex <-
function(object, layout.index)
{
  checkPtrType(object, "GtkEntry")
  layout.index <- as.integer(layout.index)

  w <- .RGtkCall("S_gtk_entry_layout_index_to_text_index", object, layout.index, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryTextIndexToLayoutIndex <-
function(object, text.index)
{
  checkPtrType(object, "GtkEntry")
  text.index <- as.integer(text.index)

  w <- .RGtkCall("S_gtk_entry_text_index_to_layout_index", object, text.index, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntrySetAlignment <-
function(object, xalign)
{
  checkPtrType(object, "GtkEntry")
  xalign <- as.numeric(xalign)

  w <- .RGtkCall("S_gtk_entry_set_alignment", object, xalign, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryGetAlignment <-
function(object)
{
  checkPtrType(object, "GtkEntry")

  w <- .RGtkCall("S_gtk_entry_get_alignment", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntrySetCompletion <-
function(object, completion)
{
  checkPtrType(object, "GtkEntry")
  checkPtrType(completion, "GtkEntryCompletion")

  w <- .RGtkCall("S_gtk_entry_set_completion", object, completion, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryGetCompletion <-
function(object)
{
  checkPtrType(object, "GtkEntry")

  w <- .RGtkCall("S_gtk_entry_get_completion", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryAppendText <-
function(object, text)
{
  if(getOption("depwarn"))
    .Deprecated("gtkEditableInsertText", "RGtk2")

  checkPtrType(object, "GtkEntry")
  text <- as.character(text)

  w <- .RGtkCall("S_gtk_entry_append_text", object, text, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryPrependText <-
function(object, text)
{
  if(getOption("depwarn"))
    .Deprecated("gtkEditableInsertText", "RGtk2")

  checkPtrType(object, "GtkEntry")
  text <- as.character(text)

  w <- .RGtkCall("S_gtk_entry_prepend_text", object, text, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntrySetPosition <-
function(object, position)
{
  checkPtrType(object, "GtkEntry")
  position <- as.integer(position)

  w <- .RGtkCall("S_gtk_entry_set_position", object, position, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntrySelectRegion <-
function(object, start, end)
{
  checkPtrType(object, "GtkEntry")
  start <- as.integer(start)
  end <- as.integer(end)

  w <- .RGtkCall("S_gtk_entry_select_region", object, start, end, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntrySetEditable <-
function(object, editable)
{
  checkPtrType(object, "GtkEntry")
  editable <- as.logical(editable)

  w <- .RGtkCall("S_gtk_entry_set_editable", object, editable, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryCompletionGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_entry_completion_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryCompletionNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_entry_completion_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryCompletionGetEntry <-
function(object)
{
  checkPtrType(object, "GtkEntryCompletion")

  w <- .RGtkCall("S_gtk_entry_completion_get_entry", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryCompletionSetModel <-
function(object, model = NULL)
{
  checkPtrType(object, "GtkEntryCompletion")
  if (!is.null( model )) checkPtrType(model, "GtkTreeModel")

  w <- .RGtkCall("S_gtk_entry_completion_set_model", object, model, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryCompletionGetModel <-
function(object)
{
  checkPtrType(object, "GtkEntryCompletion")

  w <- .RGtkCall("S_gtk_entry_completion_get_model", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryCompletionSetMatchFunc <-
function(object, func, func.data = NULL)
{
  checkPtrType(object, "GtkEntryCompletion")
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_entry_completion_set_match_func", object, func, func.data, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryCompletionSetMinimumKeyLength <-
function(object, length)
{
  checkPtrType(object, "GtkEntryCompletion")
  length <- as.integer(length)

  w <- .RGtkCall("S_gtk_entry_completion_set_minimum_key_length", object, length, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryCompletionGetMinimumKeyLength <-
function(object)
{
  checkPtrType(object, "GtkEntryCompletion")

  w <- .RGtkCall("S_gtk_entry_completion_get_minimum_key_length", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryCompletionComplete <-
function(object)
{
  checkPtrType(object, "GtkEntryCompletion")

  w <- .RGtkCall("S_gtk_entry_completion_complete", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryCompletionInsertActionText <-
function(object, index, text)
{
  checkPtrType(object, "GtkEntryCompletion")
  index <- as.integer(index)
  text <- as.character(text)

  w <- .RGtkCall("S_gtk_entry_completion_insert_action_text", object, index, text, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryCompletionInsertActionMarkup <-
function(object, index, markup)
{
  checkPtrType(object, "GtkEntryCompletion")
  index <- as.integer(index)
  markup <- as.character(markup)

  w <- .RGtkCall("S_gtk_entry_completion_insert_action_markup", object, index, markup, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryCompletionDeleteAction <-
function(object, index)
{
  checkPtrType(object, "GtkEntryCompletion")
  index <- as.integer(index)

  w <- .RGtkCall("S_gtk_entry_completion_delete_action", object, index, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryCompletionSetTextColumn <-
function(object, column)
{
  checkPtrType(object, "GtkEntryCompletion")
  column <- as.integer(column)

  w <- .RGtkCall("S_gtk_entry_completion_set_text_column", object, column, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryCompletionGetTextColumn <-
function(object)
{
  checkPtrType(object, "GtkEntryCompletion")

  w <- .RGtkCall("S_gtk_entry_completion_get_text_column", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryCompletionInsertPrefix <-
function(object)
{
  checkPtrType(object, "GtkEntryCompletion")

  w <- .RGtkCall("S_gtk_entry_completion_insert_prefix", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryCompletionSetInlineCompletion <-
function(object, inline.completion)
{
  checkPtrType(object, "GtkEntryCompletion")
  inline.completion <- as.logical(inline.completion)

  w <- .RGtkCall("S_gtk_entry_completion_set_inline_completion", object, inline.completion, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryCompletionGetInlineCompletion <-
function(object)
{
  checkPtrType(object, "GtkEntryCompletion")

  w <- .RGtkCall("S_gtk_entry_completion_get_inline_completion", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryCompletionSetPopupCompletion <-
function(object, popup.completion)
{
  checkPtrType(object, "GtkEntryCompletion")
  popup.completion <- as.logical(popup.completion)

  w <- .RGtkCall("S_gtk_entry_completion_set_popup_completion", object, popup.completion, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryCompletionGetPopupCompletion <-
function(object)
{
  checkPtrType(object, "GtkEntryCompletion")

  w <- .RGtkCall("S_gtk_entry_completion_get_popup_completion", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryCompletionSetPopupSetWidth <-
function(object, popup.set.width)
{
  checkPtrType(object, "GtkEntryCompletion")
  popup.set.width <- as.logical(popup.set.width)

  w <- .RGtkCall("S_gtk_entry_completion_set_popup_set_width", object, popup.set.width, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryCompletionGetPopupSetWidth <-
function(object)
{
  checkPtrType(object, "GtkEntryCompletion")

  w <- .RGtkCall("S_gtk_entry_completion_get_popup_set_width", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryCompletionSetPopupSingleMatch <-
function(object, popup.single.match)
{
  checkPtrType(object, "GtkEntryCompletion")
  popup.single.match <- as.logical(popup.single.match)

  w <- .RGtkCall("S_gtk_entry_completion_set_popup_single_match", object, popup.single.match, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryCompletionGetPopupSingleMatch <-
function(object)
{
  checkPtrType(object, "GtkEntryCompletion")

  w <- .RGtkCall("S_gtk_entry_completion_get_popup_single_match", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEventBoxGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_event_box_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkEventBoxNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_event_box_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkEventBoxGetVisibleWindow <-
function(object)
{
  checkPtrType(object, "GtkEventBox")

  w <- .RGtkCall("S_gtk_event_box_get_visible_window", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEventBoxSetVisibleWindow <-
function(object, visible.window)
{
  checkPtrType(object, "GtkEventBox")
  visible.window <- as.logical(visible.window)

  w <- .RGtkCall("S_gtk_event_box_set_visible_window", object, visible.window, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEventBoxGetAboveChild <-
function(object)
{
  checkPtrType(object, "GtkEventBox")

  w <- .RGtkCall("S_gtk_event_box_get_above_child", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEventBoxSetAboveChild <-
function(object, above.child)
{
  checkPtrType(object, "GtkEventBox")
  above.child <- as.logical(above.child)

  w <- .RGtkCall("S_gtk_event_box_set_above_child", object, above.child, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkExpanderGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_expander_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkExpanderNew <-
function(label = NULL, show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_expander_new", label, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkExpanderNewWithMnemonic <-
function(label = NULL)
{
  if (!is.null( label )) label <- as.character(label)

  w <- .RGtkCall("S_gtk_expander_new_with_mnemonic", label, PACKAGE = "RGtk2")

  return(w)
} 


gtkExpanderSetExpanded <-
function(object, expanded)
{
  checkPtrType(object, "GtkExpander")
  expanded <- as.logical(expanded)

  w <- .RGtkCall("S_gtk_expander_set_expanded", object, expanded, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkExpanderGetExpanded <-
function(object)
{
  checkPtrType(object, "GtkExpander")

  w <- .RGtkCall("S_gtk_expander_get_expanded", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkExpanderSetSpacing <-
function(object, spacing)
{
  checkPtrType(object, "GtkExpander")
  spacing <- as.integer(spacing)

  w <- .RGtkCall("S_gtk_expander_set_spacing", object, spacing, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkExpanderGetSpacing <-
function(object)
{
  checkPtrType(object, "GtkExpander")

  w <- .RGtkCall("S_gtk_expander_get_spacing", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkExpanderSetLabel <-
function(object, label = NULL)
{
  checkPtrType(object, "GtkExpander")
  if (!is.null( label )) label <- as.character(label)

  w <- .RGtkCall("S_gtk_expander_set_label", object, label, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkExpanderGetLabel <-
function(object)
{
  checkPtrType(object, "GtkExpander")

  w <- .RGtkCall("S_gtk_expander_get_label", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkExpanderSetUseUnderline <-
function(object, use.underline)
{
  checkPtrType(object, "GtkExpander")
  use.underline <- as.logical(use.underline)

  w <- .RGtkCall("S_gtk_expander_set_use_underline", object, use.underline, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkExpanderGetUseUnderline <-
function(object)
{
  checkPtrType(object, "GtkExpander")

  w <- .RGtkCall("S_gtk_expander_get_use_underline", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkExpanderSetUseMarkup <-
function(object, use.markup)
{
  checkPtrType(object, "GtkExpander")
  use.markup <- as.logical(use.markup)

  w <- .RGtkCall("S_gtk_expander_set_use_markup", object, use.markup, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkExpanderGetUseMarkup <-
function(object)
{
  checkPtrType(object, "GtkExpander")

  w <- .RGtkCall("S_gtk_expander_get_use_markup", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkExpanderSetLabelWidget <-
function(object, label.widget = NULL)
{
  checkPtrType(object, "GtkExpander")
  if (!is.null( label.widget )) checkPtrType(label.widget, "GtkWidget")

  w <- .RGtkCall("S_gtk_expander_set_label_widget", object, label.widget, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkExpanderGetLabelWidget <-
function(object)
{
  checkPtrType(object, "GtkExpander")

  w <- .RGtkCall("S_gtk_expander_get_label_widget", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_file_chooser_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserErrorQuark <-
function()
{
  

  w <- .RGtkCall("S_gtk_file_chooser_error_quark", PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserSetAction <-
function(object, action)
{
  checkPtrType(object, "GtkFileChooser")
  

  w <- .RGtkCall("S_gtk_file_chooser_set_action", object, action, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFileChooserGetAction <-
function(object)
{
  checkPtrType(object, "GtkFileChooser")

  w <- .RGtkCall("S_gtk_file_chooser_get_action", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserSetLocalOnly <-
function(object, local.only)
{
  checkPtrType(object, "GtkFileChooser")
  local.only <- as.logical(local.only)

  w <- .RGtkCall("S_gtk_file_chooser_set_local_only", object, local.only, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFileChooserGetLocalOnly <-
function(object)
{
  checkPtrType(object, "GtkFileChooser")

  w <- .RGtkCall("S_gtk_file_chooser_get_local_only", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserSetSelectMultiple <-
function(object, select.multiple)
{
  checkPtrType(object, "GtkFileChooser")
  select.multiple <- as.logical(select.multiple)

  w <- .RGtkCall("S_gtk_file_chooser_set_select_multiple", object, select.multiple, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFileChooserGetSelectMultiple <-
function(object)
{
  checkPtrType(object, "GtkFileChooser")

  w <- .RGtkCall("S_gtk_file_chooser_get_select_multiple", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserSetShowHidden <-
function(object, show.hidden)
{
  checkPtrType(object, "GtkFileChooser")
  show.hidden <- as.logical(show.hidden)

  w <- .RGtkCall("S_gtk_file_chooser_set_show_hidden", object, show.hidden, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFileChooserGetShowHidden <-
function(object)
{
  checkPtrType(object, "GtkFileChooser")

  w <- .RGtkCall("S_gtk_file_chooser_get_show_hidden", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserSetDoOverwriteConfirmation <-
function(object, do.overwrite.confirmation)
{
  checkPtrType(object, "GtkFileChooser")
  do.overwrite.confirmation <- as.logical(do.overwrite.confirmation)

  w <- .RGtkCall("S_gtk_file_chooser_set_do_overwrite_confirmation", object, do.overwrite.confirmation, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFileChooserGetDoOverwriteConfirmation <-
function(object)
{
  checkPtrType(object, "GtkFileChooser")

  w <- .RGtkCall("S_gtk_file_chooser_get_do_overwrite_confirmation", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserSetCurrentName <-
function(object, name)
{
  checkPtrType(object, "GtkFileChooser")
  name <- as.character(name)

  w <- .RGtkCall("S_gtk_file_chooser_set_current_name", object, name, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFileChooserGetFilename <-
function(object)
{
  checkPtrType(object, "GtkFileChooser")

  w <- .RGtkCall("S_gtk_file_chooser_get_filename", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserSetFilename <-
function(object, filename)
{
  checkPtrType(object, "GtkFileChooser")
  filename <- as.character(filename)

  w <- .RGtkCall("S_gtk_file_chooser_set_filename", object, filename, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserSelectFilename <-
function(object, filename)
{
  checkPtrType(object, "GtkFileChooser")
  filename <- as.character(filename)

  w <- .RGtkCall("S_gtk_file_chooser_select_filename", object, filename, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserUnselectFilename <-
function(object, filename)
{
  checkPtrType(object, "GtkFileChooser")
  filename <- as.character(filename)

  w <- .RGtkCall("S_gtk_file_chooser_unselect_filename", object, filename, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFileChooserSelectAll <-
function(object)
{
  checkPtrType(object, "GtkFileChooser")

  w <- .RGtkCall("S_gtk_file_chooser_select_all", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFileChooserUnselectAll <-
function(object)
{
  checkPtrType(object, "GtkFileChooser")

  w <- .RGtkCall("S_gtk_file_chooser_unselect_all", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFileChooserGetFilenames <-
function(object)
{
  checkPtrType(object, "GtkFileChooser")

  w <- .RGtkCall("S_gtk_file_chooser_get_filenames", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserSetCurrentFolder <-
function(object, filename)
{
  checkPtrType(object, "GtkFileChooser")
  filename <- as.character(filename)

  w <- .RGtkCall("S_gtk_file_chooser_set_current_folder", object, filename, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserGetCurrentFolder <-
function(object)
{
  checkPtrType(object, "GtkFileChooser")

  w <- .RGtkCall("S_gtk_file_chooser_get_current_folder", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserGetUri <-
function(object)
{
  checkPtrType(object, "GtkFileChooser")

  w <- .RGtkCall("S_gtk_file_chooser_get_uri", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserSetUri <-
function(object, uri)
{
  checkPtrType(object, "GtkFileChooser")
  uri <- as.character(uri)

  w <- .RGtkCall("S_gtk_file_chooser_set_uri", object, uri, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserSelectUri <-
function(object, uri)
{
  checkPtrType(object, "GtkFileChooser")
  uri <- as.character(uri)

  w <- .RGtkCall("S_gtk_file_chooser_select_uri", object, uri, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserUnselectUri <-
function(object, uri)
{
  checkPtrType(object, "GtkFileChooser")
  uri <- as.character(uri)

  w <- .RGtkCall("S_gtk_file_chooser_unselect_uri", object, uri, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFileChooserGetUris <-
function(object)
{
  checkPtrType(object, "GtkFileChooser")

  w <- .RGtkCall("S_gtk_file_chooser_get_uris", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserSetCurrentFolderUri <-
function(object, uri)
{
  checkPtrType(object, "GtkFileChooser")
  uri <- as.character(uri)

  w <- .RGtkCall("S_gtk_file_chooser_set_current_folder_uri", object, uri, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserGetCurrentFolderUri <-
function(object)
{
  checkPtrType(object, "GtkFileChooser")

  w <- .RGtkCall("S_gtk_file_chooser_get_current_folder_uri", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserSetPreviewWidget <-
function(object, preview.widget)
{
  checkPtrType(object, "GtkFileChooser")
  checkPtrType(preview.widget, "GtkWidget")

  w <- .RGtkCall("S_gtk_file_chooser_set_preview_widget", object, preview.widget, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFileChooserGetPreviewWidget <-
function(object)
{
  checkPtrType(object, "GtkFileChooser")

  w <- .RGtkCall("S_gtk_file_chooser_get_preview_widget", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserSetPreviewWidgetActive <-
function(object, active)
{
  checkPtrType(object, "GtkFileChooser")
  active <- as.logical(active)

  w <- .RGtkCall("S_gtk_file_chooser_set_preview_widget_active", object, active, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFileChooserGetPreviewWidgetActive <-
function(object)
{
  checkPtrType(object, "GtkFileChooser")

  w <- .RGtkCall("S_gtk_file_chooser_get_preview_widget_active", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserSetUsePreviewLabel <-
function(object, use.label)
{
  checkPtrType(object, "GtkFileChooser")
  use.label <- as.logical(use.label)

  w <- .RGtkCall("S_gtk_file_chooser_set_use_preview_label", object, use.label, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFileChooserGetUsePreviewLabel <-
function(object)
{
  checkPtrType(object, "GtkFileChooser")

  w <- .RGtkCall("S_gtk_file_chooser_get_use_preview_label", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserGetPreviewFilename <-
function(object)
{
  checkPtrType(object, "GtkFileChooser")

  w <- .RGtkCall("S_gtk_file_chooser_get_preview_filename", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserGetPreviewUri <-
function(object)
{
  checkPtrType(object, "GtkFileChooser")

  w <- .RGtkCall("S_gtk_file_chooser_get_preview_uri", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserSetExtraWidget <-
function(object, extra.widget)
{
  checkPtrType(object, "GtkFileChooser")
  checkPtrType(extra.widget, "GtkWidget")

  w <- .RGtkCall("S_gtk_file_chooser_set_extra_widget", object, extra.widget, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFileChooserGetExtraWidget <-
function(object)
{
  checkPtrType(object, "GtkFileChooser")

  w <- .RGtkCall("S_gtk_file_chooser_get_extra_widget", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserAddFilter <-
function(object, filter)
{
  checkPtrType(object, "GtkFileChooser")
  checkPtrType(filter, "GtkFileFilter")

  w <- .RGtkCall("S_gtk_file_chooser_add_filter", object, filter, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFileChooserRemoveFilter <-
function(object, filter)
{
  checkPtrType(object, "GtkFileChooser")
  checkPtrType(filter, "GtkFileFilter")

  w <- .RGtkCall("S_gtk_file_chooser_remove_filter", object, filter, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFileChooserListFilters <-
function(object)
{
  checkPtrType(object, "GtkFileChooser")

  w <- .RGtkCall("S_gtk_file_chooser_list_filters", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserSetFilter <-
function(object, filter)
{
  checkPtrType(object, "GtkFileChooser")
  checkPtrType(filter, "GtkFileFilter")

  w <- .RGtkCall("S_gtk_file_chooser_set_filter", object, filter, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFileChooserGetFilter <-
function(object)
{
  checkPtrType(object, "GtkFileChooser")

  w <- .RGtkCall("S_gtk_file_chooser_get_filter", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserAddShortcutFolder <-
function(object, folder, .errwarn = TRUE)
{
  checkPtrType(object, "GtkFileChooser")
  folder <- as.character(folder)

  w <- .RGtkCall("S_gtk_file_chooser_add_shortcut_folder", object, folder, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkFileChooserRemoveShortcutFolder <-
function(object, folder, .errwarn = TRUE)
{
  checkPtrType(object, "GtkFileChooser")
  folder <- as.character(folder)

  w <- .RGtkCall("S_gtk_file_chooser_remove_shortcut_folder", object, folder, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkFileChooserListShortcutFolders <-
function(object)
{
  checkPtrType(object, "GtkFileChooser")

  w <- .RGtkCall("S_gtk_file_chooser_list_shortcut_folders", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserAddShortcutFolderUri <-
function(object, uri, .errwarn = TRUE)
{
  checkPtrType(object, "GtkFileChooser")
  uri <- as.character(uri)

  w <- .RGtkCall("S_gtk_file_chooser_add_shortcut_folder_uri", object, uri, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkFileChooserRemoveShortcutFolderUri <-
function(object, uri, .errwarn = TRUE)
{
  checkPtrType(object, "GtkFileChooser")
  uri <- as.character(uri)

  w <- .RGtkCall("S_gtk_file_chooser_remove_shortcut_folder_uri", object, uri, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkFileChooserListShortcutFolderUris <-
function(object)
{
  checkPtrType(object, "GtkFileChooser")

  w <- .RGtkCall("S_gtk_file_chooser_list_shortcut_folder_uris", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserButtonGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_file_chooser_button_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserButtonNew <-
function(title, action, show = TRUE)
{
  title <- as.character(title)
  

  w <- .RGtkCall("S_gtk_file_chooser_button_new", title, action, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkFileChooserButtonNewWithBackend <-
function(title, action, backend, show = TRUE)
{
  if(getOption("depwarn"))
    .Deprecated("gtkFileChooserButtonNew", "RGtk2")

  title <- as.character(title)
  
  backend <- as.character(backend)

  w <- .RGtkCall("S_gtk_file_chooser_button_new_with_backend", title, action, backend, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkFileChooserButtonNewWithDialog <-
function(dialog)
{
  checkPtrType(dialog, "GtkWidget")

  w <- .RGtkCall("S_gtk_file_chooser_button_new_with_dialog", dialog, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserButtonGetTitle <-
function(object)
{
  checkPtrType(object, "GtkFileChooserButton")

  w <- .RGtkCall("S_gtk_file_chooser_button_get_title", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserButtonSetTitle <-
function(object, title)
{
  checkPtrType(object, "GtkFileChooserButton")
  title <- as.character(title)

  w <- .RGtkCall("S_gtk_file_chooser_button_set_title", object, title, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFileChooserButtonGetWidthChars <-
function(object)
{
  checkPtrType(object, "GtkFileChooserButton")

  w <- .RGtkCall("S_gtk_file_chooser_button_get_width_chars", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserButtonSetWidthChars <-
function(object, n.chars)
{
  checkPtrType(object, "GtkFileChooserButton")
  n.chars <- as.integer(n.chars)

  w <- .RGtkCall("S_gtk_file_chooser_button_set_width_chars", object, n.chars, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFileChooserDialogGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_file_chooser_dialog_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserWidgetGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_file_chooser_widget_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserWidgetNew <-
function(action, show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_file_chooser_widget_new", action, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkFileChooserWidgetNewWithBackend <-
function(action, backend, show = TRUE)
{
  if(getOption("depwarn"))
    .Deprecated("gtkFileChooserWidgetNew", "RGtk2")

  
  backend <- as.character(backend)

  w <- .RGtkCall("S_gtk_file_chooser_widget_new_with_backend", action, backend, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkFileFilterGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_file_filter_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkFileFilterNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_file_filter_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkFileFilterSetName <-
function(object, name)
{
  checkPtrType(object, "GtkFileFilter")
  name <- as.character(name)

  w <- .RGtkCall("S_gtk_file_filter_set_name", object, name, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFileFilterGetName <-
function(object)
{
  checkPtrType(object, "GtkFileFilter")

  w <- .RGtkCall("S_gtk_file_filter_get_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileFilterAddMimeType <-
function(object, mime.type)
{
  checkPtrType(object, "GtkFileFilter")
  mime.type <- as.character(mime.type)

  w <- .RGtkCall("S_gtk_file_filter_add_mime_type", object, mime.type, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFileFilterAddPattern <-
function(object, pattern)
{
  checkPtrType(object, "GtkFileFilter")
  pattern <- as.character(pattern)

  w <- .RGtkCall("S_gtk_file_filter_add_pattern", object, pattern, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFileFilterAddPixbufFormats <-
function(object)
{
  checkPtrType(object, "GtkFileFilter")

  w <- .RGtkCall("S_gtk_file_filter_add_pixbuf_formats", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFileFilterAddCustom <-
function(object, needed, func, data = NULL)
{
  checkPtrType(object, "GtkFileFilter")
  
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_file_filter_add_custom", object, needed, func, data, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileFilterGetNeeded <-
function(object)
{
  checkPtrType(object, "GtkFileFilter")

  w <- .RGtkCall("S_gtk_file_filter_get_needed", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileFilterFilter <-
function(object, filter.info)
{
  checkPtrType(object, "GtkFileFilter")
  filter.info <- as.GtkFileFilterInfo(filter.info)

  w <- .RGtkCall("S_gtk_file_filter_filter", object, filter.info, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileSelectionGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_file_selection_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkFileSelectionNew <-
function(title = NULL, show = TRUE)
{
  if(getOption("depwarn"))
    .Deprecated("GtkFileChooser", "RGtk2")

  

  w <- .RGtkCall("S_gtk_file_selection_new", title, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkFileSelectionSetFilename <-
function(object, filename)
{
  if(getOption("depwarn"))
    .Deprecated("GtkFileChooser", "RGtk2")

  checkPtrType(object, "GtkFileSelection")
  filename <- as.character(filename)

  w <- .RGtkCall("S_gtk_file_selection_set_filename", object, filename, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFileSelectionGetFilename <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("GtkFileChooser", "RGtk2")

  checkPtrType(object, "GtkFileSelection")

  w <- .RGtkCall("S_gtk_file_selection_get_filename", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileSelectionComplete <-
function(object, pattern)
{
  if(getOption("depwarn"))
    .Deprecated("GtkFileChooser", "RGtk2")

  checkPtrType(object, "GtkFileSelection")
  pattern <- as.character(pattern)

  w <- .RGtkCall("S_gtk_file_selection_complete", object, pattern, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFileSelectionShowFileopButtons <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("GtkFileChooser", "RGtk2")

  checkPtrType(object, "GtkFileSelection")

  w <- .RGtkCall("S_gtk_file_selection_show_fileop_buttons", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFileSelectionHideFileopButtons <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("GtkFileChooser", "RGtk2")

  checkPtrType(object, "GtkFileSelection")

  w <- .RGtkCall("S_gtk_file_selection_hide_fileop_buttons", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFileSelectionGetSelections <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("GtkFileChooser", "RGtk2")

  checkPtrType(object, "GtkFileSelection")

  w <- .RGtkCall("S_gtk_file_selection_get_selections", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileSelectionSetSelectMultiple <-
function(object, select.multiple)
{
  if(getOption("depwarn"))
    .Deprecated("GtkFileChooser", "RGtk2")

  checkPtrType(object, "GtkFileSelection")
  select.multiple <- as.logical(select.multiple)

  w <- .RGtkCall("S_gtk_file_selection_set_select_multiple", object, select.multiple, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFileSelectionGetSelectMultiple <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("GtkFileChooser", "RGtk2")

  checkPtrType(object, "GtkFileSelection")

  w <- .RGtkCall("S_gtk_file_selection_get_select_multiple", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFixedGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_fixed_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkFixedNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_fixed_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkFixedPut <-
function(object, widget, x, y)
{
  checkPtrType(object, "GtkFixed")
  checkPtrType(widget, "GtkWidget")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_gtk_fixed_put", object, widget, x, y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFixedMove <-
function(object, widget, x, y)
{
  checkPtrType(object, "GtkFixed")
  checkPtrType(widget, "GtkWidget")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_gtk_fixed_move", object, widget, x, y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFixedSetHasWindow <-
function(object, has.window)
{
  checkPtrType(object, "GtkFixed")
  has.window <- as.logical(has.window)

  w <- .RGtkCall("S_gtk_fixed_set_has_window", object, has.window, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFixedGetHasWindow <-
function(object)
{
  checkPtrType(object, "GtkFixed")

  w <- .RGtkCall("S_gtk_fixed_get_has_window", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFontButtonGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_font_button_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkFontButtonNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_font_button_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkFontButtonNewWithFont <-
function(fontname)
{
  fontname <- as.character(fontname)

  w <- .RGtkCall("S_gtk_font_button_new_with_font", fontname, PACKAGE = "RGtk2")

  return(w)
} 


gtkFontButtonGetTitle <-
function(object)
{
  checkPtrType(object, "GtkFontButton")

  w <- .RGtkCall("S_gtk_font_button_get_title", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFontButtonSetTitle <-
function(object, title)
{
  checkPtrType(object, "GtkFontButton")
  title <- as.character(title)

  w <- .RGtkCall("S_gtk_font_button_set_title", object, title, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFontButtonGetUseFont <-
function(object)
{
  checkPtrType(object, "GtkFontButton")

  w <- .RGtkCall("S_gtk_font_button_get_use_font", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFontButtonSetUseFont <-
function(object, use.font)
{
  checkPtrType(object, "GtkFontButton")
  use.font <- as.logical(use.font)

  w <- .RGtkCall("S_gtk_font_button_set_use_font", object, use.font, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFontButtonGetUseSize <-
function(object)
{
  checkPtrType(object, "GtkFontButton")

  w <- .RGtkCall("S_gtk_font_button_get_use_size", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFontButtonSetUseSize <-
function(object, use.size)
{
  checkPtrType(object, "GtkFontButton")
  use.size <- as.logical(use.size)

  w <- .RGtkCall("S_gtk_font_button_set_use_size", object, use.size, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFontButtonGetFontName <-
function(object)
{
  checkPtrType(object, "GtkFontButton")

  w <- .RGtkCall("S_gtk_font_button_get_font_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFontButtonSetFontName <-
function(object, fontname)
{
  checkPtrType(object, "GtkFontButton")
  fontname <- as.character(fontname)

  w <- .RGtkCall("S_gtk_font_button_set_font_name", object, fontname, PACKAGE = "RGtk2")

  return(w)
} 


gtkFontButtonGetShowStyle <-
function(object)
{
  checkPtrType(object, "GtkFontButton")

  w <- .RGtkCall("S_gtk_font_button_get_show_style", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFontButtonSetShowStyle <-
function(object, show.style)
{
  checkPtrType(object, "GtkFontButton")
  show.style <- as.logical(show.style)

  w <- .RGtkCall("S_gtk_font_button_set_show_style", object, show.style, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFontButtonGetShowSize <-
function(object)
{
  checkPtrType(object, "GtkFontButton")

  w <- .RGtkCall("S_gtk_font_button_get_show_size", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFontButtonSetShowSize <-
function(object, show.size)
{
  checkPtrType(object, "GtkFontButton")
  show.size <- as.logical(show.size)

  w <- .RGtkCall("S_gtk_font_button_set_show_size", object, show.size, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFontSelectionGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_font_selection_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkFontSelectionNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_font_selection_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkFontSelectionGetFontName <-
function(object)
{
  checkPtrType(object, "GtkFontSelection")

  w <- .RGtkCall("S_gtk_font_selection_get_font_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFontSelectionGetFont <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("gtkFontSelectionGetFontName", "RGtk2")

  checkPtrType(object, "GtkFontSelection")

  w <- .RGtkCall("S_gtk_font_selection_get_font", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFontSelectionSetFontName <-
function(object, fontname)
{
  checkPtrType(object, "GtkFontSelection")
  fontname <- as.character(fontname)

  w <- .RGtkCall("S_gtk_font_selection_set_font_name", object, fontname, PACKAGE = "RGtk2")

  return(w)
} 


gtkFontSelectionGetPreviewText <-
function(object)
{
  checkPtrType(object, "GtkFontSelection")

  w <- .RGtkCall("S_gtk_font_selection_get_preview_text", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFontSelectionSetPreviewText <-
function(object, text)
{
  checkPtrType(object, "GtkFontSelection")
  text <- as.character(text)

  w <- .RGtkCall("S_gtk_font_selection_set_preview_text", object, text, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFontSelectionDialogGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_font_selection_dialog_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkFontSelectionDialogNew <-
function(title = NULL, show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_font_selection_dialog_new", title, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkFontSelectionDialogGetFontName <-
function(object)
{
  checkPtrType(object, "GtkFontSelectionDialog")

  w <- .RGtkCall("S_gtk_font_selection_dialog_get_font_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFontSelectionDialogGetFont <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("gtkFontSelectionGetFontName", "RGtk2")

  checkPtrType(object, "GtkFontSelectionDialog")

  w <- .RGtkCall("S_gtk_font_selection_dialog_get_font", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFontSelectionDialogSetFontName <-
function(object, fontname)
{
  checkPtrType(object, "GtkFontSelectionDialog")
  fontname <- as.character(fontname)

  w <- .RGtkCall("S_gtk_font_selection_dialog_set_font_name", object, fontname, PACKAGE = "RGtk2")

  return(w)
} 


gtkFontSelectionDialogGetPreviewText <-
function(object)
{
  checkPtrType(object, "GtkFontSelectionDialog")

  w <- .RGtkCall("S_gtk_font_selection_dialog_get_preview_text", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFontSelectionDialogSetPreviewText <-
function(object, text)
{
  checkPtrType(object, "GtkFontSelectionDialog")
  text <- as.character(text)

  w <- .RGtkCall("S_gtk_font_selection_dialog_set_preview_text", object, text, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFrameGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_frame_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkFrameNew <-
function(label = NULL, show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_frame_new", label, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkFrameSetLabel <-
function(object, label = NULL)
{
  checkPtrType(object, "GtkFrame")
  if (!is.null( label )) label <- as.character(label)

  w <- .RGtkCall("S_gtk_frame_set_label", object, label, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFrameGetLabel <-
function(object)
{
  checkPtrType(object, "GtkFrame")

  w <- .RGtkCall("S_gtk_frame_get_label", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFrameSetLabelWidget <-
function(object, label.widget)
{
  checkPtrType(object, "GtkFrame")
  checkPtrType(label.widget, "GtkWidget")

  w <- .RGtkCall("S_gtk_frame_set_label_widget", object, label.widget, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFrameGetLabelWidget <-
function(object)
{
  checkPtrType(object, "GtkFrame")

  w <- .RGtkCall("S_gtk_frame_get_label_widget", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFrameSetLabelAlign <-
function(object, xalign, yalign)
{
  checkPtrType(object, "GtkFrame")
  xalign <- as.numeric(xalign)
  yalign <- as.numeric(yalign)

  w <- .RGtkCall("S_gtk_frame_set_label_align", object, xalign, yalign, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFrameGetLabelAlign <-
function(object)
{
  checkPtrType(object, "GtkFrame")

  w <- .RGtkCall("S_gtk_frame_get_label_align", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFrameSetShadowType <-
function(object, type)
{
  checkPtrType(object, "GtkFrame")
  

  w <- .RGtkCall("S_gtk_frame_set_shadow_type", object, type, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFrameGetShadowType <-
function(object)
{
  checkPtrType(object, "GtkFrame")

  w <- .RGtkCall("S_gtk_frame_get_shadow_type", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkGammaCurveGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_gamma_curve_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkGammaCurveNew <-
function(show = TRUE)
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  

  w <- .RGtkCall("S_gtk_gamma_curve_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkGcGet <-
function(depth, colormap, values)
{
  depth <- as.integer(depth)
  checkPtrType(colormap, "GdkColormap")
  values <- as.GdkGCValues(values)

  w <- .RGtkCall("S_gtk_gc_get", depth, colormap, values, PACKAGE = "RGtk2")

  return(w)
} 


gtkGcRelease <-
function(gc)
{
  checkPtrType(gc, "GdkGC")

  w <- .RGtkCall("S_gtk_gc_release", gc, PACKAGE = "RGtk2")

  return(w)
} 


gtkHandleBoxGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_handle_box_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkHandleBoxNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_handle_box_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkHandleBoxSetShadowType <-
function(object, type)
{
  checkPtrType(object, "GtkHandleBox")
  

  w <- .RGtkCall("S_gtk_handle_box_set_shadow_type", object, type, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkHandleBoxGetShadowType <-
function(object)
{
  checkPtrType(object, "GtkHandleBox")

  w <- .RGtkCall("S_gtk_handle_box_get_shadow_type", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkHandleBoxSetHandlePosition <-
function(object, position)
{
  checkPtrType(object, "GtkHandleBox")
  

  w <- .RGtkCall("S_gtk_handle_box_set_handle_position", object, position, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkHandleBoxGetHandlePosition <-
function(object)
{
  checkPtrType(object, "GtkHandleBox")

  w <- .RGtkCall("S_gtk_handle_box_get_handle_position", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkHandleBoxSetSnapEdge <-
function(object, edge)
{
  checkPtrType(object, "GtkHandleBox")
  

  w <- .RGtkCall("S_gtk_handle_box_set_snap_edge", object, edge, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkHandleBoxGetSnapEdge <-
function(object)
{
  checkPtrType(object, "GtkHandleBox")

  w <- .RGtkCall("S_gtk_handle_box_get_snap_edge", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkHButtonBoxGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_hbutton_box_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkHButtonBoxNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_hbutton_box_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkHButtonBoxGetSpacingDefault <-
function()
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  

  w <- .RGtkCall("S_gtk_hbutton_box_get_spacing_default", PACKAGE = "RGtk2")

  return(w)
} 


gtkHButtonBoxGetLayoutDefault <-
function()
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  

  w <- .RGtkCall("S_gtk_hbutton_box_get_layout_default", PACKAGE = "RGtk2")

  return(w)
} 


gtkHButtonBoxSetSpacingDefault <-
function(spacing)
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  spacing <- as.integer(spacing)

  w <- .RGtkCall("S_gtk_hbutton_box_set_spacing_default", spacing, PACKAGE = "RGtk2")

  return(w)
} 


gtkHButtonBoxSetLayoutDefault <-
function(layout)
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  

  w <- .RGtkCall("S_gtk_hbutton_box_set_layout_default", layout, PACKAGE = "RGtk2")

  return(w)
} 


gtkHBoxGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_hbox_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkHBoxNew <-
function(homogeneous = NULL, spacing = NULL, show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_hbox_new", homogeneous, spacing, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkHPanedGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_hpaned_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkHPanedNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_hpaned_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkHRulerGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_hruler_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkHRulerNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_hruler_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkHScaleGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_hscale_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkHScaleNew <-
function(adjustment = NULL, show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_hscale_new", adjustment, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkHScaleNewWithRange <-
function(min, max, step, show = TRUE)
{
  min <- as.numeric(min)
  max <- as.numeric(max)
  step <- as.numeric(step)

  w <- .RGtkCall("S_gtk_hscale_new_with_range", min, max, step, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkHScrollbarGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_hscrollbar_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkHScrollbarNew <-
function(adjustment = NULL, show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_hscrollbar_new", adjustment, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkHSeparatorGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_hseparator_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkHSeparatorNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_hseparator_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkIconFactoryGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_icon_factory_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkIconFactoryNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_icon_factory_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkIconFactoryAdd <-
function(object, stock.id, icon.set)
{
  checkPtrType(object, "GtkIconFactory")
  stock.id <- as.character(stock.id)
  checkPtrType(icon.set, "GtkIconSet")

  w <- .RGtkCall("S_gtk_icon_factory_add", object, stock.id, icon.set, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconFactoryLookup <-
function(object, stock.id)
{
  checkPtrType(object, "GtkIconFactory")
  stock.id <- as.character(stock.id)

  w <- .RGtkCall("S_gtk_icon_factory_lookup", object, stock.id, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconFactoryAddDefault <-
function(object)
{
  checkPtrType(object, "GtkIconFactory")

  w <- .RGtkCall("S_gtk_icon_factory_add_default", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconFactoryRemoveDefault <-
function(object)
{
  checkPtrType(object, "GtkIconFactory")

  w <- .RGtkCall("S_gtk_icon_factory_remove_default", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconFactoryLookupDefault <-
function(stock.id)
{
  stock.id <- as.character(stock.id)

  w <- .RGtkCall("S_gtk_icon_factory_lookup_default", stock.id, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconSizeLookup <-
function(size)
{
  

  w <- .RGtkCall("S_gtk_icon_size_lookup", size, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconSizeLookupForSettings <-
function(settings, size)
{
  checkPtrType(settings, "GtkSettings")
  

  w <- .RGtkCall("S_gtk_icon_size_lookup_for_settings", settings, size, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconSizeRegister <-
function(name, width, height)
{
  name <- as.character(name)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_icon_size_register", name, width, height, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconSizeRegisterAlias <-
function(alias, target)
{
  alias <- as.character(alias)
  

  w <- .RGtkCall("S_gtk_icon_size_register_alias", alias, target, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconSizeFromName <-
function(name)
{
  name <- as.character(name)

  w <- .RGtkCall("S_gtk_icon_size_from_name", name, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconSizeGetName <-
function(size)
{
  

  w <- .RGtkCall("S_gtk_icon_size_get_name", size, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconSetGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_icon_set_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkIconSetNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_icon_set_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkIconSetNewFromPixbuf <-
function(pixbuf)
{
  checkPtrType(pixbuf, "GdkPixbuf")

  w <- .RGtkCall("S_gtk_icon_set_new_from_pixbuf", pixbuf, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconSetCopy <-
function(object)
{
  checkPtrType(object, "GtkIconSet")

  w <- .RGtkCall("S_gtk_icon_set_copy", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconSetRenderIcon <-
function(object, style, direction, state, size, widget = NULL, detail = NULL)
{
  checkPtrType(object, "GtkIconSet")
  checkPtrType(style, "GtkStyle")
  
  
  
  if (!is.null( widget )) checkPtrType(widget, "GtkWidget")
  if (!is.null( detail )) detail <- as.character(detail)

  w <- .RGtkCall("S_gtk_icon_set_render_icon", object, style, direction, state, size, widget, detail, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconSetAddSource <-
function(object, source)
{
  checkPtrType(object, "GtkIconSet")
  checkPtrType(source, "GtkIconSource")

  w <- .RGtkCall("S_gtk_icon_set_add_source", object, source, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconSetGetSizes <-
function(object)
{
  checkPtrType(object, "GtkIconSet")

  w <- .RGtkCall("S_gtk_icon_set_get_sizes", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconSourceGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_icon_source_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkIconSourceNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_icon_source_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkIconSourceCopy <-
function(object)
{
  checkPtrType(object, "GtkIconSource")

  w <- .RGtkCall("S_gtk_icon_source_copy", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconSourceSetFilename <-
function(object, filename)
{
  checkPtrType(object, "GtkIconSource")
  filename <- as.character(filename)

  w <- .RGtkCall("S_gtk_icon_source_set_filename", object, filename, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconSourceSetIconName <-
function(object, icon.name)
{
  checkPtrType(object, "GtkIconSource")
  icon.name <- as.character(icon.name)

  w <- .RGtkCall("S_gtk_icon_source_set_icon_name", object, icon.name, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconSourceSetPixbuf <-
function(object, pixbuf)
{
  checkPtrType(object, "GtkIconSource")
  checkPtrType(pixbuf, "GdkPixbuf")

  w <- .RGtkCall("S_gtk_icon_source_set_pixbuf", object, pixbuf, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconSourceGetFilename <-
function(object)
{
  checkPtrType(object, "GtkIconSource")

  w <- .RGtkCall("S_gtk_icon_source_get_filename", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconSourceGetIconName <-
function(object)
{
  checkPtrType(object, "GtkIconSource")

  w <- .RGtkCall("S_gtk_icon_source_get_icon_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconSourceGetPixbuf <-
function(object)
{
  checkPtrType(object, "GtkIconSource")

  w <- .RGtkCall("S_gtk_icon_source_get_pixbuf", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconSourceSetDirectionWildcarded <-
function(object, setting)
{
  checkPtrType(object, "GtkIconSource")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_icon_source_set_direction_wildcarded", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconSourceSetStateWildcarded <-
function(object, setting)
{
  checkPtrType(object, "GtkIconSource")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_icon_source_set_state_wildcarded", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconSourceSetSizeWildcarded <-
function(object, setting)
{
  checkPtrType(object, "GtkIconSource")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_icon_source_set_size_wildcarded", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconSourceGetSizeWildcarded <-
function(object)
{
  checkPtrType(object, "GtkIconSource")

  w <- .RGtkCall("S_gtk_icon_source_get_size_wildcarded", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconSourceGetStateWildcarded <-
function(object)
{
  checkPtrType(object, "GtkIconSource")

  w <- .RGtkCall("S_gtk_icon_source_get_state_wildcarded", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconSourceGetDirectionWildcarded <-
function(object)
{
  checkPtrType(object, "GtkIconSource")

  w <- .RGtkCall("S_gtk_icon_source_get_direction_wildcarded", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconSourceSetDirection <-
function(object, direction)
{
  checkPtrType(object, "GtkIconSource")
  

  w <- .RGtkCall("S_gtk_icon_source_set_direction", object, direction, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconSourceSetState <-
function(object, state)
{
  checkPtrType(object, "GtkIconSource")
  

  w <- .RGtkCall("S_gtk_icon_source_set_state", object, state, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconSourceSetSize <-
function(object, size)
{
  checkPtrType(object, "GtkIconSource")
  

  w <- .RGtkCall("S_gtk_icon_source_set_size", object, size, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconSourceGetDirection <-
function(object)
{
  checkPtrType(object, "GtkIconSource")

  w <- .RGtkCall("S_gtk_icon_source_get_direction", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconSourceGetState <-
function(object)
{
  checkPtrType(object, "GtkIconSource")

  w <- .RGtkCall("S_gtk_icon_source_get_state", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconSourceGetSize <-
function(object)
{
  checkPtrType(object, "GtkIconSource")

  w <- .RGtkCall("S_gtk_icon_source_get_size", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconThemeErrorQuark <-
function()
{
  

  w <- .RGtkCall("S_gtk_icon_theme_error_quark", PACKAGE = "RGtk2")

  return(w)
} 


gtkIconThemeGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_icon_theme_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkIconThemeNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_icon_theme_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkIconThemeGetDefault <-
function()
{
  

  w <- .RGtkCall("S_gtk_icon_theme_get_default", PACKAGE = "RGtk2")

  return(w)
} 


gtkIconThemeGetForScreen <-
function(screen)
{
  checkPtrType(screen, "GdkScreen")

  w <- .RGtkCall("S_gtk_icon_theme_get_for_screen", screen, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconThemeSetScreen <-
function(object, screen)
{
  checkPtrType(object, "GtkIconTheme")
  checkPtrType(screen, "GdkScreen")

  w <- .RGtkCall("S_gtk_icon_theme_set_screen", object, screen, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconThemeGetSearchPath <-
function(object)
{
  checkPtrType(object, "GtkIconTheme")

  w <- .RGtkCall("S_gtk_icon_theme_get_search_path", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconThemeAppendSearchPath <-
function(object, path)
{
  checkPtrType(object, "GtkIconTheme")
  path <- as.character(path)

  w <- .RGtkCall("S_gtk_icon_theme_append_search_path", object, path, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconThemePrependSearchPath <-
function(object, path)
{
  checkPtrType(object, "GtkIconTheme")
  path <- as.character(path)

  w <- .RGtkCall("S_gtk_icon_theme_prepend_search_path", object, path, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconThemeSetCustomTheme <-
function(object, theme.name)
{
  checkPtrType(object, "GtkIconTheme")
  theme.name <- as.character(theme.name)

  w <- .RGtkCall("S_gtk_icon_theme_set_custom_theme", object, theme.name, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconThemeHasIcon <-
function(object, icon.name)
{
  checkPtrType(object, "GtkIconTheme")
  icon.name <- as.character(icon.name)

  w <- .RGtkCall("S_gtk_icon_theme_has_icon", object, icon.name, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconThemeLookupIcon <-
function(object, icon.name, size, flags)
{
  checkPtrType(object, "GtkIconTheme")
  icon.name <- as.character(icon.name)
  size <- as.integer(size)
  

  w <- .RGtkCall("S_gtk_icon_theme_lookup_icon", object, icon.name, size, flags, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconThemeLoadIcon <-
function(object, icon.name, size, flags, .errwarn = TRUE)
{
  checkPtrType(object, "GtkIconTheme")
  icon.name <- as.character(icon.name)
  size <- as.integer(size)
  

  w <- .RGtkCall("S_gtk_icon_theme_load_icon", object, icon.name, size, flags, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkIconThemeListIcons <-
function(object, context = NULL)
{
  checkPtrType(object, "GtkIconTheme")
  if (!is.null( context )) context <- as.character(context)

  w <- .RGtkCall("S_gtk_icon_theme_list_icons", object, context, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconThemeGetExampleIconName <-
function(object)
{
  checkPtrType(object, "GtkIconTheme")

  w <- .RGtkCall("S_gtk_icon_theme_get_example_icon_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconThemeRescanIfNeeded <-
function(object)
{
  checkPtrType(object, "GtkIconTheme")

  w <- .RGtkCall("S_gtk_icon_theme_rescan_if_needed", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconThemeAddBuiltinIcon <-
function(icon.name, size, pixbuf)
{
  icon.name <- as.character(icon.name)
  size <- as.integer(size)
  checkPtrType(pixbuf, "GdkPixbuf")

  w <- .RGtkCall("S_gtk_icon_theme_add_builtin_icon", icon.name, size, pixbuf, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconInfoGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_icon_info_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkIconInfoCopy <-
function(object)
{
  checkPtrType(object, "GtkIconInfo")

  w <- .RGtkCall("S_gtk_icon_info_copy", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconInfoGetBaseSize <-
function(object)
{
  checkPtrType(object, "GtkIconInfo")

  w <- .RGtkCall("S_gtk_icon_info_get_base_size", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconInfoGetFilename <-
function(object)
{
  checkPtrType(object, "GtkIconInfo")

  w <- .RGtkCall("S_gtk_icon_info_get_filename", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconInfoGetBuiltinPixbuf <-
function(object)
{
  checkPtrType(object, "GtkIconInfo")

  w <- .RGtkCall("S_gtk_icon_info_get_builtin_pixbuf", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconInfoLoadIcon <-
function(object, .errwarn = TRUE)
{
  checkPtrType(object, "GtkIconInfo")

  w <- .RGtkCall("S_gtk_icon_info_load_icon", object, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkIconInfoSetRawCoordinates <-
function(object, raw.coordinates)
{
  checkPtrType(object, "GtkIconInfo")
  raw.coordinates <- as.logical(raw.coordinates)

  w <- .RGtkCall("S_gtk_icon_info_set_raw_coordinates", object, raw.coordinates, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconInfoGetEmbeddedRect <-
function(object)
{
  checkPtrType(object, "GtkIconInfo")

  w <- .RGtkCall("S_gtk_icon_info_get_embedded_rect", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconInfoGetAttachPoints <-
function(object)
{
  checkPtrType(object, "GtkIconInfo")

  w <- .RGtkCall("S_gtk_icon_info_get_attach_points", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconInfoGetDisplayName <-
function(object)
{
  checkPtrType(object, "GtkIconInfo")

  w <- .RGtkCall("S_gtk_icon_info_get_display_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconThemeGetIconSizes <-
function(object, icon.name)
{
  checkPtrType(object, "GtkIconTheme")
  icon.name <- as.character(icon.name)

  w <- .RGtkCall("S_gtk_icon_theme_get_icon_sizes", object, icon.name, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconViewGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_icon_view_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkIconViewNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_icon_view_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkIconViewNewWithModel <-
function(model = NULL, show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_icon_view_new_with_model", model, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkIconViewSetModel <-
function(object, model = NULL)
{
  checkPtrType(object, "GtkIconView")
  if (!is.null( model )) checkPtrType(model, "GtkTreeModel")

  w <- .RGtkCall("S_gtk_icon_view_set_model", object, model, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconViewGetModel <-
function(object)
{
  checkPtrType(object, "GtkIconView")

  w <- .RGtkCall("S_gtk_icon_view_get_model", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconViewSetTextColumn <-
function(object, column)
{
  checkPtrType(object, "GtkIconView")
  column <- as.integer(column)

  w <- .RGtkCall("S_gtk_icon_view_set_text_column", object, column, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconViewGetTextColumn <-
function(object)
{
  checkPtrType(object, "GtkIconView")

  w <- .RGtkCall("S_gtk_icon_view_get_text_column", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconViewSetMarkupColumn <-
function(object, column)
{
  checkPtrType(object, "GtkIconView")
  column <- as.integer(column)

  w <- .RGtkCall("S_gtk_icon_view_set_markup_column", object, column, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconViewGetMarkupColumn <-
function(object)
{
  checkPtrType(object, "GtkIconView")

  w <- .RGtkCall("S_gtk_icon_view_get_markup_column", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconViewSetPixbufColumn <-
function(object, column)
{
  checkPtrType(object, "GtkIconView")
  column <- as.integer(column)

  w <- .RGtkCall("S_gtk_icon_view_set_pixbuf_column", object, column, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconViewGetPixbufColumn <-
function(object)
{
  checkPtrType(object, "GtkIconView")

  w <- .RGtkCall("S_gtk_icon_view_get_pixbuf_column", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconViewSetOrientation <-
function(object, orientation)
{
  checkPtrType(object, "GtkIconView")
  

  w <- .RGtkCall("S_gtk_icon_view_set_orientation", object, orientation, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconViewGetOrientation <-
function(object)
{
  checkPtrType(object, "GtkIconView")

  w <- .RGtkCall("S_gtk_icon_view_get_orientation", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconViewSetColumns <-
function(object, columns)
{
  checkPtrType(object, "GtkIconView")
  columns <- as.integer(columns)

  w <- .RGtkCall("S_gtk_icon_view_set_columns", object, columns, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconViewGetColumns <-
function(object)
{
  checkPtrType(object, "GtkIconView")

  w <- .RGtkCall("S_gtk_icon_view_get_columns", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconViewSetItemWidth <-
function(object, item.width)
{
  checkPtrType(object, "GtkIconView")
  item.width <- as.integer(item.width)

  w <- .RGtkCall("S_gtk_icon_view_set_item_width", object, item.width, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconViewGetItemWidth <-
function(object)
{
  checkPtrType(object, "GtkIconView")

  w <- .RGtkCall("S_gtk_icon_view_get_item_width", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconViewSetSpacing <-
function(object, spacing)
{
  checkPtrType(object, "GtkIconView")
  spacing <- as.integer(spacing)

  w <- .RGtkCall("S_gtk_icon_view_set_spacing", object, spacing, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconViewGetSpacing <-
function(object)
{
  checkPtrType(object, "GtkIconView")

  w <- .RGtkCall("S_gtk_icon_view_get_spacing", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconViewSetRowSpacing <-
function(object, row.spacing)
{
  checkPtrType(object, "GtkIconView")
  row.spacing <- as.integer(row.spacing)

  w <- .RGtkCall("S_gtk_icon_view_set_row_spacing", object, row.spacing, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconViewGetRowSpacing <-
function(object)
{
  checkPtrType(object, "GtkIconView")

  w <- .RGtkCall("S_gtk_icon_view_get_row_spacing", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconViewSetColumnSpacing <-
function(object, column.spacing)
{
  checkPtrType(object, "GtkIconView")
  column.spacing <- as.integer(column.spacing)

  w <- .RGtkCall("S_gtk_icon_view_set_column_spacing", object, column.spacing, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconViewGetColumnSpacing <-
function(object)
{
  checkPtrType(object, "GtkIconView")

  w <- .RGtkCall("S_gtk_icon_view_get_column_spacing", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconViewSetMargin <-
function(object, margin)
{
  checkPtrType(object, "GtkIconView")
  margin <- as.integer(margin)

  w <- .RGtkCall("S_gtk_icon_view_set_margin", object, margin, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconViewGetMargin <-
function(object)
{
  checkPtrType(object, "GtkIconView")

  w <- .RGtkCall("S_gtk_icon_view_get_margin", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconViewGetPathAtPos <-
function(object, x, y)
{
  checkPtrType(object, "GtkIconView")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_gtk_icon_view_get_path_at_pos", object, x, y, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconViewGetItemAtPos <-
function(object, x, y)
{
  checkPtrType(object, "GtkIconView")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_gtk_icon_view_get_item_at_pos", object, x, y, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconViewGetVisibleRange <-
function(object)
{
  checkPtrType(object, "GtkIconView")

  w <- .RGtkCall("S_gtk_icon_view_get_visible_range", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconViewSelectedForeach <-
function(object, func, data = NULL)
{
  checkPtrType(object, "GtkIconView")
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_icon_view_selected_foreach", object, func, data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconViewSetSelectionMode <-
function(object, mode)
{
  checkPtrType(object, "GtkIconView")
  

  w <- .RGtkCall("S_gtk_icon_view_set_selection_mode", object, mode, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconViewGetSelectionMode <-
function(object)
{
  checkPtrType(object, "GtkIconView")

  w <- .RGtkCall("S_gtk_icon_view_get_selection_mode", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconViewSelectPath <-
function(object, path)
{
  checkPtrType(object, "GtkIconView")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_icon_view_select_path", object, path, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconViewUnselectPath <-
function(object, path)
{
  checkPtrType(object, "GtkIconView")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_icon_view_unselect_path", object, path, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconViewPathIsSelected <-
function(object, path)
{
  checkPtrType(object, "GtkIconView")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_icon_view_path_is_selected", object, path, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconViewGetSelectedItems <-
function(object)
{
  checkPtrType(object, "GtkIconView")

  w <- .RGtkCall("S_gtk_icon_view_get_selected_items", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconViewSelectAll <-
function(object)
{
  checkPtrType(object, "GtkIconView")

  w <- .RGtkCall("S_gtk_icon_view_select_all", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconViewUnselectAll <-
function(object)
{
  checkPtrType(object, "GtkIconView")

  w <- .RGtkCall("S_gtk_icon_view_unselect_all", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconViewItemActivated <-
function(object, path)
{
  checkPtrType(object, "GtkIconView")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_icon_view_item_activated", object, path, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconViewSetCursor <-
function(object, path, cell, start.editing)
{
  checkPtrType(object, "GtkIconView")
  checkPtrType(path, "GtkTreePath")
  checkPtrType(cell, "GtkCellRenderer")
  start.editing <- as.logical(start.editing)

  w <- .RGtkCall("S_gtk_icon_view_set_cursor", object, path, cell, start.editing, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconViewGetCursor <-
function(object)
{
  checkPtrType(object, "GtkIconView")

  w <- .RGtkCall("S_gtk_icon_view_get_cursor", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconViewScrollToPath <-
function(object, path, use.align, row.align, col.align)
{
  checkPtrType(object, "GtkIconView")
  checkPtrType(path, "GtkTreePath")
  use.align <- as.logical(use.align)
  row.align <- as.numeric(row.align)
  col.align <- as.numeric(col.align)

  w <- .RGtkCall("S_gtk_icon_view_scroll_to_path", object, path, use.align, row.align, col.align, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconViewEnableModelDragSource <-
function(object, start.button.mask, targets, actions)
{
  checkPtrType(object, "GtkIconView")
  
  targets <- lapply(targets, function(x) { x <- as.GtkTargetEntry(x); x })
  

  w <- .RGtkCall("S_gtk_icon_view_enable_model_drag_source", object, start.button.mask, targets, actions, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconViewEnableModelDragDest <-
function(object, targets, actions)
{
  checkPtrType(object, "GtkIconView")
  targets <- lapply(targets, function(x) { x <- as.GtkTargetEntry(x); x })
  

  w <- .RGtkCall("S_gtk_icon_view_enable_model_drag_dest", object, targets, actions, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconViewUnsetModelDragSource <-
function(object)
{
  checkPtrType(object, "GtkIconView")

  w <- .RGtkCall("S_gtk_icon_view_unset_model_drag_source", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconViewUnsetModelDragDest <-
function(object)
{
  checkPtrType(object, "GtkIconView")

  w <- .RGtkCall("S_gtk_icon_view_unset_model_drag_dest", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconViewSetReorderable <-
function(object, reorderable)
{
  checkPtrType(object, "GtkIconView")
  reorderable <- as.logical(reorderable)

  w <- .RGtkCall("S_gtk_icon_view_set_reorderable", object, reorderable, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconViewGetReorderable <-
function(object)
{
  checkPtrType(object, "GtkIconView")

  w <- .RGtkCall("S_gtk_icon_view_get_reorderable", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconViewSetDragDestItem <-
function(object, path, pos)
{
  checkPtrType(object, "GtkIconView")
  checkPtrType(path, "GtkTreePath")
  

  w <- .RGtkCall("S_gtk_icon_view_set_drag_dest_item", object, path, pos, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconViewGetDragDestItem <-
function(object)
{
  checkPtrType(object, "GtkIconView")

  w <- .RGtkCall("S_gtk_icon_view_get_drag_dest_item", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconViewGetDestItemAtPos <-
function(object, drag.x, drag.y)
{
  checkPtrType(object, "GtkIconView")
  drag.x <- as.integer(drag.x)
  drag.y <- as.integer(drag.y)

  w <- .RGtkCall("S_gtk_icon_view_get_dest_item_at_pos", object, drag.x, drag.y, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconViewCreateDragIcon <-
function(object, path)
{
  checkPtrType(object, "GtkIconView")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_icon_view_create_drag_icon", object, path, PACKAGE = "RGtk2")

  return(w)
} 


gtkImageGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_image_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkImageNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_image_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkImageNewFromPixmap <-
function(pixmap = NULL, mask = NULL, show = TRUE)
{
  if (!is.null( pixmap )) checkPtrType(pixmap, "GdkPixmap")
  if (!is.null( mask )) checkPtrType(mask, "GdkBitmap")

  w <- .RGtkCall("S_gtk_image_new_from_pixmap", pixmap, mask, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkImageNewFromImage <-
function(image = NULL, mask = NULL, show = TRUE)
{
  if (!is.null( image )) checkPtrType(image, "GdkImage")
  if (!is.null( mask )) checkPtrType(mask, "GdkBitmap")

  w <- .RGtkCall("S_gtk_image_new_from_image", image, mask, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkImageNewFromFile <-
function(filename, show = TRUE)
{
  filename <- as.character(filename)

  w <- .RGtkCall("S_gtk_image_new_from_file", filename, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkImageNewFromPixbuf <-
function(pixbuf = NULL, show = TRUE)
{
  if (!is.null( pixbuf )) checkPtrType(pixbuf, "GdkPixbuf")

  w <- .RGtkCall("S_gtk_image_new_from_pixbuf", pixbuf, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkImageNewFromStock <-
function(stock.id, size, show = TRUE)
{
  stock.id <- as.character(stock.id)
  

  w <- .RGtkCall("S_gtk_image_new_from_stock", stock.id, size, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkImageNewFromIconSet <-
function(icon.set, size, show = TRUE)
{
  checkPtrType(icon.set, "GtkIconSet")
  

  w <- .RGtkCall("S_gtk_image_new_from_icon_set", icon.set, size, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkImageNewFromAnimation <-
function(animation, show = TRUE)
{
  checkPtrType(animation, "GdkPixbufAnimation")

  w <- .RGtkCall("S_gtk_image_new_from_animation", animation, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkImageClear <-
function(object)
{
  checkPtrType(object, "GtkImage")

  w <- .RGtkCall("S_gtk_image_clear", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkImageSetFromPixmap <-
function(object, pixmap, mask = NULL)
{
  checkPtrType(object, "GtkImage")
  checkPtrType(pixmap, "GdkPixmap")
  if (!is.null( mask )) checkPtrType(mask, "GdkBitmap")

  w <- .RGtkCall("S_gtk_image_set_from_pixmap", object, pixmap, mask, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkImageSetFromImage <-
function(object, gdk.image = NULL, mask = NULL)
{
  checkPtrType(object, "GtkImage")
  if (!is.null( gdk.image )) checkPtrType(gdk.image, "GdkImage")
  if (!is.null( mask )) checkPtrType(mask, "GdkBitmap")

  w <- .RGtkCall("S_gtk_image_set_from_image", object, gdk.image, mask, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkImageSetFromFile <-
function(object, filename = NULL)
{
  checkPtrType(object, "GtkImage")
  if (!is.null( filename )) filename <- as.character(filename)

  w <- .RGtkCall("S_gtk_image_set_from_file", object, filename, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkImageSetFromPixbuf <-
function(object, pixbuf = NULL)
{
  checkPtrType(object, "GtkImage")
  if (!is.null( pixbuf )) checkPtrType(pixbuf, "GdkPixbuf")

  w <- .RGtkCall("S_gtk_image_set_from_pixbuf", object, pixbuf, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkImageSetFromStock <-
function(object, stock.id, size)
{
  checkPtrType(object, "GtkImage")
  stock.id <- as.character(stock.id)
  

  w <- .RGtkCall("S_gtk_image_set_from_stock", object, stock.id, size, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkImageSetFromIconSet <-
function(object, icon.set, size)
{
  checkPtrType(object, "GtkImage")
  checkPtrType(icon.set, "GtkIconSet")
  

  w <- .RGtkCall("S_gtk_image_set_from_icon_set", object, icon.set, size, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkImageSetFromAnimation <-
function(object, animation)
{
  checkPtrType(object, "GtkImage")
  checkPtrType(animation, "GdkPixbufAnimation")

  w <- .RGtkCall("S_gtk_image_set_from_animation", object, animation, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkImageGetStorageType <-
function(object)
{
  checkPtrType(object, "GtkImage")

  w <- .RGtkCall("S_gtk_image_get_storage_type", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkImageGetPixmap <-
function(object)
{
  checkPtrType(object, "GtkImage")

  w <- .RGtkCall("S_gtk_image_get_pixmap", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkImageGetImage <-
function(object)
{
  checkPtrType(object, "GtkImage")

  w <- .RGtkCall("S_gtk_image_get_image", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkImageGetPixbuf <-
function(object)
{
  checkPtrType(object, "GtkImage")

  w <- .RGtkCall("S_gtk_image_get_pixbuf", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkImageGetStock <-
function(object)
{
  checkPtrType(object, "GtkImage")

  w <- .RGtkCall("S_gtk_image_get_stock", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkImageGetIconSet <-
function(object)
{
  checkPtrType(object, "GtkImage")

  w <- .RGtkCall("S_gtk_image_get_icon_set", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkImageGetAnimation <-
function(object)
{
  checkPtrType(object, "GtkImage")

  w <- .RGtkCall("S_gtk_image_get_animation", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkImageSet <-
function(object, val, mask)
{
  checkPtrType(object, "GtkImage")
  checkPtrType(val, "GdkImage")
  checkPtrType(mask, "GdkBitmap")

  w <- .RGtkCall("S_gtk_image_set", object, val, mask, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkImageGet <-
function(object)
{
  checkPtrType(object, "GtkImage")

  w <- .RGtkCall("S_gtk_image_get", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkImageNewFromIconName <-
function(icon.name, size)
{
  icon.name <- as.character(icon.name)
  

  w <- .RGtkCall("S_gtk_image_new_from_icon_name", icon.name, size, PACKAGE = "RGtk2")

  return(w)
} 


gtkImageSetFromIconName <-
function(object, icon.name, size)
{
  checkPtrType(object, "GtkImage")
  icon.name <- as.character(icon.name)
  

  w <- .RGtkCall("S_gtk_image_set_from_icon_name", object, icon.name, size, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkImageSetPixelSize <-
function(object, pixel.size)
{
  checkPtrType(object, "GtkImage")
  pixel.size <- as.integer(pixel.size)

  w <- .RGtkCall("S_gtk_image_set_pixel_size", object, pixel.size, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkImageGetIconName <-
function(object)
{
  checkPtrType(object, "GtkImage")

  w <- .RGtkCall("S_gtk_image_get_icon_name", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkImageGetPixelSize <-
function(object)
{
  checkPtrType(object, "GtkImage")

  w <- .RGtkCall("S_gtk_image_get_pixel_size", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkImageMenuItemGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_image_menu_item_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkImageMenuItemNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_image_menu_item_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkImageMenuItemNewWithLabel <-
function(label, show = TRUE)
{
  label <- as.character(label)

  w <- .RGtkCall("S_gtk_image_menu_item_new_with_label", label, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkImageMenuItemNewWithMnemonic <-
function(label, show = TRUE)
{
  label <- as.character(label)

  w <- .RGtkCall("S_gtk_image_menu_item_new_with_mnemonic", label, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkImageMenuItemNewFromStock <-
function(stock.id, accel.group, show = TRUE)
{
  stock.id <- as.character(stock.id)
  checkPtrType(accel.group, "GtkAccelGroup")

  w <- .RGtkCall("S_gtk_image_menu_item_new_from_stock", stock.id, accel.group, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkImageMenuItemSetImage <-
function(object, image = NULL)
{
  checkPtrType(object, "GtkImageMenuItem")
  if (!is.null( image )) checkPtrType(image, "GtkWidget")

  w <- .RGtkCall("S_gtk_image_menu_item_set_image", object, image, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkImageMenuItemGetImage <-
function(object)
{
  checkPtrType(object, "GtkImageMenuItem")

  w <- .RGtkCall("S_gtk_image_menu_item_get_image", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIMContextGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_im_context_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkIMContextSetClientWindow <-
function(object, window)
{
  checkPtrType(object, "GtkIMContext")
  checkPtrType(window, "GdkWindow")

  w <- .RGtkCall("S_gtk_im_context_set_client_window", object, window, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIMContextGetPreeditString <-
function(object)
{
  checkPtrType(object, "GtkIMContext")

  w <- .RGtkCall("S_gtk_im_context_get_preedit_string", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIMContextFilterKeypress <-
function(object, event)
{
  checkPtrType(object, "GtkIMContext")
  checkPtrType(event, "GdkEventKey")

  w <- .RGtkCall("S_gtk_im_context_filter_keypress", object, event, PACKAGE = "RGtk2")

  return(w)
} 


gtkIMContextFocusIn <-
function(object)
{
  checkPtrType(object, "GtkIMContext")

  w <- .RGtkCall("S_gtk_im_context_focus_in", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIMContextFocusOut <-
function(object)
{
  checkPtrType(object, "GtkIMContext")

  w <- .RGtkCall("S_gtk_im_context_focus_out", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIMContextReset <-
function(object)
{
  checkPtrType(object, "GtkIMContext")

  w <- .RGtkCall("S_gtk_im_context_reset", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIMContextSetCursorLocation <-
function(object, area)
{
  checkPtrType(object, "GtkIMContext")
  area <- as.GdkRectangle(area)

  w <- .RGtkCall("S_gtk_im_context_set_cursor_location", object, area, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIMContextSetUsePreedit <-
function(object, use.preedit)
{
  checkPtrType(object, "GtkIMContext")
  use.preedit <- as.logical(use.preedit)

  w <- .RGtkCall("S_gtk_im_context_set_use_preedit", object, use.preedit, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIMContextSetSurrounding <-
function(object, text, len, cursor.index)
{
  checkPtrType(object, "GtkIMContext")
  text <- as.character(text)
  len <- as.integer(len)
  cursor.index <- as.integer(cursor.index)

  w <- .RGtkCall("S_gtk_im_context_set_surrounding", object, text, len, cursor.index, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIMContextGetSurrounding <-
function(object)
{
  checkPtrType(object, "GtkIMContext")

  w <- .RGtkCall("S_gtk_im_context_get_surrounding", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIMContextDeleteSurrounding <-
function(object, offset, n.chars)
{
  checkPtrType(object, "GtkIMContext")
  offset <- as.integer(offset)
  n.chars <- as.integer(n.chars)

  w <- .RGtkCall("S_gtk_im_context_delete_surrounding", object, offset, n.chars, PACKAGE = "RGtk2")

  return(w)
} 


gtkIMContextSimpleGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_im_context_simple_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkIMContextSimpleNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_im_context_simple_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkIMMulticontextGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_im_multicontext_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkIMMulticontextNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_im_multicontext_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkIMMulticontextAppendMenuitems <-
function(object, menushell)
{
  checkPtrType(object, "GtkIMMulticontext")
  checkPtrType(menushell, "GtkMenuShell")

  w <- .RGtkCall("S_gtk_im_multicontext_append_menuitems", object, menushell, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkInputDialogGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_input_dialog_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkInputDialogNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_input_dialog_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkInvisibleGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_invisible_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkInvisibleNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_invisible_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkInvisibleNewForScreen <-
function(screen, show = TRUE)
{
  checkPtrType(screen, "GdkScreen")

  w <- .RGtkCall("S_gtk_invisible_new_for_screen", screen, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkInvisibleSetScreen <-
function(object, screen)
{
  checkPtrType(object, "GtkInvisible")
  checkPtrType(screen, "GdkScreen")

  w <- .RGtkCall("S_gtk_invisible_set_screen", object, screen, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkInvisibleGetScreen <-
function(object)
{
  checkPtrType(object, "GtkInvisible")

  w <- .RGtkCall("S_gtk_invisible_get_screen", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkItemGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_item_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkItemSelect <-
function(object)
{
  checkPtrType(object, "GtkItem")

  w <- .RGtkCall("S_gtk_item_select", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkItemDeselect <-
function(object)
{
  checkPtrType(object, "GtkItem")

  w <- .RGtkCall("S_gtk_item_deselect", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkItemToggle <-
function(object)
{
  checkPtrType(object, "GtkItem")

  w <- .RGtkCall("S_gtk_item_toggle", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkItemFactoryGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_item_factory_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkItemFactoryNew <-
function(container.type, path, accel.group = NULL)
{
  if(getOption("depwarn"))
    .Deprecated("GtkUIManager", "RGtk2")

  container.type <- as.GType(container.type)
  path <- as.character(path)
  if (!is.null( accel.group )) checkPtrType(accel.group, "GtkAccelGroup")

  w <- .RGtkCall("S_gtk_item_factory_new", container.type, path, accel.group, PACKAGE = "RGtk2")

  return(w)
} 


gtkItemFactoryConstruct <-
function(object, container.type, path, accel.group)
{
  checkPtrType(object, "GtkItemFactory")
  container.type <- as.GType(container.type)
  path <- as.character(path)
  checkPtrType(accel.group, "GtkAccelGroup")

  w <- .RGtkCall("S_gtk_item_factory_construct", object, container.type, path, accel.group, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkItemFactoryAddForeign <-
function(accel.widget, full.path, accel.group, keyval, modifiers)
{
  checkPtrType(accel.widget, "GtkWidget")
  full.path <- as.character(full.path)
  checkPtrType(accel.group, "GtkAccelGroup")
  keyval <- as.numeric(keyval)
  

  w <- .RGtkCall("S_gtk_item_factory_add_foreign", accel.widget, full.path, accel.group, keyval, modifiers, PACKAGE = "RGtk2")

  return(w)
} 


gtkItemFactoryFromWidget <-
function(widget)
{
  checkPtrType(widget, "GtkWidget")

  w <- .RGtkCall("S_gtk_item_factory_from_widget", widget, PACKAGE = "RGtk2")

  return(w)
} 


gtkItemFactoryPathFromWidget <-
function(widget)
{
  checkPtrType(widget, "GtkWidget")

  w <- .RGtkCall("S_gtk_item_factory_path_from_widget", widget, PACKAGE = "RGtk2")

  return(w)
} 


gtkItemFactoryGetItem <-
function(object, path)
{
  checkPtrType(object, "GtkItemFactory")
  path <- as.character(path)

  w <- .RGtkCall("S_gtk_item_factory_get_item", object, path, PACKAGE = "RGtk2")

  return(w)
} 


gtkItemFactoryGetWidget <-
function(object, path)
{
  checkPtrType(object, "GtkItemFactory")
  path <- as.character(path)

  w <- .RGtkCall("S_gtk_item_factory_get_widget", object, path, PACKAGE = "RGtk2")

  return(w)
} 


gtkItemFactoryGetWidgetByAction <-
function(object, action)
{
  checkPtrType(object, "GtkItemFactory")
  action <- as.numeric(action)

  w <- .RGtkCall("S_gtk_item_factory_get_widget_by_action", object, action, PACKAGE = "RGtk2")

  return(w)
} 


gtkItemFactoryGetItemByAction <-
function(object, action)
{
  checkPtrType(object, "GtkItemFactory")
  action <- as.numeric(action)

  w <- .RGtkCall("S_gtk_item_factory_get_item_by_action", object, action, PACKAGE = "RGtk2")

  return(w)
} 


gtkItemFactoryCreateItem <-
function(object, entry, callback.data = NULL, callback.type)
{
  checkPtrType(object, "GtkItemFactory")
  entry <- as.GtkItemFactoryEntry(entry)
  
  callback.type <- as.numeric(callback.type)

  w <- .RGtkCall("S_gtk_item_factory_create_item", object, entry, callback.data, callback.type, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkItemFactoryCreateItems <-
function(object, entries, callback.data = NULL)
{
  checkPtrType(object, "GtkItemFactory")
  entries <- lapply(entries, function(x) { x <- as.GtkItemFactoryEntry(x); x })
  

  w <- .RGtkCall("S_gtk_item_factory_create_items", object, entries, callback.data, PACKAGE = "RGtk2")

  return(w)
} 


gtkItemFactoryDeleteItem <-
function(object, path)
{
  checkPtrType(object, "GtkItemFactory")
  path <- as.character(path)

  w <- .RGtkCall("S_gtk_item_factory_delete_item", object, path, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkItemFactoryDeleteEntry <-
function(object, entry)
{
  checkPtrType(object, "GtkItemFactory")
  entry <- as.GtkItemFactoryEntry(entry)

  w <- .RGtkCall("S_gtk_item_factory_delete_entry", object, entry, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkItemFactoryDeleteEntries <-
function(object, entries)
{
  checkPtrType(object, "GtkItemFactory")
  entries <- lapply(entries, function(x) { x <- as.GtkItemFactoryEntry(x); x })

  w <- .RGtkCall("S_gtk_item_factory_delete_entries", object, entries, PACKAGE = "RGtk2")

  return(w)
} 


gtkItemFactoryPopup <-
function(object, x, y, mouse.button, time = "GDK_CURRENT_TIME")
{
  checkPtrType(object, "GtkItemFactory")
  x <- as.numeric(x)
  y <- as.numeric(y)
  mouse.button <- as.numeric(mouse.button)
  time <- as.numeric(time)

  w <- .RGtkCall("S_gtk_item_factory_popup", object, x, y, mouse.button, time, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkItemFactoryPopupWithData <-
function(object, popup.data, x, y, mouse.button, time = "GDK_CURRENT_TIME")
{
  checkPtrType(object, "GtkItemFactory")
  
  x <- as.numeric(x)
  y <- as.numeric(y)
  mouse.button <- as.numeric(mouse.button)
  time <- as.numeric(time)

  w <- .RGtkCall("S_gtk_item_factory_popup_with_data", object, popup.data, x, y, mouse.button, time, PACKAGE = "RGtk2")

  return(w)
} 


gtkItemFactoryPopupData <-
function(object)
{
  checkPtrType(object, "GtkItemFactory")

  w <- .RGtkCall("S_gtk_item_factory_popup_data", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkItemFactoryPopupDataFromWidget <-
function(widget)
{
  checkPtrType(widget, "GtkWidget")

  w <- .RGtkCall("S_gtk_item_factory_popup_data_from_widget", widget, PACKAGE = "RGtk2")

  return(w)
} 


gtkItemFactorySetTranslateFunc <-
function(object, func, data = NULL)
{
  checkPtrType(object, "GtkItemFactory")
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_item_factory_set_translate_func", object, func, data, PACKAGE = "RGtk2")

  return(w)
} 


gtkItemFactoryFromPath <-
function(path)
{
  path <- as.character(path)

  w <- .RGtkCall("S_gtk_item_factory_from_path", path, PACKAGE = "RGtk2")

  return(w)
} 


gtkItemFactoriesPathDelete <-
function(ifactory.path, path)
{
  ifactory.path <- as.character(ifactory.path)
  path <- as.character(path)

  w <- .RGtkCall("S_gtk_item_factories_path_delete", ifactory.path, path, PACKAGE = "RGtk2")

  return(w)
} 


gtkItemFactoryCreateItemsAc <-
function(object, entries, callback.data, callback.type)
{
  checkPtrType(object, "GtkItemFactory")
  entries <- lapply(entries, function(x) { x <- as.GtkItemFactoryEntry(x); x })
  
  callback.type <- as.numeric(callback.type)

  w <- .RGtkCall("S_gtk_item_factory_create_items_ac", object, entries, callback.data, callback.type, PACKAGE = "RGtk2")

  return(w)
} 


gtkLabelGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_label_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkLabelNew <-
function(str = NULL, show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_label_new", str, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkLabelNewWithMnemonic <-
function(str = NULL, show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_label_new_with_mnemonic", str, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkLabelSetText <-
function(object, str)
{
  checkPtrType(object, "GtkLabel")
  str <- as.character(str)

  w <- .RGtkCall("S_gtk_label_set_text", object, str, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkLabelGetText <-
function(object)
{
  checkPtrType(object, "GtkLabel")

  w <- .RGtkCall("S_gtk_label_get_text", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkLabelSetAttributes <-
function(object, attrs)
{
  checkPtrType(object, "GtkLabel")
  checkPtrType(attrs, "PangoAttrList")

  w <- .RGtkCall("S_gtk_label_set_attributes", object, attrs, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkLabelGetAttributes <-
function(object)
{
  checkPtrType(object, "GtkLabel")

  w <- .RGtkCall("S_gtk_label_get_attributes", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkLabelSetLabel <-
function(object, str)
{
  checkPtrType(object, "GtkLabel")
  str <- as.character(str)

  w <- .RGtkCall("S_gtk_label_set_label", object, str, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkLabelGetLabel <-
function(object)
{
  checkPtrType(object, "GtkLabel")

  w <- .RGtkCall("S_gtk_label_get_label", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkLabelSetMarkup <-
function(object, str)
{
  checkPtrType(object, "GtkLabel")
  str <- as.character(str)

  w <- .RGtkCall("S_gtk_label_set_markup", object, str, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkLabelSetUseMarkup <-
function(object, setting)
{
  checkPtrType(object, "GtkLabel")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_label_set_use_markup", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkLabelGetUseMarkup <-
function(object)
{
  checkPtrType(object, "GtkLabel")

  w <- .RGtkCall("S_gtk_label_get_use_markup", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkLabelSetUseUnderline <-
function(object, setting)
{
  checkPtrType(object, "GtkLabel")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_label_set_use_underline", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkLabelGetUseUnderline <-
function(object)
{
  checkPtrType(object, "GtkLabel")

  w <- .RGtkCall("S_gtk_label_get_use_underline", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkLabelSetMarkupWithMnemonic <-
function(object, str)
{
  checkPtrType(object, "GtkLabel")
  str <- as.character(str)

  w <- .RGtkCall("S_gtk_label_set_markup_with_mnemonic", object, str, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkLabelGetMnemonicKeyval <-
function(object)
{
  checkPtrType(object, "GtkLabel")

  w <- .RGtkCall("S_gtk_label_get_mnemonic_keyval", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkLabelSetMnemonicWidget <-
function(object, widget)
{
  checkPtrType(object, "GtkLabel")
  checkPtrType(widget, "GtkWidget")

  w <- .RGtkCall("S_gtk_label_set_mnemonic_widget", object, widget, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkLabelGetMnemonicWidget <-
function(object)
{
  checkPtrType(object, "GtkLabel")

  w <- .RGtkCall("S_gtk_label_get_mnemonic_widget", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkLabelSetTextWithMnemonic <-
function(object, str)
{
  checkPtrType(object, "GtkLabel")
  str <- as.character(str)

  w <- .RGtkCall("S_gtk_label_set_text_with_mnemonic", object, str, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkLabelSetJustify <-
function(object, jtype)
{
  checkPtrType(object, "GtkLabel")
  

  w <- .RGtkCall("S_gtk_label_set_justify", object, jtype, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkLabelGetJustify <-
function(object)
{
  checkPtrType(object, "GtkLabel")

  w <- .RGtkCall("S_gtk_label_get_justify", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkLabelSetPattern <-
function(object, pattern)
{
  checkPtrType(object, "GtkLabel")
  pattern <- as.character(pattern)

  w <- .RGtkCall("S_gtk_label_set_pattern", object, pattern, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkLabelSetLineWrap <-
function(object, wrap)
{
  checkPtrType(object, "GtkLabel")
  wrap <- as.logical(wrap)

  w <- .RGtkCall("S_gtk_label_set_line_wrap", object, wrap, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkLabelGetLineWrap <-
function(object)
{
  checkPtrType(object, "GtkLabel")

  w <- .RGtkCall("S_gtk_label_get_line_wrap", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkLabelSetSelectable <-
function(object, setting)
{
  checkPtrType(object, "GtkLabel")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_label_set_selectable", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkLabelGetSelectable <-
function(object)
{
  checkPtrType(object, "GtkLabel")

  w <- .RGtkCall("S_gtk_label_get_selectable", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkLabelSelectRegion <-
function(object, start.offset, end.offset)
{
  checkPtrType(object, "GtkLabel")
  start.offset <- as.integer(start.offset)
  end.offset <- as.integer(end.offset)

  w <- .RGtkCall("S_gtk_label_select_region", object, start.offset, end.offset, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkLabelGetSelectionBounds <-
function(object)
{
  checkPtrType(object, "GtkLabel")

  w <- .RGtkCall("S_gtk_label_get_selection_bounds", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkLabelGetLayout <-
function(object)
{
  checkPtrType(object, "GtkLabel")

  w <- .RGtkCall("S_gtk_label_get_layout", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkLabelGetLayoutOffsets <-
function(object)
{
  checkPtrType(object, "GtkLabel")

  w <- .RGtkCall("S_gtk_label_get_layout_offsets", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkLabelSet <-
function(object, str)
{
  if(getOption("depwarn"))
    .Deprecated("gtkLabelSetText", "RGtk2")

  checkPtrType(object, "GtkLabel")
  str <- as.character(str)

  w <- .RGtkCall("S_gtk_label_set", object, str, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkLabelGet <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("gtkLabelGetText", "RGtk2")

  checkPtrType(object, "GtkLabel")

  w <- .RGtkCall("S_gtk_label_get", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkLabelParseUline <-
function(object, string)
{
  checkPtrType(object, "GtkLabel")
  string <- as.character(string)

  w <- .RGtkCall("S_gtk_label_parse_uline", object, string, PACKAGE = "RGtk2")

  return(w)
} 


gtkLabelSetEllipsize <-
function(object, mode)
{
  checkPtrType(object, "GtkLabel")
  

  w <- .RGtkCall("S_gtk_label_set_ellipsize", object, mode, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkLabelGetEllipsize <-
function(object)
{
  checkPtrType(object, "GtkLabel")

  w <- .RGtkCall("S_gtk_label_get_ellipsize", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkLabelSetAngle <-
function(object, angle)
{
  checkPtrType(object, "GtkLabel")
  angle <- as.integer(angle)

  w <- .RGtkCall("S_gtk_label_set_angle", object, angle, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkLabelGetAngle <-
function(object)
{
  checkPtrType(object, "GtkLabel")

  w <- .RGtkCall("S_gtk_label_get_angle", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkLabelSetWidthChars <-
function(object, n.chars)
{
  checkPtrType(object, "GtkLabel")
  n.chars <- as.integer(n.chars)

  w <- .RGtkCall("S_gtk_label_set_width_chars", object, n.chars, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkLabelGetWidthChars <-
function(object)
{
  checkPtrType(object, "GtkLabel")

  w <- .RGtkCall("S_gtk_label_get_width_chars", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkLabelSetMaxWidthChars <-
function(object, n.chars)
{
  checkPtrType(object, "GtkLabel")
  n.chars <- as.integer(n.chars)

  w <- .RGtkCall("S_gtk_label_set_max_width_chars", object, n.chars, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkLabelGetMaxWidthChars <-
function(object)
{
  checkPtrType(object, "GtkLabel")

  w <- .RGtkCall("S_gtk_label_get_max_width_chars", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkLabelSetSingleLineMode <-
function(object, single.line.mode)
{
  checkPtrType(object, "GtkLabel")
  single.line.mode <- as.logical(single.line.mode)

  w <- .RGtkCall("S_gtk_label_set_single_line_mode", object, single.line.mode, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkLabelGetSingleLineMode <-
function(object)
{
  checkPtrType(object, "GtkLabel")

  w <- .RGtkCall("S_gtk_label_get_single_line_mode", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkLayoutGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_layout_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkLayoutNew <-
function(hadjustment = NULL, vadjustment = NULL, show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_layout_new", hadjustment, vadjustment, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkLayoutPut <-
function(object, child.widget, x, y)
{
  checkPtrType(object, "GtkLayout")
  checkPtrType(child.widget, "GtkWidget")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_gtk_layout_put", object, child.widget, x, y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkLayoutMove <-
function(object, child.widget, x, y)
{
  checkPtrType(object, "GtkLayout")
  checkPtrType(child.widget, "GtkWidget")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_gtk_layout_move", object, child.widget, x, y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkLayoutSetSize <-
function(object, width, height)
{
  checkPtrType(object, "GtkLayout")
  width <- as.numeric(width)
  height <- as.numeric(height)

  w <- .RGtkCall("S_gtk_layout_set_size", object, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkLayoutGetSize <-
function(object)
{
  checkPtrType(object, "GtkLayout")

  w <- .RGtkCall("S_gtk_layout_get_size", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkLayoutGetHadjustment <-
function(object)
{
  checkPtrType(object, "GtkLayout")

  w <- .RGtkCall("S_gtk_layout_get_hadjustment", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkLayoutGetVadjustment <-
function(object)
{
  checkPtrType(object, "GtkLayout")

  w <- .RGtkCall("S_gtk_layout_get_vadjustment", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkLayoutSetHadjustment <-
function(object, adjustment = NULL)
{
  checkPtrType(object, "GtkLayout")
  if (!is.null( adjustment )) checkPtrType(adjustment, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_layout_set_hadjustment", object, adjustment, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkLayoutSetVadjustment <-
function(object, adjustment = NULL)
{
  checkPtrType(object, "GtkLayout")
  if (!is.null( adjustment )) checkPtrType(adjustment, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_layout_set_vadjustment", object, adjustment, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkLayoutFreeze <-
function(object)
{
  checkPtrType(object, "GtkLayout")

  w <- .RGtkCall("S_gtk_layout_freeze", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkLayoutThaw <-
function(object)
{
  checkPtrType(object, "GtkLayout")

  w <- .RGtkCall("S_gtk_layout_thaw", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkListGetType <-
function()
{
  if(getOption("depwarn"))
    .Deprecated("GtkListStore/GtkTreeView", "RGtk2")

  

  w <- .RGtkCall("S_gtk_list_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkListNew <-
function(show = TRUE)
{
  if(getOption("depwarn"))
    .Deprecated("GtkListStore/GtkTreeView", "RGtk2")

  

  w <- .RGtkCall("S_gtk_list_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkListInsertItems <-
function(object, items, position)
{
  checkPtrType(object, "GtkList")
  items <- as.GList(items)
  position <- as.integer(position)

  w <- .RGtkCall("S_gtk_list_insert_items", object, items, position, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkListAppendItems <-
function(object, items)
{
  checkPtrType(object, "GtkList")
  items <- as.GList(items)

  w <- .RGtkCall("S_gtk_list_append_items", object, items, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkListPrependItems <-
function(object, items)
{
  checkPtrType(object, "GtkList")
  items <- as.GList(items)

  w <- .RGtkCall("S_gtk_list_prepend_items", object, items, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkListRemoveItems <-
function(object, items)
{
  checkPtrType(object, "GtkList")
  items <- as.GList(items)

  w <- .RGtkCall("S_gtk_list_remove_items", object, items, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkListClearItems <-
function(object, start, end)
{
  checkPtrType(object, "GtkList")
  start <- as.integer(start)
  end <- as.integer(end)

  w <- .RGtkCall("S_gtk_list_clear_items", object, start, end, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkListSelectItem <-
function(object, item)
{
  checkPtrType(object, "GtkList")
  item <- as.integer(item)

  w <- .RGtkCall("S_gtk_list_select_item", object, item, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkListUnselectItem <-
function(object, item)
{
  checkPtrType(object, "GtkList")
  item <- as.integer(item)

  w <- .RGtkCall("S_gtk_list_unselect_item", object, item, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkListSelectChild <-
function(object, child)
{
  checkPtrType(object, "GtkList")
  checkPtrType(child, "GtkWidget")

  w <- .RGtkCall("S_gtk_list_select_child", object, child, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkListUnselectChild <-
function(object, child)
{
  checkPtrType(object, "GtkList")
  checkPtrType(child, "GtkWidget")

  w <- .RGtkCall("S_gtk_list_unselect_child", object, child, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkListChildPosition <-
function(object, child)
{
  checkPtrType(object, "GtkList")
  checkPtrType(child, "GtkWidget")

  w <- .RGtkCall("S_gtk_list_child_position", object, child, PACKAGE = "RGtk2")

  return(w)
} 


gtkListSetSelectionMode <-
function(object, mode)
{
  checkPtrType(object, "GtkList")
  

  w <- .RGtkCall("S_gtk_list_set_selection_mode", object, mode, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkListExtendSelection <-
function(object, scroll.type, position, auto.start.selection)
{
  checkPtrType(object, "GtkList")
  
  position <- as.numeric(position)
  auto.start.selection <- as.logical(auto.start.selection)

  w <- .RGtkCall("S_gtk_list_extend_selection", object, scroll.type, position, auto.start.selection, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkListStartSelection <-
function(object)
{
  checkPtrType(object, "GtkList")

  w <- .RGtkCall("S_gtk_list_start_selection", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkListEndSelection <-
function(object)
{
  checkPtrType(object, "GtkList")

  w <- .RGtkCall("S_gtk_list_end_selection", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkListSelectAll <-
function(object)
{
  checkPtrType(object, "GtkList")

  w <- .RGtkCall("S_gtk_list_select_all", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkListUnselectAll <-
function(object)
{
  checkPtrType(object, "GtkList")

  w <- .RGtkCall("S_gtk_list_unselect_all", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkListScrollHorizontal <-
function(object, scroll.type, position)
{
  checkPtrType(object, "GtkList")
  
  position <- as.numeric(position)

  w <- .RGtkCall("S_gtk_list_scroll_horizontal", object, scroll.type, position, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkListScrollVertical <-
function(object, scroll.type, position)
{
  checkPtrType(object, "GtkList")
  
  position <- as.numeric(position)

  w <- .RGtkCall("S_gtk_list_scroll_vertical", object, scroll.type, position, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkListToggleAddMode <-
function(object)
{
  checkPtrType(object, "GtkList")

  w <- .RGtkCall("S_gtk_list_toggle_add_mode", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkListToggleFocusRow <-
function(object)
{
  checkPtrType(object, "GtkList")

  w <- .RGtkCall("S_gtk_list_toggle_focus_row", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkListToggleRow <-
function(object, item)
{
  checkPtrType(object, "GtkList")
  checkPtrType(item, "GtkWidget")

  w <- .RGtkCall("S_gtk_list_toggle_row", object, item, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkListUndoSelection <-
function(object)
{
  checkPtrType(object, "GtkList")

  w <- .RGtkCall("S_gtk_list_undo_selection", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkListEndDragSelection <-
function(object)
{
  checkPtrType(object, "GtkList")

  w <- .RGtkCall("S_gtk_list_end_drag_selection", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkListItemGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_list_item_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkListItemNew <-
function(show = TRUE)
{
  if(getOption("depwarn"))
    .Deprecated("GtkTreeView", "RGtk2")

  

  w <- .RGtkCall("S_gtk_list_item_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkListItemNewWithLabel <-
function(label, show = TRUE)
{
  if(getOption("depwarn"))
    .Deprecated("GtkTreeView", "RGtk2")

  label <- as.character(label)

  w <- .RGtkCall("S_gtk_list_item_new_with_label", label, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkListItemSelect <-
function(object)
{
  checkPtrType(object, "GtkListItem")

  w <- .RGtkCall("S_gtk_list_item_select", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkListItemDeselect <-
function(object)
{
  checkPtrType(object, "GtkListItem")

  w <- .RGtkCall("S_gtk_list_item_deselect", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkListStoreGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_list_store_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkListStoreNewv <-
function(value)
{
  value <- as.list(as.GType(value))

  w <- .RGtkCall("S_gtk_list_store_newv", value, PACKAGE = "RGtk2")

  return(w)
} 


gtkListStoreSetColumnTypes <-
function(object, types)
{
  checkPtrType(object, "GtkListStore")
  types <- as.list(as.GType(types))

  w <- .RGtkCall("S_gtk_list_store_set_column_types", object, types, PACKAGE = "RGtk2")

  return(w)
} 


gtkListStoreSetValue <-
function(object, iter, column, value)
{
  checkPtrType(object, "GtkListStore")
  checkPtrType(iter, "GtkTreeIter")
  column <- as.integer(column)
  

  w <- .RGtkCall("S_gtk_list_store_set_value", object, iter, column, value, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkListStoreRemove <-
function(object, iter)
{
  checkPtrType(object, "GtkListStore")
  checkPtrType(iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_list_store_remove", object, iter, PACKAGE = "RGtk2")

  return(w)
} 


gtkListStoreInsert <-
function(object, position)
{
  checkPtrType(object, "GtkListStore")
  position <- as.integer(position)

  w <- .RGtkCall("S_gtk_list_store_insert", object, position, PACKAGE = "RGtk2")

  return(w)
} 


gtkListStoreInsertBefore <-
function(object, sibling)
{
  checkPtrType(object, "GtkListStore")
  checkPtrType(sibling, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_list_store_insert_before", object, sibling, PACKAGE = "RGtk2")

  return(w)
} 


gtkListStoreInsertAfter <-
function(object, sibling)
{
  checkPtrType(object, "GtkListStore")
  checkPtrType(sibling, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_list_store_insert_after", object, sibling, PACKAGE = "RGtk2")

  return(w)
} 


gtkListStorePrepend <-
function(object, iter)
{
  checkPtrType(object, "GtkListStore")
  checkPtrType(iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_list_store_prepend", object, iter, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkListStoreAppend <-
function(object)
{
  checkPtrType(object, "GtkListStore")

  w <- .RGtkCall("S_gtk_list_store_append", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkListStoreClear <-
function(object)
{
  checkPtrType(object, "GtkListStore")

  w <- .RGtkCall("S_gtk_list_store_clear", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkListStoreIterIsValid <-
function(object, iter)
{
  checkPtrType(object, "GtkListStore")
  checkPtrType(iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_list_store_iter_is_valid", object, iter, PACKAGE = "RGtk2")

  return(w)
} 


gtkListStoreReorder <-
function(object, new.order)
{
  checkPtrType(object, "GtkListStore")
  new.order <- as.list(as.integer(new.order))

  w <- .RGtkCall("S_gtk_list_store_reorder", object, new.order, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkListStoreSwap <-
function(object, a, b)
{
  checkPtrType(object, "GtkListStore")
  checkPtrType(a, "GtkTreeIter")
  checkPtrType(b, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_list_store_swap", object, a, b, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkListStoreMoveAfter <-
function(object, iter, position = NULL)
{
  checkPtrType(object, "GtkListStore")
  checkPtrType(iter, "GtkTreeIter")
  if (!is.null( position )) checkPtrType(position, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_list_store_move_after", object, iter, position, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkListStoreMoveBefore <-
function(object, iter, position = NULL)
{
  checkPtrType(object, "GtkListStore")
  checkPtrType(iter, "GtkTreeIter")
  if (!is.null( position )) checkPtrType(position, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_list_store_move_before", object, iter, position, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCheckVersion <-
function(required.major, required.minor, required.micro)
{
  required.major <- as.numeric(required.major)
  required.minor <- as.numeric(required.minor)
  required.micro <- as.numeric(required.micro)

  w <- .RGtkCall("S_gtk_check_version", required.major, required.minor, required.micro, PACKAGE = "RGtk2")

  return(w)
} 


gtkExit <-
function(error.code)
{
  error.code <- as.integer(error.code)

  w <- .RGtkCall("S_gtk_exit", error.code, PACKAGE = "RGtk2")

  return(w)
} 


gtkGetDefaultLanguage <-
function()
{
  

  w <- .RGtkCall("S_gtk_get_default_language", PACKAGE = "RGtk2")

  return(w)
} 


gtkEventsPending <-
function()
{
  

  w <- .RGtkCall("S_gtk_events_pending", PACKAGE = "RGtk2")

  return(w)
} 


gtkMainDoEvent <-
function(event)
{
  checkPtrType(event, "GdkEvent")

  w <- .RGtkCall("S_gtk_main_do_event", event, PACKAGE = "RGtk2")

  return(w)
} 


gtkMain <-
function()
{
  

  w <- .RGtkCall("S_gtk_main", PACKAGE = "RGtk2")

  return(w)
} 


gtkMainLevel <-
function()
{
  

  w <- .RGtkCall("S_gtk_main_level", PACKAGE = "RGtk2")

  return(w)
} 


gtkMainQuit <-
function()
{
  

  w <- .RGtkCall("S_gtk_main_quit", PACKAGE = "RGtk2")

  return(w)
} 


gtkMainIteration <-
function()
{
  

  w <- .RGtkCall("S_gtk_main_iteration", PACKAGE = "RGtk2")

  return(w)
} 


gtkMainIterationDo <-
function(blocking = TRUE)
{
  blocking <- as.logical(blocking)

  w <- .RGtkCall("S_gtk_main_iteration_do", blocking, PACKAGE = "RGtk2")

  return(w)
} 


gtkTrue <-
function()
{
  

  w <- .RGtkCall("S_gtk_true", PACKAGE = "RGtk2")

  return(w)
} 


gtkFalse <-
function()
{
  

  w <- .RGtkCall("S_gtk_false", PACKAGE = "RGtk2")

  return(w)
} 


gtkGrabAdd <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_grab_add", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkGrabGetCurrent <-
function()
{
  

  w <- .RGtkCall("S_gtk_grab_get_current", PACKAGE = "RGtk2")

  return(w)
} 


gtkGrabRemove <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_grab_remove", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkInitAdd <-
function(fun, data = NULL)
{
  fun <- as.function(fun)
  

  w <- .RGtkCall("S_gtk_init_add", fun, data, PACKAGE = "RGtk2")

  return(w)
} 


gtkQuitAddDestroy <-
function(main.level, object)
{
  main.level <- as.numeric(main.level)
  checkPtrType(object, "GtkObject")

  w <- .RGtkCall("S_gtk_quit_add_destroy", main.level, object, PACKAGE = "RGtk2")

  return(w)
} 


gtkQuitAdd <-
function(main.level, fun, data = NULL)
{
  main.level <- as.numeric(main.level)
  fun <- as.function(fun)
  

  w <- .RGtkCall("S_gtk_quit_add", main.level, fun, data, PACKAGE = "RGtk2")

  return(w)
} 


gtkQuitAddFull <-
function(main.level, fun, data = NULL)
{
  main.level <- as.numeric(main.level)
  fun <- as.function(fun)
  

  w <- .RGtkCall("S_gtk_quit_add_full", main.level, fun, data, PACKAGE = "RGtk2")

  return(w)
} 


gtkQuitRemove <-
function(quit.handler.id)
{
  quit.handler.id <- as.numeric(quit.handler.id)

  w <- .RGtkCall("S_gtk_quit_remove", quit.handler.id, PACKAGE = "RGtk2")

  return(w)
} 


gtkQuitRemoveByData <-
function(data)
{
  

  w <- .RGtkCall("S_gtk_quit_remove_by_data", data, PACKAGE = "RGtk2")

  return(w)
} 


gtkTimeoutAdd <-
function(interval, fun, data = NULL)
{
  interval <- as.numeric(interval)
  fun <- as.function(fun)
  

  w <- .RGtkCall("S_gtk_timeout_add", interval, fun, data, PACKAGE = "RGtk2")

  return(w)
} 


gtkTimeoutAddFull <-
function(interval, fun, data = NULL)
{
  interval <- as.numeric(interval)
  fun <- as.function(fun)
  

  w <- .RGtkCall("S_gtk_timeout_add_full", interval, fun, data, PACKAGE = "RGtk2")

  return(w)
} 


gtkTimeoutRemove <-
function(timeout.handler.id)
{
  timeout.handler.id <- as.numeric(timeout.handler.id)

  w <- .RGtkCall("S_gtk_timeout_remove", timeout.handler.id, PACKAGE = "RGtk2")

  return(w)
} 


gtkIdleAdd <-
function(fun, data = NULL)
{
  fun <- as.function(fun)
  

  w <- .RGtkCall("S_gtk_idle_add", fun, data, PACKAGE = "RGtk2")

  return(w)
} 


gtkIdleAddPriority <-
function(priority, fun, data = NULL)
{
  priority <- as.integer(priority)
  fun <- as.function(fun)
  

  w <- .RGtkCall("S_gtk_idle_add_priority", priority, fun, data, PACKAGE = "RGtk2")

  return(w)
} 


gtkIdleAddFull <-
function(priority, fun, data = NULL)
{
  priority <- as.integer(priority)
  fun <- as.function(fun)
  

  w <- .RGtkCall("S_gtk_idle_add_full", priority, fun, data, PACKAGE = "RGtk2")

  return(w)
} 


gtkIdleRemove <-
function(idle.handler.id)
{
  idle.handler.id <- as.numeric(idle.handler.id)

  w <- .RGtkCall("S_gtk_idle_remove", idle.handler.id, PACKAGE = "RGtk2")

  return(w)
} 


gtkIdleRemoveByData <-
function(data)
{
  

  w <- .RGtkCall("S_gtk_idle_remove_by_data", data, PACKAGE = "RGtk2")

  return(w)
} 


gtkInputRemove <-
function(input.handler.id)
{
  input.handler.id <- as.numeric(input.handler.id)

  w <- .RGtkCall("S_gtk_input_remove", input.handler.id, PACKAGE = "RGtk2")

  return(w)
} 


gtkKeySnooperInstall <-
function(snooper, func.data = NULL)
{
  snooper <- as.function(snooper)
  

  w <- .RGtkCall("S_gtk_key_snooper_install", snooper, func.data, PACKAGE = "RGtk2")

  return(w)
} 


gtkKeySnooperRemove <-
function(snooper.handler.id)
{
  snooper.handler.id <- as.numeric(snooper.handler.id)

  w <- .RGtkCall("S_gtk_key_snooper_remove", snooper.handler.id, PACKAGE = "RGtk2")

  return(w)
} 


gtkGetCurrentEvent <-
function()
{
  

  w <- .RGtkCall("S_gtk_get_current_event", PACKAGE = "RGtk2")

  return(w)
} 


gtkGetCurrentEventTime <-
function()
{
  

  w <- .RGtkCall("S_gtk_get_current_event_time", PACKAGE = "RGtk2")

  return(w)
} 


gtkGetCurrentEventState <-
function()
{
  

  w <- .RGtkCall("S_gtk_get_current_event_state", PACKAGE = "RGtk2")

  return(w)
} 


gtkGetEventWidget <-
function(event)
{
  checkPtrType(event, "GdkEvent")

  w <- .RGtkCall("S_gtk_get_event_widget", event, PACKAGE = "RGtk2")

  return(w)
} 


gtkPropagateEvent <-
function(object, event)
{
  checkPtrType(object, "GtkWidget")
  checkPtrType(event, "GdkEvent")

  w <- .RGtkCall("S_gtk_propagate_event", object, event, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_menu_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkMenuNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_menu_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkMenuPopup <-
function(object, parent.menu.shell = NULL, parent.menu.item = NULL, func = NULL, data = NULL, button, activate.time)
{
  checkPtrType(object, "GtkMenu")
  if (!is.null( parent.menu.shell )) checkPtrType(parent.menu.shell, "GtkWidget")
  if (!is.null( parent.menu.item )) checkPtrType(parent.menu.item, "GtkWidget")
  if (!is.null( func )) func <- as.function(func)
  
  button <- as.numeric(button)
  activate.time <- as.numeric(activate.time)

  w <- .RGtkCall("S_gtk_menu_popup", object, parent.menu.shell, parent.menu.item, func, data, button, activate.time, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuReposition <-
function(object)
{
  checkPtrType(object, "GtkMenu")

  w <- .RGtkCall("S_gtk_menu_reposition", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuPopdown <-
function(object)
{
  checkPtrType(object, "GtkMenu")

  w <- .RGtkCall("S_gtk_menu_popdown", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuGetActive <-
function(object)
{
  checkPtrType(object, "GtkMenu")

  w <- .RGtkCall("S_gtk_menu_get_active", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkMenuSetActive <-
function(object, index)
{
  checkPtrType(object, "GtkMenu")
  index <- as.numeric(index)

  w <- .RGtkCall("S_gtk_menu_set_active", object, index, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuSetAccelGroup <-
function(object, accel.group)
{
  checkPtrType(object, "GtkMenu")
  checkPtrType(accel.group, "GtkAccelGroup")

  w <- .RGtkCall("S_gtk_menu_set_accel_group", object, accel.group, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuGetAccelGroup <-
function(object)
{
  checkPtrType(object, "GtkMenu")

  w <- .RGtkCall("S_gtk_menu_get_accel_group", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkMenuSetAccelPath <-
function(object, accel.path)
{
  checkPtrType(object, "GtkMenu")
  accel.path <- as.character(accel.path)

  w <- .RGtkCall("S_gtk_menu_set_accel_path", object, accel.path, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuDetach <-
function(object)
{
  checkPtrType(object, "GtkMenu")

  w <- .RGtkCall("S_gtk_menu_detach", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuGetAttachWidget <-
function(object)
{
  checkPtrType(object, "GtkMenu")

  w <- .RGtkCall("S_gtk_menu_get_attach_widget", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkMenuSetTearoffState <-
function(object, torn.off)
{
  checkPtrType(object, "GtkMenu")
  torn.off <- as.logical(torn.off)

  w <- .RGtkCall("S_gtk_menu_set_tearoff_state", object, torn.off, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuGetTearoffState <-
function(object)
{
  checkPtrType(object, "GtkMenu")

  w <- .RGtkCall("S_gtk_menu_get_tearoff_state", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkMenuSetTitle <-
function(object, title)
{
  checkPtrType(object, "GtkMenu")
  title <- as.character(title)

  w <- .RGtkCall("S_gtk_menu_set_title", object, title, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuGetTitle <-
function(object)
{
  checkPtrType(object, "GtkMenu")

  w <- .RGtkCall("S_gtk_menu_get_title", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkMenuReorderChild <-
function(object, child, position)
{
  checkPtrType(object, "GtkMenu")
  checkPtrType(child, "GtkWidget")
  position <- as.integer(position)

  w <- .RGtkCall("S_gtk_menu_reorder_child", object, child, position, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuSetScreen <-
function(object, screen = NULL)
{
  checkPtrType(object, "GtkMenu")
  if (!is.null( screen )) checkPtrType(screen, "GdkScreen")

  w <- .RGtkCall("S_gtk_menu_set_screen", object, screen, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuAttach <-
function(object, child, left.attach, right.attach, top.attach, bottom.attach)
{
  checkPtrType(object, "GtkMenu")
  checkPtrType(child, "GtkWidget")
  left.attach <- as.numeric(left.attach)
  right.attach <- as.numeric(right.attach)
  top.attach <- as.numeric(top.attach)
  bottom.attach <- as.numeric(bottom.attach)

  w <- .RGtkCall("S_gtk_menu_attach", object, child, left.attach, right.attach, top.attach, bottom.attach, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuSetMonitor <-
function(object, monitor.num)
{
  checkPtrType(object, "GtkMenu")
  monitor.num <- as.integer(monitor.num)

  w <- .RGtkCall("S_gtk_menu_set_monitor", object, monitor.num, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuGetForAttachWidget <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_menu_get_for_attach_widget", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkMenuBarGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_menu_bar_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkMenuBarNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_menu_bar_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkMenuBarGetPackDirection <-
function(object)
{
  checkPtrType(object, "GtkMenuBar")

  w <- .RGtkCall("S_gtk_menu_bar_get_pack_direction", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkMenuBarSetPackDirection <-
function(object, pack.dir)
{
  checkPtrType(object, "GtkMenuBar")
  

  w <- .RGtkCall("S_gtk_menu_bar_set_pack_direction", object, pack.dir, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuBarGetChildPackDirection <-
function(object)
{
  checkPtrType(object, "GtkMenuBar")

  w <- .RGtkCall("S_gtk_menu_bar_get_child_pack_direction", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkMenuBarSetChildPackDirection <-
function(object, child.pack.dir)
{
  checkPtrType(object, "GtkMenuBar")
  

  w <- .RGtkCall("S_gtk_menu_bar_set_child_pack_direction", object, child.pack.dir, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuItemGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_menu_item_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkMenuItemNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_menu_item_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkMenuItemNewWithLabel <-
function(label, show = TRUE)
{
  label <- as.character(label)

  w <- .RGtkCall("S_gtk_menu_item_new_with_label", label, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkMenuItemNewWithMnemonic <-
function(label, show = TRUE)
{
  label <- as.character(label)

  w <- .RGtkCall("S_gtk_menu_item_new_with_mnemonic", label, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkMenuItemSetSubmenu <-
function(object, submenu)
{
  checkPtrType(object, "GtkMenuItem")
  checkPtrType(submenu, "GtkWidget")

  w <- .RGtkCall("S_gtk_menu_item_set_submenu", object, submenu, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuItemGetSubmenu <-
function(object)
{
  checkPtrType(object, "GtkMenuItem")

  w <- .RGtkCall("S_gtk_menu_item_get_submenu", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkMenuItemRemoveSubmenu <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("setSubmenu", "RGtk2")

  checkPtrType(object, "GtkMenuItem")

  w <- .RGtkCall("S_gtk_menu_item_remove_submenu", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuItemSelect <-
function(object)
{
  checkPtrType(object, "GtkMenuItem")

  w <- .RGtkCall("S_gtk_menu_item_select", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuItemDeselect <-
function(object)
{
  checkPtrType(object, "GtkMenuItem")

  w <- .RGtkCall("S_gtk_menu_item_deselect", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuItemActivate <-
function(object)
{
  checkPtrType(object, "GtkMenuItem")

  w <- .RGtkCall("S_gtk_menu_item_activate", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuItemToggleSizeRequest <-
function(object, requisition)
{
  checkPtrType(object, "GtkMenuItem")
  requisition <- as.list(as.integer(requisition))

  w <- .RGtkCall("S_gtk_menu_item_toggle_size_request", object, requisition, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuItemToggleSizeAllocate <-
function(object, allocation)
{
  checkPtrType(object, "GtkMenuItem")
  allocation <- as.integer(allocation)

  w <- .RGtkCall("S_gtk_menu_item_toggle_size_allocate", object, allocation, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuItemSetRightJustified <-
function(object, right.justified)
{
  checkPtrType(object, "GtkMenuItem")
  right.justified <- as.logical(right.justified)

  w <- .RGtkCall("S_gtk_menu_item_set_right_justified", object, right.justified, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuItemGetRightJustified <-
function(object)
{
  checkPtrType(object, "GtkMenuItem")

  w <- .RGtkCall("S_gtk_menu_item_get_right_justified", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkMenuItemSetAccelPath <-
function(object, accel.path)
{
  checkPtrType(object, "GtkMenuItem")
  accel.path <- as.character(accel.path)

  w <- .RGtkCall("S_gtk_menu_item_set_accel_path", object, accel.path, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuItemRightJustify <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("gtkMenuItemSetRightJustified", "RGtk2")

  checkPtrType(object, "GtkMenuItem")

  w <- .RGtkCall("S_gtk_menu_item_right_justify", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuShellGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_menu_shell_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkMenuShellAppend <-
function(object, child)
{
  checkPtrType(object, "GtkMenuShell")
  checkPtrType(child, "GtkWidget")

  w <- .RGtkCall("S_gtk_menu_shell_append", object, child, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuShellPrepend <-
function(object, child)
{
  checkPtrType(object, "GtkMenuShell")
  checkPtrType(child, "GtkWidget")

  w <- .RGtkCall("S_gtk_menu_shell_prepend", object, child, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuShellInsert <-
function(object, child, position)
{
  checkPtrType(object, "GtkMenuShell")
  checkPtrType(child, "GtkWidget")
  position <- as.integer(position)

  w <- .RGtkCall("S_gtk_menu_shell_insert", object, child, position, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuShellDeactivate <-
function(object)
{
  checkPtrType(object, "GtkMenuShell")

  w <- .RGtkCall("S_gtk_menu_shell_deactivate", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuShellSelectItem <-
function(object, menu.item)
{
  checkPtrType(object, "GtkMenuShell")
  checkPtrType(menu.item, "GtkWidget")

  w <- .RGtkCall("S_gtk_menu_shell_select_item", object, menu.item, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuShellDeselect <-
function(object)
{
  checkPtrType(object, "GtkMenuShell")

  w <- .RGtkCall("S_gtk_menu_shell_deselect", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuShellActivateItem <-
function(object, menu.item, force.deactivate)
{
  checkPtrType(object, "GtkMenuShell")
  checkPtrType(menu.item, "GtkWidget")
  force.deactivate <- as.logical(force.deactivate)

  w <- .RGtkCall("S_gtk_menu_shell_activate_item", object, menu.item, force.deactivate, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuShellSelectFirst <-
function(object, search.sensitive)
{
  checkPtrType(object, "GtkMenuShell")
  search.sensitive <- as.logical(search.sensitive)

  w <- .RGtkCall("S_gtk_menu_shell_select_first", object, search.sensitive, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuShellCancel <-
function(object)
{
  checkPtrType(object, "GtkMenuShell")

  w <- .RGtkCall("S_gtk_menu_shell_cancel", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuShellGetTakeFocus <-
function(object)
{
  checkPtrType(object, "GtkMenuShell")

  w <- .RGtkCall("S_gtk_menu_shell_get_take_focus", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkMenuShellSetTakeFocus <-
function(object, take.focus)
{
  checkPtrType(object, "GtkMenuShell")
  take.focus <- as.logical(take.focus)

  w <- .RGtkCall("S_gtk_menu_shell_set_take_focus", object, take.focus, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuToolButtonGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_menu_tool_button_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkMenuToolButtonNew <-
function(icon.widget, label, show = TRUE)
{
  checkPtrType(icon.widget, "GtkWidget")
  label <- as.character(label)

  w <- .RGtkCall("S_gtk_menu_tool_button_new", icon.widget, label, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkMenuToolButtonNewFromStock <-
function(stock.id)
{
  stock.id <- as.character(stock.id)

  w <- .RGtkCall("S_gtk_menu_tool_button_new_from_stock", stock.id, PACKAGE = "RGtk2")

  return(w)
} 


gtkMenuToolButtonSetMenu <-
function(object, menu)
{
  checkPtrType(object, "GtkMenuToolButton")
  checkPtrType(menu, "GtkWidget")

  w <- .RGtkCall("S_gtk_menu_tool_button_set_menu", object, menu, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuToolButtonGetMenu <-
function(object)
{
  checkPtrType(object, "GtkMenuToolButton")

  w <- .RGtkCall("S_gtk_menu_tool_button_get_menu", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkMenuToolButtonSetArrowTooltip <-
function(object, tooltips, tip.text = NULL, tip.private = NULL)
{
  if(getOption("depwarn"))
    .Deprecated("setArrowTooltipText", "RGtk2")

  checkPtrType(object, "GtkMenuToolButton")
  checkPtrType(tooltips, "GtkTooltips")
  if (!is.null( tip.text )) tip.text <- as.character(tip.text)
  if (!is.null( tip.private )) tip.private <- as.character(tip.private)

  w <- .RGtkCall("S_gtk_menu_tool_button_set_arrow_tooltip", object, tooltips, tip.text, tip.private, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMessageDialogGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_message_dialog_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkMessageDialogSetMarkup <-
function(object, str)
{
  checkPtrType(object, "GtkMessageDialog")
  str <- as.character(str)

  w <- .RGtkCall("S_gtk_message_dialog_set_markup", object, str, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMiscGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_misc_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkMiscSetAlignment <-
function(object, xalign, yalign)
{
  checkPtrType(object, "GtkMisc")
  xalign <- as.numeric(xalign)
  yalign <- as.numeric(yalign)

  w <- .RGtkCall("S_gtk_misc_set_alignment", object, xalign, yalign, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMiscGetAlignment <-
function(object)
{
  checkPtrType(object, "GtkMisc")

  w <- .RGtkCall("S_gtk_misc_get_alignment", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMiscSetPadding <-
function(object, xpad, ypad)
{
  checkPtrType(object, "GtkMisc")
  xpad <- as.integer(xpad)
  ypad <- as.integer(ypad)

  w <- .RGtkCall("S_gtk_misc_set_padding", object, xpad, ypad, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMiscGetPadding <-
function(object)
{
  checkPtrType(object, "GtkMisc")

  w <- .RGtkCall("S_gtk_misc_get_padding", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkNotebookGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_notebook_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkNotebookNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_notebook_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkNotebookAppendPage <-
function(object, child, tab.label = NULL)
{
  checkPtrType(object, "GtkNotebook")
  checkPtrType(child, "GtkWidget")
  if (!is.null( tab.label )) checkPtrType(tab.label, "GtkWidget")

  w <- .RGtkCall("S_gtk_notebook_append_page", object, child, tab.label, PACKAGE = "RGtk2")

  return(w)
} 


gtkNotebookAppendPageMenu <-
function(object, child, tab.label = NULL, menu.label = NULL)
{
  checkPtrType(object, "GtkNotebook")
  checkPtrType(child, "GtkWidget")
  if (!is.null( tab.label )) checkPtrType(tab.label, "GtkWidget")
  if (!is.null( menu.label )) checkPtrType(menu.label, "GtkWidget")

  w <- .RGtkCall("S_gtk_notebook_append_page_menu", object, child, tab.label, menu.label, PACKAGE = "RGtk2")

  return(w)
} 


gtkNotebookPrependPage <-
function(object, child, tab.label = NULL)
{
  checkPtrType(object, "GtkNotebook")
  checkPtrType(child, "GtkWidget")
  if (!is.null( tab.label )) checkPtrType(tab.label, "GtkWidget")

  w <- .RGtkCall("S_gtk_notebook_prepend_page", object, child, tab.label, PACKAGE = "RGtk2")

  return(w)
} 


gtkNotebookPrependPageMenu <-
function(object, child, tab.label = NULL, menu.label = NULL)
{
  checkPtrType(object, "GtkNotebook")
  checkPtrType(child, "GtkWidget")
  if (!is.null( tab.label )) checkPtrType(tab.label, "GtkWidget")
  if (!is.null( menu.label )) checkPtrType(menu.label, "GtkWidget")

  w <- .RGtkCall("S_gtk_notebook_prepend_page_menu", object, child, tab.label, menu.label, PACKAGE = "RGtk2")

  return(w)
} 


gtkNotebookInsertPage <-
function(object, child, tab.label = NULL, position = -1)
{
  checkPtrType(object, "GtkNotebook")
  checkPtrType(child, "GtkWidget")
  if (!is.null( tab.label )) checkPtrType(tab.label, "GtkWidget")
  position <- as.integer(position)

  w <- .RGtkCall("S_gtk_notebook_insert_page", object, child, tab.label, position, PACKAGE = "RGtk2")

  return(w)
} 


gtkNotebookInsertPageMenu <-
function(object, child, tab.label = NULL, menu.label = NULL, position = -1)
{
  checkPtrType(object, "GtkNotebook")
  checkPtrType(child, "GtkWidget")
  if (!is.null( tab.label )) checkPtrType(tab.label, "GtkWidget")
  if (!is.null( menu.label )) checkPtrType(menu.label, "GtkWidget")
  position <- as.integer(position)

  w <- .RGtkCall("S_gtk_notebook_insert_page_menu", object, child, tab.label, menu.label, position, PACKAGE = "RGtk2")

  return(w)
} 


gtkNotebookRemovePage <-
function(object, page.num)
{
  checkPtrType(object, "GtkNotebook")
  page.num <- as.integer(page.num)

  w <- .RGtkCall("S_gtk_notebook_remove_page", object, page.num, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkNotebookGetCurrentPage <-
function(object)
{
  checkPtrType(object, "GtkNotebook")

  w <- .RGtkCall("S_gtk_notebook_get_current_page", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkNotebookGetNthPage <-
function(object, page.num)
{
  checkPtrType(object, "GtkNotebook")
  page.num <- as.integer(page.num)

  w <- .RGtkCall("S_gtk_notebook_get_nth_page", object, page.num, PACKAGE = "RGtk2")

  return(w)
} 


gtkNotebookGetNPages <-
function(object)
{
  checkPtrType(object, "GtkNotebook")

  w <- .RGtkCall("S_gtk_notebook_get_n_pages", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkNotebookPageNum <-
function(object, child)
{
  checkPtrType(object, "GtkNotebook")
  checkPtrType(child, "GtkWidget")

  w <- .RGtkCall("S_gtk_notebook_page_num", object, child, PACKAGE = "RGtk2")

  return(w)
} 


gtkNotebookSetCurrentPage <-
function(object, page.num)
{
  checkPtrType(object, "GtkNotebook")
  page.num <- as.integer(page.num)

  w <- .RGtkCall("S_gtk_notebook_set_current_page", object, page.num, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkNotebookNextPage <-
function(object)
{
  checkPtrType(object, "GtkNotebook")

  w <- .RGtkCall("S_gtk_notebook_next_page", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkNotebookPrevPage <-
function(object)
{
  checkPtrType(object, "GtkNotebook")

  w <- .RGtkCall("S_gtk_notebook_prev_page", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkNotebookSetShowBorder <-
function(object, show.border)
{
  checkPtrType(object, "GtkNotebook")
  show.border <- as.logical(show.border)

  w <- .RGtkCall("S_gtk_notebook_set_show_border", object, show.border, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkNotebookGetShowBorder <-
function(object)
{
  checkPtrType(object, "GtkNotebook")

  w <- .RGtkCall("S_gtk_notebook_get_show_border", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkNotebookSetShowTabs <-
function(object, show.tabs)
{
  checkPtrType(object, "GtkNotebook")
  show.tabs <- as.logical(show.tabs)

  w <- .RGtkCall("S_gtk_notebook_set_show_tabs", object, show.tabs, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkNotebookGetShowTabs <-
function(object)
{
  checkPtrType(object, "GtkNotebook")

  w <- .RGtkCall("S_gtk_notebook_get_show_tabs", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkNotebookSetTabPos <-
function(object, pos)
{
  checkPtrType(object, "GtkNotebook")
  

  w <- .RGtkCall("S_gtk_notebook_set_tab_pos", object, pos, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkNotebookGetTabPos <-
function(object)
{
  checkPtrType(object, "GtkNotebook")

  w <- .RGtkCall("S_gtk_notebook_get_tab_pos", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkNotebookSetHomogeneousTabs <-
function(object, homogeneous)
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  checkPtrType(object, "GtkNotebook")
  homogeneous <- as.logical(homogeneous)

  w <- .RGtkCall("S_gtk_notebook_set_homogeneous_tabs", object, homogeneous, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkNotebookSetTabBorder <-
function(object, border.width)
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  checkPtrType(object, "GtkNotebook")
  border.width <- as.numeric(border.width)

  w <- .RGtkCall("S_gtk_notebook_set_tab_border", object, border.width, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkNotebookSetTabHborder <-
function(object, tab.hborder)
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  checkPtrType(object, "GtkNotebook")
  tab.hborder <- as.numeric(tab.hborder)

  w <- .RGtkCall("S_gtk_notebook_set_tab_hborder", object, tab.hborder, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkNotebookSetTabVborder <-
function(object, tab.vborder)
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  checkPtrType(object, "GtkNotebook")
  tab.vborder <- as.numeric(tab.vborder)

  w <- .RGtkCall("S_gtk_notebook_set_tab_vborder", object, tab.vborder, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkNotebookSetScrollable <-
function(object, scrollable)
{
  checkPtrType(object, "GtkNotebook")
  scrollable <- as.logical(scrollable)

  w <- .RGtkCall("S_gtk_notebook_set_scrollable", object, scrollable, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkNotebookGetScrollable <-
function(object)
{
  checkPtrType(object, "GtkNotebook")

  w <- .RGtkCall("S_gtk_notebook_get_scrollable", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkNotebookPopupEnable <-
function(object)
{
  checkPtrType(object, "GtkNotebook")

  w <- .RGtkCall("S_gtk_notebook_popup_enable", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkNotebookPopupDisable <-
function(object)
{
  checkPtrType(object, "GtkNotebook")

  w <- .RGtkCall("S_gtk_notebook_popup_disable", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkNotebookGetTabLabel <-
function(object, child)
{
  checkPtrType(object, "GtkNotebook")
  checkPtrType(child, "GtkWidget")

  w <- .RGtkCall("S_gtk_notebook_get_tab_label", object, child, PACKAGE = "RGtk2")

  return(w)
} 


gtkNotebookSetTabLabel <-
function(object, child, tab.label = NULL)
{
  checkPtrType(object, "GtkNotebook")
  checkPtrType(child, "GtkWidget")
  if (!is.null( tab.label )) checkPtrType(tab.label, "GtkWidget")

  w <- .RGtkCall("S_gtk_notebook_set_tab_label", object, child, tab.label, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkNotebookSetTabLabelText <-
function(object, child, tab.text)
{
  checkPtrType(object, "GtkNotebook")
  checkPtrType(child, "GtkWidget")
  tab.text <- as.character(tab.text)

  w <- .RGtkCall("S_gtk_notebook_set_tab_label_text", object, child, tab.text, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkNotebookGetTabLabelText <-
function(object, child)
{
  checkPtrType(object, "GtkNotebook")
  checkPtrType(child, "GtkWidget")

  w <- .RGtkCall("S_gtk_notebook_get_tab_label_text", object, child, PACKAGE = "RGtk2")

  return(w)
} 


gtkNotebookGetMenuLabel <-
function(object, child)
{
  checkPtrType(object, "GtkNotebook")
  checkPtrType(child, "GtkWidget")

  w <- .RGtkCall("S_gtk_notebook_get_menu_label", object, child, PACKAGE = "RGtk2")

  return(w)
} 


gtkNotebookSetMenuLabel <-
function(object, child, menu.label = NULL)
{
  checkPtrType(object, "GtkNotebook")
  checkPtrType(child, "GtkWidget")
  if (!is.null( menu.label )) checkPtrType(menu.label, "GtkWidget")

  w <- .RGtkCall("S_gtk_notebook_set_menu_label", object, child, menu.label, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkNotebookSetMenuLabelText <-
function(object, child, menu.text)
{
  checkPtrType(object, "GtkNotebook")
  checkPtrType(child, "GtkWidget")
  menu.text <- as.character(menu.text)

  w <- .RGtkCall("S_gtk_notebook_set_menu_label_text", object, child, menu.text, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkNotebookGetMenuLabelText <-
function(object, child)
{
  checkPtrType(object, "GtkNotebook")
  checkPtrType(child, "GtkWidget")

  w <- .RGtkCall("S_gtk_notebook_get_menu_label_text", object, child, PACKAGE = "RGtk2")

  return(w)
} 


gtkNotebookQueryTabLabelPacking <-
function(object, child)
{
  if(getOption("depwarn"))
    .Deprecated("'tab-expand' and 'tab-fill' child properties", "RGtk2")

  checkPtrType(object, "GtkNotebook")
  checkPtrType(child, "GtkWidget")

  w <- .RGtkCall("S_gtk_notebook_query_tab_label_packing", object, child, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkNotebookSetTabLabelPacking <-
function(object, child, expand, fill, pack.type)
{
  if(getOption("depwarn"))
    .Deprecated("'tab-expand' and 'tab-fill' child properties", "RGtk2")

  checkPtrType(object, "GtkNotebook")
  checkPtrType(child, "GtkWidget")
  expand <- as.logical(expand)
  fill <- as.logical(fill)
  

  w <- .RGtkCall("S_gtk_notebook_set_tab_label_packing", object, child, expand, fill, pack.type, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkNotebookReorderChild <-
function(object, child, position)
{
  checkPtrType(object, "GtkNotebook")
  checkPtrType(child, "GtkWidget")
  position <- as.integer(position)

  w <- .RGtkCall("S_gtk_notebook_reorder_child", object, child, position, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkNotebookCurrentPage <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("gtkNotebookGetCurrentPage", "RGtk2")

  checkPtrType(object, "GtkNotebook")

  w <- .RGtkCall("S_gtk_notebook_current_page", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkNotebookSetPage <-
function(object, page.num)
{
  if(getOption("depwarn"))
    .Deprecated("gtkNotebookSetCurrentPage", "RGtk2")

  checkPtrType(object, "GtkNotebook")
  page.num <- as.integer(page.num)

  w <- .RGtkCall("S_gtk_notebook_set_page", object, page.num, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkObjectGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_object_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkOldEditableGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_old_editable_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkOldEditableClaimSelection <-
function(object, claim, time = "GDK_CURRENT_TIME")
{
  checkPtrType(object, "GtkOldEditable")
  claim <- as.logical(claim)
  time <- as.numeric(time)

  w <- .RGtkCall("S_gtk_old_editable_claim_selection", object, claim, time, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkOldEditableChanged <-
function(object)
{
  checkPtrType(object, "GtkOldEditable")

  w <- .RGtkCall("S_gtk_old_editable_changed", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkOptionMenuGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_option_menu_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkOptionMenuNew <-
function(show = TRUE)
{
  if(getOption("depwarn"))
    .Deprecated("GtkComboBox", "RGtk2")

  

  w <- .RGtkCall("S_gtk_option_menu_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkOptionMenuGetMenu <-
function(object)
{
  checkPtrType(object, "GtkOptionMenu")

  w <- .RGtkCall("S_gtk_option_menu_get_menu", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkOptionMenuSetMenu <-
function(object, menu)
{
  checkPtrType(object, "GtkOptionMenu")
  checkPtrType(menu, "GtkWidget")

  w <- .RGtkCall("S_gtk_option_menu_set_menu", object, menu, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkOptionMenuRemoveMenu <-
function(object)
{
  checkPtrType(object, "GtkOptionMenu")

  w <- .RGtkCall("S_gtk_option_menu_remove_menu", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkOptionMenuGetHistory <-
function(object)
{
  checkPtrType(object, "GtkOptionMenu")

  w <- .RGtkCall("S_gtk_option_menu_get_history", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkOptionMenuSetHistory <-
function(object, index)
{
  checkPtrType(object, "GtkOptionMenu")
  index <- as.numeric(index)

  w <- .RGtkCall("S_gtk_option_menu_set_history", object, index, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPanedGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_paned_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkPanedAdd1 <-
function(object, child)
{
  checkPtrType(object, "GtkPaned")
  checkPtrType(child, "GtkWidget")

  w <- .RGtkCall("S_gtk_paned_add1", object, child, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPanedAdd2 <-
function(object, child)
{
  checkPtrType(object, "GtkPaned")
  checkPtrType(child, "GtkWidget")

  w <- .RGtkCall("S_gtk_paned_add2", object, child, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPanedPack1 <-
function(object, child, resize = FALSE, shrink = TRUE)
{
  checkPtrType(object, "GtkPaned")
  checkPtrType(child, "GtkWidget")
  resize <- as.logical(resize)
  shrink <- as.logical(shrink)

  w <- .RGtkCall("S_gtk_paned_pack1", object, child, resize, shrink, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPanedPack2 <-
function(object, child, resize = TRUE, shrink = TRUE)
{
  checkPtrType(object, "GtkPaned")
  checkPtrType(child, "GtkWidget")
  resize <- as.logical(resize)
  shrink <- as.logical(shrink)

  w <- .RGtkCall("S_gtk_paned_pack2", object, child, resize, shrink, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPanedGetPosition <-
function(object)
{
  checkPtrType(object, "GtkPaned")

  w <- .RGtkCall("S_gtk_paned_get_position", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPanedSetPosition <-
function(object, position)
{
  checkPtrType(object, "GtkPaned")
  position <- as.integer(position)

  w <- .RGtkCall("S_gtk_paned_set_position", object, position, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPanedGetChild1 <-
function(object)
{
  checkPtrType(object, "GtkPaned")

  w <- .RGtkCall("S_gtk_paned_get_child1", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPanedGetChild2 <-
function(object)
{
  checkPtrType(object, "GtkPaned")

  w <- .RGtkCall("S_gtk_paned_get_child2", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPixmapGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_pixmap_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkPixmapNew <-
function(pixmap, mask = NULL, show = TRUE)
{
  if(getOption("depwarn"))
    .Deprecated("GtkImage", "RGtk2")

  checkPtrType(pixmap, "GdkPixmap")
  if (!is.null( mask )) checkPtrType(mask, "GdkBitmap")

  w <- .RGtkCall("S_gtk_pixmap_new", pixmap, mask, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkPixmapSet <-
function(object, val, mask = NULL)
{
  checkPtrType(object, "GtkPixmap")
  checkPtrType(val, "GdkPixmap")
  if (!is.null( mask )) checkPtrType(mask, "GdkBitmap")

  w <- .RGtkCall("S_gtk_pixmap_set", object, val, mask, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPixmapGet <-
function(object)
{
  checkPtrType(object, "GtkPixmap")

  w <- .RGtkCall("S_gtk_pixmap_get", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPixmapSetBuildInsensitive <-
function(object, build)
{
  checkPtrType(object, "GtkPixmap")
  build <- as.logical(build)

  w <- .RGtkCall("S_gtk_pixmap_set_build_insensitive", object, build, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPlugGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_plug_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkPlugConstruct <-
function(object, socket.id)
{
  checkPtrType(object, "GtkPlug")
  socket.id <- as.GdkNativeWindow(socket.id)

  w <- .RGtkCall("S_gtk_plug_construct", object, socket.id, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPlugNew <-
function(socket.id, show = TRUE)
{
  socket.id <- as.GdkNativeWindow(socket.id)

  w <- .RGtkCall("S_gtk_plug_new", socket.id, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkPlugConstructForDisplay <-
function(object, display, socket.id)
{
  checkPtrType(object, "GtkPlug")
  checkPtrType(display, "GdkDisplay")
  socket.id <- as.GdkNativeWindow(socket.id)

  w <- .RGtkCall("S_gtk_plug_construct_for_display", object, display, socket.id, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPlugNewForDisplay <-
function(display, socket.id)
{
  checkPtrType(display, "GdkDisplay")
  socket.id <- as.GdkNativeWindow(socket.id)

  w <- .RGtkCall("S_gtk_plug_new_for_display", display, socket.id, PACKAGE = "RGtk2")

  return(w)
} 


gtkPlugGetId <-
function(object)
{
  checkPtrType(object, "GtkPlug")

  w <- .RGtkCall("S_gtk_plug_get_id", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPreviewGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_preview_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkPreviewUninit <-
function()
{
  

  w <- .RGtkCall("S_gtk_preview_uninit", PACKAGE = "RGtk2")

  return(w)
} 


gtkPreviewNew <-
function(type, show = TRUE)
{
  if(getOption("depwarn"))
    .Deprecated("GtkDrawingArea", "RGtk2")

  

  w <- .RGtkCall("S_gtk_preview_new", type, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkPreviewSize <-
function(object, width, height)
{
  checkPtrType(object, "GtkPreview")
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_preview_size", object, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPreviewPut <-
function(object, window, gc, srcx, srcy, destx, desty, width, height)
{
  checkPtrType(object, "GtkPreview")
  checkPtrType(window, "GdkWindow")
  checkPtrType(gc, "GdkGC")
  srcx <- as.integer(srcx)
  srcy <- as.integer(srcy)
  destx <- as.integer(destx)
  desty <- as.integer(desty)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_preview_put", object, window, gc, srcx, srcy, destx, desty, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPreviewDrawRow <-
function(object, data, y, w)
{
  checkPtrType(object, "GtkPreview")
  data <- as.list(as.raw(data))
  y <- as.integer(y)
  w <- as.integer(w)

  w <- .RGtkCall("S_gtk_preview_draw_row", object, data, y, w, PACKAGE = "RGtk2")

  return(w)
} 


gtkPreviewSetExpand <-
function(object, expand)
{
  checkPtrType(object, "GtkPreview")
  expand <- as.logical(expand)

  w <- .RGtkCall("S_gtk_preview_set_expand", object, expand, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPreviewSetGamma <-
function(gamma)
{
  gamma <- as.numeric(gamma)

  w <- .RGtkCall("S_gtk_preview_set_gamma", gamma, PACKAGE = "RGtk2")

  return(w)
} 


gtkPreviewSetColorCube <-
function(nred.shades, ngreen.shades, nblue.shades, ngray.shades)
{
  nred.shades <- as.numeric(nred.shades)
  ngreen.shades <- as.numeric(ngreen.shades)
  nblue.shades <- as.numeric(nblue.shades)
  ngray.shades <- as.numeric(ngray.shades)

  w <- .RGtkCall("S_gtk_preview_set_color_cube", nred.shades, ngreen.shades, nblue.shades, ngray.shades, PACKAGE = "RGtk2")

  return(w)
} 


gtkPreviewSetInstallCmap <-
function(install.cmap)
{
  install.cmap <- as.integer(install.cmap)

  w <- .RGtkCall("S_gtk_preview_set_install_cmap", install.cmap, PACKAGE = "RGtk2")

  return(w)
} 


gtkPreviewSetReserved <-
function(nreserved)
{
  nreserved <- as.integer(nreserved)

  w <- .RGtkCall("S_gtk_preview_set_reserved", nreserved, PACKAGE = "RGtk2")

  return(w)
} 


gtkPreviewSetDither <-
function(object, dither)
{
  checkPtrType(object, "GtkPreview")
  

  w <- .RGtkCall("S_gtk_preview_set_dither", object, dither, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPreviewGetVisual <-
function()
{
  

  w <- .RGtkCall("S_gtk_preview_get_visual", PACKAGE = "RGtk2")

  return(w)
} 


gtkPreviewGetCmap <-
function()
{
  

  w <- .RGtkCall("S_gtk_preview_get_cmap", PACKAGE = "RGtk2")

  return(w)
} 


gtkPreviewGetInfo <-
function()
{
  

  w <- .RGtkCall("S_gtk_preview_get_info", PACKAGE = "RGtk2")

  return(w)
} 


gtkPreviewReset <-
function()
{
  

  w <- .RGtkCall("S_gtk_preview_reset", PACKAGE = "RGtk2")

  return(w)
} 


gtkProgressGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_progress_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkProgressSetShowText <-
function(object, show.text)
{
  if(getOption("depwarn"))
    .Deprecated("GtkProgressBar methods", "RGtk2")

  checkPtrType(object, "GtkProgress")
  show.text <- as.logical(show.text)

  w <- .RGtkCall("S_gtk_progress_set_show_text", object, show.text, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkProgressSetTextAlignment <-
function(object, x.align, y.align)
{
  if(getOption("depwarn"))
    .Deprecated("GtkProgressBar methods", "RGtk2")

  checkPtrType(object, "GtkProgress")
  x.align <- as.numeric(x.align)
  y.align <- as.numeric(y.align)

  w <- .RGtkCall("S_gtk_progress_set_text_alignment", object, x.align, y.align, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkProgressSetFormatString <-
function(object, format)
{
  if(getOption("depwarn"))
    .Deprecated("GtkProgressBar methods", "RGtk2")

  checkPtrType(object, "GtkProgress")
  format <- as.character(format)

  w <- .RGtkCall("S_gtk_progress_set_format_string", object, format, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkProgressSetAdjustment <-
function(object, adjustment)
{
  if(getOption("depwarn"))
    .Deprecated("GtkProgressBar methods", "RGtk2")

  checkPtrType(object, "GtkProgress")
  checkPtrType(adjustment, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_progress_set_adjustment", object, adjustment, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkProgressConfigure <-
function(object, value, min, max)
{
  if(getOption("depwarn"))
    .Deprecated("GtkProgressBar methods", "RGtk2")

  checkPtrType(object, "GtkProgress")
  value <- as.numeric(value)
  min <- as.numeric(min)
  max <- as.numeric(max)

  w <- .RGtkCall("S_gtk_progress_configure", object, value, min, max, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkProgressSetPercentage <-
function(object, percentage)
{
  if(getOption("depwarn"))
    .Deprecated("GtkProgressBar methods", "RGtk2")

  checkPtrType(object, "GtkProgress")
  percentage <- as.numeric(percentage)

  w <- .RGtkCall("S_gtk_progress_set_percentage", object, percentage, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkProgressSetValue <-
function(object, value)
{
  if(getOption("depwarn"))
    .Deprecated("GtkProgressBar methods", "RGtk2")

  checkPtrType(object, "GtkProgress")
  value <- as.numeric(value)

  w <- .RGtkCall("S_gtk_progress_set_value", object, value, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkProgressGetValue <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("GtkProgressBar methods", "RGtk2")

  checkPtrType(object, "GtkProgress")

  w <- .RGtkCall("S_gtk_progress_get_value", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkProgressSetActivityMode <-
function(object, activity.mode)
{
  if(getOption("depwarn"))
    .Deprecated("GtkProgressBar methods", "RGtk2")

  checkPtrType(object, "GtkProgress")
  activity.mode <- as.logical(activity.mode)

  w <- .RGtkCall("S_gtk_progress_set_activity_mode", object, activity.mode, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkProgressGetCurrentText <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("GtkProgressBar methods", "RGtk2")

  checkPtrType(object, "GtkProgress")

  w <- .RGtkCall("S_gtk_progress_get_current_text", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkProgressGetTextFromValue <-
function(object, value)
{
  if(getOption("depwarn"))
    .Deprecated("GtkProgressBar methods", "RGtk2")

  checkPtrType(object, "GtkProgress")
  value <- as.numeric(value)

  w <- .RGtkCall("S_gtk_progress_get_text_from_value", object, value, PACKAGE = "RGtk2")

  return(w)
} 


gtkProgressGetCurrentPercentage <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("GtkProgressBar methods", "RGtk2")

  checkPtrType(object, "GtkProgress")

  w <- .RGtkCall("S_gtk_progress_get_current_percentage", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkProgressGetPercentageFromValue <-
function(object, value)
{
  if(getOption("depwarn"))
    .Deprecated("GtkProgressBar methods", "RGtk2")

  checkPtrType(object, "GtkProgress")
  value <- as.numeric(value)

  w <- .RGtkCall("S_gtk_progress_get_percentage_from_value", object, value, PACKAGE = "RGtk2")

  return(w)
} 


gtkProgressBarGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_progress_bar_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkProgressBarNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_progress_bar_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkProgressBarPulse <-
function(object)
{
  checkPtrType(object, "GtkProgressBar")

  w <- .RGtkCall("S_gtk_progress_bar_pulse", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkProgressBarSetText <-
function(object, text)
{
  checkPtrType(object, "GtkProgressBar")
  text <- as.character(text)

  w <- .RGtkCall("S_gtk_progress_bar_set_text", object, text, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkProgressBarSetFraction <-
function(object, fraction)
{
  checkPtrType(object, "GtkProgressBar")
  fraction <- as.numeric(fraction)

  w <- .RGtkCall("S_gtk_progress_bar_set_fraction", object, fraction, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkProgressBarSetPulseStep <-
function(object, fraction)
{
  checkPtrType(object, "GtkProgressBar")
  fraction <- as.numeric(fraction)

  w <- .RGtkCall("S_gtk_progress_bar_set_pulse_step", object, fraction, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkProgressBarSetOrientation <-
function(object, orientation)
{
  checkPtrType(object, "GtkProgressBar")
  

  w <- .RGtkCall("S_gtk_progress_bar_set_orientation", object, orientation, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkProgressBarGetText <-
function(object)
{
  checkPtrType(object, "GtkProgressBar")

  w <- .RGtkCall("S_gtk_progress_bar_get_text", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkProgressBarGetFraction <-
function(object)
{
  checkPtrType(object, "GtkProgressBar")

  w <- .RGtkCall("S_gtk_progress_bar_get_fraction", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkProgressBarGetPulseStep <-
function(object)
{
  checkPtrType(object, "GtkProgressBar")

  w <- .RGtkCall("S_gtk_progress_bar_get_pulse_step", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkProgressBarGetOrientation <-
function(object)
{
  checkPtrType(object, "GtkProgressBar")

  w <- .RGtkCall("S_gtk_progress_bar_get_orientation", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkProgressBarNewWithAdjustment <-
function(adjustment = NULL, show = TRUE)
{
  if(getOption("depwarn"))
    .Deprecated("", "RGtk2")

  if (!is.null( adjustment )) checkPtrType(adjustment, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_progress_bar_new_with_adjustment", adjustment, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkProgressBarSetBarStyle <-
function(object, style)
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  checkPtrType(object, "GtkProgressBar")
  

  w <- .RGtkCall("S_gtk_progress_bar_set_bar_style", object, style, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkProgressBarSetDiscreteBlocks <-
function(object, blocks)
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  checkPtrType(object, "GtkProgressBar")
  blocks <- as.numeric(blocks)

  w <- .RGtkCall("S_gtk_progress_bar_set_discrete_blocks", object, blocks, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkProgressBarSetActivityStep <-
function(object, step)
{
  if(getOption("depwarn"))
    .Deprecated("'pulse-step' property", "RGtk2")

  checkPtrType(object, "GtkProgressBar")
  step <- as.numeric(step)

  w <- .RGtkCall("S_gtk_progress_bar_set_activity_step", object, step, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkProgressBarSetActivityBlocks <-
function(object, blocks)
{
  if(getOption("depwarn"))
    .Deprecated("'pulse-step' property", "RGtk2")

  checkPtrType(object, "GtkProgressBar")
  blocks <- as.numeric(blocks)

  w <- .RGtkCall("S_gtk_progress_bar_set_activity_blocks", object, blocks, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkProgressBarUpdate <-
function(object, percentage)
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  checkPtrType(object, "GtkProgressBar")
  percentage <- as.numeric(percentage)

  w <- .RGtkCall("S_gtk_progress_bar_update", object, percentage, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkProgressBarSetEllipsize <-
function(object, mode)
{
  checkPtrType(object, "GtkProgressBar")
  

  w <- .RGtkCall("S_gtk_progress_bar_set_ellipsize", object, mode, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkProgressBarGetEllipsize <-
function(object)
{
  checkPtrType(object, "GtkProgressBar")

  w <- .RGtkCall("S_gtk_progress_bar_get_ellipsize", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRadioActionGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_radio_action_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkRadioActionNew <-
function(name = NULL, label = NULL, tooltip = NULL, stock.id = NULL, value = NULL)
{
  

  w <- .RGtkCall("S_gtk_radio_action_new", name, label, tooltip, stock.id, value, PACKAGE = "RGtk2")

  return(w)
} 


gtkRadioActionSetGroup <-
function(object, group)
{
  checkPtrType(object, "GtkRadioAction")
  group <- as.GSList(group)

  w <- .RGtkCall("S_gtk_radio_action_set_group", object, group, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRadioActionGetGroup <-
function(object)
{
  checkPtrType(object, "GtkRadioAction")

  w <- .RGtkCall("S_gtk_radio_action_get_group", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRadioActionGetCurrentValue <-
function(object)
{
  checkPtrType(object, "GtkRadioAction")

  w <- .RGtkCall("S_gtk_radio_action_get_current_value", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRadioButtonGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_radio_button_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkRadioButtonNewFromWidget <-
function(group = NULL, show = TRUE)
{
  if (!is.null( group )) checkPtrType(group, "GtkRadioButton")

  w <- .RGtkCall("S_gtk_radio_button_new_from_widget", group, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkRadioButtonNewWithLabelFromWidget <-
function(group = NULL, label, show = TRUE)
{
  if (!is.null( group )) checkPtrType(group, "GtkRadioButton")
  label <- as.character(label)

  w <- .RGtkCall("S_gtk_radio_button_new_with_label_from_widget", group, label, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkRadioButtonNewWithMnemonic <-
function(group, label, show = TRUE)
{
  group <- as.GSList(group)
  label <- as.character(label)

  w <- .RGtkCall("S_gtk_radio_button_new_with_mnemonic", group, label, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkRadioButtonNewWithMnemonicFromWidget <-
function(group = NULL, label, show = TRUE)
{
  if (!is.null( group )) checkPtrType(group, "GtkRadioButton")
  label <- as.character(label)

  w <- .RGtkCall("S_gtk_radio_button_new_with_mnemonic_from_widget", group, label, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkRadioButtonGetGroup <-
function(object)
{
  checkPtrType(object, "GtkRadioButton")

  w <- .RGtkCall("S_gtk_radio_button_get_group", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRadioButtonSetGroup <-
function(object, group)
{
  checkPtrType(object, "GtkRadioButton")
  group <- as.GSList(group)

  w <- .RGtkCall("S_gtk_radio_button_set_group", object, group, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRadioButtonGroup <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("gtkRadioButtonGetGroup", "RGtk2")

  checkPtrType(object, "GtkRadioButton")

  w <- .RGtkCall("S_gtk_radio_button_group", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRadioMenuItemGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_radio_menu_item_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkRadioMenuItemNewFromWidget <-
function(group = NULL, show = TRUE)
{
  if (!is.null( group )) checkPtrType(group, "GtkRadioMenuItem")

  w <- .RGtkCall("S_gtk_radio_menu_item_new_from_widget", group, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkRadioMenuItemNewWithMnemonicFromWidget <-
function(group = NULL, label, show = TRUE)
{
  if (!is.null( group )) checkPtrType(group, "GtkRadioMenuItem")
  label <- as.character(label)

  w <- .RGtkCall("S_gtk_radio_menu_item_new_with_mnemonic_from_widget", group, label, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkRadioMenuItemNewWithLabelFromWidget <-
function(group = NULL, label, show = TRUE)
{
  if (!is.null( group )) checkPtrType(group, "GtkRadioMenuItem")
  label <- as.character(label)

  w <- .RGtkCall("S_gtk_radio_menu_item_new_with_label_from_widget", group, label, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkRadioMenuItemGetGroup <-
function(object)
{
  checkPtrType(object, "GtkRadioMenuItem")

  w <- .RGtkCall("S_gtk_radio_menu_item_get_group", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRadioMenuItemSetGroup <-
function(object, group)
{
  checkPtrType(object, "GtkRadioMenuItem")
  group <- as.GSList(group)

  w <- .RGtkCall("S_gtk_radio_menu_item_set_group", object, group, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRadioMenuItemGroup <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("gtkRadioMenuItemGetGroup", "RGtk2")

  checkPtrType(object, "GtkRadioMenuItem")

  w <- .RGtkCall("S_gtk_radio_menu_item_group", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRadioToolButtonGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_radio_tool_button_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkRadioToolButtonNewFromWidget <-
function(group = NULL, show = TRUE)
{
  if (!is.null( group )) checkPtrType(group, "GtkRadioToolButton")

  w <- .RGtkCall("S_gtk_radio_tool_button_new_from_widget", group, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkRadioToolButtonNewWithStockFromWidget <-
function(group = NULL, stock.id, show = TRUE)
{
  if (!is.null( group )) checkPtrType(group, "GtkRadioToolButton")
  stock.id <- as.character(stock.id)

  w <- .RGtkCall("S_gtk_radio_tool_button_new_with_stock_from_widget", group, stock.id, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkRadioToolButtonSetGroup <-
function(object, group)
{
  checkPtrType(object, "GtkRadioToolButton")
  group <- as.GSList(group)

  w <- .RGtkCall("S_gtk_radio_tool_button_set_group", object, group, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRadioToolButtonGetGroup <-
function(object)
{
  checkPtrType(object, "GtkRadioToolButton")

  w <- .RGtkCall("S_gtk_radio_tool_button_get_group", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRangeGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_range_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkRangeSetUpdatePolicy <-
function(object, policy)
{
  checkPtrType(object, "GtkRange")
  

  w <- .RGtkCall("S_gtk_range_set_update_policy", object, policy, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRangeGetUpdatePolicy <-
function(object)
{
  checkPtrType(object, "GtkRange")

  w <- .RGtkCall("S_gtk_range_get_update_policy", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRangeSetAdjustment <-
function(object, adjustment)
{
  checkPtrType(object, "GtkRange")
  checkPtrType(adjustment, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_range_set_adjustment", object, adjustment, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRangeGetAdjustment <-
function(object)
{
  checkPtrType(object, "GtkRange")

  w <- .RGtkCall("S_gtk_range_get_adjustment", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRangeSetInverted <-
function(object, setting)
{
  checkPtrType(object, "GtkRange")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_range_set_inverted", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRangeGetInverted <-
function(object)
{
  checkPtrType(object, "GtkRange")

  w <- .RGtkCall("S_gtk_range_get_inverted", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRangeSetIncrements <-
function(object, step, page)
{
  checkPtrType(object, "GtkRange")
  step <- as.numeric(step)
  page <- as.numeric(page)

  w <- .RGtkCall("S_gtk_range_set_increments", object, step, page, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRangeSetRange <-
function(object, min, max)
{
  checkPtrType(object, "GtkRange")
  min <- as.numeric(min)
  max <- as.numeric(max)

  w <- .RGtkCall("S_gtk_range_set_range", object, min, max, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRangeSetValue <-
function(object, value)
{
  checkPtrType(object, "GtkRange")
  value <- as.numeric(value)

  w <- .RGtkCall("S_gtk_range_set_value", object, value, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRangeGetValue <-
function(object)
{
  checkPtrType(object, "GtkRange")

  w <- .RGtkCall("S_gtk_range_get_value", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRcAddDefaultFile <-
function(filename)
{
  filename <- as.character(filename)

  w <- .RGtkCall("S_gtk_rc_add_default_file", filename, PACKAGE = "RGtk2")

  return(w)
} 


gtkRcSetDefaultFiles <-
function(filenames)
{
  filenames <- as.list(as.character(filenames))

  w <- .RGtkCall("S_gtk_rc_set_default_files", filenames, PACKAGE = "RGtk2")

  return(w)
} 


gtkRcGetDefaultFiles <-
function()
{
  

  w <- .RGtkCall("S_gtk_rc_get_default_files", PACKAGE = "RGtk2")

  return(w)
} 


gtkRcGetStyle <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_rc_get_style", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRcGetStyleByPaths <-
function(settings, widget.path, class.path, type)
{
  checkPtrType(settings, "GtkSettings")
  widget.path <- as.character(widget.path)
  class.path <- as.character(class.path)
  type <- as.GType(type)

  w <- .RGtkCall("S_gtk_rc_get_style_by_paths", settings, widget.path, class.path, type, PACKAGE = "RGtk2")

  return(w)
} 


gtkRcReparseAllForSettings <-
function(settings, force.load)
{
  checkPtrType(settings, "GtkSettings")
  force.load <- as.logical(force.load)

  w <- .RGtkCall("S_gtk_rc_reparse_all_for_settings", settings, force.load, PACKAGE = "RGtk2")

  return(w)
} 


gtkRcResetStyles <-
function(settings)
{
  checkPtrType(settings, "GtkSettings")

  w <- .RGtkCall("S_gtk_rc_reset_styles", settings, PACKAGE = "RGtk2")

  return(w)
} 


gtkRcFindPixmapInPath <-
function(settings, scanner = NULL, pixmap.file)
{
  checkPtrType(settings, "GtkSettings")
  if (!is.null( scanner )) checkPtrType(scanner, "GScanner")
  pixmap.file <- as.character(pixmap.file)

  w <- .RGtkCall("S_gtk_rc_find_pixmap_in_path", settings, scanner, pixmap.file, PACKAGE = "RGtk2")

  return(w)
} 


gtkRcParse <-
function(filename)
{
  filename <- as.character(filename)

  w <- .RGtkCall("S_gtk_rc_parse", filename, PACKAGE = "RGtk2")

  return(w)
} 


gtkRcParseString <-
function(rc.string)
{
  rc.string <- as.character(rc.string)

  w <- .RGtkCall("S_gtk_rc_parse_string", rc.string, PACKAGE = "RGtk2")

  return(w)
} 


gtkRcReparseAll <-
function()
{
  

  w <- .RGtkCall("S_gtk_rc_reparse_all", PACKAGE = "RGtk2")

  return(w)
} 


gtkRcAddWidgetNameStyle <-
function(object, pattern)
{
  checkPtrType(object, "GtkRcStyle")
  pattern <- as.character(pattern)

  w <- .RGtkCall("S_gtk_rc_add_widget_name_style", object, pattern, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRcAddWidgetClassStyle <-
function(object, pattern)
{
  checkPtrType(object, "GtkRcStyle")
  pattern <- as.character(pattern)

  w <- .RGtkCall("S_gtk_rc_add_widget_class_style", object, pattern, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRcAddClassStyle <-
function(object, pattern)
{
  checkPtrType(object, "GtkRcStyle")
  pattern <- as.character(pattern)

  w <- .RGtkCall("S_gtk_rc_add_class_style", object, pattern, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRcStyleGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_rc_style_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkRcStyleNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_rc_style_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkRcStyleCopy <-
function(object)
{
  checkPtrType(object, "GtkRcStyle")

  w <- .RGtkCall("S_gtk_rc_style_copy", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRcFindModuleInPath <-
function(module.file)
{
  module.file <- as.character(module.file)

  w <- .RGtkCall("S_gtk_rc_find_module_in_path", module.file, PACKAGE = "RGtk2")

  return(w)
} 


gtkRcGetThemeDir <-
function()
{
  

  w <- .RGtkCall("S_gtk_rc_get_theme_dir", PACKAGE = "RGtk2")

  return(w)
} 


gtkRcGetModuleDir <-
function()
{
  

  w <- .RGtkCall("S_gtk_rc_get_module_dir", PACKAGE = "RGtk2")

  return(w)
} 


gtkRcGetImModulePath <-
function()
{
  

  w <- .RGtkCall("S_gtk_rc_get_im_module_path", PACKAGE = "RGtk2")

  return(w)
} 


gtkRcGetImModuleFile <-
function()
{
  

  w <- .RGtkCall("S_gtk_rc_get_im_module_file", PACKAGE = "RGtk2")

  return(w)
} 


gtkRcScannerNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_rc_scanner_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkRcParseColor <-
function(scanner, color)
{
  checkPtrType(scanner, "GScanner")
  color <- as.GdkColor(color)

  w <- .RGtkCall("S_gtk_rc_parse_color", scanner, color, PACKAGE = "RGtk2")

  return(w)
} 


gtkRcParseState <-
function(scanner)
{
  checkPtrType(scanner, "GScanner")

  w <- .RGtkCall("S_gtk_rc_parse_state", scanner, PACKAGE = "RGtk2")

  return(w)
} 


gtkRcParsePriority <-
function(scanner)
{
  checkPtrType(scanner, "GScanner")

  w <- .RGtkCall("S_gtk_rc_parse_priority", scanner, PACKAGE = "RGtk2")

  return(w)
} 


gtkRulerGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_ruler_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkRulerSetMetric <-
function(object, metric)
{
  checkPtrType(object, "GtkRuler")
  

  w <- .RGtkCall("S_gtk_ruler_set_metric", object, metric, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRulerSetRange <-
function(object, lower, upper, position, max.size)
{
  checkPtrType(object, "GtkRuler")
  lower <- as.numeric(lower)
  upper <- as.numeric(upper)
  position <- as.numeric(position)
  max.size <- as.numeric(max.size)

  w <- .RGtkCall("S_gtk_ruler_set_range", object, lower, upper, position, max.size, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRulerGetMetric <-
function(object)
{
  checkPtrType(object, "GtkRuler")

  w <- .RGtkCall("S_gtk_ruler_get_metric", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRulerGetRange <-
function(object)
{
  checkPtrType(object, "GtkRuler")

  w <- .RGtkCall("S_gtk_ruler_get_range", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkScaleGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_scale_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkScaleSetDigits <-
function(object, digits)
{
  checkPtrType(object, "GtkScale")
  digits <- as.integer(digits)

  w <- .RGtkCall("S_gtk_scale_set_digits", object, digits, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkScaleGetDigits <-
function(object)
{
  checkPtrType(object, "GtkScale")

  w <- .RGtkCall("S_gtk_scale_get_digits", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkScaleSetDrawValue <-
function(object, draw.value)
{
  checkPtrType(object, "GtkScale")
  draw.value <- as.logical(draw.value)

  w <- .RGtkCall("S_gtk_scale_set_draw_value", object, draw.value, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkScaleGetDrawValue <-
function(object)
{
  checkPtrType(object, "GtkScale")

  w <- .RGtkCall("S_gtk_scale_get_draw_value", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkScaleSetValuePos <-
function(object, pos)
{
  checkPtrType(object, "GtkScale")
  

  w <- .RGtkCall("S_gtk_scale_set_value_pos", object, pos, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkScaleGetValuePos <-
function(object)
{
  checkPtrType(object, "GtkScale")

  w <- .RGtkCall("S_gtk_scale_get_value_pos", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkScaleGetLayout <-
function(object)
{
  checkPtrType(object, "GtkScale")

  w <- .RGtkCall("S_gtk_scale_get_layout", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkScaleGetLayoutOffsets <-
function(object)
{
  checkPtrType(object, "GtkScale")

  w <- .RGtkCall("S_gtk_scale_get_layout_offsets", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkScrollbarGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_scrollbar_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkScrolledWindowGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_scrolled_window_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkScrolledWindowNew <-
function(hadjustment = NULL, vadjustment = NULL, show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_scrolled_window_new", hadjustment, vadjustment, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkScrolledWindowSetHadjustment <-
function(object, hadjustment)
{
  checkPtrType(object, "GtkScrolledWindow")
  checkPtrType(hadjustment, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_scrolled_window_set_hadjustment", object, hadjustment, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkScrolledWindowSetVadjustment <-
function(object, hadjustment)
{
  checkPtrType(object, "GtkScrolledWindow")
  checkPtrType(hadjustment, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_scrolled_window_set_vadjustment", object, hadjustment, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkScrolledWindowGetHadjustment <-
function(object)
{
  checkPtrType(object, "GtkScrolledWindow")

  w <- .RGtkCall("S_gtk_scrolled_window_get_hadjustment", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkScrolledWindowGetVadjustment <-
function(object)
{
  checkPtrType(object, "GtkScrolledWindow")

  w <- .RGtkCall("S_gtk_scrolled_window_get_vadjustment", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkScrolledWindowGetHscrollbar <-
function(object)
{
  checkPtrType(object, "GtkScrolledWindow")

  w <- .RGtkCall("S_gtk_scrolled_window_get_hscrollbar", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkScrolledWindowGetVscrollbar <-
function(object)
{
  checkPtrType(object, "GtkScrolledWindow")

  w <- .RGtkCall("S_gtk_scrolled_window_get_vscrollbar", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkScrolledWindowSetPolicy <-
function(object, hscrollbar.policy, vscrollbar.policy)
{
  checkPtrType(object, "GtkScrolledWindow")
  
  

  w <- .RGtkCall("S_gtk_scrolled_window_set_policy", object, hscrollbar.policy, vscrollbar.policy, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkScrolledWindowGetPolicy <-
function(object)
{
  checkPtrType(object, "GtkScrolledWindow")

  w <- .RGtkCall("S_gtk_scrolled_window_get_policy", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkScrolledWindowSetPlacement <-
function(object, window.placement)
{
  checkPtrType(object, "GtkScrolledWindow")
  

  w <- .RGtkCall("S_gtk_scrolled_window_set_placement", object, window.placement, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkScrolledWindowGetPlacement <-
function(object)
{
  checkPtrType(object, "GtkScrolledWindow")

  w <- .RGtkCall("S_gtk_scrolled_window_get_placement", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkScrolledWindowSetShadowType <-
function(object, type)
{
  checkPtrType(object, "GtkScrolledWindow")
  

  w <- .RGtkCall("S_gtk_scrolled_window_set_shadow_type", object, type, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkScrolledWindowGetShadowType <-
function(object)
{
  checkPtrType(object, "GtkScrolledWindow")

  w <- .RGtkCall("S_gtk_scrolled_window_get_shadow_type", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkScrolledWindowAddWithViewport <-
function(object, child)
{
  checkPtrType(object, "GtkScrolledWindow")
  checkPtrType(child, "GtkWidget")

  w <- .RGtkCall("S_gtk_scrolled_window_add_with_viewport", object, child, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTargetListNew <-
function(targets)
{
  targets <- lapply(targets, function(x) { x <- as.GtkTargetEntry(x); x })

  w <- .RGtkCall("S_gtk_target_list_new", targets, PACKAGE = "RGtk2")

  return(w)
} 


gtkTargetListAdd <-
function(object, target, flags, info)
{
  checkPtrType(object, "GtkTargetList")
  target <- as.GdkAtom(target)
  flags <- as.numeric(flags)
  info <- as.numeric(info)

  w <- .RGtkCall("S_gtk_target_list_add", object, target, flags, info, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTargetListAddTable <-
function(object, targets)
{
  checkPtrType(object, "GtkTargetList")
  targets <- lapply(targets, function(x) { x <- as.GtkTargetEntry(x); x })

  w <- .RGtkCall("S_gtk_target_list_add_table", object, targets, PACKAGE = "RGtk2")

  return(w)
} 


gtkTargetListRemove <-
function(object, target)
{
  checkPtrType(object, "GtkTargetList")
  target <- as.GdkAtom(target)

  w <- .RGtkCall("S_gtk_target_list_remove", object, target, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTargetListFind <-
function(object, target)
{
  checkPtrType(object, "GtkTargetList")
  target <- as.GdkAtom(target)

  w <- .RGtkCall("S_gtk_target_list_find", object, target, PACKAGE = "RGtk2")

  return(w)
} 


gtkSelectionOwnerSet <-
function(object, selection, time = "GDK_CURRENT_TIME")
{
  checkPtrType(object, "GtkWidget")
  selection <- as.GdkAtom(selection)
  time <- as.numeric(time)

  w <- .RGtkCall("S_gtk_selection_owner_set", object, selection, time, PACKAGE = "RGtk2")

  return(w)
} 


gtkSelectionOwnerSetForDisplay <-
function(display, widget = NULL, selection, time = "GDK_CURRENT_TIME")
{
  checkPtrType(display, "GdkDisplay")
  if (!is.null( widget )) checkPtrType(widget, "GtkWidget")
  selection <- as.GdkAtom(selection)
  time <- as.numeric(time)

  w <- .RGtkCall("S_gtk_selection_owner_set_for_display", display, widget, selection, time, PACKAGE = "RGtk2")

  return(w)
} 


gtkSelectionAddTarget <-
function(object, selection, target, info)
{
  checkPtrType(object, "GtkWidget")
  selection <- as.GdkAtom(selection)
  target <- as.GdkAtom(target)
  info <- as.numeric(info)

  w <- .RGtkCall("S_gtk_selection_add_target", object, selection, target, info, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkSelectionAddTargets <-
function(object, selection, targets)
{
  checkPtrType(object, "GtkWidget")
  selection <- as.GdkAtom(selection)
  targets <- lapply(targets, function(x) { x <- as.GtkTargetEntry(x); x })

  w <- .RGtkCall("S_gtk_selection_add_targets", object, selection, targets, PACKAGE = "RGtk2")

  return(w)
} 


gtkSelectionClearTargets <-
function(object, selection)
{
  checkPtrType(object, "GtkWidget")
  selection <- as.GdkAtom(selection)

  w <- .RGtkCall("S_gtk_selection_clear_targets", object, selection, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkSelectionConvert <-
function(object, selection, target, time = "GDK_CURRENT_TIME")
{
  checkPtrType(object, "GtkWidget")
  selection <- as.GdkAtom(selection)
  target <- as.GdkAtom(target)
  time <- as.numeric(time)

  w <- .RGtkCall("S_gtk_selection_convert", object, selection, target, time, PACKAGE = "RGtk2")

  return(w)
} 


gtkSelectionDataSetText <-
function(object, str, len = -1)
{
  checkPtrType(object, "GtkSelectionData")
  str <- as.character(str)
  len <- as.integer(len)

  w <- .RGtkCall("S_gtk_selection_data_set_text", object, str, len, PACKAGE = "RGtk2")

  return(w)
} 


gtkSelectionDataGetText <-
function(object)
{
  checkPtrType(object, "GtkSelectionData")

  w <- .RGtkCall("S_gtk_selection_data_get_text", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkSelectionDataGetTargets <-
function(object)
{
  checkPtrType(object, "GtkSelectionData")

  w <- .RGtkCall("S_gtk_selection_data_get_targets", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkSelectionDataTargetsIncludeText <-
function(object)
{
  checkPtrType(object, "GtkSelectionData")

  w <- .RGtkCall("S_gtk_selection_data_targets_include_text", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkSelectionRemoveAll <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_selection_remove_all", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkSelectionClear <-
function(object, event)
{
  if(getOption("depwarn"))
    .Deprecated("a chain up from 'selection-clear-event' handler", "RGtk2")

  checkPtrType(object, "GtkWidget")
  checkPtrType(event, "GdkEventSelection")

  w <- .RGtkCall("S_gtk_selection_clear", object, event, PACKAGE = "RGtk2")

  return(w)
} 


gtkSelectionDataGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_selection_data_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkSelectionDataCopy <-
function(object)
{
  checkPtrType(object, "GtkSelectionData")

  w <- .RGtkCall("S_gtk_selection_data_copy", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkSelectionDataSetPixbuf <-
function(object, pixbuf)
{
  checkPtrType(object, "GtkSelectionData")
  checkPtrType(pixbuf, "GdkPixbuf")

  w <- .RGtkCall("S_gtk_selection_data_set_pixbuf", object, pixbuf, PACKAGE = "RGtk2")

  return(w)
} 


gtkSelectionDataGetPixbuf <-
function(object)
{
  checkPtrType(object, "GtkSelectionData")

  w <- .RGtkCall("S_gtk_selection_data_get_pixbuf", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkSelectionDataSetUris <-
function(object, uris)
{
  checkPtrType(object, "GtkSelectionData")
  uris <- as.list(as.character(uris))

  w <- .RGtkCall("S_gtk_selection_data_set_uris", object, uris, PACKAGE = "RGtk2")

  return(w)
} 


gtkSelectionDataGetUris <-
function(object)
{
  checkPtrType(object, "GtkSelectionData")

  w <- .RGtkCall("S_gtk_selection_data_get_uris", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkSelectionDataTargetsIncludeImage <-
function(object, writable)
{
  checkPtrType(object, "GtkSelectionData")
  writable <- as.logical(writable)

  w <- .RGtkCall("S_gtk_selection_data_targets_include_image", object, writable, PACKAGE = "RGtk2")

  return(w)
} 


gtkSeparatorGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_separator_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkSeparatorMenuItemGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_separator_menu_item_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkSeparatorMenuItemNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_separator_menu_item_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkSeparatorToolItemGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_separator_tool_item_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkSeparatorToolItemNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_separator_tool_item_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkSeparatorToolItemGetDraw <-
function(object)
{
  checkPtrType(object, "GtkSeparatorToolItem")

  w <- .RGtkCall("S_gtk_separator_tool_item_get_draw", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkSeparatorToolItemSetDraw <-
function(object, draw)
{
  checkPtrType(object, "GtkSeparatorToolItem")
  draw <- as.logical(draw)

  w <- .RGtkCall("S_gtk_separator_tool_item_set_draw", object, draw, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkSettingsGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_settings_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkSettingsGetDefault <-
function()
{
  

  w <- .RGtkCall("S_gtk_settings_get_default", PACKAGE = "RGtk2")

  return(w)
} 


gtkSettingsGetForScreen <-
function(screen)
{
  checkPtrType(screen, "GdkScreen")

  w <- .RGtkCall("S_gtk_settings_get_for_screen", screen, PACKAGE = "RGtk2")

  return(w)
} 


gtkSettingsInstallProperty <-
function(pspec)
{
  pspec <- as.GParamSpec(pspec)

  w <- .RGtkCall("S_gtk_settings_install_property", pspec, PACKAGE = "RGtk2")

  return(w)
} 


gtkRcPropertyParseColor <-
function(pspec, gstring)
{
  pspec <- as.GParamSpec(pspec)
  gstring <- as.GString(gstring)

  w <- .RGtkCall("S_gtk_rc_property_parse_color", pspec, gstring, PACKAGE = "RGtk2")

  return(w)
} 


gtkRcPropertyParseEnum <-
function(pspec, gstring)
{
  pspec <- as.GParamSpec(pspec)
  gstring <- as.GString(gstring)

  w <- .RGtkCall("S_gtk_rc_property_parse_enum", pspec, gstring, PACKAGE = "RGtk2")

  return(w)
} 


gtkRcPropertyParseFlags <-
function(pspec, gstring)
{
  pspec <- as.GParamSpec(pspec)
  gstring <- as.GString(gstring)

  w <- .RGtkCall("S_gtk_rc_property_parse_flags", pspec, gstring, PACKAGE = "RGtk2")

  return(w)
} 


gtkRcPropertyParseRequisition <-
function(pspec, gstring)
{
  pspec <- as.GParamSpec(pspec)
  gstring <- as.GString(gstring)

  w <- .RGtkCall("S_gtk_rc_property_parse_requisition", pspec, gstring, PACKAGE = "RGtk2")

  return(w)
} 


gtkRcPropertyParseBorder <-
function(pspec, gstring)
{
  pspec <- as.GParamSpec(pspec)
  gstring <- as.GString(gstring)

  w <- .RGtkCall("S_gtk_rc_property_parse_border", pspec, gstring, PACKAGE = "RGtk2")

  return(w)
} 


gtkSettingsSetPropertyValue <-
function(object, name, svalue)
{
  checkPtrType(object, "GtkSettings")
  name <- as.character(name)
  svalue <- as.GtkSettingsValue(svalue)

  w <- .RGtkCall("S_gtk_settings_set_property_value", object, name, svalue, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkSettingsSetStringProperty <-
function(object, name, v.string, origin)
{
  checkPtrType(object, "GtkSettings")
  name <- as.character(name)
  v.string <- as.character(v.string)
  origin <- as.character(origin)

  w <- .RGtkCall("S_gtk_settings_set_string_property", object, name, v.string, origin, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkSettingsSetLongProperty <-
function(object, name, v.long, origin)
{
  checkPtrType(object, "GtkSettings")
  name <- as.character(name)
  v.long <- as.numeric(v.long)
  origin <- as.character(origin)

  w <- .RGtkCall("S_gtk_settings_set_long_property", object, name, v.long, origin, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkSettingsSetDoubleProperty <-
function(object, name, v.double, origin)
{
  checkPtrType(object, "GtkSettings")
  name <- as.character(name)
  v.double <- as.numeric(v.double)
  origin <- as.character(origin)

  w <- .RGtkCall("S_gtk_settings_set_double_property", object, name, v.double, origin, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkSizeGroupGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_size_group_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkSizeGroupNew <-
function(mode = NULL)
{
  

  w <- .RGtkCall("S_gtk_size_group_new", mode, PACKAGE = "RGtk2")

  return(w)
} 


gtkSizeGroupSetMode <-
function(object, mode)
{
  checkPtrType(object, "GtkSizeGroup")
  

  w <- .RGtkCall("S_gtk_size_group_set_mode", object, mode, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkSizeGroupGetMode <-
function(object)
{
  checkPtrType(object, "GtkSizeGroup")

  w <- .RGtkCall("S_gtk_size_group_get_mode", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkSizeGroupSetIgnoreHidden <-
function(object, ignore.hidden)
{
  checkPtrType(object, "GtkSizeGroup")
  ignore.hidden <- as.logical(ignore.hidden)

  w <- .RGtkCall("S_gtk_size_group_set_ignore_hidden", object, ignore.hidden, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkSizeGroupGetIgnoreHidden <-
function(object)
{
  checkPtrType(object, "GtkSizeGroup")

  w <- .RGtkCall("S_gtk_size_group_get_ignore_hidden", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkSizeGroupAddWidget <-
function(object, widget)
{
  checkPtrType(object, "GtkSizeGroup")
  checkPtrType(widget, "GtkWidget")

  w <- .RGtkCall("S_gtk_size_group_add_widget", object, widget, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkSizeGroupRemoveWidget <-
function(object, widget)
{
  checkPtrType(object, "GtkSizeGroup")
  checkPtrType(widget, "GtkWidget")

  w <- .RGtkCall("S_gtk_size_group_remove_widget", object, widget, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkSocketGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_socket_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkSocketNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_socket_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkSocketAddId <-
function(object, window.id)
{
  checkPtrType(object, "GtkSocket")
  window.id <- as.GdkNativeWindow(window.id)

  w <- .RGtkCall("S_gtk_socket_add_id", object, window.id, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkSocketGetId <-
function(object)
{
  checkPtrType(object, "GtkSocket")

  w <- .RGtkCall("S_gtk_socket_get_id", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkSocketSteal <-
function(object, wid)
{
  checkPtrType(object, "GtkSocket")
  wid <- as.GdkNativeWindow(wid)

  w <- .RGtkCall("S_gtk_socket_steal", object, wid, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkSpinButtonGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_spin_button_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkSpinButtonConfigure <-
function(object, adjustment = NULL, climb.rate, digits)
{
  checkPtrType(object, "GtkSpinButton")
  if (!is.null( adjustment )) checkPtrType(adjustment, "GtkAdjustment")
  climb.rate <- as.numeric(climb.rate)
  digits <- as.numeric(digits)

  w <- .RGtkCall("S_gtk_spin_button_configure", object, adjustment, climb.rate, digits, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkSpinButtonNew <-
function(adjustment = NULL, climb.rate = NULL, digits = NULL, show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_spin_button_new", adjustment, climb.rate, digits, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkSpinButtonNewWithRange <-
function(min, max, step, show = TRUE)
{
  min <- as.numeric(min)
  max <- as.numeric(max)
  step <- as.numeric(step)

  w <- .RGtkCall("S_gtk_spin_button_new_with_range", min, max, step, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkSpinButtonSetAdjustment <-
function(object, adjustment)
{
  checkPtrType(object, "GtkSpinButton")
  checkPtrType(adjustment, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_spin_button_set_adjustment", object, adjustment, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkSpinButtonGetAdjustment <-
function(object)
{
  checkPtrType(object, "GtkSpinButton")

  w <- .RGtkCall("S_gtk_spin_button_get_adjustment", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkSpinButtonSetDigits <-
function(object, digits)
{
  checkPtrType(object, "GtkSpinButton")
  digits <- as.numeric(digits)

  w <- .RGtkCall("S_gtk_spin_button_set_digits", object, digits, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkSpinButtonGetDigits <-
function(object)
{
  checkPtrType(object, "GtkSpinButton")

  w <- .RGtkCall("S_gtk_spin_button_get_digits", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkSpinButtonSetIncrements <-
function(object, step, page)
{
  checkPtrType(object, "GtkSpinButton")
  step <- as.numeric(step)
  page <- as.numeric(page)

  w <- .RGtkCall("S_gtk_spin_button_set_increments", object, step, page, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkSpinButtonGetIncrements <-
function(object)
{
  checkPtrType(object, "GtkSpinButton")

  w <- .RGtkCall("S_gtk_spin_button_get_increments", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkSpinButtonSetRange <-
function(object, min, max)
{
  checkPtrType(object, "GtkSpinButton")
  min <- as.numeric(min)
  max <- as.numeric(max)

  w <- .RGtkCall("S_gtk_spin_button_set_range", object, min, max, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkSpinButtonGetRange <-
function(object)
{
  checkPtrType(object, "GtkSpinButton")

  w <- .RGtkCall("S_gtk_spin_button_get_range", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkSpinButtonGetValue <-
function(object)
{
  checkPtrType(object, "GtkSpinButton")

  w <- .RGtkCall("S_gtk_spin_button_get_value", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkSpinButtonGetValueAsInt <-
function(object)
{
  checkPtrType(object, "GtkSpinButton")

  w <- .RGtkCall("S_gtk_spin_button_get_value_as_int", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkSpinButtonSetValue <-
function(object, value)
{
  checkPtrType(object, "GtkSpinButton")
  value <- as.numeric(value)

  w <- .RGtkCall("S_gtk_spin_button_set_value", object, value, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkSpinButtonSetUpdatePolicy <-
function(object, policy)
{
  checkPtrType(object, "GtkSpinButton")
  

  w <- .RGtkCall("S_gtk_spin_button_set_update_policy", object, policy, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkSpinButtonGetUpdatePolicy <-
function(object)
{
  checkPtrType(object, "GtkSpinButton")

  w <- .RGtkCall("S_gtk_spin_button_get_update_policy", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkSpinButtonSetNumeric <-
function(object, numeric)
{
  checkPtrType(object, "GtkSpinButton")
  numeric <- as.logical(numeric)

  w <- .RGtkCall("S_gtk_spin_button_set_numeric", object, numeric, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkSpinButtonGetNumeric <-
function(object)
{
  checkPtrType(object, "GtkSpinButton")

  w <- .RGtkCall("S_gtk_spin_button_get_numeric", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkSpinButtonSpin <-
function(object, direction, increment)
{
  checkPtrType(object, "GtkSpinButton")
  
  increment <- as.numeric(increment)

  w <- .RGtkCall("S_gtk_spin_button_spin", object, direction, increment, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkSpinButtonSetWrap <-
function(object, wrap)
{
  checkPtrType(object, "GtkSpinButton")
  wrap <- as.logical(wrap)

  w <- .RGtkCall("S_gtk_spin_button_set_wrap", object, wrap, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkSpinButtonGetWrap <-
function(object)
{
  checkPtrType(object, "GtkSpinButton")

  w <- .RGtkCall("S_gtk_spin_button_get_wrap", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkSpinButtonSetSnapToTicks <-
function(object, snap.to.ticks)
{
  checkPtrType(object, "GtkSpinButton")
  snap.to.ticks <- as.logical(snap.to.ticks)

  w <- .RGtkCall("S_gtk_spin_button_set_snap_to_ticks", object, snap.to.ticks, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkSpinButtonGetSnapToTicks <-
function(object)
{
  checkPtrType(object, "GtkSpinButton")

  w <- .RGtkCall("S_gtk_spin_button_get_snap_to_ticks", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkSpinButtonUpdate <-
function(object)
{
  checkPtrType(object, "GtkSpinButton")

  w <- .RGtkCall("S_gtk_spin_button_update", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkStatusbarGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_statusbar_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkStatusbarNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_statusbar_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkStatusbarGetContextId <-
function(object, context.description)
{
  checkPtrType(object, "GtkStatusbar")
  context.description <- as.character(context.description)

  w <- .RGtkCall("S_gtk_statusbar_get_context_id", object, context.description, PACKAGE = "RGtk2")

  return(w)
} 


gtkStatusbarPush <-
function(object, context.id, text)
{
  checkPtrType(object, "GtkStatusbar")
  context.id <- as.numeric(context.id)
  text <- as.character(text)

  w <- .RGtkCall("S_gtk_statusbar_push", object, context.id, text, PACKAGE = "RGtk2")

  return(w)
} 


gtkStatusbarPop <-
function(object, context.id)
{
  checkPtrType(object, "GtkStatusbar")
  context.id <- as.numeric(context.id)

  w <- .RGtkCall("S_gtk_statusbar_pop", object, context.id, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkStatusbarRemove <-
function(object, context.id, message.id)
{
  checkPtrType(object, "GtkStatusbar")
  context.id <- as.numeric(context.id)
  message.id <- as.numeric(message.id)

  w <- .RGtkCall("S_gtk_statusbar_remove", object, context.id, message.id, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkStatusbarSetHasResizeGrip <-
function(object, setting)
{
  checkPtrType(object, "GtkStatusbar")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_statusbar_set_has_resize_grip", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkStatusbarGetHasResizeGrip <-
function(object)
{
  checkPtrType(object, "GtkStatusbar")

  w <- .RGtkCall("S_gtk_statusbar_get_has_resize_grip", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkStockAdd <-
function(items)
{
  items <- lapply(items, function(x) { x <- as.GtkStockItem(x); x })

  w <- .RGtkCall("S_gtk_stock_add", items, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkStockAddStatic <-
function(items)
{
  items <- lapply(items, function(x) { x <- as.GtkStockItem(x); x })

  w <- .RGtkCall("S_gtk_stock_add_static", items, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkStockLookup <-
function(stock.id)
{
  stock.id <- as.character(stock.id)

  w <- .RGtkCall("S_gtk_stock_lookup", stock.id, PACKAGE = "RGtk2")

  return(w)
} 


gtkStockListIds <-
function()
{
  

  w <- .RGtkCall("S_gtk_stock_list_ids", PACKAGE = "RGtk2")

  return(w)
} 


gtkStockItemCopy <-
function(object)
{
  object <- as.GtkStockItem(object)

  w <- .RGtkCall("S_gtk_stock_item_copy", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkStockSetTranslateFunc <-
function(domain, func, data)
{
  domain <- as.character(domain)
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_stock_set_translate_func", domain, func, data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkStyleGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_style_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkStyleNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_style_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkStyleCopy <-
function(object)
{
  checkPtrType(object, "GtkStyle")

  w <- .RGtkCall("S_gtk_style_copy", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkStyleAttach <-
function(object, window)
{
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")

  w <- .RGtkCall("S_gtk_style_attach", object, window, PACKAGE = "RGtk2")

  return(w)
} 


gtkStyleDetach <-
function(object)
{
  checkPtrType(object, "GtkStyle")

  w <- .RGtkCall("S_gtk_style_detach", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkStyleGetFont <-
function(object)
{
  checkPtrType(object, "GtkStyle")

  w <- .RGtkCall("S_gtk_style_get_font", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkStyleSetFont <-
function(object, font)
{
  checkPtrType(object, "GtkStyle")
  checkPtrType(font, "GdkFont")

  w <- .RGtkCall("S_gtk_style_set_font", object, font, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkStyleSetBackground <-
function(object, window, state.type)
{
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  

  w <- .RGtkCall("S_gtk_style_set_background", object, window, state.type, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkStyleApplyDefaultBackground <-
function(object, window, set.bg, state.type, area = NULL, x, y, width, height)
{
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  set.bg <- as.logical(set.bg)
  
  if (!is.null( area )) area <- as.GdkRectangle(area)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_style_apply_default_background", object, window, set.bg, state.type, area, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkStyleLookupIconSet <-
function(object, stock.id)
{
  checkPtrType(object, "GtkStyle")
  stock.id <- as.character(stock.id)

  w <- .RGtkCall("S_gtk_style_lookup_icon_set", object, stock.id, PACKAGE = "RGtk2")

  return(w)
} 


gtkStyleRenderIcon <-
function(object, source, direction, state, size, widget = NULL, detail = NULL)
{
  checkPtrType(object, "GtkStyle")
  checkPtrType(source, "GtkIconSource")
  
  
  
  if (!is.null( widget )) checkPtrType(widget, "GtkWidget")
  if (!is.null( detail )) detail <- as.character(detail)

  w <- .RGtkCall("S_gtk_style_render_icon", object, source, direction, state, size, widget, detail, PACKAGE = "RGtk2")

  return(w)
} 


gtkDrawHline <-
function(object, window, state.type, x1, x2, y)
{
  if(getOption("depwarn"))
    .Deprecated("gtkPaintHline", "RGtk2")

  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  x1 <- as.integer(x1)
  x2 <- as.integer(x2)
  y <- as.integer(y)

  w <- .RGtkCall("S_gtk_draw_hline", object, window, state.type, x1, x2, y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDrawVline <-
function(object, window, state.type, y1, y2, x)
{
  if(getOption("depwarn"))
    .Deprecated("gtkPaintVline", "RGtk2")

  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  y1 <- as.integer(y1)
  y2 <- as.integer(y2)
  x <- as.integer(x)

  w <- .RGtkCall("S_gtk_draw_vline", object, window, state.type, y1, y2, x, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDrawShadow <-
function(object, window, state.type, shadow.type, x, y, width, height)
{
  if(getOption("depwarn"))
    .Deprecated("gtkPaintShadow", "RGtk2")

  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_draw_shadow", object, window, state.type, shadow.type, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDrawPolygon <-
function(object, window, state.type, shadow.type, points, fill)
{
  if(getOption("depwarn"))
    .Deprecated("gtkPaintPolygon", "RGtk2")

  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  points <- lapply(points, function(x) { x <- as.GdkPoint(x); x })
  fill <- as.logical(fill)

  w <- .RGtkCall("S_gtk_draw_polygon", object, window, state.type, shadow.type, points, fill, PACKAGE = "RGtk2")

  return(w)
} 


gtkDrawArrow <-
function(object, window, state.type, shadow.type, arrow.type, fill, x, y, width, height)
{
  if(getOption("depwarn"))
    .Deprecated("gtkPaintArrow", "RGtk2")

  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  
  fill <- as.logical(fill)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_draw_arrow", object, window, state.type, shadow.type, arrow.type, fill, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDrawDiamond <-
function(object, window, state.type, shadow.type, x, y, width, height)
{
  if(getOption("depwarn"))
    .Deprecated("gtkPaintDiamond", "RGtk2")

  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_draw_diamond", object, window, state.type, shadow.type, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDrawBox <-
function(object, window, state.type, shadow.type, x, y, width, height)
{
  if(getOption("depwarn"))
    .Deprecated("gtkPaintBox", "RGtk2")

  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_draw_box", object, window, state.type, shadow.type, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDrawFlatBox <-
function(object, window, state.type, shadow.type, x, y, width, height)
{
  if(getOption("depwarn"))
    .Deprecated("gtkPaintFlatBox", "RGtk2")

  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_draw_flat_box", object, window, state.type, shadow.type, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDrawCheck <-
function(object, window, state.type, shadow.type, x, y, width, height)
{
  if(getOption("depwarn"))
    .Deprecated("gtkPaintCheck", "RGtk2")

  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_draw_check", object, window, state.type, shadow.type, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDrawOption <-
function(object, window, state.type, shadow.type, x, y, width, height)
{
  if(getOption("depwarn"))
    .Deprecated("gtkPaintOption", "RGtk2")

  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_draw_option", object, window, state.type, shadow.type, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDrawTab <-
function(object, window, state.type, shadow.type, x, y, width, height)
{
  if(getOption("depwarn"))
    .Deprecated("gtkPaintTab", "RGtk2")

  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_draw_tab", object, window, state.type, shadow.type, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDrawShadowGap <-
function(object, window, state.type, shadow.type, x, y, width, height, gap.side, gap.x, gap.width)
{
  if(getOption("depwarn"))
    .Deprecated("gtkPaintShadowGap", "RGtk2")

  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)
  
  gap.x <- as.integer(gap.x)
  gap.width <- as.integer(gap.width)

  w <- .RGtkCall("S_gtk_draw_shadow_gap", object, window, state.type, shadow.type, x, y, width, height, gap.side, gap.x, gap.width, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDrawBoxGap <-
function(object, window, state.type, shadow.type, x, y, width, height, gap.side, gap.x, gap.width)
{
  if(getOption("depwarn"))
    .Deprecated("gtkPaintBoxGap", "RGtk2")

  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)
  
  gap.x <- as.integer(gap.x)
  gap.width <- as.integer(gap.width)

  w <- .RGtkCall("S_gtk_draw_box_gap", object, window, state.type, shadow.type, x, y, width, height, gap.side, gap.x, gap.width, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDrawExtension <-
function(object, window, state.type, shadow.type, x, y, width, height, gap.side)
{
  if(getOption("depwarn"))
    .Deprecated("gtkPaintExtension", "RGtk2")

  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)
  

  w <- .RGtkCall("S_gtk_draw_extension", object, window, state.type, shadow.type, x, y, width, height, gap.side, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDrawFocus <-
function(object, window, x, y, width, height)
{
  if(getOption("depwarn"))
    .Deprecated("gtkPaintFocus", "RGtk2")

  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_draw_focus", object, window, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDrawSlider <-
function(object, window, state.type, shadow.type, x, y, width, height, orientation)
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)
  

  w <- .RGtkCall("S_gtk_draw_slider", object, window, state.type, shadow.type, x, y, width, height, orientation, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDrawHandle <-
function(object, window, state.type, shadow.type, x, y, width, height, orientation)
{
  if(getOption("depwarn"))
    .Deprecated("gtkPaintHandle", "RGtk2")

  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)
  

  w <- .RGtkCall("S_gtk_draw_handle", object, window, state.type, shadow.type, x, y, width, height, orientation, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDrawExpander <-
function(object, window, state.type, x, y, is.open)
{
  if(getOption("depwarn"))
    .Deprecated("gtkPaintExpander", "RGtk2")

  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  x <- as.integer(x)
  y <- as.integer(y)
  is.open <- as.logical(is.open)

  w <- .RGtkCall("S_gtk_draw_expander", object, window, state.type, x, y, is.open, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDrawLayout <-
function(object, window, state.type, use.text, x, y, layout)
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  use.text <- as.logical(use.text)
  x <- as.integer(x)
  y <- as.integer(y)
  checkPtrType(layout, "PangoLayout")

  w <- .RGtkCall("S_gtk_draw_layout", object, window, state.type, use.text, x, y, layout, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDrawResizeGrip <-
function(object, window, state.type, edge, x, y, width, height)
{
  if(getOption("depwarn"))
    .Deprecated("gtkPaintResizeGrip", "RGtk2")

  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_draw_resize_grip", object, window, state.type, edge, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPaintHline <-
function(object, window, state.type, area = NULL, widget = NULL, detail = NULL, x1, x2, y)
{
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  if (!is.null( area )) area <- as.GdkRectangle(area)
  if (!is.null( widget )) checkPtrType(widget, "GtkWidget")
  if (!is.null( detail )) detail <- as.character(detail)
  x1 <- as.integer(x1)
  x2 <- as.integer(x2)
  y <- as.integer(y)

  w <- .RGtkCall("S_gtk_paint_hline", object, window, state.type, area, widget, detail, x1, x2, y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPaintVline <-
function(object, window, state.type, area = NULL, widget = NULL, detail = NULL, y1, y2, x)
{
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  if (!is.null( area )) area <- as.GdkRectangle(area)
  if (!is.null( widget )) checkPtrType(widget, "GtkWidget")
  if (!is.null( detail )) detail <- as.character(detail)
  y1 <- as.integer(y1)
  y2 <- as.integer(y2)
  x <- as.integer(x)

  w <- .RGtkCall("S_gtk_paint_vline", object, window, state.type, area, widget, detail, y1, y2, x, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPaintShadow <-
function(object, window, state.type, shadow.type, area = NULL, widget = NULL, detail = NULL, x, y, width, height)
{
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  if (!is.null( area )) area <- as.GdkRectangle(area)
  if (!is.null( widget )) checkPtrType(widget, "GtkWidget")
  if (!is.null( detail )) detail <- as.character(detail)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_paint_shadow", object, window, state.type, shadow.type, area, widget, detail, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPaintPolygon <-
function(object, window, state.type, shadow.type, area = NULL, widget = NULL, detail = NULL, points, fill)
{
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  if (!is.null( area )) area <- as.GdkRectangle(area)
  if (!is.null( widget )) checkPtrType(widget, "GtkWidget")
  if (!is.null( detail )) detail <- as.character(detail)
  points <- lapply(points, function(x) { x <- as.GdkPoint(x); x })
  fill <- as.logical(fill)

  w <- .RGtkCall("S_gtk_paint_polygon", object, window, state.type, shadow.type, area, widget, detail, points, fill, PACKAGE = "RGtk2")

  return(w)
} 


gtkPaintArrow <-
function(object, window, state.type, shadow.type, area = NULL, widget = NULL, detail = NULL, arrow.type, fill, x, y, width, height)
{
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  if (!is.null( area )) area <- as.GdkRectangle(area)
  if (!is.null( widget )) checkPtrType(widget, "GtkWidget")
  if (!is.null( detail )) detail <- as.character(detail)
  
  fill <- as.logical(fill)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_paint_arrow", object, window, state.type, shadow.type, area, widget, detail, arrow.type, fill, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPaintDiamond <-
function(object, window, state.type, shadow.type, area = NULL, widget = NULL, detail = NULL, x, y, width, height)
{
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  if (!is.null( area )) area <- as.GdkRectangle(area)
  if (!is.null( widget )) checkPtrType(widget, "GtkWidget")
  if (!is.null( detail )) detail <- as.character(detail)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_paint_diamond", object, window, state.type, shadow.type, area, widget, detail, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPaintBox <-
function(object, window, state.type, shadow.type, area = NULL, widget = NULL, detail = NULL, x, y, width, height)
{
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  if (!is.null( area )) area <- as.GdkRectangle(area)
  if (!is.null( widget )) checkPtrType(widget, "GtkWidget")
  if (!is.null( detail )) detail <- as.character(detail)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_paint_box", object, window, state.type, shadow.type, area, widget, detail, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPaintFlatBox <-
function(object, window, state.type, shadow.type, area = NULL, widget = NULL, detail = NULL, x, y, width, height)
{
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  if (!is.null( area )) area <- as.GdkRectangle(area)
  if (!is.null( widget )) checkPtrType(widget, "GtkWidget")
  if (!is.null( detail )) detail <- as.character(detail)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_paint_flat_box", object, window, state.type, shadow.type, area, widget, detail, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPaintCheck <-
function(object, window, state.type, shadow.type, area = NULL, widget = NULL, detail = NULL, x, y, width, height)
{
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  if (!is.null( area )) area <- as.GdkRectangle(area)
  if (!is.null( widget )) checkPtrType(widget, "GtkWidget")
  if (!is.null( detail )) detail <- as.character(detail)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_paint_check", object, window, state.type, shadow.type, area, widget, detail, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPaintOption <-
function(object, window, state.type, shadow.type, area = NULL, widget = NULL, detail = NULL, x, y, width, height)
{
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  if (!is.null( area )) area <- as.GdkRectangle(area)
  if (!is.null( widget )) checkPtrType(widget, "GtkWidget")
  if (!is.null( detail )) detail <- as.character(detail)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_paint_option", object, window, state.type, shadow.type, area, widget, detail, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPaintTab <-
function(object, window, state.type, shadow.type, area = NULL, widget = NULL, detail = NULL, x, y, width, height)
{
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  if (!is.null( area )) area <- as.GdkRectangle(area)
  if (!is.null( widget )) checkPtrType(widget, "GtkWidget")
  if (!is.null( detail )) detail <- as.character(detail)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_paint_tab", object, window, state.type, shadow.type, area, widget, detail, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPaintShadowGap <-
function(object, window, state.type, shadow.type, area = NULL, widget = NULL, detail = NULL, x, y, width, height, gap.side, gap.x, gap.width)
{
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  if (!is.null( area )) area <- as.GdkRectangle(area)
  if (!is.null( widget )) checkPtrType(widget, "GtkWidget")
  if (!is.null( detail )) detail <- as.character(detail)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)
  
  gap.x <- as.integer(gap.x)
  gap.width <- as.integer(gap.width)

  w <- .RGtkCall("S_gtk_paint_shadow_gap", object, window, state.type, shadow.type, area, widget, detail, x, y, width, height, gap.side, gap.x, gap.width, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPaintBoxGap <-
function(object, window, state.type, shadow.type, area = NULL, widget = NULL, detail = NULL, x, y, width, height, gap.side, gap.x, gap.width)
{
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  if (!is.null( area )) area <- as.GdkRectangle(area)
  if (!is.null( widget )) checkPtrType(widget, "GtkWidget")
  if (!is.null( detail )) detail <- as.character(detail)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)
  
  gap.x <- as.integer(gap.x)
  gap.width <- as.integer(gap.width)

  w <- .RGtkCall("S_gtk_paint_box_gap", object, window, state.type, shadow.type, area, widget, detail, x, y, width, height, gap.side, gap.x, gap.width, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPaintExtension <-
function(object, window, state.type, shadow.type, area = NULL, widget = NULL, detail = NULL, x, y, width, height, gap.side)
{
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  if (!is.null( area )) area <- as.GdkRectangle(area)
  if (!is.null( widget )) checkPtrType(widget, "GtkWidget")
  if (!is.null( detail )) detail <- as.character(detail)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)
  

  w <- .RGtkCall("S_gtk_paint_extension", object, window, state.type, shadow.type, area, widget, detail, x, y, width, height, gap.side, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPaintFocus <-
function(object, window, state.type, area = NULL, widget = NULL, detail = NULL, x, y, width, height)
{
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  if (!is.null( area )) area <- as.GdkRectangle(area)
  if (!is.null( widget )) checkPtrType(widget, "GtkWidget")
  if (!is.null( detail )) detail <- as.character(detail)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_paint_focus", object, window, state.type, area, widget, detail, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPaintSlider <-
function(object, window, state.type, shadow.type, area = NULL, widget = NULL, detail = NULL, x, y, width, height, orientation)
{
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  if (!is.null( area )) area <- as.GdkRectangle(area)
  if (!is.null( widget )) checkPtrType(widget, "GtkWidget")
  if (!is.null( detail )) detail <- as.character(detail)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)
  

  w <- .RGtkCall("S_gtk_paint_slider", object, window, state.type, shadow.type, area, widget, detail, x, y, width, height, orientation, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPaintHandle <-
function(object, window, state.type, shadow.type, area = NULL, widget = NULL, detail = NULL, x, y, width, height, orientation)
{
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  
  if (!is.null( area )) area <- as.GdkRectangle(area)
  if (!is.null( widget )) checkPtrType(widget, "GtkWidget")
  if (!is.null( detail )) detail <- as.character(detail)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)
  

  w <- .RGtkCall("S_gtk_paint_handle", object, window, state.type, shadow.type, area, widget, detail, x, y, width, height, orientation, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPaintExpander <-
function(object, window, state.type, area = NULL, widget = NULL, detail = NULL, x, y, expander.style)
{
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  if (!is.null( area )) area <- as.GdkRectangle(area)
  if (!is.null( widget )) checkPtrType(widget, "GtkWidget")
  if (!is.null( detail )) detail <- as.character(detail)
  x <- as.integer(x)
  y <- as.integer(y)
  

  w <- .RGtkCall("S_gtk_paint_expander", object, window, state.type, area, widget, detail, x, y, expander.style, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPaintLayout <-
function(object, window, state.type, use.text, area = NULL, widget = NULL, detail = NULL, x, y, layout)
{
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  use.text <- as.logical(use.text)
  if (!is.null( area )) area <- as.GdkRectangle(area)
  if (!is.null( widget )) checkPtrType(widget, "GtkWidget")
  if (!is.null( detail )) detail <- as.character(detail)
  x <- as.integer(x)
  y <- as.integer(y)
  checkPtrType(layout, "PangoLayout")

  w <- .RGtkCall("S_gtk_paint_layout", object, window, state.type, use.text, area, widget, detail, x, y, layout, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPaintResizeGrip <-
function(object, window, state.type, area = NULL, widget = NULL, detail = NULL, edge, x, y, width, height)
{
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  if (!is.null( area )) area <- as.GdkRectangle(area)
  if (!is.null( widget )) checkPtrType(widget, "GtkWidget")
  if (!is.null( detail )) detail <- as.character(detail)
  
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_paint_resize_grip", object, window, state.type, area, widget, detail, edge, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkBorderGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_border_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkBorderCopy <-
function(object)
{
  checkPtrType(object, "GtkBorder")

  w <- .RGtkCall("S_gtk_border_copy", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkStyleApplyDefaultPixmap <-
function(object, window, set.bg, area, x, y, width, height)
{
  if(getOption("depwarn"))
    .Deprecated("gtkStyleApplyDefaultBackground", "RGtk2")

  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  set.bg <- as.logical(set.bg)
  area <- as.GdkRectangle(area)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_style_apply_default_pixmap", object, window, set.bg, area, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDrawString <-
function(object, window, state.type, x, y, string)
{
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  x <- as.integer(x)
  y <- as.integer(y)
  string <- as.character(string)

  w <- .RGtkCall("S_gtk_draw_string", object, window, state.type, x, y, string, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPaintString <-
function(object, window, state.type, area = NULL, widget = NULL, detail = NULL, x, y, string)
{
  if(getOption("depwarn"))
    .Deprecated("gtkPaintLayout", "RGtk2")

  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  if (!is.null( area )) area <- as.GdkRectangle(area)
  if (!is.null( widget )) checkPtrType(widget, "GtkWidget")
  if (!is.null( detail )) detail <- as.character(detail)
  x <- as.integer(x)
  y <- as.integer(y)
  string <- as.character(string)

  w <- .RGtkCall("S_gtk_paint_string", object, window, state.type, area, widget, detail, x, y, string, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDrawInsertionCursor <-
function(widget, drawable, area = NULL, location, is.primary, direction, draw.arrow)
{
  checkPtrType(widget, "GtkWidget")
  checkPtrType(drawable, "GdkDrawable")
  if (!is.null( area )) area <- as.GdkRectangle(area)
  location <- as.GdkRectangle(location)
  is.primary <- as.logical(is.primary)
  
  draw.arrow <- as.logical(draw.arrow)

  w <- .RGtkCall("S_gtk_draw_insertion_cursor", widget, drawable, area, location, is.primary, direction, draw.arrow, PACKAGE = "RGtk2")

  return(w)
} 


gtkTableGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_table_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkTableNew <-
function(rows = NULL, columns = NULL, homogeneous = NULL, show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_table_new", rows, columns, homogeneous, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkTableResize <-
function(object, rows, columns)
{
  checkPtrType(object, "GtkTable")
  rows <- as.numeric(rows)
  columns <- as.numeric(columns)

  w <- .RGtkCall("S_gtk_table_resize", object, rows, columns, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTableAttach <-
function(object, child, left.attach, right.attach, top.attach, bottom.attach, xoptions = 5, yoptions = 5, xpadding = 0, ypadding = 0)
{
  checkPtrType(object, "GtkTable")
  checkPtrType(child, "GtkWidget")
  left.attach <- as.numeric(left.attach)
  right.attach <- as.numeric(right.attach)
  top.attach <- as.numeric(top.attach)
  bottom.attach <- as.numeric(bottom.attach)
  
  
  xpadding <- as.numeric(xpadding)
  ypadding <- as.numeric(ypadding)

  w <- .RGtkCall("S_gtk_table_attach", object, child, left.attach, right.attach, top.attach, bottom.attach, xoptions, yoptions, xpadding, ypadding, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTableAttachDefaults <-
function(object, widget, left.attach, right.attach, top.attach, bottom.attach)
{
  checkPtrType(object, "GtkTable")
  checkPtrType(widget, "GtkWidget")
  left.attach <- as.numeric(left.attach)
  right.attach <- as.numeric(right.attach)
  top.attach <- as.numeric(top.attach)
  bottom.attach <- as.numeric(bottom.attach)

  w <- .RGtkCall("S_gtk_table_attach_defaults", object, widget, left.attach, right.attach, top.attach, bottom.attach, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTableSetRowSpacing <-
function(object, row, spacing)
{
  checkPtrType(object, "GtkTable")
  row <- as.numeric(row)
  spacing <- as.numeric(spacing)

  w <- .RGtkCall("S_gtk_table_set_row_spacing", object, row, spacing, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTableGetRowSpacing <-
function(object, row)
{
  checkPtrType(object, "GtkTable")
  row <- as.numeric(row)

  w <- .RGtkCall("S_gtk_table_get_row_spacing", object, row, PACKAGE = "RGtk2")

  return(w)
} 


gtkTableSetColSpacing <-
function(object, column, spacing)
{
  checkPtrType(object, "GtkTable")
  column <- as.numeric(column)
  spacing <- as.numeric(spacing)

  w <- .RGtkCall("S_gtk_table_set_col_spacing", object, column, spacing, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTableGetColSpacing <-
function(object, column)
{
  checkPtrType(object, "GtkTable")
  column <- as.numeric(column)

  w <- .RGtkCall("S_gtk_table_get_col_spacing", object, column, PACKAGE = "RGtk2")

  return(w)
} 


gtkTableSetRowSpacings <-
function(object, spacing)
{
  checkPtrType(object, "GtkTable")
  spacing <- as.numeric(spacing)

  w <- .RGtkCall("S_gtk_table_set_row_spacings", object, spacing, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTableGetDefaultRowSpacing <-
function(object)
{
  checkPtrType(object, "GtkTable")

  w <- .RGtkCall("S_gtk_table_get_default_row_spacing", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTableSetColSpacings <-
function(object, spacing)
{
  checkPtrType(object, "GtkTable")
  spacing <- as.numeric(spacing)

  w <- .RGtkCall("S_gtk_table_set_col_spacings", object, spacing, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTableGetDefaultColSpacing <-
function(object)
{
  checkPtrType(object, "GtkTable")

  w <- .RGtkCall("S_gtk_table_get_default_col_spacing", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTableSetHomogeneous <-
function(object, homogeneous)
{
  checkPtrType(object, "GtkTable")
  homogeneous <- as.logical(homogeneous)

  w <- .RGtkCall("S_gtk_table_set_homogeneous", object, homogeneous, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTableGetHomogeneous <-
function(object)
{
  checkPtrType(object, "GtkTable")

  w <- .RGtkCall("S_gtk_table_get_homogeneous", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTearoffMenuItemGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_tearoff_menu_item_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkTearoffMenuItemNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_tearoff_menu_item_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkTextBufferGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_text_buffer_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferNew <-
function(table = NULL)
{
  

  w <- .RGtkCall("S_gtk_text_buffer_new", table, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferGetLineCount <-
function(object)
{
  checkPtrType(object, "GtkTextBuffer")

  w <- .RGtkCall("S_gtk_text_buffer_get_line_count", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferGetCharCount <-
function(object)
{
  checkPtrType(object, "GtkTextBuffer")

  w <- .RGtkCall("S_gtk_text_buffer_get_char_count", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferGetTagTable <-
function(object)
{
  checkPtrType(object, "GtkTextBuffer")

  w <- .RGtkCall("S_gtk_text_buffer_get_tag_table", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferSetText <-
function(object, text, len = -1)
{
  checkPtrType(object, "GtkTextBuffer")
  text <- as.character(text)
  len <- as.integer(len)

  w <- .RGtkCall("S_gtk_text_buffer_set_text", object, text, len, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextBufferInsert <-
function(object, iter, text, len = -1)
{
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(iter, "GtkTextIter")
  text <- as.character(text)
  len <- as.integer(len)

  w <- .RGtkCall("S_gtk_text_buffer_insert", object, iter, text, len, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextBufferInsertAtCursor <-
function(object, text, len = -1)
{
  checkPtrType(object, "GtkTextBuffer")
  text <- as.character(text)
  len <- as.integer(len)

  w <- .RGtkCall("S_gtk_text_buffer_insert_at_cursor", object, text, len, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextBufferInsertRange <-
function(object, iter, start, end)
{
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(iter, "GtkTextIter")
  checkPtrType(start, "GtkTextIter")
  checkPtrType(end, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_buffer_insert_range", object, iter, start, end, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextBufferInsertRangeInteractive <-
function(object, iter, start, end, default.editable)
{
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(iter, "GtkTextIter")
  checkPtrType(start, "GtkTextIter")
  checkPtrType(end, "GtkTextIter")
  default.editable <- as.logical(default.editable)

  w <- .RGtkCall("S_gtk_text_buffer_insert_range_interactive", object, iter, start, end, default.editable, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferDelete <-
function(object, start, end)
{
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(start, "GtkTextIter")
  checkPtrType(end, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_buffer_delete", object, start, end, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextBufferDeleteInteractive <-
function(object, start.iter, end.iter, default.editable)
{
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(start.iter, "GtkTextIter")
  checkPtrType(end.iter, "GtkTextIter")
  default.editable <- as.logical(default.editable)

  w <- .RGtkCall("S_gtk_text_buffer_delete_interactive", object, start.iter, end.iter, default.editable, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferGetText <-
function(object, start, end, include.hidden.chars = TRUE)
{
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(start, "GtkTextIter")
  checkPtrType(end, "GtkTextIter")
  include.hidden.chars <- as.logical(include.hidden.chars)

  w <- .RGtkCall("S_gtk_text_buffer_get_text", object, start, end, include.hidden.chars, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferGetSlice <-
function(object, start, end, include.hidden.chars = TRUE)
{
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(start, "GtkTextIter")
  checkPtrType(end, "GtkTextIter")
  include.hidden.chars <- as.logical(include.hidden.chars)

  w <- .RGtkCall("S_gtk_text_buffer_get_slice", object, start, end, include.hidden.chars, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferInsertPixbuf <-
function(object, iter, pixbuf)
{
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(iter, "GtkTextIter")
  checkPtrType(pixbuf, "GdkPixbuf")

  w <- .RGtkCall("S_gtk_text_buffer_insert_pixbuf", object, iter, pixbuf, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextBufferInsertChildAnchor <-
function(object, iter, anchor)
{
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(iter, "GtkTextIter")
  checkPtrType(anchor, "GtkTextChildAnchor")

  w <- .RGtkCall("S_gtk_text_buffer_insert_child_anchor", object, iter, anchor, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextBufferCreateChildAnchor <-
function(object, iter)
{
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(iter, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_buffer_create_child_anchor", object, iter, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferCreateMark <-
function(object, mark.name = NULL, where, left.gravity = FALSE)
{
  checkPtrType(object, "GtkTextBuffer")
  if (!is.null( mark.name )) mark.name <- as.character(mark.name)
  checkPtrType(where, "GtkTextIter")
  left.gravity <- as.logical(left.gravity)

  w <- .RGtkCall("S_gtk_text_buffer_create_mark", object, mark.name, where, left.gravity, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferMoveMark <-
function(object, mark, where)
{
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(mark, "GtkTextMark")
  checkPtrType(where, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_buffer_move_mark", object, mark, where, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextBufferDeleteMark <-
function(object, mark)
{
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(mark, "GtkTextMark")

  w <- .RGtkCall("S_gtk_text_buffer_delete_mark", object, mark, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextBufferGetMark <-
function(object, name)
{
  checkPtrType(object, "GtkTextBuffer")
  name <- as.character(name)

  w <- .RGtkCall("S_gtk_text_buffer_get_mark", object, name, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferMoveMarkByName <-
function(object, name, where)
{
  checkPtrType(object, "GtkTextBuffer")
  name <- as.character(name)
  checkPtrType(where, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_buffer_move_mark_by_name", object, name, where, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextBufferDeleteMarkByName <-
function(object, name)
{
  checkPtrType(object, "GtkTextBuffer")
  name <- as.character(name)

  w <- .RGtkCall("S_gtk_text_buffer_delete_mark_by_name", object, name, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextBufferGetInsert <-
function(object)
{
  checkPtrType(object, "GtkTextBuffer")

  w <- .RGtkCall("S_gtk_text_buffer_get_insert", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferGetSelectionBound <-
function(object)
{
  checkPtrType(object, "GtkTextBuffer")

  w <- .RGtkCall("S_gtk_text_buffer_get_selection_bound", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferPlaceCursor <-
function(object, where)
{
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(where, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_buffer_place_cursor", object, where, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextBufferSelectRange <-
function(object, ins, bound)
{
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(ins, "GtkTextIter")
  checkPtrType(bound, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_buffer_select_range", object, ins, bound, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextBufferApplyTag <-
function(object, tag, start, end)
{
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(tag, "GtkTextTag")
  checkPtrType(start, "GtkTextIter")
  checkPtrType(end, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_buffer_apply_tag", object, tag, start, end, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextBufferRemoveTag <-
function(object, tag, start, end)
{
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(tag, "GtkTextTag")
  checkPtrType(start, "GtkTextIter")
  checkPtrType(end, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_buffer_remove_tag", object, tag, start, end, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextBufferApplyTagByName <-
function(object, name, start, end)
{
  checkPtrType(object, "GtkTextBuffer")
  name <- as.character(name)
  checkPtrType(start, "GtkTextIter")
  checkPtrType(end, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_buffer_apply_tag_by_name", object, name, start, end, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextBufferRemoveTagByName <-
function(object, name, start, end)
{
  checkPtrType(object, "GtkTextBuffer")
  name <- as.character(name)
  checkPtrType(start, "GtkTextIter")
  checkPtrType(end, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_buffer_remove_tag_by_name", object, name, start, end, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextBufferRemoveAllTags <-
function(object, start, end)
{
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(start, "GtkTextIter")
  checkPtrType(end, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_buffer_remove_all_tags", object, start, end, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextBufferGetIterAtLineOffset <-
function(object, line.number, char.offset)
{
  checkPtrType(object, "GtkTextBuffer")
  line.number <- as.integer(line.number)
  char.offset <- as.integer(char.offset)

  w <- .RGtkCall("S_gtk_text_buffer_get_iter_at_line_offset", object, line.number, char.offset, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferGetIterAtLineIndex <-
function(object, line.number, byte.index)
{
  checkPtrType(object, "GtkTextBuffer")
  line.number <- as.integer(line.number)
  byte.index <- as.integer(byte.index)

  w <- .RGtkCall("S_gtk_text_buffer_get_iter_at_line_index", object, line.number, byte.index, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferGetIterAtOffset <-
function(object, char.offset)
{
  checkPtrType(object, "GtkTextBuffer")
  char.offset <- as.integer(char.offset)

  w <- .RGtkCall("S_gtk_text_buffer_get_iter_at_offset", object, char.offset, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferGetIterAtLine <-
function(object, line.number)
{
  checkPtrType(object, "GtkTextBuffer")
  line.number <- as.integer(line.number)

  w <- .RGtkCall("S_gtk_text_buffer_get_iter_at_line", object, line.number, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferGetStartIter <-
function(object)
{
  checkPtrType(object, "GtkTextBuffer")

  w <- .RGtkCall("S_gtk_text_buffer_get_start_iter", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferGetEndIter <-
function(object)
{
  checkPtrType(object, "GtkTextBuffer")

  w <- .RGtkCall("S_gtk_text_buffer_get_end_iter", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferGetBounds <-
function(object)
{
  checkPtrType(object, "GtkTextBuffer")

  w <- .RGtkCall("S_gtk_text_buffer_get_bounds", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextBufferGetIterAtMark <-
function(object, mark)
{
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(mark, "GtkTextMark")

  w <- .RGtkCall("S_gtk_text_buffer_get_iter_at_mark", object, mark, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferGetIterAtChildAnchor <-
function(object, anchor)
{
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(anchor, "GtkTextChildAnchor")

  w <- .RGtkCall("S_gtk_text_buffer_get_iter_at_child_anchor", object, anchor, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferGetModified <-
function(object)
{
  checkPtrType(object, "GtkTextBuffer")

  w <- .RGtkCall("S_gtk_text_buffer_get_modified", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferSetModified <-
function(object, setting)
{
  checkPtrType(object, "GtkTextBuffer")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_text_buffer_set_modified", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextBufferAddSelectionClipboard <-
function(object, clipboard)
{
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(clipboard, "GtkClipboard")

  w <- .RGtkCall("S_gtk_text_buffer_add_selection_clipboard", object, clipboard, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextBufferRemoveSelectionClipboard <-
function(object, clipboard)
{
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(clipboard, "GtkClipboard")

  w <- .RGtkCall("S_gtk_text_buffer_remove_selection_clipboard", object, clipboard, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextBufferCutClipboard <-
function(object, clipboard, default.editable)
{
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(clipboard, "GtkClipboard")
  default.editable <- as.logical(default.editable)

  w <- .RGtkCall("S_gtk_text_buffer_cut_clipboard", object, clipboard, default.editable, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextBufferCopyClipboard <-
function(object, clipboard)
{
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(clipboard, "GtkClipboard")

  w <- .RGtkCall("S_gtk_text_buffer_copy_clipboard", object, clipboard, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextBufferPasteClipboard <-
function(object, clipboard, override.location = NULL, default.editable)
{
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(clipboard, "GtkClipboard")
  if (!is.null( override.location )) checkPtrType(override.location, "GtkTextIter")
  default.editable <- as.logical(default.editable)

  w <- .RGtkCall("S_gtk_text_buffer_paste_clipboard", object, clipboard, override.location, default.editable, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextBufferGetSelectionBounds <-
function(object)
{
  checkPtrType(object, "GtkTextBuffer")

  w <- .RGtkCall("S_gtk_text_buffer_get_selection_bounds", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferDeleteSelection <-
function(object, interactive, default.editable)
{
  checkPtrType(object, "GtkTextBuffer")
  interactive <- as.logical(interactive)
  default.editable <- as.logical(default.editable)

  w <- .RGtkCall("S_gtk_text_buffer_delete_selection", object, interactive, default.editable, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferBeginUserAction <-
function(object)
{
  checkPtrType(object, "GtkTextBuffer")

  w <- .RGtkCall("S_gtk_text_buffer_begin_user_action", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextBufferEndUserAction <-
function(object)
{
  checkPtrType(object, "GtkTextBuffer")

  w <- .RGtkCall("S_gtk_text_buffer_end_user_action", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextChildAnchorGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_text_child_anchor_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkTextChildAnchorNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_text_child_anchor_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkTextChildAnchorGetWidgets <-
function(object)
{
  checkPtrType(object, "GtkTextChildAnchor")

  w <- .RGtkCall("S_gtk_text_child_anchor_get_widgets", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextChildAnchorGetDeleted <-
function(object)
{
  checkPtrType(object, "GtkTextChildAnchor")

  w <- .RGtkCall("S_gtk_text_child_anchor_get_deleted", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferBackspace <-
function(object, iter, interactive, default.editable)
{
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(iter, "GtkTextIter")
  interactive <- as.logical(interactive)
  default.editable <- as.logical(default.editable)

  w <- .RGtkCall("S_gtk_text_buffer_backspace", object, iter, interactive, default.editable, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterGetBuffer <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_get_buffer", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterCopy <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_copy", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_text_iter_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterGetOffset <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_get_offset", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterGetLine <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_get_line", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterGetLineOffset <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_get_line_offset", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterGetLineIndex <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_get_line_index", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterGetVisibleLineOffset <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_get_visible_line_offset", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterGetVisibleLineIndex <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_get_visible_line_index", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterGetChar <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_get_char", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterGetSlice <-
function(object, end)
{
  checkPtrType(object, "GtkTextIter")
  checkPtrType(end, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_get_slice", object, end, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterGetText <-
function(object, end)
{
  checkPtrType(object, "GtkTextIter")
  checkPtrType(end, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_get_text", object, end, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterGetVisibleSlice <-
function(object, end)
{
  checkPtrType(object, "GtkTextIter")
  checkPtrType(end, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_get_visible_slice", object, end, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterGetVisibleText <-
function(object, end)
{
  checkPtrType(object, "GtkTextIter")
  checkPtrType(end, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_get_visible_text", object, end, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterGetPixbuf <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_get_pixbuf", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterGetMarks <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_get_marks", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterGetChildAnchor <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_get_child_anchor", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterGetToggledTags <-
function(object, toggled.on)
{
  checkPtrType(object, "GtkTextIter")
  toggled.on <- as.logical(toggled.on)

  w <- .RGtkCall("S_gtk_text_iter_get_toggled_tags", object, toggled.on, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterBeginsTag <-
function(object, tag = NULL)
{
  checkPtrType(object, "GtkTextIter")
  if (!is.null( tag )) checkPtrType(tag, "GtkTextTag")

  w <- .RGtkCall("S_gtk_text_iter_begins_tag", object, tag, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterEndsTag <-
function(object, tag = NULL)
{
  checkPtrType(object, "GtkTextIter")
  if (!is.null( tag )) checkPtrType(tag, "GtkTextTag")

  w <- .RGtkCall("S_gtk_text_iter_ends_tag", object, tag, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterTogglesTag <-
function(object, tag = NULL)
{
  checkPtrType(object, "GtkTextIter")
  if (!is.null( tag )) checkPtrType(tag, "GtkTextTag")

  w <- .RGtkCall("S_gtk_text_iter_toggles_tag", object, tag, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterHasTag <-
function(object, tag)
{
  checkPtrType(object, "GtkTextIter")
  checkPtrType(tag, "GtkTextTag")

  w <- .RGtkCall("S_gtk_text_iter_has_tag", object, tag, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterGetTags <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_get_tags", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterEditable <-
function(object, default.setting)
{
  checkPtrType(object, "GtkTextIter")
  default.setting <- as.logical(default.setting)

  w <- .RGtkCall("S_gtk_text_iter_editable", object, default.setting, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterCanInsert <-
function(object, default.editability)
{
  checkPtrType(object, "GtkTextIter")
  default.editability <- as.logical(default.editability)

  w <- .RGtkCall("S_gtk_text_iter_can_insert", object, default.editability, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterStartsWord <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_starts_word", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterEndsWord <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_ends_word", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterInsideWord <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_inside_word", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterStartsSentence <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_starts_sentence", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterEndsSentence <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_ends_sentence", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterInsideSentence <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_inside_sentence", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterStartsLine <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_starts_line", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterEndsLine <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_ends_line", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterIsCursorPosition <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_is_cursor_position", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterGetCharsInLine <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_get_chars_in_line", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterGetBytesInLine <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_get_bytes_in_line", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterGetAttributes <-
function(object, values)
{
  checkPtrType(object, "GtkTextIter")
  checkPtrType(values, "GtkTextAttributes")

  w <- .RGtkCall("S_gtk_text_iter_get_attributes", object, values, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterGetLanguage <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_get_language", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterIsEnd <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_is_end", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterIsStart <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_is_start", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterForwardChar <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_forward_char", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterBackwardChar <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_backward_char", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterForwardChars <-
function(object, count)
{
  checkPtrType(object, "GtkTextIter")
  count <- as.integer(count)

  w <- .RGtkCall("S_gtk_text_iter_forward_chars", object, count, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterBackwardChars <-
function(object, count)
{
  checkPtrType(object, "GtkTextIter")
  count <- as.integer(count)

  w <- .RGtkCall("S_gtk_text_iter_backward_chars", object, count, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterForwardLine <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_forward_line", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterBackwardLine <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_backward_line", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterForwardLines <-
function(object, count)
{
  checkPtrType(object, "GtkTextIter")
  count <- as.integer(count)

  w <- .RGtkCall("S_gtk_text_iter_forward_lines", object, count, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterBackwardLines <-
function(object, count)
{
  checkPtrType(object, "GtkTextIter")
  count <- as.integer(count)

  w <- .RGtkCall("S_gtk_text_iter_backward_lines", object, count, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterForwardWordEnd <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_forward_word_end", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterBackwardWordStart <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_backward_word_start", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterForwardWordEnds <-
function(object, count)
{
  checkPtrType(object, "GtkTextIter")
  count <- as.integer(count)

  w <- .RGtkCall("S_gtk_text_iter_forward_word_ends", object, count, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterBackwardWordStarts <-
function(object, count)
{
  checkPtrType(object, "GtkTextIter")
  count <- as.integer(count)

  w <- .RGtkCall("S_gtk_text_iter_backward_word_starts", object, count, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterForwardVisibleLine <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_forward_visible_line", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterBackwardVisibleLine <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_backward_visible_line", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterForwardVisibleLines <-
function(object, count)
{
  checkPtrType(object, "GtkTextIter")
  count <- as.integer(count)

  w <- .RGtkCall("S_gtk_text_iter_forward_visible_lines", object, count, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterBackwardVisibleLines <-
function(object, count)
{
  checkPtrType(object, "GtkTextIter")
  count <- as.integer(count)

  w <- .RGtkCall("S_gtk_text_iter_backward_visible_lines", object, count, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterForwardVisibleWordEnd <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_forward_visible_word_end", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterBackwardVisibleWordStart <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_backward_visible_word_start", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterForwardVisibleWordEnds <-
function(object, count)
{
  checkPtrType(object, "GtkTextIter")
  count <- as.integer(count)

  w <- .RGtkCall("S_gtk_text_iter_forward_visible_word_ends", object, count, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterBackwardVisibleWordStarts <-
function(object, count)
{
  checkPtrType(object, "GtkTextIter")
  count <- as.integer(count)

  w <- .RGtkCall("S_gtk_text_iter_backward_visible_word_starts", object, count, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterForwardSentenceEnd <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_forward_sentence_end", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterBackwardSentenceStart <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_backward_sentence_start", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterForwardSentenceEnds <-
function(object, count)
{
  checkPtrType(object, "GtkTextIter")
  count <- as.integer(count)

  w <- .RGtkCall("S_gtk_text_iter_forward_sentence_ends", object, count, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterBackwardSentenceStarts <-
function(object, count)
{
  checkPtrType(object, "GtkTextIter")
  count <- as.integer(count)

  w <- .RGtkCall("S_gtk_text_iter_backward_sentence_starts", object, count, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterForwardCursorPosition <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_forward_cursor_position", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterBackwardCursorPosition <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_backward_cursor_position", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterForwardCursorPositions <-
function(object, count)
{
  checkPtrType(object, "GtkTextIter")
  count <- as.integer(count)

  w <- .RGtkCall("S_gtk_text_iter_forward_cursor_positions", object, count, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterBackwardCursorPositions <-
function(object, count)
{
  checkPtrType(object, "GtkTextIter")
  count <- as.integer(count)

  w <- .RGtkCall("S_gtk_text_iter_backward_cursor_positions", object, count, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterForwardVisibleCursorPosition <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_forward_visible_cursor_position", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterBackwardVisibleCursorPosition <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_backward_visible_cursor_position", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterForwardVisibleCursorPositions <-
function(object, count)
{
  checkPtrType(object, "GtkTextIter")
  count <- as.integer(count)

  w <- .RGtkCall("S_gtk_text_iter_forward_visible_cursor_positions", object, count, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterBackwardVisibleCursorPositions <-
function(object, count)
{
  checkPtrType(object, "GtkTextIter")
  count <- as.integer(count)

  w <- .RGtkCall("S_gtk_text_iter_backward_visible_cursor_positions", object, count, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterSetOffset <-
function(object, char.offset)
{
  checkPtrType(object, "GtkTextIter")
  char.offset <- as.integer(char.offset)

  w <- .RGtkCall("S_gtk_text_iter_set_offset", object, char.offset, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextIterSetLine <-
function(object, line.number)
{
  checkPtrType(object, "GtkTextIter")
  line.number <- as.integer(line.number)

  w <- .RGtkCall("S_gtk_text_iter_set_line", object, line.number, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextIterSetLineOffset <-
function(object, char.on.line)
{
  checkPtrType(object, "GtkTextIter")
  char.on.line <- as.integer(char.on.line)

  w <- .RGtkCall("S_gtk_text_iter_set_line_offset", object, char.on.line, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextIterSetLineIndex <-
function(object, byte.on.line)
{
  checkPtrType(object, "GtkTextIter")
  byte.on.line <- as.integer(byte.on.line)

  w <- .RGtkCall("S_gtk_text_iter_set_line_index", object, byte.on.line, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextIterForwardToEnd <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_forward_to_end", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextIterForwardToLineEnd <-
function(object)
{
  checkPtrType(object, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_forward_to_line_end", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterSetVisibleLineOffset <-
function(object, char.on.line)
{
  checkPtrType(object, "GtkTextIter")
  char.on.line <- as.integer(char.on.line)

  w <- .RGtkCall("S_gtk_text_iter_set_visible_line_offset", object, char.on.line, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextIterSetVisibleLineIndex <-
function(object, byte.on.line)
{
  checkPtrType(object, "GtkTextIter")
  byte.on.line <- as.integer(byte.on.line)

  w <- .RGtkCall("S_gtk_text_iter_set_visible_line_index", object, byte.on.line, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextIterForwardToTagToggle <-
function(object, tag = NULL)
{
  checkPtrType(object, "GtkTextIter")
  if (!is.null( tag )) checkPtrType(tag, "GtkTextTag")

  w <- .RGtkCall("S_gtk_text_iter_forward_to_tag_toggle", object, tag, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterBackwardToTagToggle <-
function(object, tag = NULL)
{
  checkPtrType(object, "GtkTextIter")
  if (!is.null( tag )) checkPtrType(tag, "GtkTextTag")

  w <- .RGtkCall("S_gtk_text_iter_backward_to_tag_toggle", object, tag, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterForwardFindChar <-
function(object, pred, user.data = NULL, limit)
{
  checkPtrType(object, "GtkTextIter")
  pred <- as.function(pred)
  
  checkPtrType(limit, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_forward_find_char", object, pred, user.data, limit, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterBackwardFindChar <-
function(object, pred, user.data = NULL, limit)
{
  checkPtrType(object, "GtkTextIter")
  pred <- as.function(pred)
  
  checkPtrType(limit, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_backward_find_char", object, pred, user.data, limit, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterForwardSearch <-
function(object, str, flags, limit = NULL)
{
  checkPtrType(object, "GtkTextIter")
  str <- as.character(str)
  
  if (!is.null( limit )) checkPtrType(limit, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_forward_search", object, str, flags, limit, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterBackwardSearch <-
function(object, str, flags, limit = NULL)
{
  checkPtrType(object, "GtkTextIter")
  str <- as.character(str)
  
  if (!is.null( limit )) checkPtrType(limit, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_backward_search", object, str, flags, limit, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterEqual <-
function(object, rhs)
{
  checkPtrType(object, "GtkTextIter")
  checkPtrType(rhs, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_equal", object, rhs, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterCompare <-
function(object, rhs)
{
  checkPtrType(object, "GtkTextIter")
  checkPtrType(rhs, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_compare", object, rhs, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterInRange <-
function(object, start, end)
{
  checkPtrType(object, "GtkTextIter")
  checkPtrType(start, "GtkTextIter")
  checkPtrType(end, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_in_range", object, start, end, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextIterOrder <-
function(object, second)
{
  checkPtrType(object, "GtkTextIter")
  checkPtrType(second, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_iter_order", object, second, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextMarkGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_text_mark_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkTextMarkSetVisible <-
function(object, setting)
{
  checkPtrType(object, "GtkTextMark")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_text_mark_set_visible", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextMarkGetVisible <-
function(object)
{
  checkPtrType(object, "GtkTextMark")

  w <- .RGtkCall("S_gtk_text_mark_get_visible", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextMarkGetName <-
function(object)
{
  checkPtrType(object, "GtkTextMark")

  w <- .RGtkCall("S_gtk_text_mark_get_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextMarkGetDeleted <-
function(object)
{
  checkPtrType(object, "GtkTextMark")

  w <- .RGtkCall("S_gtk_text_mark_get_deleted", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextMarkGetBuffer <-
function(object)
{
  checkPtrType(object, "GtkTextMark")

  w <- .RGtkCall("S_gtk_text_mark_get_buffer", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextMarkGetLeftGravity <-
function(object)
{
  checkPtrType(object, "GtkTextMark")

  w <- .RGtkCall("S_gtk_text_mark_get_left_gravity", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextTagGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_text_tag_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkTextTagNew <-
function(name = NULL)
{
  

  w <- .RGtkCall("S_gtk_text_tag_new", name, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextTagGetPriority <-
function(object)
{
  checkPtrType(object, "GtkTextTag")

  w <- .RGtkCall("S_gtk_text_tag_get_priority", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextTagSetPriority <-
function(object, priority)
{
  checkPtrType(object, "GtkTextTag")
  priority <- as.integer(priority)

  w <- .RGtkCall("S_gtk_text_tag_set_priority", object, priority, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextTagEvent <-
function(object, event.object, event, iter)
{
  checkPtrType(object, "GtkTextTag")
  checkPtrType(event.object, "GObject")
  checkPtrType(event, "GdkEvent")
  checkPtrType(iter, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_tag_event", object, event.object, event, iter, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextAttributesNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_text_attributes_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkTextAttributesCopy <-
function(object)
{
  checkPtrType(object, "GtkTextAttributes")

  w <- .RGtkCall("S_gtk_text_attributes_copy", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextAttributesCopyValues <-
function(object, dest)
{
  checkPtrType(object, "GtkTextAttributes")
  checkPtrType(dest, "GtkTextAttributes")

  w <- .RGtkCall("S_gtk_text_attributes_copy_values", object, dest, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextAttributesGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_text_attributes_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkTextTagTableGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_text_tag_table_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkTextTagTableNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_text_tag_table_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkTextTagTableAdd <-
function(object, tag)
{
  checkPtrType(object, "GtkTextTagTable")
  checkPtrType(tag, "GtkTextTag")

  w <- .RGtkCall("S_gtk_text_tag_table_add", object, tag, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextTagTableRemove <-
function(object, tag)
{
  checkPtrType(object, "GtkTextTagTable")
  checkPtrType(tag, "GtkTextTag")

  w <- .RGtkCall("S_gtk_text_tag_table_remove", object, tag, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextTagTableLookup <-
function(object, name)
{
  checkPtrType(object, "GtkTextTagTable")
  name <- as.character(name)

  w <- .RGtkCall("S_gtk_text_tag_table_lookup", object, name, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextTagTableForeach <-
function(object, func, data = NULL)
{
  checkPtrType(object, "GtkTextTagTable")
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_text_tag_table_foreach", object, func, data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextTagTableGetSize <-
function(object)
{
  checkPtrType(object, "GtkTextTagTable")

  w <- .RGtkCall("S_gtk_text_tag_table_get_size", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextViewGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_text_view_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkTextViewNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_text_view_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkTextViewNewWithBuffer <-
function(buffer = NULL, show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_text_view_new_with_buffer", buffer, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkTextViewSetBuffer <-
function(object, buffer)
{
  checkPtrType(object, "GtkTextView")
  checkPtrType(buffer, "GtkTextBuffer")

  w <- .RGtkCall("S_gtk_text_view_set_buffer", object, buffer, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextViewGetBuffer <-
function(object)
{
  checkPtrType(object, "GtkTextView")

  w <- .RGtkCall("S_gtk_text_view_get_buffer", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextViewScrollToIter <-
function(object, iter, within.margin, use.align = FALSE, xalign = 0.5, yalign = 0.5)
{
  checkPtrType(object, "GtkTextView")
  checkPtrType(iter, "GtkTextIter")
  within.margin <- as.numeric(within.margin)
  use.align <- as.logical(use.align)
  xalign <- as.numeric(xalign)
  yalign <- as.numeric(yalign)

  w <- .RGtkCall("S_gtk_text_view_scroll_to_iter", object, iter, within.margin, use.align, xalign, yalign, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextViewScrollToMark <-
function(object, mark, within.margin, use.align = FALSE, xalign = 0.5, yalign = 0.5)
{
  checkPtrType(object, "GtkTextView")
  checkPtrType(mark, "GtkTextMark")
  within.margin <- as.numeric(within.margin)
  use.align <- as.logical(use.align)
  xalign <- as.numeric(xalign)
  yalign <- as.numeric(yalign)

  w <- .RGtkCall("S_gtk_text_view_scroll_to_mark", object, mark, within.margin, use.align, xalign, yalign, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextViewScrollMarkOnscreen <-
function(object, mark)
{
  checkPtrType(object, "GtkTextView")
  checkPtrType(mark, "GtkTextMark")

  w <- .RGtkCall("S_gtk_text_view_scroll_mark_onscreen", object, mark, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextViewMoveMarkOnscreen <-
function(object, mark)
{
  checkPtrType(object, "GtkTextView")
  checkPtrType(mark, "GtkTextMark")

  w <- .RGtkCall("S_gtk_text_view_move_mark_onscreen", object, mark, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextViewPlaceCursorOnscreen <-
function(object)
{
  checkPtrType(object, "GtkTextView")

  w <- .RGtkCall("S_gtk_text_view_place_cursor_onscreen", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextViewGetVisibleRect <-
function(object)
{
  checkPtrType(object, "GtkTextView")

  w <- .RGtkCall("S_gtk_text_view_get_visible_rect", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextViewSetCursorVisible <-
function(object, setting)
{
  checkPtrType(object, "GtkTextView")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_text_view_set_cursor_visible", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextViewGetCursorVisible <-
function(object)
{
  checkPtrType(object, "GtkTextView")

  w <- .RGtkCall("S_gtk_text_view_get_cursor_visible", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextViewGetIterLocation <-
function(object, iter)
{
  checkPtrType(object, "GtkTextView")
  checkPtrType(iter, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_view_get_iter_location", object, iter, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextViewGetIterAtLocation <-
function(object, x, y)
{
  checkPtrType(object, "GtkTextView")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_gtk_text_view_get_iter_at_location", object, x, y, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextViewGetIterAtPosition <-
function(object, x, y)
{
  checkPtrType(object, "GtkTextView")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_gtk_text_view_get_iter_at_position", object, x, y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextViewGetLineYrange <-
function(object, iter)
{
  checkPtrType(object, "GtkTextView")
  checkPtrType(iter, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_view_get_line_yrange", object, iter, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextViewGetLineAtY <-
function(object, y)
{
  checkPtrType(object, "GtkTextView")
  y <- as.integer(y)

  w <- .RGtkCall("S_gtk_text_view_get_line_at_y", object, y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextViewBufferToWindowCoords <-
function(object, win, buffer.x, buffer.y)
{
  checkPtrType(object, "GtkTextView")
  
  buffer.x <- as.integer(buffer.x)
  buffer.y <- as.integer(buffer.y)

  w <- .RGtkCall("S_gtk_text_view_buffer_to_window_coords", object, win, buffer.x, buffer.y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextViewWindowToBufferCoords <-
function(object, win, window.x, window.y)
{
  checkPtrType(object, "GtkTextView")
  
  window.x <- as.integer(window.x)
  window.y <- as.integer(window.y)

  w <- .RGtkCall("S_gtk_text_view_window_to_buffer_coords", object, win, window.x, window.y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextViewGetWindow <-
function(object, win)
{
  checkPtrType(object, "GtkTextView")
  

  w <- .RGtkCall("S_gtk_text_view_get_window", object, win, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextViewGetWindowType <-
function(object, window)
{
  checkPtrType(object, "GtkTextView")
  checkPtrType(window, "GdkWindow")

  w <- .RGtkCall("S_gtk_text_view_get_window_type", object, window, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextViewSetBorderWindowSize <-
function(object, type, size)
{
  checkPtrType(object, "GtkTextView")
  
  size <- as.integer(size)

  w <- .RGtkCall("S_gtk_text_view_set_border_window_size", object, type, size, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextViewGetBorderWindowSize <-
function(object, type)
{
  checkPtrType(object, "GtkTextView")
  

  w <- .RGtkCall("S_gtk_text_view_get_border_window_size", object, type, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextViewForwardDisplayLine <-
function(object, iter)
{
  checkPtrType(object, "GtkTextView")
  checkPtrType(iter, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_view_forward_display_line", object, iter, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextViewBackwardDisplayLine <-
function(object, iter)
{
  checkPtrType(object, "GtkTextView")
  checkPtrType(iter, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_view_backward_display_line", object, iter, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextViewForwardDisplayLineEnd <-
function(object, iter)
{
  checkPtrType(object, "GtkTextView")
  checkPtrType(iter, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_view_forward_display_line_end", object, iter, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextViewBackwardDisplayLineStart <-
function(object, iter)
{
  checkPtrType(object, "GtkTextView")
  checkPtrType(iter, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_view_backward_display_line_start", object, iter, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextViewStartsDisplayLine <-
function(object, iter)
{
  checkPtrType(object, "GtkTextView")
  checkPtrType(iter, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_view_starts_display_line", object, iter, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextViewMoveVisually <-
function(object, iter, count)
{
  checkPtrType(object, "GtkTextView")
  checkPtrType(iter, "GtkTextIter")
  count <- as.integer(count)

  w <- .RGtkCall("S_gtk_text_view_move_visually", object, iter, count, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextViewAddChildAtAnchor <-
function(object, child, anchor)
{
  checkPtrType(object, "GtkTextView")
  checkPtrType(child, "GtkWidget")
  checkPtrType(anchor, "GtkTextChildAnchor")

  w <- .RGtkCall("S_gtk_text_view_add_child_at_anchor", object, child, anchor, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextViewAddChildInWindow <-
function(object, child, which.window, xpos, ypos)
{
  checkPtrType(object, "GtkTextView")
  checkPtrType(child, "GtkWidget")
  
  xpos <- as.integer(xpos)
  ypos <- as.integer(ypos)

  w <- .RGtkCall("S_gtk_text_view_add_child_in_window", object, child, which.window, xpos, ypos, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextViewMoveChild <-
function(object, child, xpos, ypos)
{
  checkPtrType(object, "GtkTextView")
  checkPtrType(child, "GtkWidget")
  xpos <- as.integer(xpos)
  ypos <- as.integer(ypos)

  w <- .RGtkCall("S_gtk_text_view_move_child", object, child, xpos, ypos, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextViewSetWrapMode <-
function(object, wrap.mode)
{
  checkPtrType(object, "GtkTextView")
  

  w <- .RGtkCall("S_gtk_text_view_set_wrap_mode", object, wrap.mode, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextViewGetWrapMode <-
function(object)
{
  checkPtrType(object, "GtkTextView")

  w <- .RGtkCall("S_gtk_text_view_get_wrap_mode", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextViewSetEditable <-
function(object, setting)
{
  checkPtrType(object, "GtkTextView")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_text_view_set_editable", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextViewGetEditable <-
function(object)
{
  checkPtrType(object, "GtkTextView")

  w <- .RGtkCall("S_gtk_text_view_get_editable", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextViewSetOverwrite <-
function(object, overwrite)
{
  checkPtrType(object, "GtkTextView")
  overwrite <- as.logical(overwrite)

  w <- .RGtkCall("S_gtk_text_view_set_overwrite", object, overwrite, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextViewGetOverwrite <-
function(object)
{
  checkPtrType(object, "GtkTextView")

  w <- .RGtkCall("S_gtk_text_view_get_overwrite", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextViewSetAcceptsTab <-
function(object, accepts.tab)
{
  checkPtrType(object, "GtkTextView")
  accepts.tab <- as.logical(accepts.tab)

  w <- .RGtkCall("S_gtk_text_view_set_accepts_tab", object, accepts.tab, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextViewGetAcceptsTab <-
function(object)
{
  checkPtrType(object, "GtkTextView")

  w <- .RGtkCall("S_gtk_text_view_get_accepts_tab", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextViewSetPixelsAboveLines <-
function(object, pixels.above.lines)
{
  checkPtrType(object, "GtkTextView")
  pixels.above.lines <- as.integer(pixels.above.lines)

  w <- .RGtkCall("S_gtk_text_view_set_pixels_above_lines", object, pixels.above.lines, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextViewGetPixelsAboveLines <-
function(object)
{
  checkPtrType(object, "GtkTextView")

  w <- .RGtkCall("S_gtk_text_view_get_pixels_above_lines", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextViewSetPixelsBelowLines <-
function(object, pixels.below.lines)
{
  checkPtrType(object, "GtkTextView")
  pixels.below.lines <- as.integer(pixels.below.lines)

  w <- .RGtkCall("S_gtk_text_view_set_pixels_below_lines", object, pixels.below.lines, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextViewGetPixelsBelowLines <-
function(object)
{
  checkPtrType(object, "GtkTextView")

  w <- .RGtkCall("S_gtk_text_view_get_pixels_below_lines", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextViewSetPixelsInsideWrap <-
function(object, pixels.inside.wrap)
{
  checkPtrType(object, "GtkTextView")
  pixels.inside.wrap <- as.integer(pixels.inside.wrap)

  w <- .RGtkCall("S_gtk_text_view_set_pixels_inside_wrap", object, pixels.inside.wrap, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextViewGetPixelsInsideWrap <-
function(object)
{
  checkPtrType(object, "GtkTextView")

  w <- .RGtkCall("S_gtk_text_view_get_pixels_inside_wrap", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextViewSetJustification <-
function(object, justification)
{
  checkPtrType(object, "GtkTextView")
  

  w <- .RGtkCall("S_gtk_text_view_set_justification", object, justification, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextViewGetJustification <-
function(object)
{
  checkPtrType(object, "GtkTextView")

  w <- .RGtkCall("S_gtk_text_view_get_justification", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextViewSetLeftMargin <-
function(object, left.margin)
{
  checkPtrType(object, "GtkTextView")
  left.margin <- as.integer(left.margin)

  w <- .RGtkCall("S_gtk_text_view_set_left_margin", object, left.margin, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextViewGetLeftMargin <-
function(object)
{
  checkPtrType(object, "GtkTextView")

  w <- .RGtkCall("S_gtk_text_view_get_left_margin", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextViewSetRightMargin <-
function(object, right.margin)
{
  checkPtrType(object, "GtkTextView")
  right.margin <- as.integer(right.margin)

  w <- .RGtkCall("S_gtk_text_view_set_right_margin", object, right.margin, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextViewGetRightMargin <-
function(object)
{
  checkPtrType(object, "GtkTextView")

  w <- .RGtkCall("S_gtk_text_view_get_right_margin", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextViewSetIndent <-
function(object, indent)
{
  checkPtrType(object, "GtkTextView")
  indent <- as.integer(indent)

  w <- .RGtkCall("S_gtk_text_view_set_indent", object, indent, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextViewGetIndent <-
function(object)
{
  checkPtrType(object, "GtkTextView")

  w <- .RGtkCall("S_gtk_text_view_get_indent", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextViewSetTabs <-
function(object, tabs)
{
  checkPtrType(object, "GtkTextView")
  checkPtrType(tabs, "PangoTabArray")

  w <- .RGtkCall("S_gtk_text_view_set_tabs", object, tabs, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextViewGetTabs <-
function(object)
{
  checkPtrType(object, "GtkTextView")

  w <- .RGtkCall("S_gtk_text_view_get_tabs", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextViewGetDefaultAttributes <-
function(object)
{
  checkPtrType(object, "GtkTextView")

  w <- .RGtkCall("S_gtk_text_view_get_default_attributes", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTipsQueryGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_tips_query_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkTipsQueryNew <-
function(show = TRUE)
{
  if(getOption("depwarn"))
    .Deprecated("GtkTooltips", "RGtk2")

  

  w <- .RGtkCall("S_gtk_tips_query_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkTipsQueryStartQuery <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("GtkTooltips", "RGtk2")

  checkPtrType(object, "GtkTipsQuery")

  w <- .RGtkCall("S_gtk_tips_query_start_query", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTipsQueryStopQuery <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("GtkTooltips", "RGtk2")

  checkPtrType(object, "GtkTipsQuery")

  w <- .RGtkCall("S_gtk_tips_query_stop_query", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTipsQuerySetCaller <-
function(object, caller)
{
  if(getOption("depwarn"))
    .Deprecated("GtkTooltips", "RGtk2")

  checkPtrType(object, "GtkTipsQuery")
  checkPtrType(caller, "GtkWidget")

  w <- .RGtkCall("S_gtk_tips_query_set_caller", object, caller, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTipsQuerySetLabels <-
function(object, label.inactive, label.no.tip)
{
  if(getOption("depwarn"))
    .Deprecated("GtkTooltips", "RGtk2")

  checkPtrType(object, "GtkTipsQuery")
  label.inactive <- as.character(label.inactive)
  label.no.tip <- as.character(label.no.tip)

  w <- .RGtkCall("S_gtk_tips_query_set_labels", object, label.inactive, label.no.tip, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToggleActionGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_toggle_action_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkToggleActionNew <-
function(name = NULL, label = NULL, tooltip = NULL, stock.id = NULL)
{
  

  w <- .RGtkCall("S_gtk_toggle_action_new", name, label, tooltip, stock.id, PACKAGE = "RGtk2")

  return(w)
} 


gtkToggleActionToggled <-
function(object)
{
  checkPtrType(object, "GtkToggleAction")

  w <- .RGtkCall("S_gtk_toggle_action_toggled", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToggleActionSetActive <-
function(object, is.active)
{
  checkPtrType(object, "GtkToggleAction")
  is.active <- as.logical(is.active)

  w <- .RGtkCall("S_gtk_toggle_action_set_active", object, is.active, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToggleActionGetActive <-
function(object)
{
  checkPtrType(object, "GtkToggleAction")

  w <- .RGtkCall("S_gtk_toggle_action_get_active", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToggleActionSetDrawAsRadio <-
function(object, draw.as.radio)
{
  checkPtrType(object, "GtkToggleAction")
  draw.as.radio <- as.logical(draw.as.radio)

  w <- .RGtkCall("S_gtk_toggle_action_set_draw_as_radio", object, draw.as.radio, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToggleActionGetDrawAsRadio <-
function(object)
{
  checkPtrType(object, "GtkToggleAction")

  w <- .RGtkCall("S_gtk_toggle_action_get_draw_as_radio", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToggleButtonGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_toggle_button_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkToggleButtonNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_toggle_button_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkToggleButtonNewWithLabel <-
function(label, show = TRUE)
{
  label <- as.character(label)

  w <- .RGtkCall("S_gtk_toggle_button_new_with_label", label, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkToggleButtonNewWithMnemonic <-
function(label, show = TRUE)
{
  label <- as.character(label)

  w <- .RGtkCall("S_gtk_toggle_button_new_with_mnemonic", label, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkToggleButtonSetMode <-
function(object, draw.indicator)
{
  checkPtrType(object, "GtkToggleButton")
  draw.indicator <- as.logical(draw.indicator)

  w <- .RGtkCall("S_gtk_toggle_button_set_mode", object, draw.indicator, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToggleButtonGetMode <-
function(object)
{
  checkPtrType(object, "GtkToggleButton")

  w <- .RGtkCall("S_gtk_toggle_button_get_mode", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToggleButtonSetActive <-
function(object, is.active)
{
  checkPtrType(object, "GtkToggleButton")
  is.active <- as.logical(is.active)

  w <- .RGtkCall("S_gtk_toggle_button_set_active", object, is.active, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToggleButtonGetActive <-
function(object)
{
  checkPtrType(object, "GtkToggleButton")

  w <- .RGtkCall("S_gtk_toggle_button_get_active", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToggleButtonToggled <-
function(object)
{
  checkPtrType(object, "GtkToggleButton")

  w <- .RGtkCall("S_gtk_toggle_button_toggled", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToggleButtonSetInconsistent <-
function(object, setting)
{
  checkPtrType(object, "GtkToggleButton")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_toggle_button_set_inconsistent", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToggleButtonGetInconsistent <-
function(object)
{
  checkPtrType(object, "GtkToggleButton")

  w <- .RGtkCall("S_gtk_toggle_button_get_inconsistent", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToggleButtonSetState <-
function(object, is.active)
{
  if(getOption("depwarn"))
    .Deprecated("gtkToggleButtonSetActive", "RGtk2")

  checkPtrType(object, "GtkToggleButton")
  is.active <- as.logical(is.active)

  w <- .RGtkCall("S_gtk_toggle_button_set_state", object, is.active, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToggleToolButtonGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_toggle_tool_button_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkToggleToolButtonNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_toggle_tool_button_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkToggleToolButtonNewFromStock <-
function(stock.id)
{
  stock.id <- as.character(stock.id)

  w <- .RGtkCall("S_gtk_toggle_tool_button_new_from_stock", stock.id, PACKAGE = "RGtk2")

  return(w)
} 


gtkToggleToolButtonSetActive <-
function(object, is.active)
{
  checkPtrType(object, "GtkToggleToolButton")
  is.active <- as.logical(is.active)

  w <- .RGtkCall("S_gtk_toggle_tool_button_set_active", object, is.active, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToggleToolButtonGetActive <-
function(object)
{
  checkPtrType(object, "GtkToggleToolButton")

  w <- .RGtkCall("S_gtk_toggle_tool_button_get_active", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolbarGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_toolbar_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkToolbarNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_toolbar_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkToolbarInsert <-
function(object, item, pos)
{
  checkPtrType(object, "GtkToolbar")
  checkPtrType(item, "GtkToolItem")
  pos <- as.integer(pos)

  w <- .RGtkCall("S_gtk_toolbar_insert", object, item, pos, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolbarGetItemIndex <-
function(object, item)
{
  checkPtrType(object, "GtkToolbar")
  checkPtrType(item, "GtkToolItem")

  w <- .RGtkCall("S_gtk_toolbar_get_item_index", object, item, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolbarGetNItems <-
function(object)
{
  checkPtrType(object, "GtkToolbar")

  w <- .RGtkCall("S_gtk_toolbar_get_n_items", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolbarGetNthItem <-
function(object, n)
{
  checkPtrType(object, "GtkToolbar")
  n <- as.integer(n)

  w <- .RGtkCall("S_gtk_toolbar_get_nth_item", object, n, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolbarGetDropIndex <-
function(object, x, y)
{
  checkPtrType(object, "GtkToolbar")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_gtk_toolbar_get_drop_index", object, x, y, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolbarSetDropHighlightItem <-
function(object, tool.item = NULL, index)
{
  checkPtrType(object, "GtkToolbar")
  if (!is.null( tool.item )) checkPtrType(tool.item, "GtkToolItem")
  index <- as.integer(index)

  w <- .RGtkCall("S_gtk_toolbar_set_drop_highlight_item", object, tool.item, index, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolbarSetShowArrow <-
function(object, show.arrow)
{
  checkPtrType(object, "GtkToolbar")
  show.arrow <- as.logical(show.arrow)

  w <- .RGtkCall("S_gtk_toolbar_set_show_arrow", object, show.arrow, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolbarGetShowArrow <-
function(object)
{
  checkPtrType(object, "GtkToolbar")

  w <- .RGtkCall("S_gtk_toolbar_get_show_arrow", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolbarGetReliefStyle <-
function(object)
{
  checkPtrType(object, "GtkToolbar")

  w <- .RGtkCall("S_gtk_toolbar_get_relief_style", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolbarAppendItem <-
function(object, text, tooltip.text, tooltip.private.text, icon, callback, user.data = NULL)
{
  if(getOption("depwarn"))
    .Deprecated("gtkToolbarInsert", "RGtk2")

  checkPtrType(object, "GtkToolbar")
  text <- as.character(text)
  tooltip.text <- as.character(tooltip.text)
  tooltip.private.text <- as.character(tooltip.private.text)
  checkPtrType(icon, "GtkWidget")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gtk_toolbar_append_item", object, text, tooltip.text, tooltip.private.text, icon, callback, user.data, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolbarPrependItem <-
function(object, text, tooltip.text, tooltip.private.text, icon, callback, user.data)
{
  if(getOption("depwarn"))
    .Deprecated("gtkToolbarInsert", "RGtk2")

  checkPtrType(object, "GtkToolbar")
  text <- as.character(text)
  tooltip.text <- as.character(tooltip.text)
  tooltip.private.text <- as.character(tooltip.private.text)
  checkPtrType(icon, "GtkWidget")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gtk_toolbar_prepend_item", object, text, tooltip.text, tooltip.private.text, icon, callback, user.data, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolbarInsertItem <-
function(object, text, tooltip.text, tooltip.private.text, icon, callback, user.data, position)
{
  if(getOption("depwarn"))
    .Deprecated("gtkToolbarInsert", "RGtk2")

  checkPtrType(object, "GtkToolbar")
  text <- as.character(text)
  tooltip.text <- as.character(tooltip.text)
  tooltip.private.text <- as.character(tooltip.private.text)
  checkPtrType(icon, "GtkWidget")
  callback <- as.function(callback)
  
  position <- as.integer(position)

  w <- .RGtkCall("S_gtk_toolbar_insert_item", object, text, tooltip.text, tooltip.private.text, icon, callback, user.data, position, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolbarInsertStock <-
function(object, stock.id, tooltip.text, tooltip.private.text, callback, user.data, position)
{
  if(getOption("depwarn"))
    .Deprecated("gtkToolbarInsert", "RGtk2")

  checkPtrType(object, "GtkToolbar")
  stock.id <- as.character(stock.id)
  tooltip.text <- as.character(tooltip.text)
  tooltip.private.text <- as.character(tooltip.private.text)
  callback <- as.function(callback)
  
  position <- as.integer(position)

  w <- .RGtkCall("S_gtk_toolbar_insert_stock", object, stock.id, tooltip.text, tooltip.private.text, callback, user.data, position, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolbarAppendSpace <-
function(object)
{
  checkPtrType(object, "GtkToolbar")

  w <- .RGtkCall("S_gtk_toolbar_append_space", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolbarPrependSpace <-
function(object)
{
  checkPtrType(object, "GtkToolbar")

  w <- .RGtkCall("S_gtk_toolbar_prepend_space", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolbarInsertSpace <-
function(object, position)
{
  checkPtrType(object, "GtkToolbar")
  position <- as.integer(position)

  w <- .RGtkCall("S_gtk_toolbar_insert_space", object, position, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolbarRemoveSpace <-
function(object, position)
{
  checkPtrType(object, "GtkToolbar")
  position <- as.integer(position)

  w <- .RGtkCall("S_gtk_toolbar_remove_space", object, position, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolbarAppendElement <-
function(object, type, widget, text, tooltip.text, tooltip.private.text, icon, callback, user.data = NULL)
{
  if(getOption("depwarn"))
    .Deprecated("gtkToolbarInsert", "RGtk2")

  checkPtrType(object, "GtkToolbar")
  
  checkPtrType(widget, "GtkWidget")
  text <- as.character(text)
  tooltip.text <- as.character(tooltip.text)
  tooltip.private.text <- as.character(tooltip.private.text)
  checkPtrType(icon, "GtkWidget")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gtk_toolbar_append_element", object, type, widget, text, tooltip.text, tooltip.private.text, icon, callback, user.data, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolbarPrependElement <-
function(object, type, widget, text, tooltip.text, tooltip.private.text, icon, callback, user.data = NULL)
{
  if(getOption("depwarn"))
    .Deprecated("gtkToolbarInsert", "RGtk2")

  checkPtrType(object, "GtkToolbar")
  
  checkPtrType(widget, "GtkWidget")
  text <- as.character(text)
  tooltip.text <- as.character(tooltip.text)
  tooltip.private.text <- as.character(tooltip.private.text)
  checkPtrType(icon, "GtkWidget")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gtk_toolbar_prepend_element", object, type, widget, text, tooltip.text, tooltip.private.text, icon, callback, user.data, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolbarInsertElement <-
function(object, type, widget, text, tooltip.text, tooltip.private.text, icon, callback, user.data = NULL, position)
{
  if(getOption("depwarn"))
    .Deprecated("gtkToolbarInsert", "RGtk2")

  checkPtrType(object, "GtkToolbar")
  
  checkPtrType(widget, "GtkWidget")
  text <- as.character(text)
  tooltip.text <- as.character(tooltip.text)
  tooltip.private.text <- as.character(tooltip.private.text)
  checkPtrType(icon, "GtkWidget")
  callback <- as.function(callback)
  
  position <- as.integer(position)

  w <- .RGtkCall("S_gtk_toolbar_insert_element", object, type, widget, text, tooltip.text, tooltip.private.text, icon, callback, user.data, position, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolbarAppendWidget <-
function(object, widget, tooltip.text = NULL, tooltip.private.text = NULL)
{
  checkPtrType(object, "GtkToolbar")
  checkPtrType(widget, "GtkWidget")
  if (!is.null( tooltip.text )) tooltip.text <- as.character(tooltip.text)
  if (!is.null( tooltip.private.text )) tooltip.private.text <- as.character(tooltip.private.text)

  w <- .RGtkCall("S_gtk_toolbar_append_widget", object, widget, tooltip.text, tooltip.private.text, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolbarPrependWidget <-
function(object, widget, tooltip.text = NULL, tooltip.private.text = NULL)
{
  checkPtrType(object, "GtkToolbar")
  checkPtrType(widget, "GtkWidget")
  if (!is.null( tooltip.text )) tooltip.text <- as.character(tooltip.text)
  if (!is.null( tooltip.private.text )) tooltip.private.text <- as.character(tooltip.private.text)

  w <- .RGtkCall("S_gtk_toolbar_prepend_widget", object, widget, tooltip.text, tooltip.private.text, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolbarInsertWidget <-
function(object, widget, tooltip.text = NULL, tooltip.private.text = NULL, position)
{
  checkPtrType(object, "GtkToolbar")
  checkPtrType(widget, "GtkWidget")
  if (!is.null( tooltip.text )) tooltip.text <- as.character(tooltip.text)
  if (!is.null( tooltip.private.text )) tooltip.private.text <- as.character(tooltip.private.text)
  position <- as.integer(position)

  w <- .RGtkCall("S_gtk_toolbar_insert_widget", object, widget, tooltip.text, tooltip.private.text, position, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolbarSetOrientation <-
function(object, orientation)
{
  checkPtrType(object, "GtkToolbar")
  

  w <- .RGtkCall("S_gtk_toolbar_set_orientation", object, orientation, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolbarSetStyle <-
function(object, style)
{
  checkPtrType(object, "GtkToolbar")
  

  w <- .RGtkCall("S_gtk_toolbar_set_style", object, style, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolbarSetIconSize <-
function(object, icon.size)
{
  checkPtrType(object, "GtkToolbar")
  

  w <- .RGtkCall("S_gtk_toolbar_set_icon_size", object, icon.size, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolbarSetTooltips <-
function(object, enable)
{
  if(getOption("depwarn"))
    .Deprecated("'gtk-enable-tooltips' property on GtkSettings", "RGtk2")

  checkPtrType(object, "GtkToolbar")
  enable <- as.logical(enable)

  w <- .RGtkCall("S_gtk_toolbar_set_tooltips", object, enable, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolbarUnsetStyle <-
function(object)
{
  checkPtrType(object, "GtkToolbar")

  w <- .RGtkCall("S_gtk_toolbar_unset_style", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolbarUnsetIconSize <-
function(object)
{
  checkPtrType(object, "GtkToolbar")

  w <- .RGtkCall("S_gtk_toolbar_unset_icon_size", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolbarGetOrientation <-
function(object)
{
  checkPtrType(object, "GtkToolbar")

  w <- .RGtkCall("S_gtk_toolbar_get_orientation", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolbarGetStyle <-
function(object)
{
  checkPtrType(object, "GtkToolbar")

  w <- .RGtkCall("S_gtk_toolbar_get_style", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolbarGetIconSize <-
function(object)
{
  checkPtrType(object, "GtkToolbar")

  w <- .RGtkCall("S_gtk_toolbar_get_icon_size", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolbarGetTooltips <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("'gtk-enable-tooltips' property on GtkSettings", "RGtk2")

  checkPtrType(object, "GtkToolbar")

  w <- .RGtkCall("S_gtk_toolbar_get_tooltips", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolButtonGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_tool_button_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkToolButtonNew <-
function(icon.widget = NULL, label = NULL, show = TRUE)
{
  if (!is.null( icon.widget )) checkPtrType(icon.widget, "GtkWidget")
  if (!is.null( label )) label <- as.character(label)

  w <- .RGtkCall("S_gtk_tool_button_new", icon.widget, label, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkToolButtonNewFromStock <-
function(stock.id, show = TRUE)
{
  stock.id <- as.character(stock.id)

  w <- .RGtkCall("S_gtk_tool_button_new_from_stock", stock.id, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkToolButtonSetLabel <-
function(object, label = NULL)
{
  checkPtrType(object, "GtkToolButton")
  if (!is.null( label )) label <- as.character(label)

  w <- .RGtkCall("S_gtk_tool_button_set_label", object, label, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolButtonGetLabel <-
function(object)
{
  checkPtrType(object, "GtkToolButton")

  w <- .RGtkCall("S_gtk_tool_button_get_label", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolButtonSetUseUnderline <-
function(object, use.underline)
{
  checkPtrType(object, "GtkToolButton")
  use.underline <- as.logical(use.underline)

  w <- .RGtkCall("S_gtk_tool_button_set_use_underline", object, use.underline, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolButtonGetUseUnderline <-
function(object)
{
  checkPtrType(object, "GtkToolButton")

  w <- .RGtkCall("S_gtk_tool_button_get_use_underline", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolButtonSetStockId <-
function(object, stock.id = NULL)
{
  checkPtrType(object, "GtkToolButton")
  if (!is.null( stock.id )) stock.id <- as.character(stock.id)

  w <- .RGtkCall("S_gtk_tool_button_set_stock_id", object, stock.id, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolButtonSetIconName <-
function(object, icon.name)
{
  checkPtrType(object, "GtkToolButton")
  icon.name <- as.character(icon.name)

  w <- .RGtkCall("S_gtk_tool_button_set_icon_name", object, icon.name, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolButtonGetIconName <-
function(object)
{
  checkPtrType(object, "GtkToolButton")

  w <- .RGtkCall("S_gtk_tool_button_get_icon_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolButtonGetStockId <-
function(object)
{
  checkPtrType(object, "GtkToolButton")

  w <- .RGtkCall("S_gtk_tool_button_get_stock_id", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolButtonSetIconWidget <-
function(object, icon.widget = NULL)
{
  checkPtrType(object, "GtkToolButton")
  if (!is.null( icon.widget )) checkPtrType(icon.widget, "GtkWidget")

  w <- .RGtkCall("S_gtk_tool_button_set_icon_widget", object, icon.widget, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolButtonGetIconWidget <-
function(object)
{
  checkPtrType(object, "GtkToolButton")

  w <- .RGtkCall("S_gtk_tool_button_get_icon_widget", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolButtonSetLabelWidget <-
function(object, label.widget = NULL)
{
  checkPtrType(object, "GtkToolButton")
  if (!is.null( label.widget )) checkPtrType(label.widget, "GtkWidget")

  w <- .RGtkCall("S_gtk_tool_button_set_label_widget", object, label.widget, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolButtonGetLabelWidget <-
function(object)
{
  checkPtrType(object, "GtkToolButton")

  w <- .RGtkCall("S_gtk_tool_button_get_label_widget", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolItemGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_tool_item_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkToolItemNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_tool_item_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkToolItemSetHomogeneous <-
function(object, homogeneous)
{
  checkPtrType(object, "GtkToolItem")
  homogeneous <- as.logical(homogeneous)

  w <- .RGtkCall("S_gtk_tool_item_set_homogeneous", object, homogeneous, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolItemGetHomogeneous <-
function(object)
{
  checkPtrType(object, "GtkToolItem")

  w <- .RGtkCall("S_gtk_tool_item_get_homogeneous", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolItemSetExpand <-
function(object, expand)
{
  checkPtrType(object, "GtkToolItem")
  expand <- as.logical(expand)

  w <- .RGtkCall("S_gtk_tool_item_set_expand", object, expand, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolItemGetExpand <-
function(object)
{
  checkPtrType(object, "GtkToolItem")

  w <- .RGtkCall("S_gtk_tool_item_get_expand", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolItemSetTooltip <-
function(object, tooltips, tip.text = NULL, tip.private = NULL)
{
  checkPtrType(object, "GtkToolItem")
  checkPtrType(tooltips, "GtkTooltips")
  if (!is.null( tip.text )) tip.text <- as.character(tip.text)
  if (!is.null( tip.private )) tip.private <- as.character(tip.private)

  w <- .RGtkCall("S_gtk_tool_item_set_tooltip", object, tooltips, tip.text, tip.private, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolItemSetUseDragWindow <-
function(object, use.drag.window)
{
  checkPtrType(object, "GtkToolItem")
  use.drag.window <- as.logical(use.drag.window)

  w <- .RGtkCall("S_gtk_tool_item_set_use_drag_window", object, use.drag.window, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolItemGetUseDragWindow <-
function(object)
{
  checkPtrType(object, "GtkToolItem")

  w <- .RGtkCall("S_gtk_tool_item_get_use_drag_window", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolItemSetVisibleHorizontal <-
function(object, visible.horizontal)
{
  checkPtrType(object, "GtkToolItem")
  visible.horizontal <- as.logical(visible.horizontal)

  w <- .RGtkCall("S_gtk_tool_item_set_visible_horizontal", object, visible.horizontal, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolItemGetVisibleHorizontal <-
function(object)
{
  checkPtrType(object, "GtkToolItem")

  w <- .RGtkCall("S_gtk_tool_item_get_visible_horizontal", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolItemSetVisibleVertical <-
function(object, visible.vertical)
{
  checkPtrType(object, "GtkToolItem")
  visible.vertical <- as.logical(visible.vertical)

  w <- .RGtkCall("S_gtk_tool_item_set_visible_vertical", object, visible.vertical, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolItemGetVisibleVertical <-
function(object)
{
  checkPtrType(object, "GtkToolItem")

  w <- .RGtkCall("S_gtk_tool_item_get_visible_vertical", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolItemSetIsImportant <-
function(object, is.important)
{
  checkPtrType(object, "GtkToolItem")
  is.important <- as.logical(is.important)

  w <- .RGtkCall("S_gtk_tool_item_set_is_important", object, is.important, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolItemGetIsImportant <-
function(object)
{
  checkPtrType(object, "GtkToolItem")

  w <- .RGtkCall("S_gtk_tool_item_get_is_important", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolItemGetIconSize <-
function(object)
{
  checkPtrType(object, "GtkToolItem")

  w <- .RGtkCall("S_gtk_tool_item_get_icon_size", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolItemGetOrientation <-
function(object)
{
  checkPtrType(object, "GtkToolItem")

  w <- .RGtkCall("S_gtk_tool_item_get_orientation", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolItemGetToolbarStyle <-
function(object)
{
  checkPtrType(object, "GtkToolItem")

  w <- .RGtkCall("S_gtk_tool_item_get_toolbar_style", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolItemGetReliefStyle <-
function(object)
{
  checkPtrType(object, "GtkToolItem")

  w <- .RGtkCall("S_gtk_tool_item_get_relief_style", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolItemRetrieveProxyMenuItem <-
function(object)
{
  checkPtrType(object, "GtkToolItem")

  w <- .RGtkCall("S_gtk_tool_item_retrieve_proxy_menu_item", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolItemSetProxyMenuItem <-
function(object, menu.item.id, menu.item = NULL)
{
  checkPtrType(object, "GtkToolItem")
  menu.item.id <- as.character(menu.item.id)
  if (!is.null( menu.item )) checkPtrType(menu.item, "GtkWidget")

  w <- .RGtkCall("S_gtk_tool_item_set_proxy_menu_item", object, menu.item.id, menu.item, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolItemGetProxyMenuItem <-
function(object, menu.item.id)
{
  checkPtrType(object, "GtkToolItem")
  menu.item.id <- as.character(menu.item.id)

  w <- .RGtkCall("S_gtk_tool_item_get_proxy_menu_item", object, menu.item.id, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolItemRebuildMenu <-
function(object)
{
  checkPtrType(object, "GtkToolItem")

  w <- .RGtkCall("S_gtk_tool_item_rebuild_menu", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTooltipsGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_tooltips_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkTooltipsNew <-
function()
{
  if(getOption("depwarn"))
    .Deprecated("GtkTooltip", "RGtk2")

  

  w <- .RGtkCall("S_gtk_tooltips_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkTooltipsEnable <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("GtkTooltip", "RGtk2")

  checkPtrType(object, "GtkTooltips")

  w <- .RGtkCall("S_gtk_tooltips_enable", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTooltipsDisable <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("GtkTooltip", "RGtk2")

  checkPtrType(object, "GtkTooltips")

  w <- .RGtkCall("S_gtk_tooltips_disable", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTooltipsSetDelay <-
function(object, delay)
{
  if(getOption("depwarn"))
    .Deprecated("GtkTooltip", "RGtk2")

  checkPtrType(object, "GtkTooltips")
  delay <- as.numeric(delay)

  w <- .RGtkCall("S_gtk_tooltips_set_delay", object, delay, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTooltipsSetTip <-
function(object, widget, tip.text = NULL, tip.private = NULL)
{
  if(getOption("depwarn"))
    .Deprecated("GtkTooltip", "RGtk2")

  checkPtrType(object, "GtkTooltips")
  checkPtrType(widget, "GtkWidget")
  if (!is.null( tip.text )) tip.text <- as.character(tip.text)
  if (!is.null( tip.private )) tip.private <- as.character(tip.private)

  w <- .RGtkCall("S_gtk_tooltips_set_tip", object, widget, tip.text, tip.private, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTooltipsDataGet <-
function(widget)
{
  if(getOption("depwarn"))
    .Deprecated("GtkTooltip", "RGtk2")

  checkPtrType(widget, "GtkWidget")

  w <- .RGtkCall("S_gtk_tooltips_data_get", widget, PACKAGE = "RGtk2")

  return(w)
} 


gtkTooltipsForceWindow <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("GtkTooltip", "RGtk2")

  checkPtrType(object, "GtkTooltips")

  w <- .RGtkCall("S_gtk_tooltips_force_window", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTooltipsGetInfoFromTipWindow <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("GtkTooltip", "RGtk2")

  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_tooltips_get_info_from_tip_window", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeDragSourceGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_tree_drag_source_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeDragSourceRowDraggable <-
function(object, path)
{
  checkPtrType(object, "GtkTreeDragSource")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_drag_source_row_draggable", object, path, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeDragSourceDragDataDelete <-
function(object, path)
{
  checkPtrType(object, "GtkTreeDragSource")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_drag_source_drag_data_delete", object, path, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeDragSourceDragDataGet <-
function(object, path)
{
  checkPtrType(object, "GtkTreeDragSource")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_drag_source_drag_data_get", object, path, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeDragDestGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_tree_drag_dest_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeDragDestDragDataReceived <-
function(object, dest, selection.data)
{
  checkPtrType(object, "GtkTreeDragDest")
  checkPtrType(dest, "GtkTreePath")
  checkPtrType(selection.data, "GtkSelectionData")

  w <- .RGtkCall("S_gtk_tree_drag_dest_drag_data_received", object, dest, selection.data, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeDragDestRowDropPossible <-
function(object, dest.path, selection.data)
{
  checkPtrType(object, "GtkTreeDragDest")
  checkPtrType(dest.path, "GtkTreePath")
  checkPtrType(selection.data, "GtkSelectionData")

  w <- .RGtkCall("S_gtk_tree_drag_dest_row_drop_possible", object, dest.path, selection.data, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeSetRowDragData <-
function(object, tree.model, path)
{
  checkPtrType(object, "GtkSelectionData")
  checkPtrType(tree.model, "GtkTreeModel")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_set_row_drag_data", object, tree.model, path, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeGetRowDragData <-
function(object)
{
  checkPtrType(object, "GtkSelectionData")

  w <- .RGtkCall("S_gtk_tree_get_row_drag_data", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreePathNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_tree_path_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkTreePathNewFromString <-
function(path)
{
  path <- as.character(path)

  w <- .RGtkCall("S_gtk_tree_path_new_from_string", path, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreePathToString <-
function(object)
{
  checkPtrType(object, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_path_to_string", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreePathNewFirst <-
function()
{
  

  w <- .RGtkCall("S_gtk_tree_path_new_first", PACKAGE = "RGtk2")

  return(w)
} 


gtkTreePathAppendIndex <-
function(object, index)
{
  checkPtrType(object, "GtkTreePath")
  index <- as.integer(index)

  w <- .RGtkCall("S_gtk_tree_path_append_index", object, index, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreePathPrependIndex <-
function(object, index)
{
  checkPtrType(object, "GtkTreePath")
  index <- as.integer(index)

  w <- .RGtkCall("S_gtk_tree_path_prepend_index", object, index, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreePathGetDepth <-
function(object)
{
  checkPtrType(object, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_path_get_depth", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreePathGetIndices <-
function(object)
{
  checkPtrType(object, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_path_get_indices", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreePathCopy <-
function(object)
{
  checkPtrType(object, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_path_copy", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreePathCompare <-
function(object, b)
{
  checkPtrType(object, "GtkTreePath")
  checkPtrType(b, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_path_compare", object, b, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreePathNext <-
function(object)
{
  checkPtrType(object, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_path_next", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreePathPrev <-
function(object)
{
  checkPtrType(object, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_path_prev", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreePathUp <-
function(object)
{
  checkPtrType(object, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_path_up", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreePathDown <-
function(object)
{
  checkPtrType(object, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_path_down", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreePathIsAncestor <-
function(object, descendant)
{
  checkPtrType(object, "GtkTreePath")
  checkPtrType(descendant, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_path_is_ancestor", object, descendant, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreePathIsDescendant <-
function(object, ancestor)
{
  checkPtrType(object, "GtkTreePath")
  checkPtrType(ancestor, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_path_is_descendant", object, ancestor, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeRowReferenceGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_tree_row_reference_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeRowReferenceNew <-
function(model, path)
{
  checkPtrType(model, "GtkTreeModel")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_row_reference_new", model, path, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeRowReferenceNewProxy <-
function(proxy, model, path)
{
  checkPtrType(proxy, "GObject")
  checkPtrType(model, "GtkTreeModel")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_row_reference_new_proxy", proxy, model, path, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeRowReferenceGetPath <-
function(object)
{
  checkPtrType(object, "GtkTreeRowReference")

  w <- .RGtkCall("S_gtk_tree_row_reference_get_path", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeRowReferenceGetModel <-
function(object)
{
  checkPtrType(object, "GtkTreeRowReference")

  w <- .RGtkCall("S_gtk_tree_row_reference_get_model", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeRowReferenceValid <-
function(object)
{
  checkPtrType(object, "GtkTreeRowReference")

  w <- .RGtkCall("S_gtk_tree_row_reference_valid", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeRowReferenceCopy <-
function(object)
{
  checkPtrType(object, "GtkTreeRowReference")

  w <- .RGtkCall("S_gtk_tree_row_reference_copy", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeRowReferenceInserted <-
function(proxy, path)
{
  checkPtrType(proxy, "GObject")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_row_reference_inserted", proxy, path, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeRowReferenceDeleted <-
function(proxy, path)
{
  checkPtrType(proxy, "GObject")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_row_reference_deleted", proxy, path, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeRowReferenceReordered <-
function(proxy, path, iter, new.order)
{
  checkPtrType(proxy, "GObject")
  checkPtrType(path, "GtkTreePath")
  checkPtrType(iter, "GtkTreeIter")
  new.order <- as.list(as.integer(new.order))

  w <- .RGtkCall("S_gtk_tree_row_reference_reordered", proxy, path, iter, new.order, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeIterCopy <-
function(object)
{
  checkPtrType(object, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_iter_copy", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeIterGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_tree_iter_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeModelGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_tree_model_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeModelGetFlags <-
function(object)
{
  checkPtrType(object, "GtkTreeModel")

  w <- .RGtkCall("S_gtk_tree_model_get_flags", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeModelGetNColumns <-
function(object)
{
  checkPtrType(object, "GtkTreeModel")

  w <- .RGtkCall("S_gtk_tree_model_get_n_columns", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeModelGetColumnType <-
function(object, index)
{
  checkPtrType(object, "GtkTreeModel")
  index <- as.integer(index)

  w <- .RGtkCall("S_gtk_tree_model_get_column_type", object, index, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeModelGetIter <-
function(object, path)
{
  checkPtrType(object, "GtkTreeModel")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_model_get_iter", object, path, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeModelGetIterFromString <-
function(object, path.string)
{
  checkPtrType(object, "GtkTreeModel")
  path.string <- as.character(path.string)

  w <- .RGtkCall("S_gtk_tree_model_get_iter_from_string", object, path.string, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeModelGetStringFromIter <-
function(object, iter)
{
  checkPtrType(object, "GtkTreeModel")
  checkPtrType(iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_model_get_string_from_iter", object, iter, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeModelGetIterRoot <-
function(object)
{
  checkPtrType(object, "GtkTreeModel")

  w <- .RGtkCall("S_gtk_tree_model_get_iter_root", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeModelGetIterFirst <-
function(object)
{
  checkPtrType(object, "GtkTreeModel")

  w <- .RGtkCall("S_gtk_tree_model_get_iter_first", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeModelGetPath <-
function(object, iter)
{
  checkPtrType(object, "GtkTreeModel")
  checkPtrType(iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_model_get_path", object, iter, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeModelGetValue <-
function(object, iter, column)
{
  checkPtrType(object, "GtkTreeModel")
  checkPtrType(iter, "GtkTreeIter")
  column <- as.integer(column)

  w <- .RGtkCall("S_gtk_tree_model_get_value", object, iter, column, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeModelIterNext <-
function(object, iter)
{
  checkPtrType(object, "GtkTreeModel")
  checkPtrType(iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_model_iter_next", object, iter, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeModelIterChildren <-
function(object, parent = NULL)
{
  checkPtrType(object, "GtkTreeModel")
  if (!is.null( parent )) checkPtrType(parent, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_model_iter_children", object, parent, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeModelIterHasChild <-
function(object, iter)
{
  checkPtrType(object, "GtkTreeModel")
  checkPtrType(iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_model_iter_has_child", object, iter, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeModelIterNChildren <-
function(object, iter = NULL)
{
  checkPtrType(object, "GtkTreeModel")
  if (!is.null( iter )) checkPtrType(iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_model_iter_n_children", object, iter, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeModelIterNthChild <-
function(object, parent = NULL, n)
{
  checkPtrType(object, "GtkTreeModel")
  if (!is.null( parent )) checkPtrType(parent, "GtkTreeIter")
  n <- as.integer(n)

  w <- .RGtkCall("S_gtk_tree_model_iter_nth_child", object, parent, n, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeModelIterParent <-
function(object, child)
{
  checkPtrType(object, "GtkTreeModel")
  checkPtrType(child, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_model_iter_parent", object, child, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeModelRefNode <-
function(object, iter)
{
  checkPtrType(object, "GtkTreeModel")
  checkPtrType(iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_model_ref_node", object, iter, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeModelUnrefNode <-
function(object, iter)
{
  checkPtrType(object, "GtkTreeModel")
  checkPtrType(iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_model_unref_node", object, iter, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeModelForeach <-
function(object, func, user.data = NULL)
{
  checkPtrType(object, "GtkTreeModel")
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_tree_model_foreach", object, func, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeModelRowChanged <-
function(object, path, iter)
{
  checkPtrType(object, "GtkTreeModel")
  checkPtrType(path, "GtkTreePath")
  checkPtrType(iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_model_row_changed", object, path, iter, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeModelRowInserted <-
function(object, path, iter)
{
  checkPtrType(object, "GtkTreeModel")
  checkPtrType(path, "GtkTreePath")
  checkPtrType(iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_model_row_inserted", object, path, iter, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeModelRowHasChildToggled <-
function(object, path, iter)
{
  checkPtrType(object, "GtkTreeModel")
  checkPtrType(path, "GtkTreePath")
  checkPtrType(iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_model_row_has_child_toggled", object, path, iter, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeModelRowDeleted <-
function(object, path)
{
  checkPtrType(object, "GtkTreeModel")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_model_row_deleted", object, path, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeModelRowsReordered <-
function(object, path, iter, new.order)
{
  checkPtrType(object, "GtkTreeModel")
  checkPtrType(path, "GtkTreePath")
  checkPtrType(iter, "GtkTreeIter")
  new.order <- as.list(as.integer(new.order))

  w <- .RGtkCall("S_gtk_tree_model_rows_reordered", object, path, iter, new.order, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeModelFilterGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_tree_model_filter_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeModelFilterNew <-
function(child.model, root = NULL)
{
  checkPtrType(child.model, "GtkTreeModel")
  if (!is.null( root )) checkPtrType(root, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_model_filter_new", child.model, root, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeModelFilterSetVisibleFunc <-
function(object, func, data = NULL)
{
  checkPtrType(object, "GtkTreeModelFilter")
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_tree_model_filter_set_visible_func", object, func, data, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeModelFilterSetModifyFunc <-
function(object, types, func, data = NULL)
{
  checkPtrType(object, "GtkTreeModelFilter")
  types <- as.list(as.GType(types))
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_tree_model_filter_set_modify_func", object, types, func, data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeModelFilterSetVisibleColumn <-
function(object, column)
{
  checkPtrType(object, "GtkTreeModelFilter")
  column <- as.integer(column)

  w <- .RGtkCall("S_gtk_tree_model_filter_set_visible_column", object, column, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeModelFilterGetModel <-
function(object)
{
  checkPtrType(object, "GtkTreeModelFilter")

  w <- .RGtkCall("S_gtk_tree_model_filter_get_model", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeModelFilterConvertChildIterToIter <-
function(object, child.iter)
{
  checkPtrType(object, "GtkTreeModelFilter")
  checkPtrType(child.iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_model_filter_convert_child_iter_to_iter", object, child.iter, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeModelFilterConvertIterToChildIter <-
function(object, filter.iter)
{
  checkPtrType(object, "GtkTreeModelFilter")
  checkPtrType(filter.iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_model_filter_convert_iter_to_child_iter", object, filter.iter, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeModelFilterConvertChildPathToPath <-
function(object, child.path)
{
  checkPtrType(object, "GtkTreeModelFilter")
  checkPtrType(child.path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_model_filter_convert_child_path_to_path", object, child.path, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeModelFilterConvertPathToChildPath <-
function(object, filter.path)
{
  checkPtrType(object, "GtkTreeModelFilter")
  checkPtrType(filter.path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_model_filter_convert_path_to_child_path", object, filter.path, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeModelFilterRefilter <-
function(object)
{
  checkPtrType(object, "GtkTreeModelFilter")

  w <- .RGtkCall("S_gtk_tree_model_filter_refilter", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeModelFilterClearCache <-
function(object)
{
  checkPtrType(object, "GtkTreeModelFilter")

  w <- .RGtkCall("S_gtk_tree_model_filter_clear_cache", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeModelSortGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_tree_model_sort_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeModelSortNewWithModel <-
function(child.model = NULL)
{
  

  w <- .RGtkCall("S_gtk_tree_model_sort_new_with_model", child.model, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeModelSortGetModel <-
function(object)
{
  checkPtrType(object, "GtkTreeModelSort")

  w <- .RGtkCall("S_gtk_tree_model_sort_get_model", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeModelSortConvertChildPathToPath <-
function(object, child.path)
{
  checkPtrType(object, "GtkTreeModelSort")
  checkPtrType(child.path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_model_sort_convert_child_path_to_path", object, child.path, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeModelSortConvertChildIterToIter <-
function(object, child.iter)
{
  checkPtrType(object, "GtkTreeModelSort")
  checkPtrType(child.iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_model_sort_convert_child_iter_to_iter", object, child.iter, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeModelSortConvertPathToChildPath <-
function(object, sorted.path)
{
  checkPtrType(object, "GtkTreeModelSort")
  checkPtrType(sorted.path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_model_sort_convert_path_to_child_path", object, sorted.path, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeModelSortConvertIterToChildIter <-
function(object, sorted.iter)
{
  checkPtrType(object, "GtkTreeModelSort")
  checkPtrType(sorted.iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_model_sort_convert_iter_to_child_iter", object, sorted.iter, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeModelSortResetDefaultSortFunc <-
function(object)
{
  checkPtrType(object, "GtkTreeModelSort")

  w <- .RGtkCall("S_gtk_tree_model_sort_reset_default_sort_func", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeModelSortClearCache <-
function(object)
{
  checkPtrType(object, "GtkTreeModelSort")

  w <- .RGtkCall("S_gtk_tree_model_sort_clear_cache", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeModelSortIterIsValid <-
function(object, iter)
{
  checkPtrType(object, "GtkTreeModelSort")
  checkPtrType(iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_model_sort_iter_is_valid", object, iter, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeSelectionGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_tree_selection_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeSelectionSetMode <-
function(object, type)
{
  checkPtrType(object, "GtkTreeSelection")
  

  w <- .RGtkCall("S_gtk_tree_selection_set_mode", object, type, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeSelectionGetMode <-
function(object)
{
  checkPtrType(object, "GtkTreeSelection")

  w <- .RGtkCall("S_gtk_tree_selection_get_mode", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeSelectionSetSelectFunction <-
function(object, func, data = NULL)
{
  checkPtrType(object, "GtkTreeSelection")
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_tree_selection_set_select_function", object, func, data, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeSelectionGetUserData <-
function(object)
{
  checkPtrType(object, "GtkTreeSelection")

  w <- .RGtkCall("S_gtk_tree_selection_get_user_data", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeSelectionGetTreeView <-
function(object)
{
  checkPtrType(object, "GtkTreeSelection")

  w <- .RGtkCall("S_gtk_tree_selection_get_tree_view", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeSelectionGetSelected <-
function(object)
{
  checkPtrType(object, "GtkTreeSelection")

  w <- .RGtkCall("S_gtk_tree_selection_get_selected", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeSelectionGetSelectedRows <-
function(object)
{
  checkPtrType(object, "GtkTreeSelection")

  w <- .RGtkCall("S_gtk_tree_selection_get_selected_rows", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeSelectionCountSelectedRows <-
function(object)
{
  checkPtrType(object, "GtkTreeSelection")

  w <- .RGtkCall("S_gtk_tree_selection_count_selected_rows", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeSelectionSelectedForeach <-
function(object, func, data = NULL)
{
  checkPtrType(object, "GtkTreeSelection")
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_tree_selection_selected_foreach", object, func, data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeSelectionSelectPath <-
function(object, path)
{
  checkPtrType(object, "GtkTreeSelection")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_selection_select_path", object, path, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeSelectionUnselectPath <-
function(object, path)
{
  checkPtrType(object, "GtkTreeSelection")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_selection_unselect_path", object, path, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeSelectionSelectIter <-
function(object, iter)
{
  checkPtrType(object, "GtkTreeSelection")
  checkPtrType(iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_selection_select_iter", object, iter, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeSelectionUnselectIter <-
function(object, iter)
{
  checkPtrType(object, "GtkTreeSelection")
  checkPtrType(iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_selection_unselect_iter", object, iter, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeSelectionPathIsSelected <-
function(object, path)
{
  checkPtrType(object, "GtkTreeSelection")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_selection_path_is_selected", object, path, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeSelectionIterIsSelected <-
function(object, iter)
{
  checkPtrType(object, "GtkTreeSelection")
  checkPtrType(iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_selection_iter_is_selected", object, iter, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeSelectionSelectAll <-
function(object)
{
  checkPtrType(object, "GtkTreeSelection")

  w <- .RGtkCall("S_gtk_tree_selection_select_all", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeSelectionUnselectAll <-
function(object)
{
  checkPtrType(object, "GtkTreeSelection")

  w <- .RGtkCall("S_gtk_tree_selection_unselect_all", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeSelectionSelectRange <-
function(object, start.path, end.path)
{
  checkPtrType(object, "GtkTreeSelection")
  checkPtrType(start.path, "GtkTreePath")
  checkPtrType(end.path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_selection_select_range", object, start.path, end.path, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeSelectionUnselectRange <-
function(object, start.path, end.path)
{
  checkPtrType(object, "GtkTreeSelection")
  checkPtrType(start.path, "GtkTreePath")
  checkPtrType(end.path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_selection_unselect_range", object, start.path, end.path, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeSortableGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_tree_sortable_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeSortableSortColumnChanged <-
function(object)
{
  checkPtrType(object, "GtkTreeSortable")

  w <- .RGtkCall("S_gtk_tree_sortable_sort_column_changed", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeSortableGetSortColumnId <-
function(object)
{
  checkPtrType(object, "GtkTreeSortable")

  w <- .RGtkCall("S_gtk_tree_sortable_get_sort_column_id", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeSortableSetSortColumnId <-
function(object, sort.column.id, order)
{
  checkPtrType(object, "GtkTreeSortable")
  sort.column.id <- as.integer(sort.column.id)
  

  w <- .RGtkCall("S_gtk_tree_sortable_set_sort_column_id", object, sort.column.id, order, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeSortableSetSortFunc <-
function(object, sort.column.id, sort.func, user.data = NULL)
{
  checkPtrType(object, "GtkTreeSortable")
  sort.column.id <- as.integer(sort.column.id)
  sort.func <- as.function(sort.func)
  

  w <- .RGtkCall("S_gtk_tree_sortable_set_sort_func", object, sort.column.id, sort.func, user.data, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeSortableSetDefaultSortFunc <-
function(object, sort.func, user.data = NULL)
{
  checkPtrType(object, "GtkTreeSortable")
  sort.func <- as.function(sort.func)
  

  w <- .RGtkCall("S_gtk_tree_sortable_set_default_sort_func", object, sort.func, user.data, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeSortableHasDefaultSortFunc <-
function(object)
{
  checkPtrType(object, "GtkTreeSortable")

  w <- .RGtkCall("S_gtk_tree_sortable_has_default_sort_func", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeStoreGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_tree_store_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeStoreNewv <-
function(types)
{
  types <- as.list(as.GType(types))

  w <- .RGtkCall("S_gtk_tree_store_newv", types, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeStoreSetColumnTypes <-
function(object, types)
{
  checkPtrType(object, "GtkTreeStore")
  types <- as.list(as.GType(types))

  w <- .RGtkCall("S_gtk_tree_store_set_column_types", object, types, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeStoreSetValue <-
function(object, iter, column, value)
{
  checkPtrType(object, "GtkTreeStore")
  checkPtrType(iter, "GtkTreeIter")
  column <- as.integer(column)
  

  w <- .RGtkCall("S_gtk_tree_store_set_value", object, iter, column, value, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeStoreRemove <-
function(object, iter)
{
  checkPtrType(object, "GtkTreeStore")
  checkPtrType(iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_store_remove", object, iter, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeStoreInsert <-
function(object, parent = NULL, position)
{
  checkPtrType(object, "GtkTreeStore")
  if (!is.null( parent )) checkPtrType(parent, "GtkTreeIter")
  position <- as.integer(position)

  w <- .RGtkCall("S_gtk_tree_store_insert", object, parent, position, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeStoreInsertBefore <-
function(object, parent, sibling)
{
  checkPtrType(object, "GtkTreeStore")
  checkPtrType(parent, "GtkTreeIter")
  checkPtrType(sibling, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_store_insert_before", object, parent, sibling, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeStoreInsertAfter <-
function(object, parent, sibling)
{
  checkPtrType(object, "GtkTreeStore")
  checkPtrType(parent, "GtkTreeIter")
  checkPtrType(sibling, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_store_insert_after", object, parent, sibling, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeStorePrepend <-
function(object, parent = NULL)
{
  checkPtrType(object, "GtkTreeStore")
  if (!is.null( parent )) checkPtrType(parent, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_store_prepend", object, parent, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeStoreAppend <-
function(object, parent = NULL)
{
  checkPtrType(object, "GtkTreeStore")
  if (!is.null( parent )) checkPtrType(parent, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_store_append", object, parent, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeStoreIsAncestor <-
function(object, iter, descendant)
{
  checkPtrType(object, "GtkTreeStore")
  checkPtrType(iter, "GtkTreeIter")
  checkPtrType(descendant, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_store_is_ancestor", object, iter, descendant, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeStoreIterDepth <-
function(object, iter)
{
  checkPtrType(object, "GtkTreeStore")
  checkPtrType(iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_store_iter_depth", object, iter, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeStoreClear <-
function(object)
{
  checkPtrType(object, "GtkTreeStore")

  w <- .RGtkCall("S_gtk_tree_store_clear", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeStoreIterIsValid <-
function(object, iter)
{
  checkPtrType(object, "GtkTreeStore")
  checkPtrType(iter, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_store_iter_is_valid", object, iter, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeStoreReorder <-
function(object, parent, new.order)
{
  checkPtrType(object, "GtkTreeStore")
  checkPtrType(parent, "GtkTreeIter")
  new.order <- as.list(as.integer(new.order))

  w <- .RGtkCall("S_gtk_tree_store_reorder", object, parent, new.order, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeStoreSwap <-
function(object, a, b)
{
  checkPtrType(object, "GtkTreeStore")
  checkPtrType(a, "GtkTreeIter")
  checkPtrType(b, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_store_swap", object, a, b, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeStoreMoveAfter <-
function(object, iter, position = NULL)
{
  checkPtrType(object, "GtkTreeStore")
  checkPtrType(iter, "GtkTreeIter")
  if (!is.null( position )) checkPtrType(position, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_store_move_after", object, iter, position, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeStoreMoveBefore <-
function(object, iter, position = NULL)
{
  checkPtrType(object, "GtkTreeStore")
  checkPtrType(iter, "GtkTreeIter")
  if (!is.null( position )) checkPtrType(position, "GtkTreeIter")

  w <- .RGtkCall("S_gtk_tree_store_move_before", object, iter, position, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewColumnQueueResize <-
function(object)
{
  checkPtrType(object, "GtkTreeViewColumn")

  w <- .RGtkCall("S_gtk_tree_view_column_queue_resize", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_tree_view_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_tree_view_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkTreeViewNewWithModel <-
function(model = NULL, show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_tree_view_new_with_model", model, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkTreeViewGetModel <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_get_model", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewSetModel <-
function(object, model = NULL)
{
  checkPtrType(object, "GtkTreeView")
  if (!is.null( model )) checkPtrType(model, "GtkTreeModel")

  w <- .RGtkCall("S_gtk_tree_view_set_model", object, model, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewGetSelection <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_get_selection", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewGetHadjustment <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_get_hadjustment", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewSetHadjustment <-
function(object, adjustment)
{
  checkPtrType(object, "GtkTreeView")
  checkPtrType(adjustment, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_tree_view_set_hadjustment", object, adjustment, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewGetVadjustment <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_get_vadjustment", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewSetVadjustment <-
function(object, adjustment)
{
  checkPtrType(object, "GtkTreeView")
  checkPtrType(adjustment, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_tree_view_set_vadjustment", object, adjustment, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewGetHeadersVisible <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_get_headers_visible", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewSetHeadersVisible <-
function(object, headers.visible)
{
  checkPtrType(object, "GtkTreeView")
  headers.visible <- as.logical(headers.visible)

  w <- .RGtkCall("S_gtk_tree_view_set_headers_visible", object, headers.visible, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewColumnsAutosize <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_columns_autosize", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewSetHeadersClickable <-
function(object, active)
{
  checkPtrType(object, "GtkTreeView")
  active <- as.logical(active)

  w <- .RGtkCall("S_gtk_tree_view_set_headers_clickable", object, active, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewSetRulesHint <-
function(object, setting)
{
  checkPtrType(object, "GtkTreeView")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_tree_view_set_rules_hint", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewGetRulesHint <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_get_rules_hint", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewAppendColumn <-
function(object, column)
{
  checkPtrType(object, "GtkTreeView")
  checkPtrType(column, "GtkTreeViewColumn")

  w <- .RGtkCall("S_gtk_tree_view_append_column", object, column, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewRemoveColumn <-
function(object, column)
{
  checkPtrType(object, "GtkTreeView")
  checkPtrType(column, "GtkTreeViewColumn")

  w <- .RGtkCall("S_gtk_tree_view_remove_column", object, column, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewInsertColumn <-
function(object, column, position)
{
  checkPtrType(object, "GtkTreeView")
  checkPtrType(column, "GtkTreeViewColumn")
  position <- as.integer(position)

  w <- .RGtkCall("S_gtk_tree_view_insert_column", object, column, position, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewInsertColumnWithDataFunc <-
function(object, position, title, cell, func, data = NULL)
{
  checkPtrType(object, "GtkTreeView")
  position <- as.integer(position)
  title <- as.character(title)
  checkPtrType(cell, "GtkCellRenderer")
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_tree_view_insert_column_with_data_func", object, position, title, cell, func, data, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewGetColumn <-
function(object, n)
{
  checkPtrType(object, "GtkTreeView")
  n <- as.integer(n)

  w <- .RGtkCall("S_gtk_tree_view_get_column", object, n, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewGetColumns <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_get_columns", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewMoveColumnAfter <-
function(object, column, base.column = NULL)
{
  checkPtrType(object, "GtkTreeView")
  checkPtrType(column, "GtkTreeViewColumn")
  if (!is.null( base.column )) checkPtrType(base.column, "GtkTreeViewColumn")

  w <- .RGtkCall("S_gtk_tree_view_move_column_after", object, column, base.column, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewSetExpanderColumn <-
function(object, column)
{
  checkPtrType(object, "GtkTreeView")
  checkPtrType(column, "GtkTreeViewColumn")

  w <- .RGtkCall("S_gtk_tree_view_set_expander_column", object, column, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewGetExpanderColumn <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_get_expander_column", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewSetColumnDragFunction <-
function(object, func, user.data = NULL)
{
  checkPtrType(object, "GtkTreeView")
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_tree_view_set_column_drag_function", object, func, user.data, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewScrollToPoint <-
function(object, tree.x, tree.y)
{
  checkPtrType(object, "GtkTreeView")
  tree.x <- as.integer(tree.x)
  tree.y <- as.integer(tree.y)

  w <- .RGtkCall("S_gtk_tree_view_scroll_to_point", object, tree.x, tree.y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewScrollToCell <-
function(object, path, column = NULL, use.align = FALSE, row.align = 0, col.align = 0)
{
  checkPtrType(object, "GtkTreeView")
  checkPtrType(path, "GtkTreePath")
  if (!is.null( column )) checkPtrType(column, "GtkTreeViewColumn")
  use.align <- as.logical(use.align)
  row.align <- as.numeric(row.align)
  col.align <- as.numeric(col.align)

  w <- .RGtkCall("S_gtk_tree_view_scroll_to_cell", object, path, column, use.align, row.align, col.align, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewRowActivated <-
function(object, path, column)
{
  checkPtrType(object, "GtkTreeView")
  checkPtrType(path, "GtkTreePath")
  checkPtrType(column, "GtkTreeViewColumn")

  w <- .RGtkCall("S_gtk_tree_view_row_activated", object, path, column, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewExpandAll <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_expand_all", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewCollapseAll <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_collapse_all", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewExpandToPath <-
function(object, path)
{
  checkPtrType(object, "GtkTreeView")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_view_expand_to_path", object, path, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewExpandRow <-
function(object, path, open.all)
{
  checkPtrType(object, "GtkTreeView")
  checkPtrType(path, "GtkTreePath")
  open.all <- as.logical(open.all)

  w <- .RGtkCall("S_gtk_tree_view_expand_row", object, path, open.all, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewCollapseRow <-
function(object, path)
{
  checkPtrType(object, "GtkTreeView")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_view_collapse_row", object, path, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewMapExpandedRows <-
function(object, func, data = NULL)
{
  checkPtrType(object, "GtkTreeView")
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_tree_view_map_expanded_rows", object, func, data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewRowExpanded <-
function(object, path)
{
  checkPtrType(object, "GtkTreeView")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_view_row_expanded", object, path, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewSetReorderable <-
function(object, reorderable)
{
  checkPtrType(object, "GtkTreeView")
  reorderable <- as.logical(reorderable)

  w <- .RGtkCall("S_gtk_tree_view_set_reorderable", object, reorderable, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewGetReorderable <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_get_reorderable", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewSetCursor <-
function(object, path, focus.column = NULL, start.editing = FALSE)
{
  checkPtrType(object, "GtkTreeView")
  checkPtrType(path, "GtkTreePath")
  if (!is.null( focus.column )) checkPtrType(focus.column, "GtkTreeViewColumn")
  start.editing <- as.logical(start.editing)

  w <- .RGtkCall("S_gtk_tree_view_set_cursor", object, path, focus.column, start.editing, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewSetCursorOnCell <-
function(object, path, focus.column = NULL, focus.cell = NULL, start.editing = FALSE)
{
  checkPtrType(object, "GtkTreeView")
  checkPtrType(path, "GtkTreePath")
  if (!is.null( focus.column )) checkPtrType(focus.column, "GtkTreeViewColumn")
  if (!is.null( focus.cell )) checkPtrType(focus.cell, "GtkCellRenderer")
  start.editing <- as.logical(start.editing)

  w <- .RGtkCall("S_gtk_tree_view_set_cursor_on_cell", object, path, focus.column, focus.cell, start.editing, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewGetCursor <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_get_cursor", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewGetBinWindow <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_get_bin_window", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewGetPathAtPos <-
function(object, x, y)
{
  checkPtrType(object, "GtkTreeView")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_gtk_tree_view_get_path_at_pos", object, x, y, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewGetCellArea <-
function(object, path, column)
{
  checkPtrType(object, "GtkTreeView")
  checkPtrType(path, "GtkTreePath")
  checkPtrType(column, "GtkTreeViewColumn")

  w <- .RGtkCall("S_gtk_tree_view_get_cell_area", object, path, column, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewGetBackgroundArea <-
function(object, path, column)
{
  checkPtrType(object, "GtkTreeView")
  checkPtrType(path, "GtkTreePath")
  checkPtrType(column, "GtkTreeViewColumn")

  w <- .RGtkCall("S_gtk_tree_view_get_background_area", object, path, column, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewGetVisibleRect <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_get_visible_rect", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewGetVisibleRange <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_get_visible_range", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewWidgetToTreeCoords <-
function(object, wx, wy)
{
  checkPtrType(object, "GtkTreeView")
  wx <- as.integer(wx)
  wy <- as.integer(wy)

  w <- .RGtkCall("S_gtk_tree_view_widget_to_tree_coords", object, wx, wy, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewTreeToWidgetCoords <-
function(object, tx, ty)
{
  checkPtrType(object, "GtkTreeView")
  tx <- as.integer(tx)
  ty <- as.integer(ty)

  w <- .RGtkCall("S_gtk_tree_view_tree_to_widget_coords", object, tx, ty, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewEnableModelDragSource <-
function(object, start.button.mask, targets, actions)
{
  checkPtrType(object, "GtkTreeView")
  
  targets <- lapply(targets, function(x) { x <- as.GtkTargetEntry(x); x })
  

  w <- .RGtkCall("S_gtk_tree_view_enable_model_drag_source", object, start.button.mask, targets, actions, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewEnableModelDragDest <-
function(object, targets, actions)
{
  checkPtrType(object, "GtkTreeView")
  targets <- lapply(targets, function(x) { x <- as.GtkTargetEntry(x); x })
  

  w <- .RGtkCall("S_gtk_tree_view_enable_model_drag_dest", object, targets, actions, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewUnsetRowsDragSource <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_unset_rows_drag_source", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewUnsetRowsDragDest <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_unset_rows_drag_dest", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewSetDragDestRow <-
function(object, path, pos)
{
  checkPtrType(object, "GtkTreeView")
  checkPtrType(path, "GtkTreePath")
  

  w <- .RGtkCall("S_gtk_tree_view_set_drag_dest_row", object, path, pos, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewGetDragDestRow <-
function(object, path)
{
  checkPtrType(object, "GtkTreeView")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_view_get_drag_dest_row", object, path, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewGetDestRowAtPos <-
function(object, drag.x, drag.y)
{
  checkPtrType(object, "GtkTreeView")
  drag.x <- as.integer(drag.x)
  drag.y <- as.integer(drag.y)

  w <- .RGtkCall("S_gtk_tree_view_get_dest_row_at_pos", object, drag.x, drag.y, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewCreateRowDragIcon <-
function(object, path)
{
  checkPtrType(object, "GtkTreeView")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_view_create_row_drag_icon", object, path, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewSetEnableSearch <-
function(object, enable.search)
{
  checkPtrType(object, "GtkTreeView")
  enable.search <- as.logical(enable.search)

  w <- .RGtkCall("S_gtk_tree_view_set_enable_search", object, enable.search, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewGetEnableSearch <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_get_enable_search", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewGetSearchColumn <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_get_search_column", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewSetSearchColumn <-
function(object, column)
{
  checkPtrType(object, "GtkTreeView")
  column <- as.integer(column)

  w <- .RGtkCall("S_gtk_tree_view_set_search_column", object, column, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewGetSearchEqualFunc <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_get_search_equal_func", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewSetSearchEqualFunc <-
function(object, search.equal.func, search.user.data = NULL)
{
  checkPtrType(object, "GtkTreeView")
  search.equal.func <- as.function(search.equal.func)
  

  w <- .RGtkCall("S_gtk_tree_view_set_search_equal_func", object, search.equal.func, search.user.data, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewSetDestroyCountFunc <-
function(object, func, data = NULL)
{
  checkPtrType(object, "GtkTreeView")
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_tree_view_set_destroy_count_func", object, func, data, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewSetFixedHeightMode <-
function(object, enable)
{
  checkPtrType(object, "GtkTreeView")
  enable <- as.logical(enable)

  w <- .RGtkCall("S_gtk_tree_view_set_fixed_height_mode", object, enable, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewGetFixedHeightMode <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_get_fixed_height_mode", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewSetHoverSelection <-
function(object, hover)
{
  checkPtrType(object, "GtkTreeView")
  hover <- as.logical(hover)

  w <- .RGtkCall("S_gtk_tree_view_set_hover_selection", object, hover, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewGetHoverSelection <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_get_hover_selection", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewSetHoverExpand <-
function(object, expand)
{
  checkPtrType(object, "GtkTreeView")
  expand <- as.logical(expand)

  w <- .RGtkCall("S_gtk_tree_view_set_hover_expand", object, expand, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewGetHoverExpand <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_get_hover_expand", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewGetRowSeparatorFunc <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_get_row_separator_func", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewSetRowSeparatorFunc <-
function(object, func, data = NULL)
{
  checkPtrType(object, "GtkTreeView")
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_tree_view_set_row_separator_func", object, func, data, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewColumnGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_tree_view_column_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewColumnNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_tree_view_column_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewColumnPackStart <-
function(object, cell, expand = TRUE)
{
  checkPtrType(object, "GtkTreeViewColumn")
  checkPtrType(cell, "GtkCellRenderer")
  expand <- as.logical(expand)

  w <- .RGtkCall("S_gtk_tree_view_column_pack_start", object, cell, expand, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewColumnPackEnd <-
function(object, cell, expand = TRUE)
{
  checkPtrType(object, "GtkTreeViewColumn")
  checkPtrType(cell, "GtkCellRenderer")
  expand <- as.logical(expand)

  w <- .RGtkCall("S_gtk_tree_view_column_pack_end", object, cell, expand, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewColumnClear <-
function(object)
{
  checkPtrType(object, "GtkTreeViewColumn")

  w <- .RGtkCall("S_gtk_tree_view_column_clear", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewColumnGetCellRenderers <-
function(object)
{
  checkPtrType(object, "GtkTreeViewColumn")

  w <- .RGtkCall("S_gtk_tree_view_column_get_cell_renderers", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewColumnAddAttribute <-
function(object, cell.renderer, attribute, column)
{
  checkPtrType(object, "GtkTreeViewColumn")
  checkPtrType(cell.renderer, "GtkCellRenderer")
  attribute <- as.character(attribute)
  column <- as.integer(column)

  w <- .RGtkCall("S_gtk_tree_view_column_add_attribute", object, cell.renderer, attribute, column, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewColumnSetCellDataFunc <-
function(object, cell.renderer, func, func.data = NULL)
{
  checkPtrType(object, "GtkTreeViewColumn")
  checkPtrType(cell.renderer, "GtkCellRenderer")
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_tree_view_column_set_cell_data_func", object, cell.renderer, func, func.data, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewColumnClearAttributes <-
function(object, cell.renderer)
{
  checkPtrType(object, "GtkTreeViewColumn")
  checkPtrType(cell.renderer, "GtkCellRenderer")

  w <- .RGtkCall("S_gtk_tree_view_column_clear_attributes", object, cell.renderer, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewColumnSetSpacing <-
function(object, spacing)
{
  checkPtrType(object, "GtkTreeViewColumn")
  spacing <- as.integer(spacing)

  w <- .RGtkCall("S_gtk_tree_view_column_set_spacing", object, spacing, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewColumnGetSpacing <-
function(object)
{
  checkPtrType(object, "GtkTreeViewColumn")

  w <- .RGtkCall("S_gtk_tree_view_column_get_spacing", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewColumnSetVisible <-
function(object, visible)
{
  checkPtrType(object, "GtkTreeViewColumn")
  visible <- as.logical(visible)

  w <- .RGtkCall("S_gtk_tree_view_column_set_visible", object, visible, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewColumnGetVisible <-
function(object)
{
  checkPtrType(object, "GtkTreeViewColumn")

  w <- .RGtkCall("S_gtk_tree_view_column_get_visible", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewColumnSetResizable <-
function(object, resizable)
{
  checkPtrType(object, "GtkTreeViewColumn")
  resizable <- as.logical(resizable)

  w <- .RGtkCall("S_gtk_tree_view_column_set_resizable", object, resizable, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewColumnGetResizable <-
function(object)
{
  checkPtrType(object, "GtkTreeViewColumn")

  w <- .RGtkCall("S_gtk_tree_view_column_get_resizable", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewColumnSetSizing <-
function(object, type)
{
  checkPtrType(object, "GtkTreeViewColumn")
  

  w <- .RGtkCall("S_gtk_tree_view_column_set_sizing", object, type, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewColumnGetSizing <-
function(object)
{
  checkPtrType(object, "GtkTreeViewColumn")

  w <- .RGtkCall("S_gtk_tree_view_column_get_sizing", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewColumnGetWidth <-
function(object)
{
  checkPtrType(object, "GtkTreeViewColumn")

  w <- .RGtkCall("S_gtk_tree_view_column_get_width", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewColumnGetFixedWidth <-
function(object)
{
  checkPtrType(object, "GtkTreeViewColumn")

  w <- .RGtkCall("S_gtk_tree_view_column_get_fixed_width", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewColumnSetFixedWidth <-
function(object, fixed.width)
{
  checkPtrType(object, "GtkTreeViewColumn")
  fixed.width <- as.integer(fixed.width)

  w <- .RGtkCall("S_gtk_tree_view_column_set_fixed_width", object, fixed.width, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewColumnSetMinWidth <-
function(object, min.width)
{
  checkPtrType(object, "GtkTreeViewColumn")
  min.width <- as.integer(min.width)

  w <- .RGtkCall("S_gtk_tree_view_column_set_min_width", object, min.width, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewColumnGetMinWidth <-
function(object)
{
  checkPtrType(object, "GtkTreeViewColumn")

  w <- .RGtkCall("S_gtk_tree_view_column_get_min_width", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewColumnSetMaxWidth <-
function(object, max.width)
{
  checkPtrType(object, "GtkTreeViewColumn")
  max.width <- as.integer(max.width)

  w <- .RGtkCall("S_gtk_tree_view_column_set_max_width", object, max.width, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewColumnGetMaxWidth <-
function(object)
{
  checkPtrType(object, "GtkTreeViewColumn")

  w <- .RGtkCall("S_gtk_tree_view_column_get_max_width", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewColumnClicked <-
function(object)
{
  checkPtrType(object, "GtkTreeViewColumn")

  w <- .RGtkCall("S_gtk_tree_view_column_clicked", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewColumnSetTitle <-
function(object, title)
{
  checkPtrType(object, "GtkTreeViewColumn")
  title <- as.character(title)

  w <- .RGtkCall("S_gtk_tree_view_column_set_title", object, title, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewColumnGetTitle <-
function(object)
{
  checkPtrType(object, "GtkTreeViewColumn")

  w <- .RGtkCall("S_gtk_tree_view_column_get_title", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewColumnSetExpand <-
function(object, expand)
{
  checkPtrType(object, "GtkTreeViewColumn")
  expand <- as.logical(expand)

  w <- .RGtkCall("S_gtk_tree_view_column_set_expand", object, expand, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewColumnGetExpand <-
function(object)
{
  checkPtrType(object, "GtkTreeViewColumn")

  w <- .RGtkCall("S_gtk_tree_view_column_get_expand", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewColumnSetClickable <-
function(object, active)
{
  checkPtrType(object, "GtkTreeViewColumn")
  active <- as.logical(active)

  w <- .RGtkCall("S_gtk_tree_view_column_set_clickable", object, active, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewColumnGetClickable <-
function(object)
{
  checkPtrType(object, "GtkTreeViewColumn")

  w <- .RGtkCall("S_gtk_tree_view_column_get_clickable", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewColumnSetWidget <-
function(object, widget = NULL)
{
  checkPtrType(object, "GtkTreeViewColumn")
  if (!is.null( widget )) checkPtrType(widget, "GtkWidget")

  w <- .RGtkCall("S_gtk_tree_view_column_set_widget", object, widget, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewColumnGetWidget <-
function(object)
{
  checkPtrType(object, "GtkTreeViewColumn")

  w <- .RGtkCall("S_gtk_tree_view_column_get_widget", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewColumnSetAlignment <-
function(object, xalign)
{
  checkPtrType(object, "GtkTreeViewColumn")
  xalign <- as.numeric(xalign)

  w <- .RGtkCall("S_gtk_tree_view_column_set_alignment", object, xalign, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewColumnGetAlignment <-
function(object)
{
  checkPtrType(object, "GtkTreeViewColumn")

  w <- .RGtkCall("S_gtk_tree_view_column_get_alignment", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewColumnSetReorderable <-
function(object, reorderable)
{
  checkPtrType(object, "GtkTreeViewColumn")
  reorderable <- as.logical(reorderable)

  w <- .RGtkCall("S_gtk_tree_view_column_set_reorderable", object, reorderable, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewColumnGetReorderable <-
function(object)
{
  checkPtrType(object, "GtkTreeViewColumn")

  w <- .RGtkCall("S_gtk_tree_view_column_get_reorderable", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewColumnSetSortColumnId <-
function(object, sort.column.id)
{
  checkPtrType(object, "GtkTreeViewColumn")
  sort.column.id <- as.integer(sort.column.id)

  w <- .RGtkCall("S_gtk_tree_view_column_set_sort_column_id", object, sort.column.id, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewColumnGetSortColumnId <-
function(object)
{
  checkPtrType(object, "GtkTreeViewColumn")

  w <- .RGtkCall("S_gtk_tree_view_column_get_sort_column_id", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewColumnSetSortIndicator <-
function(object, setting)
{
  checkPtrType(object, "GtkTreeViewColumn")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_tree_view_column_set_sort_indicator", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewColumnGetSortIndicator <-
function(object)
{
  checkPtrType(object, "GtkTreeViewColumn")

  w <- .RGtkCall("S_gtk_tree_view_column_get_sort_indicator", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewColumnSetSortOrder <-
function(object, order)
{
  checkPtrType(object, "GtkTreeViewColumn")
  

  w <- .RGtkCall("S_gtk_tree_view_column_set_sort_order", object, order, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewColumnGetSortOrder <-
function(object)
{
  checkPtrType(object, "GtkTreeViewColumn")

  w <- .RGtkCall("S_gtk_tree_view_column_get_sort_order", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewColumnCellSetCellData <-
function(object, tree.model, iter, is.expander, is.expanded)
{
  checkPtrType(object, "GtkTreeViewColumn")
  checkPtrType(tree.model, "GtkTreeModel")
  checkPtrType(iter, "GtkTreeIter")
  is.expander <- as.logical(is.expander)
  is.expanded <- as.logical(is.expanded)

  w <- .RGtkCall("S_gtk_tree_view_column_cell_set_cell_data", object, tree.model, iter, is.expander, is.expanded, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewColumnCellGetSize <-
function(object)
{
  checkPtrType(object, "GtkTreeViewColumn")

  w <- .RGtkCall("S_gtk_tree_view_column_cell_get_size", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewColumnCellIsVisible <-
function(object)
{
  checkPtrType(object, "GtkTreeViewColumn")

  w <- .RGtkCall("S_gtk_tree_view_column_cell_is_visible", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewColumnFocusCell <-
function(object, cell)
{
  checkPtrType(object, "GtkTreeViewColumn")
  checkPtrType(cell, "GtkCellRenderer")

  w <- .RGtkCall("S_gtk_tree_view_column_focus_cell", object, cell, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewColumnCellGetPosition <-
function(object, cell.renderer)
{
  checkPtrType(object, "GtkTreeViewColumn")
  checkPtrType(cell.renderer, "GtkCellRenderer")

  w <- .RGtkCall("S_gtk_tree_view_column_cell_get_position", object, cell.renderer, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkUIManagerGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_ui_manager_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkUIManagerNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_ui_manager_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkUIManagerSetAddTearoffs <-
function(object, add.tearoffs)
{
  checkPtrType(object, "GtkUIManager")
  add.tearoffs <- as.logical(add.tearoffs)

  w <- .RGtkCall("S_gtk_ui_manager_set_add_tearoffs", object, add.tearoffs, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkUIManagerGetAddTearoffs <-
function(object)
{
  checkPtrType(object, "GtkUIManager")

  w <- .RGtkCall("S_gtk_ui_manager_get_add_tearoffs", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkUIManagerInsertActionGroup <-
function(object, action.group, pos)
{
  checkPtrType(object, "GtkUIManager")
  checkPtrType(action.group, "GtkActionGroup")
  pos <- as.integer(pos)

  w <- .RGtkCall("S_gtk_ui_manager_insert_action_group", object, action.group, pos, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkUIManagerRemoveActionGroup <-
function(object, action.group)
{
  checkPtrType(object, "GtkUIManager")
  checkPtrType(action.group, "GtkActionGroup")

  w <- .RGtkCall("S_gtk_ui_manager_remove_action_group", object, action.group, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkUIManagerGetActionGroups <-
function(object)
{
  checkPtrType(object, "GtkUIManager")

  w <- .RGtkCall("S_gtk_ui_manager_get_action_groups", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkUIManagerGetAccelGroup <-
function(object)
{
  checkPtrType(object, "GtkUIManager")

  w <- .RGtkCall("S_gtk_ui_manager_get_accel_group", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkUIManagerGetWidget <-
function(object, path)
{
  checkPtrType(object, "GtkUIManager")
  path <- as.character(path)

  w <- .RGtkCall("S_gtk_ui_manager_get_widget", object, path, PACKAGE = "RGtk2")

  return(w)
} 


gtkUIManagerGetToplevels <-
function(object, types)
{
  checkPtrType(object, "GtkUIManager")
  

  w <- .RGtkCall("S_gtk_ui_manager_get_toplevels", object, types, PACKAGE = "RGtk2")

  return(w)
} 


gtkUIManagerGetAction <-
function(object, path)
{
  checkPtrType(object, "GtkUIManager")
  path <- as.character(path)

  w <- .RGtkCall("S_gtk_ui_manager_get_action", object, path, PACKAGE = "RGtk2")

  return(w)
} 


gtkUIManagerAddUiFromString <-
function(object, buffer, length = -1, .errwarn = TRUE)
{
  checkPtrType(object, "GtkUIManager")
  buffer <- as.character(buffer)
  length <- as.integer(length)

  w <- .RGtkCall("S_gtk_ui_manager_add_ui_from_string", object, buffer, length, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkUIManagerAddUiFromFile <-
function(object, filename, .errwarn = TRUE)
{
  checkPtrType(object, "GtkUIManager")
  filename <- as.character(filename)

  w <- .RGtkCall("S_gtk_ui_manager_add_ui_from_file", object, filename, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkUIManagerAddUi <-
function(object, merge.id, path, name, action = NULL, type, top)
{
  checkPtrType(object, "GtkUIManager")
  merge.id <- as.numeric(merge.id)
  path <- as.character(path)
  name <- as.character(name)
  if (!is.null( action )) action <- as.character(action)
  
  top <- as.logical(top)

  w <- .RGtkCall("S_gtk_ui_manager_add_ui", object, merge.id, path, name, action, type, top, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkUIManagerRemoveUi <-
function(object, merge.id)
{
  checkPtrType(object, "GtkUIManager")
  merge.id <- as.numeric(merge.id)

  w <- .RGtkCall("S_gtk_ui_manager_remove_ui", object, merge.id, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkUIManagerGetUi <-
function(object)
{
  checkPtrType(object, "GtkUIManager")

  w <- .RGtkCall("S_gtk_ui_manager_get_ui", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkUIManagerEnsureUpdate <-
function(object)
{
  checkPtrType(object, "GtkUIManager")

  w <- .RGtkCall("S_gtk_ui_manager_ensure_update", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkUIManagerNewMergeId <-
function(object)
{
  checkPtrType(object, "GtkUIManager")

  w <- .RGtkCall("S_gtk_ui_manager_new_merge_id", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkVButtonBoxGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_vbutton_box_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkVButtonBoxNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_vbutton_box_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkVButtonBoxGetSpacingDefault <-
function()
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  

  w <- .RGtkCall("S_gtk_vbutton_box_get_spacing_default", PACKAGE = "RGtk2")

  return(w)
} 


gtkVButtonBoxSetSpacingDefault <-
function(spacing)
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  spacing <- as.integer(spacing)

  w <- .RGtkCall("S_gtk_vbutton_box_set_spacing_default", spacing, PACKAGE = "RGtk2")

  return(w)
} 


gtkVButtonBoxGetLayoutDefault <-
function()
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  

  w <- .RGtkCall("S_gtk_vbutton_box_get_layout_default", PACKAGE = "RGtk2")

  return(w)
} 


gtkVButtonBoxSetLayoutDefault <-
function(layout)
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  

  w <- .RGtkCall("S_gtk_vbutton_box_set_layout_default", layout, PACKAGE = "RGtk2")

  return(w)
} 


gtkVBoxGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_vbox_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkVBoxNew <-
function(homogeneous = NULL, spacing = NULL, show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_vbox_new", homogeneous, spacing, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkViewportGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_viewport_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkViewportNew <-
function(hadjustment = NULL, vadjustment = NULL, show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_viewport_new", hadjustment, vadjustment, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkViewportGetHadjustment <-
function(object)
{
  checkPtrType(object, "GtkViewport")

  w <- .RGtkCall("S_gtk_viewport_get_hadjustment", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkViewportGetVadjustment <-
function(object)
{
  checkPtrType(object, "GtkViewport")

  w <- .RGtkCall("S_gtk_viewport_get_vadjustment", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkViewportSetHadjustment <-
function(object, adjustment = NULL)
{
  checkPtrType(object, "GtkViewport")
  if (!is.null( adjustment )) checkPtrType(adjustment, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_viewport_set_hadjustment", object, adjustment, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkViewportSetVadjustment <-
function(object, adjustment = NULL)
{
  checkPtrType(object, "GtkViewport")
  if (!is.null( adjustment )) checkPtrType(adjustment, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_viewport_set_vadjustment", object, adjustment, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkViewportSetShadowType <-
function(object, type)
{
  checkPtrType(object, "GtkViewport")
  

  w <- .RGtkCall("S_gtk_viewport_set_shadow_type", object, type, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkViewportGetShadowType <-
function(object)
{
  checkPtrType(object, "GtkViewport")

  w <- .RGtkCall("S_gtk_viewport_get_shadow_type", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkVPanedGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_vpaned_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkVPanedNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_vpaned_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkVRulerGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_vruler_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkVRulerNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_vruler_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkVScaleGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_vscale_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkVScaleNew <-
function(adjustment = NULL, show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_vscale_new", adjustment, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkVScaleNewWithRange <-
function(min, max, step, show = TRUE)
{
  min <- as.numeric(min)
  max <- as.numeric(max)
  step <- as.numeric(step)

  w <- .RGtkCall("S_gtk_vscale_new_with_range", min, max, step, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkVScrollbarGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_vscrollbar_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkVScrollbarNew <-
function(adjustment = NULL, show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_vscrollbar_new", adjustment, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkVSeparatorGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_vseparator_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkVSeparatorNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_vseparator_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkWidgetGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_widget_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetUnparent <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_unparent", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetShow <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_show", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetShowNow <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_show_now", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetHide <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_hide", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetShowAll <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_show_all", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetHideAll <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_hide_all", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetSetNoShowAll <-
function(object, no.show.all)
{
  checkPtrType(object, "GtkWidget")
  no.show.all <- as.logical(no.show.all)

  w <- .RGtkCall("S_gtk_widget_set_no_show_all", object, no.show.all, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetGetNoShowAll <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_no_show_all", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetMap <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_map", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetUnmap <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_unmap", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetRealize <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_realize", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetUnrealize <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_unrealize", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetQueueDraw <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_queue_draw", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetQueueDrawArea <-
function(object, x, y, width, height)
{
  checkPtrType(object, "GtkWidget")
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_widget_queue_draw_area", object, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetQueueClear <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("gtkWidgetQueueDraw", "RGtk2")

  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_queue_clear", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetQueueClearArea <-
function(object, x, y, width, height)
{
  if(getOption("depwarn"))
    .Deprecated("gtkWidgetQueueDrawArea", "RGtk2")

  checkPtrType(object, "GtkWidget")
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_widget_queue_clear_area", object, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetQueueResize <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_queue_resize", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetQueueResizeNoRedraw <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_queue_resize_no_redraw", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetDraw <-
function(object, area)
{
  if(getOption("depwarn"))
    .Deprecated("gtkWidgetQueueDrawArea", "RGtk2")

  checkPtrType(object, "GtkWidget")
  area <- as.GdkRectangle(area)

  w <- .RGtkCall("S_gtk_widget_draw", object, area, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetSizeRequest <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_size_request", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetSizeAllocate <-
function(object, allocation)
{
  checkPtrType(object, "GtkWidget")
  allocation <- as.GtkAllocation(allocation)

  w <- .RGtkCall("S_gtk_widget_size_allocate", object, allocation, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetGetChildRequisition <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_child_requisition", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetAddAccelerator <-
function(object, accel.signal, accel.group, accel.key, accel.mods, accel.flags)
{
  checkPtrType(object, "GtkWidget")
  accel.signal <- as.character(accel.signal)
  checkPtrType(accel.group, "GtkAccelGroup")
  accel.key <- as.numeric(accel.key)
  
  

  w <- .RGtkCall("S_gtk_widget_add_accelerator", object, accel.signal, accel.group, accel.key, accel.mods, accel.flags, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetRemoveAccelerator <-
function(object, accel.group, accel.key, accel.mods)
{
  checkPtrType(object, "GtkWidget")
  checkPtrType(accel.group, "GtkAccelGroup")
  accel.key <- as.numeric(accel.key)
  

  w <- .RGtkCall("S_gtk_widget_remove_accelerator", object, accel.group, accel.key, accel.mods, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetSetAccelPath <-
function(object, accel.path, accel.group)
{
  checkPtrType(object, "GtkWidget")
  accel.path <- as.character(accel.path)
  checkPtrType(accel.group, "GtkAccelGroup")

  w <- .RGtkCall("S_gtk_widget_set_accel_path", object, accel.path, accel.group, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetListAccelClosures <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_list_accel_closures", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetCanActivateAccel <-
function(object, signal.id)
{
  checkPtrType(object, "GtkWidget")
  signal.id <- as.numeric(signal.id)

  w <- .RGtkCall("S_gtk_widget_can_activate_accel", object, signal.id, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetMnemonicActivate <-
function(object, group.cycling)
{
  checkPtrType(object, "GtkWidget")
  group.cycling <- as.logical(group.cycling)

  w <- .RGtkCall("S_gtk_widget_mnemonic_activate", object, group.cycling, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetEvent <-
function(object, event)
{
  checkPtrType(object, "GtkWidget")
  checkPtrType(event, "GdkEvent")

  w <- .RGtkCall("S_gtk_widget_event", object, event, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetSendExpose <-
function(object, event)
{
  checkPtrType(object, "GtkWidget")
  checkPtrType(event, "GdkEvent")

  w <- .RGtkCall("S_gtk_widget_send_expose", object, event, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetActivate <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_activate", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetSetScrollAdjustments <-
function(object, hadjustment = NULL, vadjustment = NULL)
{
  checkPtrType(object, "GtkWidget")
  if (!is.null( hadjustment )) checkPtrType(hadjustment, "GtkAdjustment")
  if (!is.null( vadjustment )) checkPtrType(vadjustment, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_widget_set_scroll_adjustments", object, hadjustment, vadjustment, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetReparent <-
function(object, new.parent)
{
  checkPtrType(object, "GtkWidget")
  checkPtrType(new.parent, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_reparent", object, new.parent, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetIntersect <-
function(object, area, intersection)
{
  checkPtrType(object, "GtkWidget")
  area <- as.GdkRectangle(area)
  intersection <- as.GdkRectangle(intersection)

  w <- .RGtkCall("S_gtk_widget_intersect", object, area, intersection, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetRegionIntersect <-
function(object, region)
{
  checkPtrType(object, "GtkWidget")
  checkPtrType(region, "GdkRegion")

  w <- .RGtkCall("S_gtk_widget_region_intersect", object, region, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetFreezeChildNotify <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_freeze_child_notify", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetChildNotify <-
function(object, child.property)
{
  checkPtrType(object, "GtkWidget")
  child.property <- as.character(child.property)

  w <- .RGtkCall("S_gtk_widget_child_notify", object, child.property, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetThawChildNotify <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_thaw_child_notify", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetIsFocus <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_is_focus", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetGrabFocus <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_grab_focus", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetGrabDefault <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_grab_default", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetSetName <-
function(object, name)
{
  checkPtrType(object, "GtkWidget")
  name <- as.character(name)

  w <- .RGtkCall("S_gtk_widget_set_name", object, name, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetGetName <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetSetState <-
function(object, state)
{
  checkPtrType(object, "GtkWidget")
  

  w <- .RGtkCall("S_gtk_widget_set_state", object, state, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetSetSensitive <-
function(object, sensitive)
{
  checkPtrType(object, "GtkWidget")
  sensitive <- as.logical(sensitive)

  w <- .RGtkCall("S_gtk_widget_set_sensitive", object, sensitive, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetSetAppPaintable <-
function(object, app.paintable)
{
  checkPtrType(object, "GtkWidget")
  app.paintable <- as.logical(app.paintable)

  w <- .RGtkCall("S_gtk_widget_set_app_paintable", object, app.paintable, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetSetDoubleBuffered <-
function(object, double.buffered)
{
  checkPtrType(object, "GtkWidget")
  double.buffered <- as.logical(double.buffered)

  w <- .RGtkCall("S_gtk_widget_set_double_buffered", object, double.buffered, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetSetRedrawOnAllocate <-
function(object, redraw.on.allocate)
{
  checkPtrType(object, "GtkWidget")
  redraw.on.allocate <- as.logical(redraw.on.allocate)

  w <- .RGtkCall("S_gtk_widget_set_redraw_on_allocate", object, redraw.on.allocate, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetSetParent <-
function(object, parent)
{
  checkPtrType(object, "GtkWidget")
  checkPtrType(parent, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_set_parent", object, parent, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetSetParentWindow <-
function(object, parent.window)
{
  checkPtrType(object, "GtkWidget")
  checkPtrType(parent.window, "GdkWindow")

  w <- .RGtkCall("S_gtk_widget_set_parent_window", object, parent.window, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetSetChildVisible <-
function(object, is.visible)
{
  checkPtrType(object, "GtkWidget")
  is.visible <- as.logical(is.visible)

  w <- .RGtkCall("S_gtk_widget_set_child_visible", object, is.visible, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetGetChildVisible <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_child_visible", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetGetParent <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_parent", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetGetParentWindow <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_parent_window", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetChildFocus <-
function(object, direction)
{
  checkPtrType(object, "GtkWidget")
  

  w <- .RGtkCall("S_gtk_widget_child_focus", object, direction, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetSetSizeRequest <-
function(object, width, height)
{
  checkPtrType(object, "GtkWidget")
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_widget_set_size_request", object, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetGetSizeRequest <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_size_request", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetSetUposition <-
function(object, x, y)
{
  checkPtrType(object, "GtkWidget")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_gtk_widget_set_uposition", object, x, y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetSetUsize <-
function(object, width, height)
{
  if(getOption("depwarn"))
    .Deprecated("gtkWidgetSetSizeRequest", "RGtk2")

  checkPtrType(object, "GtkWidget")
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_widget_set_usize", object, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetSetEvents <-
function(object, events)
{
  checkPtrType(object, "GtkWidget")
  events <- as.integer(events)

  w <- .RGtkCall("S_gtk_widget_set_events", object, events, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetAddEvents <-
function(object, events)
{
  checkPtrType(object, "GtkWidget")
  events <- as.integer(events)

  w <- .RGtkCall("S_gtk_widget_add_events", object, events, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetSetExtensionEvents <-
function(object, mode)
{
  checkPtrType(object, "GtkWidget")
  

  w <- .RGtkCall("S_gtk_widget_set_extension_events", object, mode, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetGetExtensionEvents <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_extension_events", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetGetToplevel <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_toplevel", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetGetAncestor <-
function(object, widget.type)
{
  checkPtrType(object, "GtkWidget")
  widget.type <- as.GType(widget.type)

  w <- .RGtkCall("S_gtk_widget_get_ancestor", object, widget.type, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetGetColormap <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_colormap", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetGetVisual <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_visual", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetGetScreen <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_screen", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetHasScreen <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_has_screen", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetGetDisplay <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_display", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetGetRootWindow <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_root_window", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetGetSettings <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_settings", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetGetClipboard <-
function(object, selection)
{
  checkPtrType(object, "GtkWidget")
  selection <- as.GdkAtom(selection)

  w <- .RGtkCall("S_gtk_widget_get_clipboard", object, selection, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetGetAccessible <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_accessible", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetSetColormap <-
function(object, colormap)
{
  checkPtrType(object, "GtkWidget")
  checkPtrType(colormap, "GdkColormap")

  w <- .RGtkCall("S_gtk_widget_set_colormap", object, colormap, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetGetEvents <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_events", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetGetPointer <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_pointer", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetIsAncestor <-
function(object, ancestor)
{
  checkPtrType(object, "GtkWidget")
  checkPtrType(ancestor, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_is_ancestor", object, ancestor, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetTranslateCoordinates <-
function(object, dest.widget, src.x, src.y)
{
  checkPtrType(object, "GtkWidget")
  checkPtrType(dest.widget, "GtkWidget")
  src.x <- as.integer(src.x)
  src.y <- as.integer(src.y)

  w <- .RGtkCall("S_gtk_widget_translate_coordinates", object, dest.widget, src.x, src.y, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetHideOnDelete <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_hide_on_delete", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetSetStyle <-
function(object, style = NULL)
{
  checkPtrType(object, "GtkWidget")
  if (!is.null( style )) checkPtrType(style, "GtkStyle")

  w <- .RGtkCall("S_gtk_widget_set_style", object, style, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetEnsureStyle <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_ensure_style", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetGetStyle <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_style", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetModifyStyle <-
function(object, style)
{
  checkPtrType(object, "GtkWidget")
  checkPtrType(style, "GtkRcStyle")

  w <- .RGtkCall("S_gtk_widget_modify_style", object, style, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetGetModifierStyle <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_modifier_style", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetModifyFg <-
function(object, state, color = NULL)
{
  checkPtrType(object, "GtkWidget")
  
  if (!is.null( color )) color <- as.GdkColor(color)

  w <- .RGtkCall("S_gtk_widget_modify_fg", object, state, color, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetModifyBg <-
function(object, state, color = NULL)
{
  checkPtrType(object, "GtkWidget")
  
  if (!is.null( color )) color <- as.GdkColor(color)

  w <- .RGtkCall("S_gtk_widget_modify_bg", object, state, color, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetModifyText <-
function(object, state, color = NULL)
{
  checkPtrType(object, "GtkWidget")
  
  if (!is.null( color )) color <- as.GdkColor(color)

  w <- .RGtkCall("S_gtk_widget_modify_text", object, state, color, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetModifyBase <-
function(object, state, color = NULL)
{
  checkPtrType(object, "GtkWidget")
  
  if (!is.null( color )) color <- as.GdkColor(color)

  w <- .RGtkCall("S_gtk_widget_modify_base", object, state, color, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetModifyFont <-
function(object, font.desc = NULL)
{
  checkPtrType(object, "GtkWidget")
  if (!is.null( font.desc )) checkPtrType(font.desc, "PangoFontDescription")

  w <- .RGtkCall("S_gtk_widget_modify_font", object, font.desc, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetCreatePangoContext <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_create_pango_context", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetGetPangoContext <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_pango_context", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetCreatePangoLayout <-
function(object, text)
{
  checkPtrType(object, "GtkWidget")
  text <- as.character(text)

  w <- .RGtkCall("S_gtk_widget_create_pango_layout", object, text, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetRenderIcon <-
function(object, stock.id, size, detail = NULL)
{
  checkPtrType(object, "GtkWidget")
  stock.id <- as.character(stock.id)
  
  if (!is.null( detail )) detail <- as.character(detail)

  w <- .RGtkCall("S_gtk_widget_render_icon", object, stock.id, size, detail, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetSetCompositeName <-
function(object, name)
{
  checkPtrType(object, "GtkWidget")
  name <- as.character(name)

  w <- .RGtkCall("S_gtk_widget_set_composite_name", object, name, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetGetCompositeName <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_composite_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetResetRcStyles <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_reset_rc_styles", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetPushColormap <-
function(cmap)
{
  checkPtrType(cmap, "GdkColormap")

  w <- .RGtkCall("S_gtk_widget_push_colormap", cmap, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetPushCompositeChild <-
function()
{
  

  w <- .RGtkCall("S_gtk_widget_push_composite_child", PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetPopCompositeChild <-
function()
{
  

  w <- .RGtkCall("S_gtk_widget_pop_composite_child", PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetPopColormap <-
function()
{
  

  w <- .RGtkCall("S_gtk_widget_pop_colormap", PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetClassInstallStyleProperty <-
function(klass, pspec)
{
  checkPtrType(klass, "GtkWidgetClass")
  pspec <- as.GParamSpec(pspec)

  w <- .RGtkCall("S_gtk_widget_class_install_style_property", klass, pspec, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetClassFindStyleProperty <-
function(klass, property.name)
{
  checkPtrType(klass, "GtkWidgetClass")
  property.name <- as.character(property.name)

  w <- .RGtkCall("S_gtk_widget_class_find_style_property", klass, property.name, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetClassListStyleProperties <-
function(klass)
{
  checkPtrType(klass, "GtkWidgetClass")

  w <- .RGtkCall("S_gtk_widget_class_list_style_properties", klass, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetStyleGetProperty <-
function(object, property.name)
{
  checkPtrType(object, "GtkWidget")
  property.name <- as.character(property.name)

  w <- .RGtkCall("S_gtk_widget_style_get_property", object, property.name, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetGetDefaultStyle <-
function()
{
  

  w <- .RGtkCall("S_gtk_widget_get_default_style", PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetSetDefaultColormap <-
function(colormap)
{
  checkPtrType(colormap, "GdkColormap")

  w <- .RGtkCall("S_gtk_widget_set_default_colormap", colormap, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetGetDefaultColormap <-
function()
{
  

  w <- .RGtkCall("S_gtk_widget_get_default_colormap", PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetGetDefaultVisual <-
function()
{
  

  w <- .RGtkCall("S_gtk_widget_get_default_visual", PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetSetDirection <-
function(object, dir)
{
  checkPtrType(object, "GtkWidget")
  

  w <- .RGtkCall("S_gtk_widget_set_direction", object, dir, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetGetDirection <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_direction", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetSetDefaultDirection <-
function(dir)
{
  

  w <- .RGtkCall("S_gtk_widget_set_default_direction", dir, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetGetDefaultDirection <-
function()
{
  

  w <- .RGtkCall("S_gtk_widget_get_default_direction", PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetShapeCombineMask <-
function(object, shape.mask, offset.x, offset.y)
{
  checkPtrType(object, "GtkWidget")
  checkPtrType(shape.mask, "GdkBitmap")
  offset.x <- as.integer(offset.x)
  offset.y <- as.integer(offset.y)

  w <- .RGtkCall("S_gtk_widget_shape_combine_mask", object, shape.mask, offset.x, offset.y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetResetShapes <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_reset_shapes", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetPath <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_path", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetClassPath <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_class_path", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetListMnemonicLabels <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_list_mnemonic_labels", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetAddMnemonicLabel <-
function(object, label)
{
  checkPtrType(object, "GtkWidget")
  checkPtrType(label, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_add_mnemonic_label", object, label, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetRemoveMnemonicLabel <-
function(object, label)
{
  checkPtrType(object, "GtkWidget")
  checkPtrType(label, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_remove_mnemonic_label", object, label, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRequisitionGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_requisition_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkRequisitionCopy <-
function(object)
{
  checkPtrType(object, "GtkRequisition")

  w <- .RGtkCall("S_gtk_requisition_copy", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_window_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowNew <-
function(type = NULL, show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_window_new", type, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkWindowSetTitle <-
function(object, title)
{
  checkPtrType(object, "GtkWindow")
  title <- as.character(title)

  w <- .RGtkCall("S_gtk_window_set_title", object, title, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowGetTitle <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_get_title", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowSetWmclass <-
function(object, wmclass.name, wmclass.class)
{
  checkPtrType(object, "GtkWindow")
  wmclass.name <- as.character(wmclass.name)
  wmclass.class <- as.character(wmclass.class)

  w <- .RGtkCall("S_gtk_window_set_wmclass", object, wmclass.name, wmclass.class, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowSetRole <-
function(object, role)
{
  checkPtrType(object, "GtkWindow")
  role <- as.character(role)

  w <- .RGtkCall("S_gtk_window_set_role", object, role, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowGetRole <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_get_role", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowAddAccelGroup <-
function(object, accel.group)
{
  checkPtrType(object, "GtkWindow")
  checkPtrType(accel.group, "GtkAccelGroup")

  w <- .RGtkCall("S_gtk_window_add_accel_group", object, accel.group, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowRemoveAccelGroup <-
function(object, accel.group)
{
  checkPtrType(object, "GtkWindow")
  checkPtrType(accel.group, "GtkAccelGroup")

  w <- .RGtkCall("S_gtk_window_remove_accel_group", object, accel.group, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowSetPosition <-
function(object, position)
{
  checkPtrType(object, "GtkWindow")
  

  w <- .RGtkCall("S_gtk_window_set_position", object, position, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowActivateFocus <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_activate_focus", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowSetFocus <-
function(object, focus = NULL)
{
  checkPtrType(object, "GtkWindow")
  if (!is.null( focus )) checkPtrType(focus, "GtkWidget")

  w <- .RGtkCall("S_gtk_window_set_focus", object, focus, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowGetFocus <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_get_focus", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowSetDefault <-
function(object, default.widget = NULL)
{
  checkPtrType(object, "GtkWindow")
  if (!is.null( default.widget )) checkPtrType(default.widget, "GtkWidget")

  w <- .RGtkCall("S_gtk_window_set_default", object, default.widget, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowActivateDefault <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_activate_default", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowSetTransientFor <-
function(object, parent = NULL)
{
  checkPtrType(object, "GtkWindow")
  if (!is.null( parent )) checkPtrType(parent, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_set_transient_for", object, parent, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowGetTransientFor <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_get_transient_for", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowSetTypeHint <-
function(object, hint)
{
  checkPtrType(object, "GtkWindow")
  

  w <- .RGtkCall("S_gtk_window_set_type_hint", object, hint, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowGetTypeHint <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_get_type_hint", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowSetSkipTaskbarHint <-
function(object, setting)
{
  checkPtrType(object, "GtkWindow")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_window_set_skip_taskbar_hint", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowGetSkipTaskbarHint <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_get_skip_taskbar_hint", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowSetSkipPagerHint <-
function(object, setting)
{
  checkPtrType(object, "GtkWindow")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_window_set_skip_pager_hint", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowGetSkipPagerHint <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_get_skip_pager_hint", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowSetUrgencyHint <-
function(object, setting)
{
  checkPtrType(object, "GtkWindow")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_window_set_urgency_hint", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowGetUrgencyHint <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_get_urgency_hint", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowSetAcceptFocus <-
function(object, setting)
{
  checkPtrType(object, "GtkWindow")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_window_set_accept_focus", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowGetAcceptFocus <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_get_accept_focus", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowSetDestroyWithParent <-
function(object, setting)
{
  checkPtrType(object, "GtkWindow")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_window_set_destroy_with_parent", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowGetDestroyWithParent <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_get_destroy_with_parent", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowSetResizable <-
function(object, resizable)
{
  checkPtrType(object, "GtkWindow")
  resizable <- as.logical(resizable)

  w <- .RGtkCall("S_gtk_window_set_resizable", object, resizable, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowGetResizable <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_get_resizable", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowSetGravity <-
function(object, gravity)
{
  checkPtrType(object, "GtkWindow")
  

  w <- .RGtkCall("S_gtk_window_set_gravity", object, gravity, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowGetGravity <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_get_gravity", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowSetGeometryHints <-
function(object, geometry.widget, geometry)
{
  checkPtrType(object, "GtkWindow")
  checkPtrType(geometry.widget, "GtkWidget")
  geometry <- as.GdkGeometry(geometry)

  w <- .RGtkCall("S_gtk_window_set_geometry_hints", object, geometry.widget, geometry, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowSetScreen <-
function(object, screen)
{
  checkPtrType(object, "GtkWindow")
  checkPtrType(screen, "GdkScreen")

  w <- .RGtkCall("S_gtk_window_set_screen", object, screen, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowGetScreen <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_get_screen", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowIsActive <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_is_active", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowHasToplevelFocus <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_has_toplevel_focus", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowSetHasFrame <-
function(object, setting)
{
  checkPtrType(object, "GtkWindow")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_window_set_has_frame", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowGetHasFrame <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_get_has_frame", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowSetFrameDimensions <-
function(object, left, top, right, bottom)
{
  checkPtrType(object, "GtkWindow")
  left <- as.integer(left)
  top <- as.integer(top)
  right <- as.integer(right)
  bottom <- as.integer(bottom)

  w <- .RGtkCall("S_gtk_window_set_frame_dimensions", object, left, top, right, bottom, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowGetFrameDimensions <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_get_frame_dimensions", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowSetDecorated <-
function(object, setting)
{
  checkPtrType(object, "GtkWindow")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_window_set_decorated", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowGetDecorated <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_get_decorated", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowSetIconList <-
function(object, list)
{
  checkPtrType(object, "GtkWindow")
  list <- as.GList(list)

  w <- .RGtkCall("S_gtk_window_set_icon_list", object, list, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowGetIconList <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_get_icon_list", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowSetIcon <-
function(object, icon = NULL)
{
  checkPtrType(object, "GtkWindow")
  if (!is.null( icon )) checkPtrType(icon, "GdkPixbuf")

  w <- .RGtkCall("S_gtk_window_set_icon", object, icon, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowSetIconFromFile <-
function(object, filename, .errwarn = TRUE)
{
  checkPtrType(object, "GtkWindow")
  filename <- as.character(filename)

  w <- .RGtkCall("S_gtk_window_set_icon_from_file", object, filename, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkWindowGetIcon <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_get_icon", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowSetDefaultIconList <-
function(list)
{
  list <- as.GList(list)

  w <- .RGtkCall("S_gtk_window_set_default_icon_list", list, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowGetDefaultIconList <-
function()
{
  

  w <- .RGtkCall("S_gtk_window_get_default_icon_list", PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowSetDefaultIcon <-
function(icon)
{
  checkPtrType(icon, "GdkPixbuf")

  w <- .RGtkCall("S_gtk_window_set_default_icon", icon, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowSetDefaultIconFromFile <-
function(filename, .errwarn = TRUE)
{
  filename <- as.character(filename)

  w <- .RGtkCall("S_gtk_window_set_default_icon_from_file", filename, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(invisible(w))
} 


gtkWindowSetAutoStartupNotification <-
function(setting)
{
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_window_set_auto_startup_notification", setting, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowSetModal <-
function(object, modal)
{
  checkPtrType(object, "GtkWindow")
  modal <- as.logical(modal)

  w <- .RGtkCall("S_gtk_window_set_modal", object, modal, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowGetModal <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_get_modal", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowListToplevels <-
function()
{
  

  w <- .RGtkCall("S_gtk_window_list_toplevels", PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowAddMnemonic <-
function(object, keyval, target)
{
  checkPtrType(object, "GtkWindow")
  keyval <- as.numeric(keyval)
  checkPtrType(target, "GtkWidget")

  w <- .RGtkCall("S_gtk_window_add_mnemonic", object, keyval, target, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowRemoveMnemonic <-
function(object, keyval, target)
{
  checkPtrType(object, "GtkWindow")
  keyval <- as.numeric(keyval)
  checkPtrType(target, "GtkWidget")

  w <- .RGtkCall("S_gtk_window_remove_mnemonic", object, keyval, target, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowMnemonicActivate <-
function(object, keyval, modifier)
{
  checkPtrType(object, "GtkWindow")
  keyval <- as.numeric(keyval)
  

  w <- .RGtkCall("S_gtk_window_mnemonic_activate", object, keyval, modifier, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowSetMnemonicModifier <-
function(object, modifier)
{
  checkPtrType(object, "GtkWindow")
  

  w <- .RGtkCall("S_gtk_window_set_mnemonic_modifier", object, modifier, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowGetMnemonicModifier <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_get_mnemonic_modifier", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowActivateKey <-
function(object, event)
{
  checkPtrType(object, "GtkWindow")
  checkPtrType(event, "GdkEventKey")

  w <- .RGtkCall("S_gtk_window_activate_key", object, event, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowPropagateKeyEvent <-
function(object, event)
{
  checkPtrType(object, "GtkWindow")
  checkPtrType(event, "GdkEventKey")

  w <- .RGtkCall("S_gtk_window_propagate_key_event", object, event, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowPresent <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_present", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowPresentWithTime <-
function(object, timestamp)
{
  checkPtrType(object, "GtkWindow")
  timestamp <- as.numeric(timestamp)

  w <- .RGtkCall("S_gtk_window_present_with_time", object, timestamp, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowIconify <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_iconify", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowDeiconify <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_deiconify", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowStick <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_stick", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowUnstick <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_unstick", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowMaximize <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_maximize", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowUnmaximize <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_unmaximize", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowFullscreen <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_fullscreen", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowUnfullscreen <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_unfullscreen", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowSetKeepAbove <-
function(object, setting)
{
  checkPtrType(object, "GtkWindow")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_window_set_keep_above", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowSetKeepBelow <-
function(object, setting)
{
  checkPtrType(object, "GtkWindow")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_window_set_keep_below", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowBeginResizeDrag <-
function(object, edge, button, root.x, root.y, timestamp)
{
  checkPtrType(object, "GtkWindow")
  
  button <- as.integer(button)
  root.x <- as.integer(root.x)
  root.y <- as.integer(root.y)
  timestamp <- as.numeric(timestamp)

  w <- .RGtkCall("S_gtk_window_begin_resize_drag", object, edge, button, root.x, root.y, timestamp, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowBeginMoveDrag <-
function(object, button, root.x, root.y, timestamp)
{
  checkPtrType(object, "GtkWindow")
  button <- as.integer(button)
  root.x <- as.integer(root.x)
  root.y <- as.integer(root.y)
  timestamp <- as.numeric(timestamp)

  w <- .RGtkCall("S_gtk_window_begin_move_drag", object, button, root.x, root.y, timestamp, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowSetPolicy <-
function(object, allow.shrink, allow.grow, auto.shrink)
{
  if(getOption("depwarn"))
    .Deprecated("gtkWindowSetResizable", "RGtk2")

  checkPtrType(object, "GtkWindow")
  allow.shrink <- as.integer(allow.shrink)
  allow.grow <- as.integer(allow.grow)
  auto.shrink <- as.integer(auto.shrink)

  w <- .RGtkCall("S_gtk_window_set_policy", object, allow.shrink, allow.grow, auto.shrink, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowSetDefaultSize <-
function(object, width, height)
{
  checkPtrType(object, "GtkWindow")
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_window_set_default_size", object, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowGetDefaultSize <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_get_default_size", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowResize <-
function(object, width, height)
{
  checkPtrType(object, "GtkWindow")
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_window_resize", object, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowGetSize <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_get_size", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowMove <-
function(object, x, y)
{
  checkPtrType(object, "GtkWindow")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_gtk_window_move", object, x, y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowGetPosition <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_get_position", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowReshowWithInitialSize <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_reshow_with_initial_size", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowGroupGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_window_group_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowGroupNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_window_group_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowGroupAddWindow <-
function(object, window)
{
  checkPtrType(object, "GtkWindowGroup")
  checkPtrType(window, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_group_add_window", object, window, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowGroupRemoveWindow <-
function(object, window)
{
  checkPtrType(object, "GtkWindowGroup")
  checkPtrType(window, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_group_remove_window", object, window, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowSetFocusOnMap <-
function(object, setting)
{
  checkPtrType(object, "GtkWindow")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_window_set_focus_on_map", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowGetFocusOnMap <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_get_focus_on_map", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowSetIconName <-
function(object, name = NULL)
{
  checkPtrType(object, "GtkWindow")
  if (!is.null( name )) name <- as.character(name)

  w <- .RGtkCall("S_gtk_window_set_icon_name", object, name, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowGetIconName <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_get_icon_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowSetDefaultIconName <-
function(name)
{
  name <- as.character(name)

  w <- .RGtkCall("S_gtk_window_set_default_icon_name", name, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetGetAction <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_action", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkAssistantGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_assistant_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkAssistantNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_assistant_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkAssistantGetCurrentPage <-
function(object)
{
  checkPtrType(object, "GtkAssistant")

  w <- .RGtkCall("S_gtk_assistant_get_current_page", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkAssistantSetCurrentPage <-
function(object, page.num)
{
  checkPtrType(object, "GtkAssistant")
  page.num <- as.integer(page.num)

  w <- .RGtkCall("S_gtk_assistant_set_current_page", object, page.num, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAssistantGetNPages <-
function(object)
{
  checkPtrType(object, "GtkAssistant")

  w <- .RGtkCall("S_gtk_assistant_get_n_pages", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkAssistantGetNthPage <-
function(object, page.num)
{
  checkPtrType(object, "GtkAssistant")
  page.num <- as.integer(page.num)

  w <- .RGtkCall("S_gtk_assistant_get_nth_page", object, page.num, PACKAGE = "RGtk2")

  return(w)
} 


gtkAssistantPrependPage <-
function(object, page)
{
  checkPtrType(object, "GtkAssistant")
  checkPtrType(page, "GtkWidget")

  w <- .RGtkCall("S_gtk_assistant_prepend_page", object, page, PACKAGE = "RGtk2")

  return(w)
} 


gtkAssistantAppendPage <-
function(object, page)
{
  checkPtrType(object, "GtkAssistant")
  checkPtrType(page, "GtkWidget")

  w <- .RGtkCall("S_gtk_assistant_append_page", object, page, PACKAGE = "RGtk2")

  return(w)
} 


gtkAssistantInsertPage <-
function(object, page, position)
{
  checkPtrType(object, "GtkAssistant")
  checkPtrType(page, "GtkWidget")
  position <- as.integer(position)

  w <- .RGtkCall("S_gtk_assistant_insert_page", object, page, position, PACKAGE = "RGtk2")

  return(w)
} 


gtkAssistantSetForwardPageFunc <-
function(object, page.func, data)
{
  checkPtrType(object, "GtkAssistant")
  page.func <- as.function(page.func)
  

  w <- .RGtkCall("S_gtk_assistant_set_forward_page_func", object, page.func, data, PACKAGE = "RGtk2")

  return(w)
} 


gtkAssistantSetPageType <-
function(object, page, type)
{
  checkPtrType(object, "GtkAssistant")
  checkPtrType(page, "GtkWidget")
  

  w <- .RGtkCall("S_gtk_assistant_set_page_type", object, page, type, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAssistantGetPageType <-
function(object, page)
{
  checkPtrType(object, "GtkAssistant")
  checkPtrType(page, "GtkWidget")

  w <- .RGtkCall("S_gtk_assistant_get_page_type", object, page, PACKAGE = "RGtk2")

  return(w)
} 


gtkAssistantSetPageTitle <-
function(object, page, title)
{
  checkPtrType(object, "GtkAssistant")
  checkPtrType(page, "GtkWidget")
  title <- as.character(title)

  w <- .RGtkCall("S_gtk_assistant_set_page_title", object, page, title, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAssistantGetPageTitle <-
function(object, page)
{
  checkPtrType(object, "GtkAssistant")
  checkPtrType(page, "GtkWidget")

  w <- .RGtkCall("S_gtk_assistant_get_page_title", object, page, PACKAGE = "RGtk2")

  return(w)
} 


gtkAssistantSetPageHeaderImage <-
function(object, page, pixbuf = NULL)
{
  checkPtrType(object, "GtkAssistant")
  checkPtrType(page, "GtkWidget")
  if (!is.null( pixbuf )) checkPtrType(pixbuf, "GdkPixbuf")

  w <- .RGtkCall("S_gtk_assistant_set_page_header_image", object, page, pixbuf, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAssistantGetPageHeaderImage <-
function(object, page)
{
  checkPtrType(object, "GtkAssistant")
  checkPtrType(page, "GtkWidget")

  w <- .RGtkCall("S_gtk_assistant_get_page_header_image", object, page, PACKAGE = "RGtk2")

  return(w)
} 


gtkAssistantSetPageSideImage <-
function(object, page, pixbuf = NULL)
{
  checkPtrType(object, "GtkAssistant")
  checkPtrType(page, "GtkWidget")
  if (!is.null( pixbuf )) checkPtrType(pixbuf, "GdkPixbuf")

  w <- .RGtkCall("S_gtk_assistant_set_page_side_image", object, page, pixbuf, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAssistantGetPageSideImage <-
function(object, page)
{
  checkPtrType(object, "GtkAssistant")
  checkPtrType(page, "GtkWidget")

  w <- .RGtkCall("S_gtk_assistant_get_page_side_image", object, page, PACKAGE = "RGtk2")

  return(w)
} 


gtkAssistantSetPageComplete <-
function(object, page, complete)
{
  checkPtrType(object, "GtkAssistant")
  checkPtrType(page, "GtkWidget")
  complete <- as.logical(complete)

  w <- .RGtkCall("S_gtk_assistant_set_page_complete", object, page, complete, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAssistantGetPageComplete <-
function(object, page)
{
  checkPtrType(object, "GtkAssistant")
  checkPtrType(page, "GtkWidget")

  w <- .RGtkCall("S_gtk_assistant_get_page_complete", object, page, PACKAGE = "RGtk2")

  return(w)
} 


gtkAssistantAddActionWidget <-
function(object, child)
{
  checkPtrType(object, "GtkAssistant")
  checkPtrType(child, "GtkWidget")

  w <- .RGtkCall("S_gtk_assistant_add_action_widget", object, child, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAssistantRemoveActionWidget <-
function(object, child)
{
  checkPtrType(object, "GtkAssistant")
  checkPtrType(child, "GtkWidget")

  w <- .RGtkCall("S_gtk_assistant_remove_action_widget", object, child, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAssistantUpdateButtonsState <-
function(object)
{
  checkPtrType(object, "GtkAssistant")

  w <- .RGtkCall("S_gtk_assistant_update_buttons_state", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkButtonSetImagePosition <-
function(object, position)
{
  checkPtrType(object, "GtkButton")
  

  w <- .RGtkCall("S_gtk_button_set_image_position", object, position, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkButtonGetImagePosition <-
function(object)
{
  checkPtrType(object, "GtkButton")

  w <- .RGtkCall("S_gtk_button_get_image_position", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkCellRendererAccelGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_cell_renderer_accel_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkCellRendererAccelNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_cell_renderer_accel_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkCellRendererSpinGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_cell_renderer_spin_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkCellRendererSpinNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_cell_renderer_spin_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkClipboardRequestRichText <-
function(object, buffer, callback, user.data)
{
  checkPtrType(object, "GtkClipboard")
  checkPtrType(buffer, "GtkTextBuffer")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gtk_clipboard_request_rich_text", object, buffer, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkClipboardWaitForRichText <-
function(object, buffer)
{
  checkPtrType(object, "GtkClipboard")
  checkPtrType(buffer, "GtkTextBuffer")

  w <- .RGtkCall("S_gtk_clipboard_wait_for_rich_text", object, buffer, PACKAGE = "RGtk2")

  return(w)
} 


gtkClipboardWaitIsRichTextAvailable <-
function(object, buffer)
{
  checkPtrType(object, "GtkClipboard")
  checkPtrType(buffer, "GtkTextBuffer")

  w <- .RGtkCall("S_gtk_clipboard_wait_is_rich_text_available", object, buffer, PACKAGE = "RGtk2")

  return(w)
} 


gtkComboBoxGetTitle <-
function(object)
{
  checkPtrType(object, "GtkComboBox")

  w <- .RGtkCall("S_gtk_combo_box_get_title", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkComboBoxSetTitle <-
function(object, title)
{
  checkPtrType(object, "GtkComboBox")
  title <- as.character(title)

  w <- .RGtkCall("S_gtk_combo_box_set_title", object, title, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDragDestSetTrackMotion <-
function(object, track.motion)
{
  checkPtrType(object, "GtkWidget")
  track.motion <- as.logical(track.motion)

  w <- .RGtkCall("S_gtk_drag_dest_set_track_motion", object, track.motion, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDragDestGetTrackMotion <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_drag_dest_get_track_motion", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntrySetInnerBorder <-
function(object, border = NULL)
{
  checkPtrType(object, "GtkEntry")
  if (!is.null( border )) checkPtrType(border, "GtkBorder")

  w <- .RGtkCall("S_gtk_entry_set_inner_border", object, border, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryGetInnerBorder <-
function(object)
{
  checkPtrType(object, "GtkEntry")

  w <- .RGtkCall("S_gtk_entry_get_inner_border", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserButtonGetFocusOnClick <-
function(object)
{
  checkPtrType(object, "GtkFileChooserButton")

  w <- .RGtkCall("S_gtk_file_chooser_button_get_focus_on_click", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserButtonSetFocusOnClick <-
function(object, focus.on.click)
{
  checkPtrType(object, "GtkFileChooserButton")
  focus.on.click <- as.logical(focus.on.click)

  w <- .RGtkCall("S_gtk_file_chooser_button_set_focus_on_click", object, focus.on.click, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkLabelGetLineWrapMode <-
function(object)
{
  checkPtrType(object, "GtkLabel")

  w <- .RGtkCall("S_gtk_label_get_line_wrap_mode", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkLabelSetLineWrapMode <-
function(object, wrap.mode)
{
  checkPtrType(object, "GtkLabel")
  

  w <- .RGtkCall("S_gtk_label_set_line_wrap_mode", object, wrap.mode, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkLinkButtonGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_link_button_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkLinkButtonNew <-
function(uri)
{
  uri <- as.character(uri)

  w <- .RGtkCall("S_gtk_link_button_new", uri, PACKAGE = "RGtk2")

  return(w)
} 


gtkLinkButtonNewWithLabel <-
function(uri, label = NULL, show = TRUE)
{
  uri <- as.character(uri)
  if (!is.null( label )) label <- as.character(label)

  w <- .RGtkCall("S_gtk_link_button_new_with_label", uri, label, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkLinkButtonGetUri <-
function(object)
{
  checkPtrType(object, "GtkLinkButton")

  w <- .RGtkCall("S_gtk_link_button_get_uri", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkLinkButtonSetUri <-
function(object, uri)
{
  checkPtrType(object, "GtkLinkButton")
  uri <- as.character(uri)

  w <- .RGtkCall("S_gtk_link_button_set_uri", object, uri, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkLinkButtonSetUriHook <-
function(func, data)
{
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_link_button_set_uri_hook", func, data, PACKAGE = "RGtk2")

  return(w)
} 


gtkMessageDialogSetImage <-
function(object, image)
{
  checkPtrType(object, "GtkMessageDialog")
  checkPtrType(image, "GtkWidget")

  w <- .RGtkCall("S_gtk_message_dialog_set_image", object, image, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkNotebookSetWindowCreationHook <-
function(func, data)
{
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_notebook_set_window_creation_hook", func, data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkNotebookSetGroupId <-
function(object, group.id)
{
  if(getOption("depwarn"))
    .Deprecated("setGroup", "RGtk2")

  checkPtrType(object, "GtkNotebook")
  group.id <- as.integer(group.id)

  w <- .RGtkCall("S_gtk_notebook_set_group_id", object, group.id, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkNotebookGetGroupId <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("getGroup", "RGtk2")

  checkPtrType(object, "GtkNotebook")

  w <- .RGtkCall("S_gtk_notebook_get_group_id", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkNotebookGetTabReorderable <-
function(object, child)
{
  checkPtrType(object, "GtkNotebook")
  checkPtrType(child, "GtkWidget")

  w <- .RGtkCall("S_gtk_notebook_get_tab_reorderable", object, child, PACKAGE = "RGtk2")

  return(w)
} 


gtkNotebookSetTabReorderable <-
function(object, child, reorderable)
{
  checkPtrType(object, "GtkNotebook")
  checkPtrType(child, "GtkWidget")
  reorderable <- as.logical(reorderable)

  w <- .RGtkCall("S_gtk_notebook_set_tab_reorderable", object, child, reorderable, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkNotebookGetTabDetachable <-
function(object, child)
{
  checkPtrType(object, "GtkNotebook")
  checkPtrType(child, "GtkWidget")

  w <- .RGtkCall("S_gtk_notebook_get_tab_detachable", object, child, PACKAGE = "RGtk2")

  return(w)
} 


gtkNotebookSetTabDetachable <-
function(object, child, detachable)
{
  checkPtrType(object, "GtkNotebook")
  checkPtrType(child, "GtkWidget")
  detachable <- as.logical(detachable)

  w <- .RGtkCall("S_gtk_notebook_set_tab_detachable", object, child, detachable, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPageSetupGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_page_setup_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkPageSetupNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_page_setup_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkPageSetupCopy <-
function(object)
{
  checkPtrType(object, "GtkPageSetup")

  w <- .RGtkCall("S_gtk_page_setup_copy", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPageSetupGetOrientation <-
function(object)
{
  checkPtrType(object, "GtkPageSetup")

  w <- .RGtkCall("S_gtk_page_setup_get_orientation", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPageSetupSetOrientation <-
function(object, orientation)
{
  checkPtrType(object, "GtkPageSetup")
  

  w <- .RGtkCall("S_gtk_page_setup_set_orientation", object, orientation, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPageSetupGetPaperSize <-
function(object)
{
  checkPtrType(object, "GtkPageSetup")

  w <- .RGtkCall("S_gtk_page_setup_get_paper_size", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPageSetupSetPaperSize <-
function(object, size)
{
  checkPtrType(object, "GtkPageSetup")
  checkPtrType(size, "GtkPaperSize")

  w <- .RGtkCall("S_gtk_page_setup_set_paper_size", object, size, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPageSetupGetTopMargin <-
function(object, unit)
{
  checkPtrType(object, "GtkPageSetup")
  

  w <- .RGtkCall("S_gtk_page_setup_get_top_margin", object, unit, PACKAGE = "RGtk2")

  return(w)
} 


gtkPageSetupSetTopMargin <-
function(object, margin, unit)
{
  checkPtrType(object, "GtkPageSetup")
  margin <- as.numeric(margin)
  

  w <- .RGtkCall("S_gtk_page_setup_set_top_margin", object, margin, unit, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPageSetupGetBottomMargin <-
function(object, unit)
{
  checkPtrType(object, "GtkPageSetup")
  

  w <- .RGtkCall("S_gtk_page_setup_get_bottom_margin", object, unit, PACKAGE = "RGtk2")

  return(w)
} 


gtkPageSetupSetBottomMargin <-
function(object, margin, unit)
{
  checkPtrType(object, "GtkPageSetup")
  margin <- as.numeric(margin)
  

  w <- .RGtkCall("S_gtk_page_setup_set_bottom_margin", object, margin, unit, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPageSetupGetLeftMargin <-
function(object, unit)
{
  checkPtrType(object, "GtkPageSetup")
  

  w <- .RGtkCall("S_gtk_page_setup_get_left_margin", object, unit, PACKAGE = "RGtk2")

  return(w)
} 


gtkPageSetupSetLeftMargin <-
function(object, margin, unit)
{
  checkPtrType(object, "GtkPageSetup")
  margin <- as.numeric(margin)
  

  w <- .RGtkCall("S_gtk_page_setup_set_left_margin", object, margin, unit, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPageSetupGetRightMargin <-
function(object, unit)
{
  checkPtrType(object, "GtkPageSetup")
  

  w <- .RGtkCall("S_gtk_page_setup_get_right_margin", object, unit, PACKAGE = "RGtk2")

  return(w)
} 


gtkPageSetupSetRightMargin <-
function(object, margin, unit)
{
  checkPtrType(object, "GtkPageSetup")
  margin <- as.numeric(margin)
  

  w <- .RGtkCall("S_gtk_page_setup_set_right_margin", object, margin, unit, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPageSetupSetPaperSizeAndDefaultMargins <-
function(object, size)
{
  checkPtrType(object, "GtkPageSetup")
  checkPtrType(size, "GtkPaperSize")

  w <- .RGtkCall("S_gtk_page_setup_set_paper_size_and_default_margins", object, size, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPageSetupGetPaperWidth <-
function(object, unit)
{
  checkPtrType(object, "GtkPageSetup")
  

  w <- .RGtkCall("S_gtk_page_setup_get_paper_width", object, unit, PACKAGE = "RGtk2")

  return(w)
} 


gtkPageSetupGetPaperHeight <-
function(object, unit)
{
  checkPtrType(object, "GtkPageSetup")
  

  w <- .RGtkCall("S_gtk_page_setup_get_paper_height", object, unit, PACKAGE = "RGtk2")

  return(w)
} 


gtkPageSetupGetPageWidth <-
function(object, unit)
{
  checkPtrType(object, "GtkPageSetup")
  

  w <- .RGtkCall("S_gtk_page_setup_get_page_width", object, unit, PACKAGE = "RGtk2")

  return(w)
} 


gtkPageSetupGetPageHeight <-
function(object, unit)
{
  checkPtrType(object, "GtkPageSetup")
  

  w <- .RGtkCall("S_gtk_page_setup_get_page_height", object, unit, PACKAGE = "RGtk2")

  return(w)
} 


gtkPaperSizeGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_paper_size_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkPaperSizeNew <-
function(name = NULL)
{
  if (!is.null( name )) name <- as.character(name)

  w <- .RGtkCall("S_gtk_paper_size_new", name, PACKAGE = "RGtk2")

  return(w)
} 


gtkPaperSizeNewFromPpd <-
function(ppd.name, ppd.display.name, width, height)
{
  ppd.name <- as.character(ppd.name)
  ppd.display.name <- as.character(ppd.display.name)
  width <- as.numeric(width)
  height <- as.numeric(height)

  w <- .RGtkCall("S_gtk_paper_size_new_from_ppd", ppd.name, ppd.display.name, width, height, PACKAGE = "RGtk2")

  return(w)
} 


gtkPaperSizeNewCustom <-
function(name, display.name, width, height, unit)
{
  name <- as.character(name)
  display.name <- as.character(display.name)
  width <- as.numeric(width)
  height <- as.numeric(height)
  

  w <- .RGtkCall("S_gtk_paper_size_new_custom", name, display.name, width, height, unit, PACKAGE = "RGtk2")

  return(w)
} 


gtkPaperSizeCopy <-
function(object)
{
  checkPtrType(object, "GtkPaperSize")

  w <- .RGtkCall("S_gtk_paper_size_copy", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPaperSizeIsEqual <-
function(object, size2)
{
  checkPtrType(object, "GtkPaperSize")
  checkPtrType(size2, "GtkPaperSize")

  w <- .RGtkCall("S_gtk_paper_size_is_equal", object, size2, PACKAGE = "RGtk2")

  return(w)
} 


gtkPaperSizeGetName <-
function(object)
{
  checkPtrType(object, "GtkPaperSize")

  w <- .RGtkCall("S_gtk_paper_size_get_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPaperSizeGetDisplayName <-
function(object)
{
  checkPtrType(object, "GtkPaperSize")

  w <- .RGtkCall("S_gtk_paper_size_get_display_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPaperSizeGetPpdName <-
function(object)
{
  checkPtrType(object, "GtkPaperSize")

  w <- .RGtkCall("S_gtk_paper_size_get_ppd_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPaperSizeGetWidth <-
function(object, unit)
{
  checkPtrType(object, "GtkPaperSize")
  

  w <- .RGtkCall("S_gtk_paper_size_get_width", object, unit, PACKAGE = "RGtk2")

  return(w)
} 


gtkPaperSizeGetHeight <-
function(object, unit)
{
  checkPtrType(object, "GtkPaperSize")
  

  w <- .RGtkCall("S_gtk_paper_size_get_height", object, unit, PACKAGE = "RGtk2")

  return(w)
} 


gtkPaperSizeIsCustom <-
function(object)
{
  checkPtrType(object, "GtkPaperSize")

  w <- .RGtkCall("S_gtk_paper_size_is_custom", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPaperSizeSetSize <-
function(object, width, height, unit)
{
  checkPtrType(object, "GtkPaperSize")
  width <- as.numeric(width)
  height <- as.numeric(height)
  

  w <- .RGtkCall("S_gtk_paper_size_set_size", object, width, height, unit, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPaperSizeGetDefaultTopMargin <-
function(object, unit)
{
  checkPtrType(object, "GtkPaperSize")
  

  w <- .RGtkCall("S_gtk_paper_size_get_default_top_margin", object, unit, PACKAGE = "RGtk2")

  return(w)
} 


gtkPaperSizeGetDefaultBottomMargin <-
function(object, unit)
{
  checkPtrType(object, "GtkPaperSize")
  

  w <- .RGtkCall("S_gtk_paper_size_get_default_bottom_margin", object, unit, PACKAGE = "RGtk2")

  return(w)
} 


gtkPaperSizeGetDefaultLeftMargin <-
function(object, unit)
{
  checkPtrType(object, "GtkPaperSize")
  

  w <- .RGtkCall("S_gtk_paper_size_get_default_left_margin", object, unit, PACKAGE = "RGtk2")

  return(w)
} 


gtkPaperSizeGetDefaultRightMargin <-
function(object, unit)
{
  checkPtrType(object, "GtkPaperSize")
  

  w <- .RGtkCall("S_gtk_paper_size_get_default_right_margin", object, unit, PACKAGE = "RGtk2")

  return(w)
} 


gtkPaperSizeGetDefault <-
function()
{
  

  w <- .RGtkCall("S_gtk_paper_size_get_default", PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintContextGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_print_context_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintContextGetCairoContext <-
function(object)
{
  checkPtrType(object, "GtkPrintContext")

  w <- .RGtkCall("S_gtk_print_context_get_cairo_context", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintContextGetPageSetup <-
function(object)
{
  checkPtrType(object, "GtkPrintContext")

  w <- .RGtkCall("S_gtk_print_context_get_page_setup", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintContextGetWidth <-
function(object)
{
  checkPtrType(object, "GtkPrintContext")

  w <- .RGtkCall("S_gtk_print_context_get_width", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintContextGetHeight <-
function(object)
{
  checkPtrType(object, "GtkPrintContext")

  w <- .RGtkCall("S_gtk_print_context_get_height", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintContextGetDpiX <-
function(object)
{
  checkPtrType(object, "GtkPrintContext")

  w <- .RGtkCall("S_gtk_print_context_get_dpi_x", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintContextGetDpiY <-
function(object)
{
  checkPtrType(object, "GtkPrintContext")

  w <- .RGtkCall("S_gtk_print_context_get_dpi_y", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintContextGetPangoFontmap <-
function(object)
{
  checkPtrType(object, "GtkPrintContext")

  w <- .RGtkCall("S_gtk_print_context_get_pango_fontmap", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintContextCreatePangoContext <-
function(object)
{
  checkPtrType(object, "GtkPrintContext")

  w <- .RGtkCall("S_gtk_print_context_create_pango_context", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintContextCreatePangoLayout <-
function(object)
{
  checkPtrType(object, "GtkPrintContext")

  w <- .RGtkCall("S_gtk_print_context_create_pango_layout", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintContextSetCairoContext <-
function(object, cr, dpi.x, dpi.y)
{
  checkPtrType(object, "GtkPrintContext")
  checkPtrType(cr, "Cairo")
  dpi.x <- as.numeric(dpi.x)
  dpi.y <- as.numeric(dpi.y)

  w <- .RGtkCall("S_gtk_print_context_set_cairo_context", object, cr, dpi.x, dpi.y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintErrorQuark <-
function()
{
  

  w <- .RGtkCall("S_gtk_print_error_quark", PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintOperationGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_print_operation_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintOperationNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_print_operation_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintOperationSetDefaultPageSetup <-
function(object, default.page.setup = NULL)
{
  checkPtrType(object, "GtkPrintOperation")
  if (!is.null( default.page.setup )) checkPtrType(default.page.setup, "GtkPageSetup")

  w <- .RGtkCall("S_gtk_print_operation_set_default_page_setup", object, default.page.setup, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintOperationGetDefaultPageSetup <-
function(object)
{
  checkPtrType(object, "GtkPrintOperation")

  w <- .RGtkCall("S_gtk_print_operation_get_default_page_setup", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintOperationSetPrintSettings <-
function(object, print.settings = NULL)
{
  checkPtrType(object, "GtkPrintOperation")
  if (!is.null( print.settings )) checkPtrType(print.settings, "GtkPrintSettings")

  w <- .RGtkCall("S_gtk_print_operation_set_print_settings", object, print.settings, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintOperationGetPrintSettings <-
function(object)
{
  checkPtrType(object, "GtkPrintOperation")

  w <- .RGtkCall("S_gtk_print_operation_get_print_settings", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintOperationSetJobName <-
function(object, job.name)
{
  checkPtrType(object, "GtkPrintOperation")
  job.name <- as.character(job.name)

  w <- .RGtkCall("S_gtk_print_operation_set_job_name", object, job.name, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintOperationSetNPages <-
function(object, n.pages)
{
  checkPtrType(object, "GtkPrintOperation")
  n.pages <- as.integer(n.pages)

  w <- .RGtkCall("S_gtk_print_operation_set_n_pages", object, n.pages, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintOperationSetCurrentPage <-
function(object, current.page)
{
  checkPtrType(object, "GtkPrintOperation")
  current.page <- as.integer(current.page)

  w <- .RGtkCall("S_gtk_print_operation_set_current_page", object, current.page, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintOperationSetUseFullPage <-
function(object, full.page)
{
  checkPtrType(object, "GtkPrintOperation")
  full.page <- as.logical(full.page)

  w <- .RGtkCall("S_gtk_print_operation_set_use_full_page", object, full.page, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintOperationSetUnit <-
function(object, unit)
{
  checkPtrType(object, "GtkPrintOperation")
  

  w <- .RGtkCall("S_gtk_print_operation_set_unit", object, unit, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintOperationSetExportFilename <-
function(object, filename)
{
  checkPtrType(object, "GtkPrintOperation")
  filename <- as.character(filename)

  w <- .RGtkCall("S_gtk_print_operation_set_export_filename", object, filename, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintOperationSetTrackPrintStatus <-
function(object, track.status)
{
  checkPtrType(object, "GtkPrintOperation")
  track.status <- as.logical(track.status)

  w <- .RGtkCall("S_gtk_print_operation_set_track_print_status", object, track.status, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintOperationSetShowProgress <-
function(object, show.progress)
{
  checkPtrType(object, "GtkPrintOperation")
  show.progress <- as.logical(show.progress)

  w <- .RGtkCall("S_gtk_print_operation_set_show_progress", object, show.progress, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintOperationSetAllowAsync <-
function(object, allow.async)
{
  checkPtrType(object, "GtkPrintOperation")
  allow.async <- as.logical(allow.async)

  w <- .RGtkCall("S_gtk_print_operation_set_allow_async", object, allow.async, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintOperationSetCustomTabLabel <-
function(object, label)
{
  checkPtrType(object, "GtkPrintOperation")
  label <- as.character(label)

  w <- .RGtkCall("S_gtk_print_operation_set_custom_tab_label", object, label, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintOperationRun <-
function(object, action, parent = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GtkPrintOperation")
  
  if (!is.null( parent )) checkPtrType(parent, "GtkWindow")

  w <- .RGtkCall("S_gtk_print_operation_run", object, action, parent, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkPrintOperationGetError <-
function(object, .errwarn = TRUE)
{
  checkPtrType(object, "GtkPrintOperation")

  w <- .RGtkCall("S_gtk_print_operation_get_error", object, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkPrintOperationGetStatus <-
function(object)
{
  checkPtrType(object, "GtkPrintOperation")

  w <- .RGtkCall("S_gtk_print_operation_get_status", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintOperationGetStatusString <-
function(object)
{
  checkPtrType(object, "GtkPrintOperation")

  w <- .RGtkCall("S_gtk_print_operation_get_status_string", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintOperationIsFinished <-
function(object)
{
  checkPtrType(object, "GtkPrintOperation")

  w <- .RGtkCall("S_gtk_print_operation_is_finished", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintOperationCancel <-
function(object)
{
  checkPtrType(object, "GtkPrintOperation")

  w <- .RGtkCall("S_gtk_print_operation_cancel", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintRunPageSetupDialog <-
function(parent, page.setup = NULL, settings)
{
  checkPtrType(parent, "GtkWindow")
  if (!is.null( page.setup )) checkPtrType(page.setup, "GtkPageSetup")
  checkPtrType(settings, "GtkPrintSettings")

  w <- .RGtkCall("S_gtk_print_run_page_setup_dialog", parent, page.setup, settings, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintRunPageSetupDialogAsync <-
function(parent, page.setup, settings, done.cb, data)
{
  checkPtrType(parent, "GtkWindow")
  checkPtrType(page.setup, "GtkPageSetup")
  checkPtrType(settings, "GtkPrintSettings")
  done.cb <- as.function(done.cb)
  

  w <- .RGtkCall("S_gtk_print_run_page_setup_dialog_async", parent, page.setup, settings, done.cb, data, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintOperationPreviewGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_print_operation_preview_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintOperationPreviewRenderPage <-
function(object, page.nr)
{
  checkPtrType(object, "GtkPrintOperationPreview")
  page.nr <- as.integer(page.nr)

  w <- .RGtkCall("S_gtk_print_operation_preview_render_page", object, page.nr, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintOperationPreviewEndPreview <-
function(object)
{
  checkPtrType(object, "GtkPrintOperationPreview")

  w <- .RGtkCall("S_gtk_print_operation_preview_end_preview", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintOperationPreviewIsSelected <-
function(object, page.nr)
{
  checkPtrType(object, "GtkPrintOperationPreview")
  page.nr <- as.integer(page.nr)

  w <- .RGtkCall("S_gtk_print_operation_preview_is_selected", object, page.nr, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_print_settings_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_print_settings_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsCopy <-
function(object)
{
  checkPtrType(object, "GtkPrintSettings")

  w <- .RGtkCall("S_gtk_print_settings_copy", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsHasKey <-
function(object, key)
{
  checkPtrType(object, "GtkPrintSettings")
  key <- as.character(key)

  w <- .RGtkCall("S_gtk_print_settings_has_key", object, key, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsGet <-
function(object, key)
{
  checkPtrType(object, "GtkPrintSettings")
  key <- as.character(key)

  w <- .RGtkCall("S_gtk_print_settings_get", object, key, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsSet <-
function(object, key, value)
{
  checkPtrType(object, "GtkPrintSettings")
  key <- as.character(key)
  value <- as.character(value)

  w <- .RGtkCall("S_gtk_print_settings_set", object, key, value, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintSettingsUnset <-
function(object, key)
{
  checkPtrType(object, "GtkPrintSettings")
  key <- as.character(key)

  w <- .RGtkCall("S_gtk_print_settings_unset", object, key, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintSettingsForeach <-
function(object, func, user.data = NULL)
{
  checkPtrType(object, "GtkPrintSettings")
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_print_settings_foreach", object, func, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintSettingsGetBool <-
function(object, key)
{
  checkPtrType(object, "GtkPrintSettings")
  key <- as.character(key)

  w <- .RGtkCall("S_gtk_print_settings_get_bool", object, key, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsSetBool <-
function(object, key, value)
{
  checkPtrType(object, "GtkPrintSettings")
  key <- as.character(key)
  value <- as.logical(value)

  w <- .RGtkCall("S_gtk_print_settings_set_bool", object, key, value, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintSettingsGetDouble <-
function(object, key)
{
  checkPtrType(object, "GtkPrintSettings")
  key <- as.character(key)

  w <- .RGtkCall("S_gtk_print_settings_get_double", object, key, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsGetDoubleWithDefault <-
function(object, key, def)
{
  checkPtrType(object, "GtkPrintSettings")
  key <- as.character(key)
  def <- as.numeric(def)

  w <- .RGtkCall("S_gtk_print_settings_get_double_with_default", object, key, def, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsSetDouble <-
function(object, key, value)
{
  checkPtrType(object, "GtkPrintSettings")
  key <- as.character(key)
  value <- as.numeric(value)

  w <- .RGtkCall("S_gtk_print_settings_set_double", object, key, value, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintSettingsGetLength <-
function(object, key, unit)
{
  checkPtrType(object, "GtkPrintSettings")
  key <- as.character(key)
  

  w <- .RGtkCall("S_gtk_print_settings_get_length", object, key, unit, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsSetLength <-
function(object, key, value, unit)
{
  checkPtrType(object, "GtkPrintSettings")
  key <- as.character(key)
  value <- as.numeric(value)
  

  w <- .RGtkCall("S_gtk_print_settings_set_length", object, key, value, unit, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintSettingsGetInt <-
function(object, key)
{
  checkPtrType(object, "GtkPrintSettings")
  key <- as.character(key)

  w <- .RGtkCall("S_gtk_print_settings_get_int", object, key, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsGetIntWithDefault <-
function(object, key, def)
{
  checkPtrType(object, "GtkPrintSettings")
  key <- as.character(key)
  def <- as.integer(def)

  w <- .RGtkCall("S_gtk_print_settings_get_int_with_default", object, key, def, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsSetInt <-
function(object, key, value)
{
  checkPtrType(object, "GtkPrintSettings")
  key <- as.character(key)
  value <- as.integer(value)

  w <- .RGtkCall("S_gtk_print_settings_set_int", object, key, value, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintSettingsGetPrinter <-
function(object)
{
  checkPtrType(object, "GtkPrintSettings")

  w <- .RGtkCall("S_gtk_print_settings_get_printer", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsSetPrinter <-
function(object, printer)
{
  checkPtrType(object, "GtkPrintSettings")
  printer <- as.character(printer)

  w <- .RGtkCall("S_gtk_print_settings_set_printer", object, printer, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintSettingsGetOrientation <-
function(object)
{
  checkPtrType(object, "GtkPrintSettings")

  w <- .RGtkCall("S_gtk_print_settings_get_orientation", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsSetOrientation <-
function(object, orientation)
{
  checkPtrType(object, "GtkPrintSettings")
  

  w <- .RGtkCall("S_gtk_print_settings_set_orientation", object, orientation, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintSettingsGetPaperSize <-
function(object)
{
  checkPtrType(object, "GtkPrintSettings")

  w <- .RGtkCall("S_gtk_print_settings_get_paper_size", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsSetPaperSize <-
function(object, paper.size)
{
  checkPtrType(object, "GtkPrintSettings")
  checkPtrType(paper.size, "GtkPaperSize")

  w <- .RGtkCall("S_gtk_print_settings_set_paper_size", object, paper.size, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintSettingsGetPaperWidth <-
function(object, unit)
{
  checkPtrType(object, "GtkPrintSettings")
  

  w <- .RGtkCall("S_gtk_print_settings_get_paper_width", object, unit, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsSetPaperWidth <-
function(object, width, unit)
{
  checkPtrType(object, "GtkPrintSettings")
  width <- as.numeric(width)
  

  w <- .RGtkCall("S_gtk_print_settings_set_paper_width", object, width, unit, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintSettingsGetPaperHeight <-
function(object, unit)
{
  checkPtrType(object, "GtkPrintSettings")
  

  w <- .RGtkCall("S_gtk_print_settings_get_paper_height", object, unit, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsSetPaperHeight <-
function(object, height, unit)
{
  checkPtrType(object, "GtkPrintSettings")
  height <- as.numeric(height)
  

  w <- .RGtkCall("S_gtk_print_settings_set_paper_height", object, height, unit, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintSettingsGetUseColor <-
function(object)
{
  checkPtrType(object, "GtkPrintSettings")

  w <- .RGtkCall("S_gtk_print_settings_get_use_color", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsSetUseColor <-
function(object, use.color)
{
  checkPtrType(object, "GtkPrintSettings")
  use.color <- as.logical(use.color)

  w <- .RGtkCall("S_gtk_print_settings_set_use_color", object, use.color, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintSettingsGetCollate <-
function(object)
{
  checkPtrType(object, "GtkPrintSettings")

  w <- .RGtkCall("S_gtk_print_settings_get_collate", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsSetCollate <-
function(object, collate)
{
  checkPtrType(object, "GtkPrintSettings")
  collate <- as.logical(collate)

  w <- .RGtkCall("S_gtk_print_settings_set_collate", object, collate, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintSettingsGetReverse <-
function(object)
{
  checkPtrType(object, "GtkPrintSettings")

  w <- .RGtkCall("S_gtk_print_settings_get_reverse", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsSetReverse <-
function(object, reverse)
{
  checkPtrType(object, "GtkPrintSettings")
  reverse <- as.logical(reverse)

  w <- .RGtkCall("S_gtk_print_settings_set_reverse", object, reverse, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintSettingsGetDuplex <-
function(object)
{
  checkPtrType(object, "GtkPrintSettings")

  w <- .RGtkCall("S_gtk_print_settings_get_duplex", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsSetDuplex <-
function(object, duplex)
{
  checkPtrType(object, "GtkPrintSettings")
  

  w <- .RGtkCall("S_gtk_print_settings_set_duplex", object, duplex, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintSettingsGetQuality <-
function(object)
{
  checkPtrType(object, "GtkPrintSettings")

  w <- .RGtkCall("S_gtk_print_settings_get_quality", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsSetQuality <-
function(object, quality)
{
  checkPtrType(object, "GtkPrintSettings")
  

  w <- .RGtkCall("S_gtk_print_settings_set_quality", object, quality, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintSettingsGetNCopies <-
function(object)
{
  checkPtrType(object, "GtkPrintSettings")

  w <- .RGtkCall("S_gtk_print_settings_get_n_copies", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsSetNCopies <-
function(object, num.copies)
{
  checkPtrType(object, "GtkPrintSettings")
  num.copies <- as.integer(num.copies)

  w <- .RGtkCall("S_gtk_print_settings_set_n_copies", object, num.copies, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintSettingsGetNumberUp <-
function(object)
{
  checkPtrType(object, "GtkPrintSettings")

  w <- .RGtkCall("S_gtk_print_settings_get_number_up", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsSetNumberUp <-
function(object, number.up)
{
  checkPtrType(object, "GtkPrintSettings")
  number.up <- as.integer(number.up)

  w <- .RGtkCall("S_gtk_print_settings_set_number_up", object, number.up, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintSettingsGetResolution <-
function(object)
{
  checkPtrType(object, "GtkPrintSettings")

  w <- .RGtkCall("S_gtk_print_settings_get_resolution", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsSetResolution <-
function(object, resolution)
{
  checkPtrType(object, "GtkPrintSettings")
  resolution <- as.integer(resolution)

  w <- .RGtkCall("S_gtk_print_settings_set_resolution", object, resolution, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintSettingsGetScale <-
function(object)
{
  checkPtrType(object, "GtkPrintSettings")

  w <- .RGtkCall("S_gtk_print_settings_get_scale", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsSetScale <-
function(object, scale)
{
  checkPtrType(object, "GtkPrintSettings")
  scale <- as.numeric(scale)

  w <- .RGtkCall("S_gtk_print_settings_set_scale", object, scale, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintSettingsGetPrintPages <-
function(object)
{
  checkPtrType(object, "GtkPrintSettings")

  w <- .RGtkCall("S_gtk_print_settings_get_print_pages", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsSetPrintPages <-
function(object, pages)
{
  checkPtrType(object, "GtkPrintSettings")
  

  w <- .RGtkCall("S_gtk_print_settings_set_print_pages", object, pages, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintSettingsGetPageRanges <-
function(object, num.ranges)
{
  checkPtrType(object, "GtkPrintSettings")
  num.ranges <- as.list(as.integer(num.ranges))

  w <- .RGtkCall("S_gtk_print_settings_get_page_ranges", object, num.ranges, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsSetPageRanges <-
function(object, page.ranges, num.ranges)
{
  checkPtrType(object, "GtkPrintSettings")
  page.ranges <- as.GtkPageRange(page.ranges)
  num.ranges <- as.integer(num.ranges)

  w <- .RGtkCall("S_gtk_print_settings_set_page_ranges", object, page.ranges, num.ranges, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintSettingsGetPageSet <-
function(object)
{
  checkPtrType(object, "GtkPrintSettings")

  w <- .RGtkCall("S_gtk_print_settings_get_page_set", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsSetPageSet <-
function(object, page.set)
{
  checkPtrType(object, "GtkPrintSettings")
  

  w <- .RGtkCall("S_gtk_print_settings_set_page_set", object, page.set, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintSettingsGetDefaultSource <-
function(object)
{
  checkPtrType(object, "GtkPrintSettings")

  w <- .RGtkCall("S_gtk_print_settings_get_default_source", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsSetDefaultSource <-
function(object, default.source)
{
  checkPtrType(object, "GtkPrintSettings")
  default.source <- as.character(default.source)

  w <- .RGtkCall("S_gtk_print_settings_set_default_source", object, default.source, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintSettingsGetMediaType <-
function(object)
{
  checkPtrType(object, "GtkPrintSettings")

  w <- .RGtkCall("S_gtk_print_settings_get_media_type", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsSetMediaType <-
function(object, media.type)
{
  checkPtrType(object, "GtkPrintSettings")
  media.type <- as.character(media.type)

  w <- .RGtkCall("S_gtk_print_settings_set_media_type", object, media.type, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintSettingsGetDither <-
function(object)
{
  checkPtrType(object, "GtkPrintSettings")

  w <- .RGtkCall("S_gtk_print_settings_get_dither", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsSetDither <-
function(object, dither)
{
  checkPtrType(object, "GtkPrintSettings")
  dither <- as.character(dither)

  w <- .RGtkCall("S_gtk_print_settings_set_dither", object, dither, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintSettingsGetFinishings <-
function(object)
{
  checkPtrType(object, "GtkPrintSettings")

  w <- .RGtkCall("S_gtk_print_settings_get_finishings", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsSetFinishings <-
function(object, finishings)
{
  checkPtrType(object, "GtkPrintSettings")
  finishings <- as.character(finishings)

  w <- .RGtkCall("S_gtk_print_settings_set_finishings", object, finishings, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintSettingsGetOutputBin <-
function(object)
{
  checkPtrType(object, "GtkPrintSettings")

  w <- .RGtkCall("S_gtk_print_settings_get_output_bin", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsSetOutputBin <-
function(object, output.bin)
{
  checkPtrType(object, "GtkPrintSettings")
  output.bin <- as.character(output.bin)

  w <- .RGtkCall("S_gtk_print_settings_set_output_bin", object, output.bin, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRadioActionSetCurrentValue <-
function(object, current.value)
{
  checkPtrType(object, "GtkRadioAction")
  current.value <- as.integer(current.value)

  w <- .RGtkCall("S_gtk_radio_action_set_current_value", object, current.value, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRangeSetLowerStepperSensitivity <-
function(object, sensitivity)
{
  checkPtrType(object, "GtkRange")
  

  w <- .RGtkCall("S_gtk_range_set_lower_stepper_sensitivity", object, sensitivity, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRangeGetLowerStepperSensitivity <-
function(object)
{
  checkPtrType(object, "GtkRange")

  w <- .RGtkCall("S_gtk_range_get_lower_stepper_sensitivity", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRangeSetUpperStepperSensitivity <-
function(object, sensitivity)
{
  checkPtrType(object, "GtkRange")
  

  w <- .RGtkCall("S_gtk_range_set_upper_stepper_sensitivity", object, sensitivity, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRangeGetUpperStepperSensitivity <-
function(object)
{
  checkPtrType(object, "GtkRange")

  w <- .RGtkCall("S_gtk_range_get_upper_stepper_sensitivity", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentChooserDialogGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_recent_chooser_dialog_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentChooserErrorQuark <-
function()
{
  

  w <- .RGtkCall("S_gtk_recent_chooser_error_quark", PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentChooserGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_recent_chooser_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentChooserSetShowPrivate <-
function(object, show.private)
{
  checkPtrType(object, "GtkRecentChooser")
  show.private <- as.logical(show.private)

  w <- .RGtkCall("S_gtk_recent_chooser_set_show_private", object, show.private, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRecentChooserGetShowPrivate <-
function(object)
{
  checkPtrType(object, "GtkRecentChooser")

  w <- .RGtkCall("S_gtk_recent_chooser_get_show_private", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentChooserSetShowNotFound <-
function(object, show.not.found)
{
  checkPtrType(object, "GtkRecentChooser")
  show.not.found <- as.logical(show.not.found)

  w <- .RGtkCall("S_gtk_recent_chooser_set_show_not_found", object, show.not.found, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRecentChooserGetShowNotFound <-
function(object)
{
  checkPtrType(object, "GtkRecentChooser")

  w <- .RGtkCall("S_gtk_recent_chooser_get_show_not_found", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentChooserSetSelectMultiple <-
function(object, select.multiple)
{
  checkPtrType(object, "GtkRecentChooser")
  select.multiple <- as.logical(select.multiple)

  w <- .RGtkCall("S_gtk_recent_chooser_set_select_multiple", object, select.multiple, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRecentChooserGetSelectMultiple <-
function(object)
{
  checkPtrType(object, "GtkRecentChooser")

  w <- .RGtkCall("S_gtk_recent_chooser_get_select_multiple", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentChooserSetLimit <-
function(object, limit)
{
  checkPtrType(object, "GtkRecentChooser")
  limit <- as.integer(limit)

  w <- .RGtkCall("S_gtk_recent_chooser_set_limit", object, limit, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRecentChooserGetLimit <-
function(object)
{
  checkPtrType(object, "GtkRecentChooser")

  w <- .RGtkCall("S_gtk_recent_chooser_get_limit", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentChooserSetLocalOnly <-
function(object, local.only)
{
  checkPtrType(object, "GtkRecentChooser")
  local.only <- as.logical(local.only)

  w <- .RGtkCall("S_gtk_recent_chooser_set_local_only", object, local.only, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRecentChooserGetLocalOnly <-
function(object)
{
  checkPtrType(object, "GtkRecentChooser")

  w <- .RGtkCall("S_gtk_recent_chooser_get_local_only", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentChooserSetShowTips <-
function(object, show.tips)
{
  checkPtrType(object, "GtkRecentChooser")
  show.tips <- as.logical(show.tips)

  w <- .RGtkCall("S_gtk_recent_chooser_set_show_tips", object, show.tips, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRecentChooserGetShowTips <-
function(object)
{
  checkPtrType(object, "GtkRecentChooser")

  w <- .RGtkCall("S_gtk_recent_chooser_get_show_tips", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentChooserSetShowIcons <-
function(object, show.icons)
{
  checkPtrType(object, "GtkRecentChooser")
  show.icons <- as.logical(show.icons)

  w <- .RGtkCall("S_gtk_recent_chooser_set_show_icons", object, show.icons, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRecentChooserGetShowIcons <-
function(object)
{
  checkPtrType(object, "GtkRecentChooser")

  w <- .RGtkCall("S_gtk_recent_chooser_get_show_icons", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentChooserSetSortType <-
function(object, sort.type)
{
  checkPtrType(object, "GtkRecentChooser")
  

  w <- .RGtkCall("S_gtk_recent_chooser_set_sort_type", object, sort.type, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRecentChooserGetSortType <-
function(object)
{
  checkPtrType(object, "GtkRecentChooser")

  w <- .RGtkCall("S_gtk_recent_chooser_get_sort_type", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentChooserSetSortFunc <-
function(object, sort.func, sort.data)
{
  checkPtrType(object, "GtkRecentChooser")
  sort.func <- as.function(sort.func)
  

  w <- .RGtkCall("S_gtk_recent_chooser_set_sort_func", object, sort.func, sort.data, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentChooserSetCurrentUri <-
function(object, uri, .errwarn = TRUE)
{
  checkPtrType(object, "GtkRecentChooser")
  uri <- as.character(uri)

  w <- .RGtkCall("S_gtk_recent_chooser_set_current_uri", object, uri, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkRecentChooserGetCurrentUri <-
function(object)
{
  checkPtrType(object, "GtkRecentChooser")

  w <- .RGtkCall("S_gtk_recent_chooser_get_current_uri", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentChooserGetCurrentItem <-
function(object)
{
  checkPtrType(object, "GtkRecentChooser")

  w <- .RGtkCall("S_gtk_recent_chooser_get_current_item", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentChooserSelectUri <-
function(object, uri, .errwarn = TRUE)
{
  checkPtrType(object, "GtkRecentChooser")
  uri <- as.character(uri)

  w <- .RGtkCall("S_gtk_recent_chooser_select_uri", object, uri, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkRecentChooserUnselectUri <-
function(object, uri)
{
  checkPtrType(object, "GtkRecentChooser")
  uri <- as.character(uri)

  w <- .RGtkCall("S_gtk_recent_chooser_unselect_uri", object, uri, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRecentChooserSelectAll <-
function(object)
{
  checkPtrType(object, "GtkRecentChooser")

  w <- .RGtkCall("S_gtk_recent_chooser_select_all", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRecentChooserUnselectAll <-
function(object)
{
  checkPtrType(object, "GtkRecentChooser")

  w <- .RGtkCall("S_gtk_recent_chooser_unselect_all", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRecentChooserGetItems <-
function(object)
{
  checkPtrType(object, "GtkRecentChooser")

  w <- .RGtkCall("S_gtk_recent_chooser_get_items", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentChooserGetUris <-
function(object)
{
  checkPtrType(object, "GtkRecentChooser")

  w <- .RGtkCall("S_gtk_recent_chooser_get_uris", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentChooserAddFilter <-
function(object, filter)
{
  checkPtrType(object, "GtkRecentChooser")
  checkPtrType(filter, "GtkRecentFilter")

  w <- .RGtkCall("S_gtk_recent_chooser_add_filter", object, filter, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRecentChooserRemoveFilter <-
function(object, filter)
{
  checkPtrType(object, "GtkRecentChooser")
  checkPtrType(filter, "GtkRecentFilter")

  w <- .RGtkCall("S_gtk_recent_chooser_remove_filter", object, filter, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRecentChooserListFilters <-
function(object)
{
  checkPtrType(object, "GtkRecentChooser")

  w <- .RGtkCall("S_gtk_recent_chooser_list_filters", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentChooserSetFilter <-
function(object, filter)
{
  checkPtrType(object, "GtkRecentChooser")
  checkPtrType(filter, "GtkRecentFilter")

  w <- .RGtkCall("S_gtk_recent_chooser_set_filter", object, filter, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRecentChooserGetFilter <-
function(object)
{
  checkPtrType(object, "GtkRecentChooser")

  w <- .RGtkCall("S_gtk_recent_chooser_get_filter", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentChooserMenuGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_recent_chooser_menu_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentChooserMenuNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_recent_chooser_menu_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentChooserMenuNewForManager <-
function(manager = NULL, show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_recent_chooser_menu_new_for_manager", manager, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkRecentChooserMenuGetShowNumbers <-
function(object)
{
  checkPtrType(object, "GtkRecentChooserMenu")

  w <- .RGtkCall("S_gtk_recent_chooser_menu_get_show_numbers", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentChooserMenuSetShowNumbers <-
function(object, show.numbers)
{
  checkPtrType(object, "GtkRecentChooserMenu")
  show.numbers <- as.logical(show.numbers)

  w <- .RGtkCall("S_gtk_recent_chooser_menu_set_show_numbers", object, show.numbers, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRecentChooserWidgetGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_recent_chooser_widget_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentChooserWidgetNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_recent_chooser_widget_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentChooserWidgetNewForManager <-
function(manager = NULL, show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_recent_chooser_widget_new_for_manager", manager, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkRecentFilterGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_recent_filter_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentFilterNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_recent_filter_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentFilterSetName <-
function(object, name)
{
  checkPtrType(object, "GtkRecentFilter")
  name <- as.character(name)

  w <- .RGtkCall("S_gtk_recent_filter_set_name", object, name, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRecentFilterGetName <-
function(object)
{
  checkPtrType(object, "GtkRecentFilter")

  w <- .RGtkCall("S_gtk_recent_filter_get_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentFilterAddMimeType <-
function(object, mime.type)
{
  checkPtrType(object, "GtkRecentFilter")
  mime.type <- as.character(mime.type)

  w <- .RGtkCall("S_gtk_recent_filter_add_mime_type", object, mime.type, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRecentFilterAddPattern <-
function(object, pattern)
{
  checkPtrType(object, "GtkRecentFilter")
  pattern <- as.character(pattern)

  w <- .RGtkCall("S_gtk_recent_filter_add_pattern", object, pattern, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRecentFilterAddPixbufFormats <-
function(object)
{
  checkPtrType(object, "GtkRecentFilter")

  w <- .RGtkCall("S_gtk_recent_filter_add_pixbuf_formats", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRecentFilterAddApplication <-
function(object, application)
{
  checkPtrType(object, "GtkRecentFilter")
  application <- as.character(application)

  w <- .RGtkCall("S_gtk_recent_filter_add_application", object, application, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRecentFilterAddGroup <-
function(object, group)
{
  checkPtrType(object, "GtkRecentFilter")
  group <- as.character(group)

  w <- .RGtkCall("S_gtk_recent_filter_add_group", object, group, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRecentFilterAddAge <-
function(object, days)
{
  checkPtrType(object, "GtkRecentFilter")
  days <- as.integer(days)

  w <- .RGtkCall("S_gtk_recent_filter_add_age", object, days, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRecentFilterAddCustom <-
function(object, needed, func, data)
{
  checkPtrType(object, "GtkRecentFilter")
  
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_recent_filter_add_custom", object, needed, func, data, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentFilterGetNeeded <-
function(object)
{
  checkPtrType(object, "GtkRecentFilter")

  w <- .RGtkCall("S_gtk_recent_filter_get_needed", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentFilterFilter <-
function(object, filter.info)
{
  checkPtrType(object, "GtkRecentFilter")
  filter.info <- as.GtkRecentFilterInfo(filter.info)

  w <- .RGtkCall("S_gtk_recent_filter_filter", object, filter.info, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentManagerErrorQuark <-
function()
{
  

  w <- .RGtkCall("S_gtk_recent_manager_error_quark", PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentManagerGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_recent_manager_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentManagerNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_recent_manager_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentManagerGetDefault <-
function()
{
  

  w <- .RGtkCall("S_gtk_recent_manager_get_default", PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentManagerGetForScreen <-
function(screen)
{
  if(getOption("depwarn"))
    .Deprecated("gtkRecentManagerGetDefault", "RGtk2")

  checkPtrType(screen, "GdkScreen")

  w <- .RGtkCall("S_gtk_recent_manager_get_for_screen", screen, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentManagerSetScreen <-
function(object, screen)
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  checkPtrType(object, "GtkRecentManager")
  checkPtrType(screen, "GdkScreen")

  w <- .RGtkCall("S_gtk_recent_manager_set_screen", object, screen, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRecentManagerAddItem <-
function(object, uri)
{
  checkPtrType(object, "GtkRecentManager")
  uri <- as.character(uri)

  w <- .RGtkCall("S_gtk_recent_manager_add_item", object, uri, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentManagerAddFull <-
function(object, uri, recent.data)
{
  checkPtrType(object, "GtkRecentManager")
  uri <- as.character(uri)
  recent.data <- as.GtkRecentData(recent.data)

  w <- .RGtkCall("S_gtk_recent_manager_add_full", object, uri, recent.data, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentManagerRemoveItem <-
function(object, uri, .errwarn = TRUE)
{
  checkPtrType(object, "GtkRecentManager")
  uri <- as.character(uri)

  w <- .RGtkCall("S_gtk_recent_manager_remove_item", object, uri, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkRecentManagerLookupItem <-
function(object, uri, .errwarn = TRUE)
{
  checkPtrType(object, "GtkRecentManager")
  uri <- as.character(uri)

  w <- .RGtkCall("S_gtk_recent_manager_lookup_item", object, uri, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkRecentManagerHasItem <-
function(object, uri)
{
  checkPtrType(object, "GtkRecentManager")
  uri <- as.character(uri)

  w <- .RGtkCall("S_gtk_recent_manager_has_item", object, uri, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentManagerMoveItem <-
function(object, uri, new.uri, .errwarn = TRUE)
{
  checkPtrType(object, "GtkRecentManager")
  uri <- as.character(uri)
  new.uri <- as.character(new.uri)

  w <- .RGtkCall("S_gtk_recent_manager_move_item", object, uri, new.uri, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkRecentManagerSetLimit <-
function(object, limit)
{
  checkPtrType(object, "GtkRecentManager")
  limit <- as.integer(limit)

  w <- .RGtkCall("S_gtk_recent_manager_set_limit", object, limit, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRecentManagerGetLimit <-
function(object)
{
  checkPtrType(object, "GtkRecentManager")

  w <- .RGtkCall("S_gtk_recent_manager_get_limit", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentManagerGetItems <-
function(object)
{
  checkPtrType(object, "GtkRecentManager")

  w <- .RGtkCall("S_gtk_recent_manager_get_items", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentManagerPurgeItems <-
function(object, .errwarn = TRUE)
{
  checkPtrType(object, "GtkRecentManager")

  w <- .RGtkCall("S_gtk_recent_manager_purge_items", object, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkRecentInfoGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_recent_info_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentInfoRef <-
function(object)
{
  checkPtrType(object, "GtkRecentInfo")

  w <- .RGtkCall("S_gtk_recent_info_ref", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentInfoUnref <-
function(object)
{
  checkPtrType(object, "GtkRecentInfo")

  w <- .RGtkCall("S_gtk_recent_info_unref", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRecentInfoGetUri <-
function(object)
{
  checkPtrType(object, "GtkRecentInfo")

  w <- .RGtkCall("S_gtk_recent_info_get_uri", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentInfoGetDisplayName <-
function(object)
{
  checkPtrType(object, "GtkRecentInfo")

  w <- .RGtkCall("S_gtk_recent_info_get_display_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentInfoGetDescription <-
function(object)
{
  checkPtrType(object, "GtkRecentInfo")

  w <- .RGtkCall("S_gtk_recent_info_get_description", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentInfoGetMimeType <-
function(object)
{
  checkPtrType(object, "GtkRecentInfo")

  w <- .RGtkCall("S_gtk_recent_info_get_mime_type", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentInfoGetAdded <-
function(object)
{
  checkPtrType(object, "GtkRecentInfo")

  w <- .RGtkCall("S_gtk_recent_info_get_added", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentInfoGetModified <-
function(object)
{
  checkPtrType(object, "GtkRecentInfo")

  w <- .RGtkCall("S_gtk_recent_info_get_modified", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentInfoGetVisited <-
function(object)
{
  checkPtrType(object, "GtkRecentInfo")

  w <- .RGtkCall("S_gtk_recent_info_get_visited", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentInfoGetPrivateHint <-
function(object)
{
  checkPtrType(object, "GtkRecentInfo")

  w <- .RGtkCall("S_gtk_recent_info_get_private_hint", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentInfoGetApplicationInfo <-
function(object, app.name)
{
  checkPtrType(object, "GtkRecentInfo")
  app.name <- as.character(app.name)

  w <- .RGtkCall("S_gtk_recent_info_get_application_info", object, app.name, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentInfoGetApplications <-
function(object, length)
{
  checkPtrType(object, "GtkRecentInfo")
  length <- as.list(as.numeric(length))

  w <- .RGtkCall("S_gtk_recent_info_get_applications", object, length, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentInfoLastApplication <-
function(object)
{
  checkPtrType(object, "GtkRecentInfo")

  w <- .RGtkCall("S_gtk_recent_info_last_application", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentInfoHasApplication <-
function(object, app.name)
{
  checkPtrType(object, "GtkRecentInfo")
  app.name <- as.character(app.name)

  w <- .RGtkCall("S_gtk_recent_info_has_application", object, app.name, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentInfoGetGroups <-
function(object)
{
  checkPtrType(object, "GtkRecentInfo")

  w <- .RGtkCall("S_gtk_recent_info_get_groups", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentInfoHasGroup <-
function(object, group.name)
{
  checkPtrType(object, "GtkRecentInfo")
  group.name <- as.character(group.name)

  w <- .RGtkCall("S_gtk_recent_info_has_group", object, group.name, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentInfoGetIcon <-
function(object, size)
{
  checkPtrType(object, "GtkRecentInfo")
  size <- as.integer(size)

  w <- .RGtkCall("S_gtk_recent_info_get_icon", object, size, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentInfoGetShortName <-
function(object)
{
  checkPtrType(object, "GtkRecentInfo")

  w <- .RGtkCall("S_gtk_recent_info_get_short_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentInfoGetUriDisplay <-
function(object)
{
  checkPtrType(object, "GtkRecentInfo")

  w <- .RGtkCall("S_gtk_recent_info_get_uri_display", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentInfoGetAge <-
function(object)
{
  checkPtrType(object, "GtkRecentInfo")

  w <- .RGtkCall("S_gtk_recent_info_get_age", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentInfoIsLocal <-
function(object)
{
  checkPtrType(object, "GtkRecentInfo")

  w <- .RGtkCall("S_gtk_recent_info_is_local", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentInfoExists <-
function(object)
{
  checkPtrType(object, "GtkRecentInfo")

  w <- .RGtkCall("S_gtk_recent_info_exists", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentInfoMatch <-
function(object, info.b)
{
  checkPtrType(object, "GtkRecentInfo")
  checkPtrType(info.b, "GtkRecentInfo")

  w <- .RGtkCall("S_gtk_recent_info_match", object, info.b, PACKAGE = "RGtk2")

  return(w)
} 


gtkScrolledWindowUnsetPlacement <-
function(object)
{
  checkPtrType(object, "GtkScrolledWindow")

  w <- .RGtkCall("S_gtk_scrolled_window_unset_placement", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTargetListAddRichTextTargets <-
function(list, info, deserializable, buffer)
{
  checkPtrType(list, "GtkTargetList")
  info <- as.numeric(info)
  deserializable <- as.logical(deserializable)
  checkPtrType(buffer, "GtkTextBuffer")

  w <- .RGtkCall("S_gtk_target_list_add_rich_text_targets", list, info, deserializable, buffer, PACKAGE = "RGtk2")

  return(w)
} 


gtkTargetTableNewFromList <-
function(list)
{
  checkPtrType(list, "GtkTargetList")

  w <- .RGtkCall("S_gtk_target_table_new_from_list", list, PACKAGE = "RGtk2")

  return(w)
} 


gtkSelectionDataTargetsIncludeRichText <-
function(object, buffer)
{
  checkPtrType(object, "GtkSelectionData")
  checkPtrType(buffer, "GtkTextBuffer")

  w <- .RGtkCall("S_gtk_selection_data_targets_include_rich_text", object, buffer, PACKAGE = "RGtk2")

  return(w)
} 


gtkSelectionDataTargetsIncludeUri <-
function(object)
{
  checkPtrType(object, "GtkSelectionData")

  w <- .RGtkCall("S_gtk_selection_data_targets_include_uri", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTargetsIncludeText <-
function(targets)
{
  targets <- lapply(targets, function(x) { x <- as.GdkAtom(x); x })

  w <- .RGtkCall("S_gtk_targets_include_text", targets, PACKAGE = "RGtk2")

  return(w)
} 


gtkTargetsIncludeRichText <-
function(targets, buffer)
{
  targets <- lapply(targets, function(x) { x <- as.GdkAtom(x); x })
  checkPtrType(buffer, "GtkTextBuffer")

  w <- .RGtkCall("S_gtk_targets_include_rich_text", targets, buffer, PACKAGE = "RGtk2")

  return(w)
} 


gtkTargetsIncludeImage <-
function(targets, writable)
{
  targets <- lapply(targets, function(x) { x <- as.GdkAtom(x); x })
  writable <- as.logical(writable)

  w <- .RGtkCall("S_gtk_targets_include_image", targets, writable, PACKAGE = "RGtk2")

  return(w)
} 


gtkTargetsIncludeUri <-
function(targets)
{
  targets <- lapply(targets, function(x) { x <- as.GdkAtom(x); x })

  w <- .RGtkCall("S_gtk_targets_include_uri", targets, PACKAGE = "RGtk2")

  return(w)
} 


gtkTargetListGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_target_list_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkSizeGroupGetWidgets <-
function(object)
{
  checkPtrType(object, "GtkSizeGroup")

  w <- .RGtkCall("S_gtk_size_group_get_widgets", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkStatusIconGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_status_icon_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkStatusIconNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_status_icon_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkStatusIconNewFromPixbuf <-
function(pixbuf)
{
  checkPtrType(pixbuf, "GdkPixbuf")

  w <- .RGtkCall("S_gtk_status_icon_new_from_pixbuf", pixbuf, PACKAGE = "RGtk2")

  return(w)
} 


gtkStatusIconNewFromFile <-
function(filename)
{
  filename <- as.character(filename)

  w <- .RGtkCall("S_gtk_status_icon_new_from_file", filename, PACKAGE = "RGtk2")

  return(w)
} 


gtkStatusIconNewFromStock <-
function(stock.id)
{
  stock.id <- as.character(stock.id)

  w <- .RGtkCall("S_gtk_status_icon_new_from_stock", stock.id, PACKAGE = "RGtk2")

  return(w)
} 


gtkStatusIconNewFromIconName <-
function(icon.name)
{
  icon.name <- as.character(icon.name)

  w <- .RGtkCall("S_gtk_status_icon_new_from_icon_name", icon.name, PACKAGE = "RGtk2")

  return(w)
} 


gtkStatusIconSetFromPixbuf <-
function(object, pixbuf)
{
  checkPtrType(object, "GtkStatusIcon")
  checkPtrType(pixbuf, "GdkPixbuf")

  w <- .RGtkCall("S_gtk_status_icon_set_from_pixbuf", object, pixbuf, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkStatusIconSetFromFile <-
function(object, filename)
{
  checkPtrType(object, "GtkStatusIcon")
  filename <- as.character(filename)

  w <- .RGtkCall("S_gtk_status_icon_set_from_file", object, filename, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkStatusIconSetFromStock <-
function(object, stock.id)
{
  checkPtrType(object, "GtkStatusIcon")
  stock.id <- as.character(stock.id)

  w <- .RGtkCall("S_gtk_status_icon_set_from_stock", object, stock.id, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkStatusIconSetFromIconName <-
function(object, icon.name)
{
  checkPtrType(object, "GtkStatusIcon")
  icon.name <- as.character(icon.name)

  w <- .RGtkCall("S_gtk_status_icon_set_from_icon_name", object, icon.name, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkStatusIconGetStorageType <-
function(object)
{
  checkPtrType(object, "GtkStatusIcon")

  w <- .RGtkCall("S_gtk_status_icon_get_storage_type", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkStatusIconGetPixbuf <-
function(object)
{
  checkPtrType(object, "GtkStatusIcon")

  w <- .RGtkCall("S_gtk_status_icon_get_pixbuf", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkStatusIconGetStock <-
function(object)
{
  checkPtrType(object, "GtkStatusIcon")

  w <- .RGtkCall("S_gtk_status_icon_get_stock", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkStatusIconGetIconName <-
function(object)
{
  checkPtrType(object, "GtkStatusIcon")

  w <- .RGtkCall("S_gtk_status_icon_get_icon_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkStatusIconGetSize <-
function(object)
{
  checkPtrType(object, "GtkStatusIcon")

  w <- .RGtkCall("S_gtk_status_icon_get_size", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkStatusIconSetTooltip <-
function(object, tooltip.text)
{
  if(getOption("depwarn"))
    .Deprecated("setTooltipText", "RGtk2")

  checkPtrType(object, "GtkStatusIcon")
  tooltip.text <- as.character(tooltip.text)

  w <- .RGtkCall("S_gtk_status_icon_set_tooltip", object, tooltip.text, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkStatusIconSetVisible <-
function(object, visible)
{
  checkPtrType(object, "GtkStatusIcon")
  visible <- as.logical(visible)

  w <- .RGtkCall("S_gtk_status_icon_set_visible", object, visible, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkStatusIconGetVisible <-
function(object)
{
  checkPtrType(object, "GtkStatusIcon")

  w <- .RGtkCall("S_gtk_status_icon_get_visible", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkStatusIconSetBlinking <-
function(object, blinking)
{
  checkPtrType(object, "GtkStatusIcon")
  blinking <- as.logical(blinking)

  w <- .RGtkCall("S_gtk_status_icon_set_blinking", object, blinking, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkStatusIconGetBlinking <-
function(object)
{
  checkPtrType(object, "GtkStatusIcon")

  w <- .RGtkCall("S_gtk_status_icon_get_blinking", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkStatusIconIsEmbedded <-
function(object)
{
  checkPtrType(object, "GtkStatusIcon")

  w <- .RGtkCall("S_gtk_status_icon_is_embedded", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkStatusIconPositionMenu <-
function(menu, user.data)
{
  checkPtrType(menu, "GtkMenu")
  

  w <- .RGtkCall("S_gtk_status_icon_position_menu", menu, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkStatusIconGetGeometry <-
function(object)
{
  checkPtrType(object, "GtkStatusIcon")

  w <- .RGtkCall("S_gtk_status_icon_get_geometry", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkStyleLookupColor <-
function(object, color.name)
{
  checkPtrType(object, "GtkStyle")
  color.name <- as.character(color.name)

  w <- .RGtkCall("S_gtk_style_lookup_color", object, color.name, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferGetHasSelection <-
function(object)
{
  checkPtrType(object, "GtkTextBuffer")

  w <- .RGtkCall("S_gtk_text_buffer_get_has_selection", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferGetCopyTargetList <-
function(object)
{
  checkPtrType(object, "GtkTextBuffer")

  w <- .RGtkCall("S_gtk_text_buffer_get_copy_target_list", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferGetPasteTargetList <-
function(object)
{
  checkPtrType(object, "GtkTextBuffer")

  w <- .RGtkCall("S_gtk_text_buffer_get_paste_target_list", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferRegisterSerializeFormat <-
function(object, mime.type, fun, user.data)
{
  checkPtrType(object, "GtkTextBuffer")
  mime.type <- as.character(mime.type)
  fun <- as.function(fun)
  

  w <- .RGtkCall("S_gtk_text_buffer_register_serialize_format", object, mime.type, fun, user.data, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferRegisterSerializeTagset <-
function(object, tagset.name = NULL)
{
  checkPtrType(object, "GtkTextBuffer")
  if (!is.null( tagset.name )) tagset.name <- as.character(tagset.name)

  w <- .RGtkCall("S_gtk_text_buffer_register_serialize_tagset", object, tagset.name, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferRegisterDeserializeFormat <-
function(object, mime.type, fun, user.data)
{
  checkPtrType(object, "GtkTextBuffer")
  mime.type <- as.character(mime.type)
  fun <- as.function(fun)
  

  w <- .RGtkCall("S_gtk_text_buffer_register_deserialize_format", object, mime.type, fun, user.data, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferRegisterDeserializeTagset <-
function(object, tagset.name = NULL)
{
  checkPtrType(object, "GtkTextBuffer")
  if (!is.null( tagset.name )) tagset.name <- as.character(tagset.name)

  w <- .RGtkCall("S_gtk_text_buffer_register_deserialize_tagset", object, tagset.name, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferUnregisterSerializeFormat <-
function(object, format)
{
  checkPtrType(object, "GtkTextBuffer")
  format <- as.GdkAtom(format)

  w <- .RGtkCall("S_gtk_text_buffer_unregister_serialize_format", object, format, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextBufferUnregisterDeserializeFormat <-
function(object, format)
{
  checkPtrType(object, "GtkTextBuffer")
  format <- as.GdkAtom(format)

  w <- .RGtkCall("S_gtk_text_buffer_unregister_deserialize_format", object, format, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextBufferDeserializeSetCanCreateTags <-
function(object, format, can.create.tags)
{
  checkPtrType(object, "GtkTextBuffer")
  format <- as.GdkAtom(format)
  can.create.tags <- as.logical(can.create.tags)

  w <- .RGtkCall("S_gtk_text_buffer_deserialize_set_can_create_tags", object, format, can.create.tags, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextBufferDeserializeGetCanCreateTags <-
function(object, format)
{
  checkPtrType(object, "GtkTextBuffer")
  format <- as.GdkAtom(format)

  w <- .RGtkCall("S_gtk_text_buffer_deserialize_get_can_create_tags", object, format, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferGetSerializeFormats <-
function(object)
{
  checkPtrType(object, "GtkTextBuffer")

  w <- .RGtkCall("S_gtk_text_buffer_get_serialize_formats", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferGetDeserializeFormats <-
function(object)
{
  checkPtrType(object, "GtkTextBuffer")

  w <- .RGtkCall("S_gtk_text_buffer_get_deserialize_formats", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferSerialize <-
function(object, content.buffer, format, start, end)
{
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(content.buffer, "GtkTextBuffer")
  format <- as.GdkAtom(format)
  checkPtrType(start, "GtkTextIter")
  checkPtrType(end, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_buffer_serialize", object, content.buffer, format, start, end, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferDeserialize <-
function(object, content.buffer, format, iter, data, .errwarn = TRUE)
{
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(content.buffer, "GtkTextBuffer")
  format <- as.GdkAtom(format)
  checkPtrType(iter, "GtkTextIter")
  data <- as.list(as.raw(data))

  w <- .RGtkCall("S_gtk_text_buffer_deserialize", object, content.buffer, format, iter, data, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkTreeViewGetHeadersClickable <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_get_headers_clickable", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewGetSearchEntry <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_get_search_entry", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewSetSearchEntry <-
function(object, entry = NULL)
{
  checkPtrType(object, "GtkTreeView")
  if (!is.null( entry )) checkPtrType(entry, "GtkEntry")

  w <- .RGtkCall("S_gtk_tree_view_set_search_entry", object, entry, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewGetSearchPositionFunc <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_get_search_position_func", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewSetSearchPositionFunc <-
function(object, func, data)
{
  checkPtrType(object, "GtkTreeView")
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_tree_view_set_search_position_func", object, func, data, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewSetRubberBanding <-
function(object, enable)
{
  checkPtrType(object, "GtkTreeView")
  enable <- as.logical(enable)

  w <- .RGtkCall("S_gtk_tree_view_set_rubber_banding", object, enable, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewGetRubberBanding <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_get_rubber_banding", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewGetGridLines <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_get_grid_lines", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewSetGridLines <-
function(object, grid.lines)
{
  checkPtrType(object, "GtkTreeView")
  

  w <- .RGtkCall("S_gtk_tree_view_set_grid_lines", object, grid.lines, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewGetEnableTreeLines <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_get_enable_tree_lines", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewSetEnableTreeLines <-
function(object, enabled)
{
  checkPtrType(object, "GtkTreeView")
  enabled <- as.logical(enabled)

  w <- .RGtkCall("S_gtk_tree_view_set_enable_tree_lines", object, enabled, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAssistantPageTypeGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_assistant_page_type_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkCellRendererAccelModeGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_cell_renderer_accel_mode_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkSensitivityTypeGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_sensitivity_type_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintPagesGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_print_pages_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkPageSetGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_page_set_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkPageOrientationGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_page_orientation_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintQualityGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_print_quality_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintDuplexGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_print_duplex_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkUnitGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_unit_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewGridLinesGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_tree_view_grid_lines_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintOperationActionGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_print_operation_action_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentSortTypeGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_recent_sort_type_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentChooserErrorGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_recent_chooser_error_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentFilterFlagsGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_recent_filter_flags_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentManagerErrorGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_recent_manager_error_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferTargetInfoGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_text_buffer_target_info_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetIsComposited <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_is_composited", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetInputShapeCombineMask <-
function(object, shape.mask = NULL, offset.x, offset.y)
{
  checkPtrType(object, "GtkWidget")
  if (!is.null( shape.mask )) checkPtrType(shape.mask, "GdkBitmap")
  offset.x <- as.integer(offset.x)
  offset.y <- as.integer(offset.y)

  w <- .RGtkCall("S_gtk_widget_input_shape_combine_mask", object, shape.mask, offset.x, offset.y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowSetDeletable <-
function(object, setting)
{
  checkPtrType(object, "GtkWindow")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_window_set_deletable", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowGetDeletable <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_get_deletable", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowGetGroup <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_get_group", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkBuildableGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_buildable_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkBuildableSetName <-
function(object, name)
{
  checkPtrType(object, "GtkBuildable")
  name <- as.character(name)

  w <- .RGtkCall("S_gtk_buildable_set_name", object, name, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkBuildableGetName <-
function(object)
{
  checkPtrType(object, "GtkBuildable")

  w <- .RGtkCall("S_gtk_buildable_get_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkBuildableAddChild <-
function(object, builder, child, type)
{
  checkPtrType(object, "GtkBuildable")
  checkPtrType(builder, "GtkBuilder")
  checkPtrType(child, "GObject")
  type <- as.character(type)

  w <- .RGtkCall("S_gtk_buildable_add_child", object, builder, child, type, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkBuildableSetBuildableProperty <-
function(object, builder, name, value)
{
  checkPtrType(object, "GtkBuildable")
  checkPtrType(builder, "GtkBuilder")
  name <- as.character(name)
  

  w <- .RGtkCall("S_gtk_buildable_set_buildable_property", object, builder, name, value, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkBuildableConstructChild <-
function(object, builder, name)
{
  checkPtrType(object, "GtkBuildable")
  checkPtrType(builder, "GtkBuilder")
  name <- as.character(name)

  w <- .RGtkCall("S_gtk_buildable_construct_child", object, builder, name, PACKAGE = "RGtk2")

  return(w)
} 


gtkBuildableCustomTagStart <-
function(object, builder, child, tagname, parser, data)
{
  checkPtrType(object, "GtkBuildable")
  checkPtrType(builder, "GtkBuilder")
  checkPtrType(child, "GObject")
  tagname <- as.character(tagname)
  checkPtrType(parser, "GMarkupParser")
  

  w <- .RGtkCall("S_gtk_buildable_custom_tag_start", object, builder, child, tagname, parser, data, PACKAGE = "RGtk2")

  return(w)
} 


gtkBuildableCustomTagEnd <-
function(object, builder, child, tagname, data)
{
  checkPtrType(object, "GtkBuildable")
  checkPtrType(builder, "GtkBuilder")
  checkPtrType(child, "GObject")
  tagname <- as.character(tagname)
  

  w <- .RGtkCall("S_gtk_buildable_custom_tag_end", object, builder, child, tagname, data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkBuildableCustomFinished <-
function(object, builder, child, tagname, data)
{
  checkPtrType(object, "GtkBuildable")
  checkPtrType(builder, "GtkBuilder")
  checkPtrType(child, "GObject")
  tagname <- as.character(tagname)
  

  w <- .RGtkCall("S_gtk_buildable_custom_finished", object, builder, child, tagname, data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkBuildableParserFinished <-
function(object, builder)
{
  checkPtrType(object, "GtkBuildable")
  checkPtrType(builder, "GtkBuilder")

  w <- .RGtkCall("S_gtk_buildable_parser_finished", object, builder, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkBuildableGetInternalChild <-
function(object, builder, childname)
{
  checkPtrType(object, "GtkBuildable")
  checkPtrType(builder, "GtkBuilder")
  childname <- as.character(childname)

  w <- .RGtkCall("S_gtk_buildable_get_internal_child", object, builder, childname, PACKAGE = "RGtk2")

  return(w)
} 


gtkBuilderErrorQuark <-
function()
{
  

  w <- .RGtkCall("S_gtk_builder_error_quark", PACKAGE = "RGtk2")

  return(w)
} 


gtkBuilderGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_builder_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkBuilderAddFromFile <-
function(object, filename, .errwarn = TRUE)
{
  checkPtrType(object, "GtkBuilder")
  filename <- as.character(filename)

  w <- .RGtkCall("S_gtk_builder_add_from_file", object, filename, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkBuilderAddFromString <-
function(object, buffer, length, .errwarn = TRUE)
{
  checkPtrType(object, "GtkBuilder")
  buffer <- as.character(buffer)
  length <- as.numeric(length)

  w <- .RGtkCall("S_gtk_builder_add_from_string", object, buffer, length, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkBuilderGetObject <-
function(object, name)
{
  checkPtrType(object, "GtkBuilder")
  name <- as.character(name)

  w <- .RGtkCall("S_gtk_builder_get_object", object, name, PACKAGE = "RGtk2")

  return(w)
} 


gtkBuilderGetObjects <-
function(object)
{
  checkPtrType(object, "GtkBuilder")

  w <- .RGtkCall("S_gtk_builder_get_objects", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkBuilderConnectSignals <-
function(object, user.data = NULL)
{
  checkPtrType(object, "GtkBuilder")
  

  w <- .RGtkCall("S_gtk_builder_connect_signals", object, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkBuilderConnectSignalsFull <-
function(object, func, user.data)
{
  checkPtrType(object, "GtkBuilder")
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_builder_connect_signals_full", object, func, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkBuilderSetTranslationDomain <-
function(object, domain)
{
  checkPtrType(object, "GtkBuilder")
  domain <- as.character(domain)

  w <- .RGtkCall("S_gtk_builder_set_translation_domain", object, domain, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkBuilderGetTranslationDomain <-
function(object)
{
  checkPtrType(object, "GtkBuilder")

  w <- .RGtkCall("S_gtk_builder_get_translation_domain", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkBuilderGetTypeFromName <-
function(object, type.name)
{
  checkPtrType(object, "GtkBuilder")
  type.name <- as.character(type.name)

  w <- .RGtkCall("S_gtk_builder_get_type_from_name", object, type.name, PACKAGE = "RGtk2")

  return(w)
} 


gtkBuilderValueFromString <-
function(object, pspec, string, .errwarn = TRUE)
{
  checkPtrType(object, "GtkBuilder")
  pspec <- as.GParamSpec(pspec)
  string <- as.character(string)

  w <- .RGtkCall("S_gtk_builder_value_from_string", object, pspec, string, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkBuilderValueFromStringType <-
function(object, type, string, .errwarn = TRUE)
{
  checkPtrType(object, "GtkBuilder")
  type <- as.GType(type)
  string <- as.character(string)

  w <- .RGtkCall("S_gtk_builder_value_from_string_type", object, type, string, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkAboutDialogGetProgramName <-
function(object)
{
  checkPtrType(object, "GtkAboutDialog")

  w <- .RGtkCall("S_gtk_about_dialog_get_program_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkAboutDialogSetProgramName <-
function(object, name)
{
  checkPtrType(object, "GtkAboutDialog")
  name <- as.character(name)

  w <- .RGtkCall("S_gtk_about_dialog_set_program_name", object, name, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActionCreateMenu <-
function(object)
{
  checkPtrType(object, "GtkAction")

  w <- .RGtkCall("S_gtk_action_create_menu", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkCellLayoutGetCells <-
function(object)
{
  checkPtrType(object, "GtkCellLayout")

  w <- .RGtkCall("S_gtk_cell_layout_get_cells", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryCompletionGetCompletionPrefix <-
function(object)
{
  checkPtrType(object, "GtkEntryCompletion")

  w <- .RGtkCall("S_gtk_entry_completion_get_completion_prefix", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryCompletionSetInlineSelection <-
function(object, inline.selection)
{
  checkPtrType(object, "GtkEntryCompletion")
  inline.selection <- as.logical(inline.selection)

  w <- .RGtkCall("S_gtk_entry_completion_set_inline_selection", object, inline.selection, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryCompletionGetInlineSelection <-
function(object)
{
  checkPtrType(object, "GtkEntryCompletion")

  w <- .RGtkCall("S_gtk_entry_completion_get_inline_selection", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntrySetCursorHadjustment <-
function(object, adjustment)
{
  checkPtrType(object, "GtkEntry")
  checkPtrType(adjustment, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_entry_set_cursor_hadjustment", object, adjustment, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryGetCursorHadjustment <-
function(object)
{
  checkPtrType(object, "GtkEntry")

  w <- .RGtkCall("S_gtk_entry_get_cursor_hadjustment", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconThemeChooseIcon <-
function(object, icon.names, size, flags)
{
  checkPtrType(object, "GtkIconTheme")
  icon.names <- as.list(as.character(icon.names))
  size <- as.integer(size)
  

  w <- .RGtkCall("S_gtk_icon_theme_choose_icon", object, icon.names, size, flags, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconThemeListContexts <-
function(object)
{
  checkPtrType(object, "GtkIconTheme")

  w <- .RGtkCall("S_gtk_icon_theme_list_contexts", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconViewConvertWidgetToBinWindowCoords <-
function(object, wx, wy)
{
  checkPtrType(object, "GtkIconView")
  wx <- as.integer(wx)
  wy <- as.integer(wy)

  w <- .RGtkCall("S_gtk_icon_view_convert_widget_to_bin_window_coords", object, wx, wy, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconViewSetTooltipItem <-
function(object, tooltip, path)
{
  checkPtrType(object, "GtkIconView")
  checkPtrType(tooltip, "GtkTooltip")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_icon_view_set_tooltip_item", object, tooltip, path, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconViewSetTooltipCell <-
function(object, tooltip, path, cell)
{
  checkPtrType(object, "GtkIconView")
  checkPtrType(tooltip, "GtkTooltip")
  checkPtrType(path, "GtkTreePath")
  checkPtrType(cell, "GtkCellRenderer")

  w <- .RGtkCall("S_gtk_icon_view_set_tooltip_cell", object, tooltip, path, cell, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconViewSetTooltipColumn <-
function(object, column)
{
  checkPtrType(object, "GtkIconView")
  column <- as.integer(column)

  w <- .RGtkCall("S_gtk_icon_view_set_tooltip_column", object, column, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconViewGetTooltipColumn <-
function(object)
{
  checkPtrType(object, "GtkIconView")

  w <- .RGtkCall("S_gtk_icon_view_get_tooltip_column", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkListStoreSetValuesv <-
function(object, iter, columns, values)
{
  checkPtrType(object, "GtkListStore")
  checkPtrType(iter, "GtkTreeIter")
  columns <- as.list(as.integer(columns))
  values <- as.list(values)

  w <- .RGtkCall("S_gtk_list_store_set_valuesv", object, iter, columns, values, PACKAGE = "RGtk2")

  return(w)
} 


gtkMenuToolButtonSetArrowTooltipText <-
function(object, text)
{
  checkPtrType(object, "GtkMenuToolButton")
  text <- as.character(text)

  w <- .RGtkCall("S_gtk_menu_tool_button_set_arrow_tooltip_text", object, text, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuToolButtonSetArrowTooltipMarkup <-
function(object, markup)
{
  checkPtrType(object, "GtkMenuToolButton")
  markup <- as.character(markup)

  w <- .RGtkCall("S_gtk_menu_tool_button_set_arrow_tooltip_markup", object, markup, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkNotebookSetGroup <-
function(object, group)
{
  checkPtrType(object, "GtkNotebook")
  

  w <- .RGtkCall("S_gtk_notebook_set_group", object, group, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkNotebookGetGroup <-
function(object)
{
  checkPtrType(object, "GtkNotebook")

  w <- .RGtkCall("S_gtk_notebook_get_group", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPageSetupNewFromFile <-
function(file.name, .errwarn = TRUE)
{
  file.name <- as.character(file.name)

  w <- .RGtkCall("S_gtk_page_setup_new_from_file", file.name, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkPageSetupNewFromKeyFile <-
function(key.file, group.name, .errwarn = TRUE)
{
  checkPtrType(key.file, "GKeyFile")
  group.name <- as.character(group.name)

  w <- .RGtkCall("S_gtk_page_setup_new_from_key_file", key.file, group.name, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkPageSetupToFile <-
function(object, file.name, .errwarn = TRUE)
{
  checkPtrType(object, "GtkPageSetup")
  file.name <- as.character(file.name)

  w <- .RGtkCall("S_gtk_page_setup_to_file", object, file.name, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkPageSetupToKeyFile <-
function(object, key.file, group.name)
{
  checkPtrType(object, "GtkPageSetup")
  checkPtrType(key.file, "GKeyFile")
  group.name <- as.character(group.name)

  w <- .RGtkCall("S_gtk_page_setup_to_key_file", object, key.file, group.name, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPaperSizeGetPaperSizes <-
function(include.custom)
{
  include.custom <- as.logical(include.custom)

  w <- .RGtkCall("S_gtk_paper_size_get_paper_sizes", include.custom, PACKAGE = "RGtk2")

  return(w)
} 


gtkPaperSizeNewFromKeyFile <-
function(key.file, group.name, .errwarn = TRUE)
{
  checkPtrType(key.file, "GKeyFile")
  group.name <- as.character(group.name)

  w <- .RGtkCall("S_gtk_paper_size_new_from_key_file", key.file, group.name, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkPaperSizeToKeyFile <-
function(object, key.file, group.name)
{
  checkPtrType(object, "GtkPaperSize")
  checkPtrType(key.file, "GKeyFile")
  group.name <- as.character(group.name)

  w <- .RGtkCall("S_gtk_paper_size_to_key_file", object, key.file, group.name, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintSettingsNewFromFile <-
function(file.name, .errwarn = TRUE)
{
  file.name <- as.character(file.name)

  w <- .RGtkCall("S_gtk_print_settings_new_from_file", file.name, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkPrintSettingsToFile <-
function(object, file.name, .errwarn = TRUE)
{
  checkPtrType(object, "GtkPrintSettings")
  file.name <- as.character(file.name)

  w <- .RGtkCall("S_gtk_print_settings_to_file", object, file.name, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkPrintSettingsNewFromKeyFile <-
function(key.file, group.name, .errwarn = TRUE)
{
  checkPtrType(key.file, "GKeyFile")
  group.name <- as.character(group.name)

  w <- .RGtkCall("S_gtk_print_settings_new_from_key_file", key.file, group.name, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkPrintSettingsToKeyFile <-
function(object, key.file, group.name)
{
  checkPtrType(object, "GtkPrintSettings")
  checkPtrType(key.file, "GKeyFile")
  group.name <- as.character(group.name)

  w <- .RGtkCall("S_gtk_print_settings_to_key_file", object, key.file, group.name, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRangeSetShowFillLevel <-
function(object, show.fill.level)
{
  checkPtrType(object, "GtkRange")
  show.fill.level <- as.logical(show.fill.level)

  w <- .RGtkCall("S_gtk_range_set_show_fill_level", object, show.fill.level, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRangeGetShowFillLevel <-
function(object)
{
  checkPtrType(object, "GtkRange")

  w <- .RGtkCall("S_gtk_range_get_show_fill_level", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRangeSetRestrictToFillLevel <-
function(object, restrict.to.fill.level)
{
  checkPtrType(object, "GtkRange")
  restrict.to.fill.level <- as.logical(restrict.to.fill.level)

  w <- .RGtkCall("S_gtk_range_set_restrict_to_fill_level", object, restrict.to.fill.level, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRangeGetRestrictToFillLevel <-
function(object)
{
  checkPtrType(object, "GtkRange")

  w <- .RGtkCall("S_gtk_range_get_restrict_to_fill_level", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRangeSetFillLevel <-
function(object, fill.level)
{
  checkPtrType(object, "GtkRange")
  fill.level <- as.numeric(fill.level)

  w <- .RGtkCall("S_gtk_range_set_fill_level", object, fill.level, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRangeGetFillLevel <-
function(object)
{
  checkPtrType(object, "GtkRange")

  w <- .RGtkCall("S_gtk_range_get_fill_level", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRcParseColorFull <-
function(scanner, style)
{
  checkPtrType(scanner, "GScanner")
  checkPtrType(style, "GtkRcStyle")

  w <- .RGtkCall("S_gtk_rc_parse_color_full", scanner, style, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentActionGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_recent_action_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentActionNew <-
function(name, label, tooltip, stock.id)
{
  name <- as.character(name)
  label <- as.character(label)
  tooltip <- as.character(tooltip)
  stock.id <- as.character(stock.id)

  w <- .RGtkCall("S_gtk_recent_action_new", name, label, tooltip, stock.id, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentActionNewForManager <-
function(name, label, tooltip, stock.id, manager)
{
  name <- as.character(name)
  label <- as.character(label)
  tooltip <- as.character(tooltip)
  stock.id <- as.character(stock.id)
  checkPtrType(manager, "GtkRecentManager")

  w <- .RGtkCall("S_gtk_recent_action_new_for_manager", name, label, tooltip, stock.id, manager, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentActionGetShowNumbers <-
function(object)
{
  checkPtrType(object, "GtkRecentAction")

  w <- .RGtkCall("S_gtk_recent_action_get_show_numbers", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRecentActionSetShowNumbers <-
function(object, show.numbers)
{
  checkPtrType(object, "GtkRecentAction")
  show.numbers <- as.logical(show.numbers)

  w <- .RGtkCall("S_gtk_recent_action_set_show_numbers", object, show.numbers, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkScaleButtonGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_scale_button_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkScaleButtonNew <-
function(size, min, max, step, icons, show = TRUE)
{
  
  min <- as.numeric(min)
  max <- as.numeric(max)
  step <- as.numeric(step)
  icons <- as.list(as.character(icons))

  w <- .RGtkCall("S_gtk_scale_button_new", size, min, max, step, icons, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkScaleButtonSetIcons <-
function(object, icons)
{
  checkPtrType(object, "GtkScaleButton")
  icons <- as.list(as.character(icons))

  w <- .RGtkCall("S_gtk_scale_button_set_icons", object, icons, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkScaleButtonGetValue <-
function(object)
{
  checkPtrType(object, "GtkScaleButton")

  w <- .RGtkCall("S_gtk_scale_button_get_value", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkScaleButtonSetValue <-
function(object, value)
{
  checkPtrType(object, "GtkScaleButton")
  value <- as.numeric(value)

  w <- .RGtkCall("S_gtk_scale_button_set_value", object, value, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkScaleButtonGetAdjustment <-
function(object)
{
  checkPtrType(object, "GtkScaleButton")

  w <- .RGtkCall("S_gtk_scale_button_get_adjustment", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkScaleButtonSetAdjustment <-
function(object, adjustment)
{
  checkPtrType(object, "GtkScaleButton")
  checkPtrType(adjustment, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_scale_button_set_adjustment", object, adjustment, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkStatusIconSetScreen <-
function(object, screen)
{
  checkPtrType(object, "GtkStatusIcon")
  checkPtrType(screen, "GdkScreen")

  w <- .RGtkCall("S_gtk_status_icon_set_screen", object, screen, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkStatusIconGetScreen <-
function(object)
{
  checkPtrType(object, "GtkStatusIcon")

  w <- .RGtkCall("S_gtk_status_icon_get_screen", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTextBufferAddMark <-
function(object, mark, where)
{
  checkPtrType(object, "GtkTextBuffer")
  checkPtrType(mark, "GtkTextMark")
  checkPtrType(where, "GtkTextIter")

  w <- .RGtkCall("S_gtk_text_buffer_add_mark", object, mark, where, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTextMarkNew <-
function(name, left.gravity)
{
  name <- as.character(name)
  left.gravity <- as.logical(left.gravity)

  w <- .RGtkCall("S_gtk_text_mark_new", name, left.gravity, PACKAGE = "RGtk2")

  return(w)
} 


gtkTooltipGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_tooltip_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkTooltipSetMarkup <-
function(object, markup)
{
  checkPtrType(object, "GtkTooltip")
  markup <- as.character(markup)

  w <- .RGtkCall("S_gtk_tooltip_set_markup", object, markup, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTooltipSetText <-
function(object, text)
{
  checkPtrType(object, "GtkTooltip")
  text <- as.character(text)

  w <- .RGtkCall("S_gtk_tooltip_set_text", object, text, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTooltipSetIcon <-
function(object, pixbuf)
{
  checkPtrType(object, "GtkTooltip")
  checkPtrType(pixbuf, "GdkPixbuf")

  w <- .RGtkCall("S_gtk_tooltip_set_icon", object, pixbuf, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTooltipSetIconFromStock <-
function(object, stock.id, size)
{
  checkPtrType(object, "GtkTooltip")
  stock.id <- as.character(stock.id)
  

  w <- .RGtkCall("S_gtk_tooltip_set_icon_from_stock", object, stock.id, size, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTooltipSetCustom <-
function(object, custom.widget)
{
  checkPtrType(object, "GtkTooltip")
  checkPtrType(custom.widget, "GtkWidget")

  w <- .RGtkCall("S_gtk_tooltip_set_custom", object, custom.widget, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTooltipTriggerTooltipQuery <-
function(display)
{
  checkPtrType(display, "GdkDisplay")

  w <- .RGtkCall("S_gtk_tooltip_trigger_tooltip_query", display, PACKAGE = "RGtk2")

  return(w)
} 


gtkTooltipSetTipArea <-
function(object, area)
{
  checkPtrType(object, "GtkTooltip")
  area <- as.GdkRectangle(area)

  w <- .RGtkCall("S_gtk_tooltip_set_tip_area", object, area, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolItemSetTooltipText <-
function(object, text)
{
  checkPtrType(object, "GtkToolItem")
  text <- as.character(text)

  w <- .RGtkCall("S_gtk_tool_item_set_tooltip_text", object, text, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolItemSetTooltipMarkup <-
function(object, markup)
{
  checkPtrType(object, "GtkToolItem")
  markup <- as.character(markup)

  w <- .RGtkCall("S_gtk_tool_item_set_tooltip_markup", object, markup, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeStoreSetValuesv <-
function(object, iter, columns, values)
{
  checkPtrType(object, "GtkTreeStore")
  checkPtrType(iter, "GtkTreeIter")
  columns <- as.list(as.integer(columns))
  values <- as.list(values)

  w <- .RGtkCall("S_gtk_tree_store_set_valuesv", object, iter, columns, values, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewColumnGetTreeView <-
function(object)
{
  checkPtrType(object, "GtkTreeViewColumn")

  w <- .RGtkCall("S_gtk_tree_view_column_get_tree_view", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewConvertWidgetToTreeCoords <-
function(object, wx, wy)
{
  checkPtrType(object, "GtkTreeView")
  wx <- as.integer(wx)
  wy <- as.integer(wy)

  w <- .RGtkCall("S_gtk_tree_view_convert_widget_to_tree_coords", object, wx, wy, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewConvertTreeToWidgetCoords <-
function(object, tx, ty)
{
  checkPtrType(object, "GtkTreeView")
  tx <- as.integer(tx)
  ty <- as.integer(ty)

  w <- .RGtkCall("S_gtk_tree_view_convert_tree_to_widget_coords", object, tx, ty, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewConvertWidgetToBinWindowCoords <-
function(object, wx, wy)
{
  checkPtrType(object, "GtkTreeView")
  wx <- as.integer(wx)
  wy <- as.integer(wy)

  w <- .RGtkCall("S_gtk_tree_view_convert_widget_to_bin_window_coords", object, wx, wy, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewConvertBinWindowToWidgetCoords <-
function(object, bx, by)
{
  checkPtrType(object, "GtkTreeView")
  bx <- as.integer(bx)
  by <- as.integer(by)

  w <- .RGtkCall("S_gtk_tree_view_convert_bin_window_to_widget_coords", object, bx, by, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewConvertTreeToBinWindowCoords <-
function(object, tx, ty)
{
  checkPtrType(object, "GtkTreeView")
  tx <- as.integer(tx)
  ty <- as.integer(ty)

  w <- .RGtkCall("S_gtk_tree_view_convert_tree_to_bin_window_coords", object, tx, ty, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewConvertBinWindowToTreeCoords <-
function(object, bx, by)
{
  checkPtrType(object, "GtkTreeView")
  bx <- as.integer(bx)
  by <- as.integer(by)

  w <- .RGtkCall("S_gtk_tree_view_convert_bin_window_to_tree_coords", object, bx, by, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewSetShowExpanders <-
function(object, enabled)
{
  checkPtrType(object, "GtkTreeView")
  enabled <- as.logical(enabled)

  w <- .RGtkCall("S_gtk_tree_view_set_show_expanders", object, enabled, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewGetShowExpanders <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_get_show_expanders", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewSetLevelIndentation <-
function(object, indentation)
{
  checkPtrType(object, "GtkTreeView")
  indentation <- as.integer(indentation)

  w <- .RGtkCall("S_gtk_tree_view_set_level_indentation", object, indentation, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewGetLevelIndentation <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_get_level_indentation", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewIsRubberBandingActive <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_is_rubber_banding_active", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewSetTooltipColumn <-
function(object, column)
{
  checkPtrType(object, "GtkTreeView")
  column <- as.integer(column)

  w <- .RGtkCall("S_gtk_tree_view_set_tooltip_column", object, column, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewGetTooltipColumn <-
function(object)
{
  checkPtrType(object, "GtkTreeView")

  w <- .RGtkCall("S_gtk_tree_view_get_tooltip_column", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTreeViewSetTooltipRow <-
function(object, tooltip, path)
{
  checkPtrType(object, "GtkTreeView")
  checkPtrType(tooltip, "GtkTooltip")
  checkPtrType(path, "GtkTreePath")

  w <- .RGtkCall("S_gtk_tree_view_set_tooltip_row", object, tooltip, path, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeViewSetTooltipCell <-
function(object, tooltip, path, column, cell)
{
  checkPtrType(object, "GtkTreeView")
  checkPtrType(tooltip, "GtkTooltip")
  checkPtrType(path, "GtkTreePath")
  checkPtrType(column, "GtkTreeViewColumn")
  checkPtrType(cell, "GtkCellRenderer")

  w <- .RGtkCall("S_gtk_tree_view_set_tooltip_cell", object, tooltip, path, column, cell, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkVolumeButtonGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_volume_button_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkVolumeButtonNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_volume_button_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkWidgetKeynavFailed <-
function(object, direction)
{
  checkPtrType(object, "GtkWidget")
  

  w <- .RGtkCall("S_gtk_widget_keynav_failed", object, direction, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetErrorBell <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_error_bell", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetSetTooltipWindow <-
function(object, custom.window)
{
  checkPtrType(object, "GtkWidget")
  checkPtrType(custom.window, "GtkWindow")

  w <- .RGtkCall("S_gtk_widget_set_tooltip_window", object, custom.window, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetGetTooltipWindow <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_tooltip_window", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetTriggerTooltipQuery <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_trigger_tooltip_query", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetSetTooltipText <-
function(object, text)
{
  checkPtrType(object, "GtkWidget")
  text <- as.character(text)

  w <- .RGtkCall("S_gtk_widget_set_tooltip_text", object, text, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetGetTooltipText <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_tooltip_text", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetSetTooltipMarkup <-
function(object, markup)
{
  checkPtrType(object, "GtkWidget")
  markup <- as.character(markup)

  w <- .RGtkCall("S_gtk_widget_set_tooltip_markup", object, markup, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetGetTooltipMarkup <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_tooltip_markup", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetModifyCursor <-
function(object, primary, secondary)
{
  checkPtrType(object, "GtkWidget")
  primary <- as.GdkColor(primary)
  secondary <- as.GdkColor(secondary)

  w <- .RGtkCall("S_gtk_widget_modify_cursor", object, primary, secondary, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetGetHasTooltip <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_has_tooltip", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetSetHasTooltip <-
function(object, has.tooltip)
{
  checkPtrType(object, "GtkWidget")
  has.tooltip <- as.logical(has.tooltip)

  w <- .RGtkCall("S_gtk_widget_set_has_tooltip", object, has.tooltip, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowSetOpacity <-
function(object, opacity)
{
  checkPtrType(object, "GtkWindow")
  opacity <- as.numeric(opacity)

  w <- .RGtkCall("S_gtk_window_set_opacity", object, opacity, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowGetOpacity <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_get_opacity", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowSetStartupId <-
function(object, startup.id)
{
  checkPtrType(object, "GtkWindow")
  startup.id <- as.character(startup.id)

  w <- .RGtkCall("S_gtk_window_set_startup_id", object, startup.id, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAccelGroupGetIsLocked <-
function(object)
{
  checkPtrType(object, "GtkAccelGroup")

  w <- .RGtkCall("S_gtk_accel_group_get_is_locked", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkAccelGroupGetModifierMask <-
function(object)
{
  checkPtrType(object, "GtkAccelGroup")

  w <- .RGtkCall("S_gtk_accel_group_get_modifier_mask", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkAdjustmentGetLower <-
function(object)
{
  checkPtrType(object, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_adjustment_get_lower", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkAdjustmentSetLower <-
function(object, lower)
{
  checkPtrType(object, "GtkAdjustment")
  lower <- as.numeric(lower)

  w <- .RGtkCall("S_gtk_adjustment_set_lower", object, lower, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAdjustmentGetUpper <-
function(object)
{
  checkPtrType(object, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_adjustment_get_upper", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkAdjustmentSetUpper <-
function(object, upper)
{
  checkPtrType(object, "GtkAdjustment")
  upper <- as.numeric(upper)

  w <- .RGtkCall("S_gtk_adjustment_set_upper", object, upper, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAdjustmentGetStepIncrement <-
function(object)
{
  checkPtrType(object, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_adjustment_get_step_increment", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkAdjustmentSetStepIncrement <-
function(object, step.increment)
{
  checkPtrType(object, "GtkAdjustment")
  step.increment <- as.numeric(step.increment)

  w <- .RGtkCall("S_gtk_adjustment_set_step_increment", object, step.increment, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAdjustmentGetPageIncrement <-
function(object)
{
  checkPtrType(object, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_adjustment_get_page_increment", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkAdjustmentSetPageIncrement <-
function(object, page.increment)
{
  checkPtrType(object, "GtkAdjustment")
  page.increment <- as.numeric(page.increment)

  w <- .RGtkCall("S_gtk_adjustment_set_page_increment", object, page.increment, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAdjustmentGetPageSize <-
function(object)
{
  checkPtrType(object, "GtkAdjustment")

  w <- .RGtkCall("S_gtk_adjustment_get_page_size", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkAdjustmentSetPageSize <-
function(object, page.size)
{
  checkPtrType(object, "GtkAdjustment")
  page.size <- as.numeric(page.size)

  w <- .RGtkCall("S_gtk_adjustment_set_page_size", object, page.size, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkAdjustmentConfigure <-
function(object, value, lower, upper, step.increment, page.increment, page.size)
{
  checkPtrType(object, "GtkAdjustment")
  value <- as.numeric(value)
  lower <- as.numeric(lower)
  upper <- as.numeric(upper)
  step.increment <- as.numeric(step.increment)
  page.increment <- as.numeric(page.increment)
  page.size <- as.numeric(page.size)

  w <- .RGtkCall("S_gtk_adjustment_configure", object, value, lower, upper, step.increment, page.increment, page.size, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkBuilderAddObjectsFromFile <-
function(object, filename, object.ids, .errwarn = TRUE)
{
  checkPtrType(object, "GtkBuilder")
  filename <- as.character(filename)
  object.ids <- as.list(as.character(object.ids))

  w <- .RGtkCall("S_gtk_builder_add_objects_from_file", object, filename, object.ids, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkBuilderAddObjectsFromString <-
function(object, buffer, length, object.ids, .errwarn = TRUE)
{
  checkPtrType(object, "GtkBuilder")
  buffer <- as.character(buffer)
  length <- as.numeric(length)
  object.ids <- as.list(as.character(object.ids))

  w <- .RGtkCall("S_gtk_builder_add_objects_from_string", object, buffer, length, object.ids, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkCalendarSetDetailFunc <-
function(object, func, data)
{
  checkPtrType(object, "GtkCalendar")
  func <- as.function(func)
  

  w <- .RGtkCall("S_gtk_calendar_set_detail_func", object, func, data, PACKAGE = "RGtk2")

  return(w)
} 


gtkCalendarSetDetailWidthChars <-
function(object, chars)
{
  checkPtrType(object, "GtkCalendar")
  chars <- as.integer(chars)

  w <- .RGtkCall("S_gtk_calendar_set_detail_width_chars", object, chars, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCalendarSetDetailHeightRows <-
function(object, rows)
{
  checkPtrType(object, "GtkCalendar")
  rows <- as.integer(rows)

  w <- .RGtkCall("S_gtk_calendar_set_detail_height_rows", object, rows, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCalendarGetDetailWidthChars <-
function(object)
{
  checkPtrType(object, "GtkCalendar")

  w <- .RGtkCall("S_gtk_calendar_get_detail_width_chars", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkCalendarGetDetailHeightRows <-
function(object)
{
  checkPtrType(object, "GtkCalendar")

  w <- .RGtkCall("S_gtk_calendar_get_detail_height_rows", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkClipboardWaitIsUrisAvailable <-
function(object)
{
  checkPtrType(object, "GtkClipboard")

  w <- .RGtkCall("S_gtk_clipboard_wait_is_uris_available", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkClipboardWaitForUris <-
function(object)
{
  checkPtrType(object, "GtkClipboard")

  w <- .RGtkCall("S_gtk_clipboard_wait_for_uris", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkClipboardRequestUris <-
function(object, callback, user.data)
{
  checkPtrType(object, "GtkClipboard")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gtk_clipboard_request_uris", object, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkColorSelectionDialogGetColorSelection <-
function(object)
{
  checkPtrType(object, "GtkColorSelectionDialog")

  w <- .RGtkCall("S_gtk_color_selection_dialog_get_color_selection", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkComboBoxSetButtonSensitivity <-
function(object, sensitivity)
{
  checkPtrType(object, "GtkComboBox")
  

  w <- .RGtkCall("S_gtk_combo_box_set_button_sensitivity", object, sensitivity, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkComboBoxGetButtonSensitivity <-
function(object)
{
  checkPtrType(object, "GtkComboBox")

  w <- .RGtkCall("S_gtk_combo_box_get_button_sensitivity", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkContainerGetFocusChild <-
function(object)
{
  checkPtrType(object, "GtkContainer")

  w <- .RGtkCall("S_gtk_container_get_focus_child", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkDialogGetActionArea <-
function(object)
{
  checkPtrType(object, "GtkDialog")

  w <- .RGtkCall("S_gtk_dialog_get_action_area", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkDialogGetContentArea <-
function(object)
{
  checkPtrType(object, "GtkDialog")

  w <- .RGtkCall("S_gtk_dialog_get_content_area", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntrySetOverwriteMode <-
function(object, overwrite)
{
  checkPtrType(object, "GtkEntry")
  overwrite <- as.logical(overwrite)

  w <- .RGtkCall("S_gtk_entry_set_overwrite_mode", object, overwrite, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryGetOverwriteMode <-
function(object)
{
  checkPtrType(object, "GtkEntry")

  w <- .RGtkCall("S_gtk_entry_get_overwrite_mode", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryGetTextLength <-
function(object)
{
  checkPtrType(object, "GtkEntry")

  w <- .RGtkCall("S_gtk_entry_get_text_length", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserGetFile <-
function(object)
{
  checkPtrType(object, "GtkFileChooser")

  w <- .RGtkCall("S_gtk_file_chooser_get_file", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserSetFile <-
function(object, file, .errwarn = TRUE)
{
  checkPtrType(object, "GtkFileChooser")
  checkPtrType(file, "GFile")

  w <- .RGtkCall("S_gtk_file_chooser_set_file", object, file, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkFileChooserSelectFile <-
function(object, file, .errwarn = TRUE)
{
  checkPtrType(object, "GtkFileChooser")
  checkPtrType(file, "GFile")

  w <- .RGtkCall("S_gtk_file_chooser_select_file", object, file, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkFileChooserUnselectFile <-
function(object, file)
{
  checkPtrType(object, "GtkFileChooser")
  checkPtrType(file, "GFile")

  w <- .RGtkCall("S_gtk_file_chooser_unselect_file", object, file, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFileChooserGetFiles <-
function(object)
{
  checkPtrType(object, "GtkFileChooser")

  w <- .RGtkCall("S_gtk_file_chooser_get_files", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserSetCurrentFolderFile <-
function(object, file, .errwarn = TRUE)
{
  checkPtrType(object, "GtkFileChooser")
  checkPtrType(file, "GFile")

  w <- .RGtkCall("S_gtk_file_chooser_set_current_folder_file", object, file, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkFileChooserGetCurrentFolderFile <-
function(object)
{
  checkPtrType(object, "GtkFileChooser")

  w <- .RGtkCall("S_gtk_file_chooser_get_current_folder_file", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFileChooserGetPreviewFile <-
function(object)
{
  checkPtrType(object, "GtkFileChooser")

  w <- .RGtkCall("S_gtk_file_chooser_get_preview_file", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFontSelectionDialogGetOkButton <-
function(object)
{
  checkPtrType(object, "GtkFontSelectionDialog")

  w <- .RGtkCall("S_gtk_font_selection_dialog_get_ok_button", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFontSelectionDialogGetApplyButton <-
function(object)
{
  if(getOption("depwarn"))
    .Deprecated("don't use this method", "RGtk2")

  checkPtrType(object, "GtkFontSelectionDialog")

  w <- .RGtkCall("S_gtk_font_selection_dialog_get_apply_button", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFontSelectionDialogGetCancelButton <-
function(object)
{
  checkPtrType(object, "GtkFontSelectionDialog")

  w <- .RGtkCall("S_gtk_font_selection_dialog_get_cancel_button", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFontSelectionGetFamilyList <-
function(object)
{
  checkPtrType(object, "GtkFontSelection")

  w <- .RGtkCall("S_gtk_font_selection_get_family_list", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFontSelectionGetFaceList <-
function(object)
{
  checkPtrType(object, "GtkFontSelection")

  w <- .RGtkCall("S_gtk_font_selection_get_face_list", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFontSelectionGetSizeEntry <-
function(object)
{
  checkPtrType(object, "GtkFontSelection")

  w <- .RGtkCall("S_gtk_font_selection_get_size_entry", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFontSelectionGetSizeList <-
function(object)
{
  checkPtrType(object, "GtkFontSelection")

  w <- .RGtkCall("S_gtk_font_selection_get_size_list", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFontSelectionGetPreviewEntry <-
function(object)
{
  checkPtrType(object, "GtkFontSelection")

  w <- .RGtkCall("S_gtk_font_selection_get_preview_entry", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFontSelectionGetFamily <-
function(object)
{
  checkPtrType(object, "GtkFontSelection")

  w <- .RGtkCall("S_gtk_font_selection_get_family", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFontSelectionGetFace <-
function(object)
{
  checkPtrType(object, "GtkFontSelection")

  w <- .RGtkCall("S_gtk_font_selection_get_face", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkFontSelectionGetSize <-
function(object)
{
  checkPtrType(object, "GtkFontSelection")

  w <- .RGtkCall("S_gtk_font_selection_get_size", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkHandleBoxGetChildDetached <-
function(object)
{
  checkPtrType(object, "GtkHandleBox")

  w <- .RGtkCall("S_gtk_handle_box_get_child_detached", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconInfoNewForPixbuf <-
function(icon.theme, pixbuf)
{
  checkPtrType(icon.theme, "GtkIconTheme")
  checkPtrType(pixbuf, "GdkPixbuf")

  w <- .RGtkCall("S_gtk_icon_info_new_for_pixbuf", icon.theme, pixbuf, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconThemeLookupByGicon <-
function(object, icon, size, flags)
{
  checkPtrType(object, "GtkIconTheme")
  checkPtrType(icon, "GIcon")
  size <- as.integer(size)
  

  w <- .RGtkCall("S_gtk_icon_theme_lookup_by_gicon", object, icon, size, flags, PACKAGE = "RGtk2")

  return(w)
} 


gtkImageSetFromGicon <-
function(object, icon, size)
{
  checkPtrType(object, "GtkImage")
  checkPtrType(icon, "GIcon")
  

  w <- .RGtkCall("S_gtk_image_set_from_gicon", object, icon, size, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkImageNewFromGicon <-
function(icon, size, show = TRUE)
{
  checkPtrType(icon, "GIcon")
  

  w <- .RGtkCall("S_gtk_image_new_from_gicon", icon, size, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkImageGetGicon <-
function(object)
{
  checkPtrType(object, "GtkImage")

  w <- .RGtkCall("S_gtk_image_get_gicon", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkLayoutGetBinWindow <-
function(object)
{
  checkPtrType(object, "GtkLayout")

  w <- .RGtkCall("S_gtk_layout_get_bin_window", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkMenuGetAccelPath <-
function(object)
{
  checkPtrType(object, "GtkMenu")

  w <- .RGtkCall("S_gtk_menu_get_accel_path", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkMenuGetMonitor <-
function(object)
{
  checkPtrType(object, "GtkMenu")

  w <- .RGtkCall("S_gtk_menu_get_monitor", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkMenuItemGetAccelPath <-
function(object)
{
  checkPtrType(object, "GtkMenuItem")

  w <- .RGtkCall("S_gtk_menu_item_get_accel_path", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkMessageDialogGetImage <-
function(object)
{
  checkPtrType(object, "GtkMessageDialog")

  w <- .RGtkCall("S_gtk_message_dialog_get_image", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkMountOperationGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_mount_operation_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkMountOperationNew <-
function(parent = NULL)
{
  

  w <- .RGtkCall("S_gtk_mount_operation_new", parent, PACKAGE = "RGtk2")

  return(w)
} 


gtkMountOperationIsShowing <-
function(object)
{
  checkPtrType(object, "GtkMountOperation")

  w <- .RGtkCall("S_gtk_mount_operation_is_showing", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkMountOperationSetParent <-
function(object, parent)
{
  checkPtrType(object, "GtkMountOperation")
  checkPtrType(parent, "GtkWindow")

  w <- .RGtkCall("S_gtk_mount_operation_set_parent", object, parent, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMountOperationGetParent <-
function(object)
{
  checkPtrType(object, "GtkMountOperation")

  w <- .RGtkCall("S_gtk_mount_operation_get_parent", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkMountOperationSetScreen <-
function(object, screen)
{
  checkPtrType(object, "GtkMountOperation")
  checkPtrType(screen, "GdkScreen")

  w <- .RGtkCall("S_gtk_mount_operation_set_screen", object, screen, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMountOperationGetScreen <-
function(object)
{
  checkPtrType(object, "GtkMountOperation")

  w <- .RGtkCall("S_gtk_mount_operation_get_screen", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPlugGetEmbedded <-
function(object)
{
  checkPtrType(object, "GtkPlug")

  w <- .RGtkCall("S_gtk_plug_get_embedded", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPlugGetSocketWindow <-
function(object)
{
  checkPtrType(object, "GtkPlug")

  w <- .RGtkCall("S_gtk_plug_get_socket_window", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPageSetupLoadKeyFile <-
function(object, key.file, group.name, .errwarn = TRUE)
{
  checkPtrType(object, "GtkPageSetup")
  checkPtrType(key.file, "GKeyFile")
  group.name <- as.character(group.name)

  w <- .RGtkCall("S_gtk_page_setup_load_key_file", object, key.file, group.name, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkPageSetupLoadFile <-
function(object, file.name, .errwarn = TRUE)
{
  checkPtrType(object, "GtkPageSetup")
  file.name <- as.character(file.name)

  w <- .RGtkCall("S_gtk_page_setup_load_file", object, file.name, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkPrintSettingsLoadKeyFile <-
function(object, key.file, group.name, .errwarn = TRUE)
{
  checkPtrType(object, "GtkPrintSettings")
  checkPtrType(key.file, "GKeyFile")
  group.name <- as.character(group.name)

  w <- .RGtkCall("S_gtk_print_settings_load_key_file", object, key.file, group.name, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkPrintSettingsLoadFile <-
function(object, file.name, .errwarn = TRUE)
{
  checkPtrType(object, "GtkPrintSettings")
  file.name <- as.character(file.name)

  w <- .RGtkCall("S_gtk_print_settings_load_file", object, file.name, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkPrintSettingsGetNumberUpLayout <-
function(object)
{
  checkPtrType(object, "GtkPrintSettings")

  w <- .RGtkCall("S_gtk_print_settings_get_number_up_layout", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsSetNumberUpLayout <-
function(object, number.up.layout)
{
  checkPtrType(object, "GtkPrintSettings")
  

  w <- .RGtkCall("S_gtk_print_settings_set_number_up_layout", object, number.up.layout, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkScaleButtonGetOrientation <-
function(object)
{
  checkPtrType(object, "GtkScaleButton")

  w <- .RGtkCall("S_gtk_scale_button_get_orientation", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkScaleButtonSetOrientation <-
function(object, orientation)
{
  checkPtrType(object, "GtkScaleButton")
  

  w <- .RGtkCall("S_gtk_scale_button_set_orientation", object, orientation, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkScaleButtonGetPlusButton <-
function(object)
{
  checkPtrType(object, "GtkScaleButton")

  w <- .RGtkCall("S_gtk_scale_button_get_plus_button", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkScaleButtonGetMinusButton <-
function(object)
{
  checkPtrType(object, "GtkScaleButton")

  w <- .RGtkCall("S_gtk_scale_button_get_minus_button", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkScaleButtonGetPopup <-
function(object)
{
  checkPtrType(object, "GtkScaleButton")

  w <- .RGtkCall("S_gtk_scale_button_get_popup", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkSelectionDataGetTarget <-
function(object)
{
  checkPtrType(object, "GtkSelectionData")

  w <- .RGtkCall("S_gtk_selection_data_get_target", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkSelectionDataGetDataType <-
function(object)
{
  checkPtrType(object, "GtkSelectionData")

  w <- .RGtkCall("S_gtk_selection_data_get_data_type", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkSelectionDataGetFormat <-
function(object)
{
  checkPtrType(object, "GtkSelectionData")

  w <- .RGtkCall("S_gtk_selection_data_get_format", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkSelectionDataGetData <-
function(object)
{
  checkPtrType(object, "GtkSelectionData")

  w <- .RGtkCall("S_gtk_selection_data_get_data", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkSelectionDataGetLength <-
function(object)
{
  checkPtrType(object, "GtkSelectionData")

  w <- .RGtkCall("S_gtk_selection_data_get_length", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkSelectionDataGetDisplay <-
function(object)
{
  checkPtrType(object, "GtkSelectionData")

  w <- .RGtkCall("S_gtk_selection_data_get_display", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkShowUri <-
function(screen = NULL, uri, timestamp, .errwarn = TRUE)
{
  if (!is.null( screen )) checkPtrType(screen, "GdkScreen")
  uri <- as.character(uri)
  timestamp <- as.numeric(timestamp)

  w <- .RGtkCall("S_gtk_show_uri", screen, uri, timestamp, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gtkSocketGetPlugWindow <-
function(object)
{
  checkPtrType(object, "GtkSocket")

  w <- .RGtkCall("S_gtk_socket_get_plug_window", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkStatusIconNewFromGicon <-
function(icon)
{
  checkPtrType(icon, "GIcon")

  w <- .RGtkCall("S_gtk_status_icon_new_from_gicon", icon, PACKAGE = "RGtk2")

  return(w)
} 


gtkStatusIconGetX11WindowId <-
function(object)
{
  checkPtrType(object, "GtkStatusIcon")

  w <- .RGtkCall("S_gtk_status_icon_get_x11_window_id", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkStatusIconGetGicon <-
function(object)
{
  checkPtrType(object, "GtkStatusIcon")

  w <- .RGtkCall("S_gtk_status_icon_get_gicon", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkStatusIconSetFromGicon <-
function(object, icon)
{
  checkPtrType(object, "GtkStatusIcon")
  checkPtrType(icon, "GIcon")

  w <- .RGtkCall("S_gtk_status_icon_set_from_gicon", object, icon, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTooltipSetIconFromIconName <-
function(object, icon.name = NULL, size)
{
  checkPtrType(object, "GtkTooltip")
  if (!is.null( icon.name )) icon.name <- as.character(icon.name)
  

  w <- .RGtkCall("S_gtk_tooltip_set_icon_from_icon_name", object, icon.name, size, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolItemToolbarReconfigured <-
function(object)
{
  checkPtrType(object, "GtkToolItem")

  w <- .RGtkCall("S_gtk_tool_item_toolbar_reconfigured", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolShellGetIconSize <-
function(object)
{
  checkPtrType(object, "GtkToolShell")

  w <- .RGtkCall("S_gtk_tool_shell_get_icon_size", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolShellGetOrientation <-
function(object)
{
  checkPtrType(object, "GtkToolShell")

  w <- .RGtkCall("S_gtk_tool_shell_get_orientation", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolShellGetStyle <-
function(object)
{
  checkPtrType(object, "GtkToolShell")

  w <- .RGtkCall("S_gtk_tool_shell_get_style", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolShellGetReliefStyle <-
function(object)
{
  checkPtrType(object, "GtkToolShell")

  w <- .RGtkCall("S_gtk_tool_shell_get_relief_style", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolShellRebuildMenu <-
function(object)
{
  checkPtrType(object, "GtkToolShell")

  w <- .RGtkCall("S_gtk_tool_shell_rebuild_menu", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkTreeSelectionGetSelectFunction <-
function(object)
{
  checkPtrType(object, "GtkTreeSelection")

  w <- .RGtkCall("S_gtk_tree_selection_get_select_function", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetGetSnapshot <-
function(object, clip.rect = NULL)
{
  checkPtrType(object, "GtkWidget")
  if (!is.null( clip.rect )) clip.rect <- as.GdkRectangle(clip.rect)

  w <- .RGtkCall("S_gtk_widget_get_snapshot", object, clip.rect, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetGetWindow <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_window", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowGetDefaultWidget <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_get_default_widget", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowGroupListWindows <-
function(object)
{
  checkPtrType(object, "GtkWindowGroup")

  w <- .RGtkCall("S_gtk_window_group_list_windows", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkLinkButtonGetVisited <-
function(object)
{
  checkPtrType(object, "GtkLinkButton")

  w <- .RGtkCall("S_gtk_link_button_get_visited", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkLinkButtonSetVisited <-
function(object, visited)
{
  checkPtrType(object, "GtkLinkButton")
  visited <- as.logical(visited)

  w <- .RGtkCall("S_gtk_link_button_set_visited", object, visited, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkBorderNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_border_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkTestRegisterAllTypes <-
function()
{
  

  w <- .RGtkCall("S_gtk_test_register_all_types", PACKAGE = "RGtk2")

  return(w)
} 


gtkTestListAllTypes <-
function()
{
  

  w <- .RGtkCall("S_gtk_test_list_all_types", PACKAGE = "RGtk2")

  return(w)
} 


gtkTestFindWidget <-
function(widget, label.pattern, widget.type)
{
  checkPtrType(widget, "GtkWidget")
  label.pattern <- as.character(label.pattern)
  widget.type <- as.GType(widget.type)

  w <- .RGtkCall("S_gtk_test_find_widget", widget, label.pattern, widget.type, PACKAGE = "RGtk2")

  return(w)
} 


gtkTestSliderSetPerc <-
function(widget, percentage)
{
  checkPtrType(widget, "GtkWidget")
  percentage <- as.numeric(percentage)

  w <- .RGtkCall("S_gtk_test_slider_set_perc", widget, percentage, PACKAGE = "RGtk2")

  return(w)
} 


gtkTestSliderGetValue <-
function(widget)
{
  checkPtrType(widget, "GtkWidget")

  w <- .RGtkCall("S_gtk_test_slider_get_value", widget, PACKAGE = "RGtk2")

  return(w)
} 


gtkTestSpinButtonClick <-
function(spinner, button, upwards)
{
  checkPtrType(spinner, "GtkSpinButton")
  button <- as.numeric(button)
  upwards <- as.logical(upwards)

  w <- .RGtkCall("S_gtk_test_spin_button_click", spinner, button, upwards, PACKAGE = "RGtk2")

  return(w)
} 


gtkTestWidgetClick <-
function(widget, button, modifiers)
{
  checkPtrType(widget, "GtkWidget")
  button <- as.numeric(button)
  

  w <- .RGtkCall("S_gtk_test_widget_click", widget, button, modifiers, PACKAGE = "RGtk2")

  return(w)
} 


gtkTestWidgetSendKey <-
function(widget, keyval, modifiers)
{
  checkPtrType(widget, "GtkWidget")
  keyval <- as.numeric(keyval)
  

  w <- .RGtkCall("S_gtk_test_widget_send_key", widget, keyval, modifiers, PACKAGE = "RGtk2")

  return(w)
} 


gtkTestTextSet <-
function(widget, string)
{
  checkPtrType(widget, "GtkWidget")
  string <- as.character(string)

  w <- .RGtkCall("S_gtk_test_text_set", widget, string, PACKAGE = "RGtk2")

  return(w)
} 


gtkTestTextGet <-
function(widget)
{
  checkPtrType(widget, "GtkWidget")

  w <- .RGtkCall("S_gtk_test_text_get", widget, PACKAGE = "RGtk2")

  return(w)
} 


gtkTestFindSibling <-
function(base.widget, widget.type)
{
  checkPtrType(base.widget, "GtkWidget")
  widget.type <- as.GType(widget.type)

  w <- .RGtkCall("S_gtk_test_find_sibling", base.widget, widget.type, PACKAGE = "RGtk2")

  return(w)
} 


gtkTestFindLabel <-
function(widget, label.pattern)
{
  checkPtrType(widget, "GtkWidget")
  label.pattern <- as.character(label.pattern)

  w <- .RGtkCall("S_gtk_test_find_label", widget, label.pattern, PACKAGE = "RGtk2")

  return(w)
} 


gtkActionBlockActivate <-
function(object)
{
  checkPtrType(object, "GtkAction")

  w <- .RGtkCall("S_gtk_action_block_activate", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActionUnblockActivate <-
function(object)
{
  checkPtrType(object, "GtkAction")

  w <- .RGtkCall("S_gtk_action_unblock_activate", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActionSetGicon <-
function(object, icon)
{
  checkPtrType(object, "GtkAction")
  checkPtrType(icon, "GIcon")

  w <- .RGtkCall("S_gtk_action_set_gicon", object, icon, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActionGetGicon <-
function(object)
{
  checkPtrType(object, "GtkAction")

  w <- .RGtkCall("S_gtk_action_get_gicon", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkActionSetIconName <-
function(object, icon.name)
{
  checkPtrType(object, "GtkAction")
  icon.name <- as.character(icon.name)

  w <- .RGtkCall("S_gtk_action_set_icon_name", object, icon.name, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActionGetIconName <-
function(object)
{
  checkPtrType(object, "GtkAction")

  w <- .RGtkCall("S_gtk_action_get_icon_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkActionSetVisibleHorizontal <-
function(object, visible.horizontal)
{
  checkPtrType(object, "GtkAction")
  visible.horizontal <- as.logical(visible.horizontal)

  w <- .RGtkCall("S_gtk_action_set_visible_horizontal", object, visible.horizontal, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActionGetVisibleHorizontal <-
function(object)
{
  checkPtrType(object, "GtkAction")

  w <- .RGtkCall("S_gtk_action_get_visible_horizontal", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkActionSetVisibleVertical <-
function(object, visible.vertical)
{
  checkPtrType(object, "GtkAction")
  visible.vertical <- as.logical(visible.vertical)

  w <- .RGtkCall("S_gtk_action_set_visible_vertical", object, visible.vertical, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActionGetVisibleVertical <-
function(object)
{
  checkPtrType(object, "GtkAction")

  w <- .RGtkCall("S_gtk_action_get_visible_vertical", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkActionSetIsImportant <-
function(object, is.important)
{
  checkPtrType(object, "GtkAction")
  is.important <- as.logical(is.important)

  w <- .RGtkCall("S_gtk_action_set_is_important", object, is.important, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActionGetIsImportant <-
function(object)
{
  checkPtrType(object, "GtkAction")

  w <- .RGtkCall("S_gtk_action_get_is_important", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkActionSetLabel <-
function(object, label)
{
  checkPtrType(object, "GtkAction")
  label <- as.character(label)

  w <- .RGtkCall("S_gtk_action_set_label", object, label, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActionGetLabel <-
function(object)
{
  checkPtrType(object, "GtkAction")

  w <- .RGtkCall("S_gtk_action_get_label", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkActionSetShortLabel <-
function(object, short.label)
{
  checkPtrType(object, "GtkAction")
  short.label <- as.character(short.label)

  w <- .RGtkCall("S_gtk_action_set_short_label", object, short.label, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActionGetShortLabel <-
function(object)
{
  checkPtrType(object, "GtkAction")

  w <- .RGtkCall("S_gtk_action_get_short_label", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkActionSetTooltip <-
function(object, tooltip)
{
  checkPtrType(object, "GtkAction")
  tooltip <- as.character(tooltip)

  w <- .RGtkCall("S_gtk_action_set_tooltip", object, tooltip, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActionGetTooltip <-
function(object)
{
  checkPtrType(object, "GtkAction")

  w <- .RGtkCall("S_gtk_action_get_tooltip", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkActionSetStockId <-
function(object, stock.id)
{
  checkPtrType(object, "GtkAction")
  stock.id <- as.character(stock.id)

  w <- .RGtkCall("S_gtk_action_set_stock_id", object, stock.id, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActionGetStockId <-
function(object)
{
  checkPtrType(object, "GtkAction")

  w <- .RGtkCall("S_gtk_action_get_stock_id", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkActivatableGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_activatable_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkActivatableSyncActionProperties <-
function(object, action = NULL)
{
  checkPtrType(object, "GtkActivatable")
  if (!is.null( action )) checkPtrType(action, "GtkAction")

  w <- .RGtkCall("S_gtk_activatable_sync_action_properties", object, action, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActivatableSetRelatedAction <-
function(object, action)
{
  checkPtrType(object, "GtkActivatable")
  checkPtrType(action, "GtkAction")

  w <- .RGtkCall("S_gtk_activatable_set_related_action", object, action, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActivatableGetRelatedAction <-
function(object)
{
  checkPtrType(object, "GtkActivatable")

  w <- .RGtkCall("S_gtk_activatable_get_related_action", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkActivatableSetUseActionAppearance <-
function(object, use.appearance)
{
  checkPtrType(object, "GtkActivatable")
  use.appearance <- as.logical(use.appearance)

  w <- .RGtkCall("S_gtk_activatable_set_use_action_appearance", object, use.appearance, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkActivatableGetUseActionAppearance <-
function(object)
{
  checkPtrType(object, "GtkActivatable")

  w <- .RGtkCall("S_gtk_activatable_get_use_action_appearance", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkActivatableDoSetRelatedAction <-
function(object, action)
{
  checkPtrType(object, "GtkActivatable")
  checkPtrType(action, "GtkAction")

  w <- .RGtkCall("S_gtk_activatable_do_set_related_action", object, action, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCellViewGetModel <-
function(object)
{
  checkPtrType(object, "GtkCellView")

  w <- .RGtkCall("S_gtk_cell_view_get_model", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntrySetProgressFraction <-
function(object, fraction)
{
  checkPtrType(object, "GtkEntry")
  fraction <- as.numeric(fraction)

  w <- .RGtkCall("S_gtk_entry_set_progress_fraction", object, fraction, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryGetProgressFraction <-
function(object)
{
  checkPtrType(object, "GtkEntry")

  w <- .RGtkCall("S_gtk_entry_get_progress_fraction", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntrySetProgressPulseStep <-
function(object, fraction)
{
  checkPtrType(object, "GtkEntry")
  fraction <- as.numeric(fraction)

  w <- .RGtkCall("S_gtk_entry_set_progress_pulse_step", object, fraction, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryGetProgressPulseStep <-
function(object)
{
  checkPtrType(object, "GtkEntry")

  w <- .RGtkCall("S_gtk_entry_get_progress_pulse_step", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryProgressPulse <-
function(object)
{
  checkPtrType(object, "GtkEntry")

  w <- .RGtkCall("S_gtk_entry_progress_pulse", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntrySetIconFromPixbuf <-
function(object, icon.pos, pixbuf = NULL)
{
  checkPtrType(object, "GtkEntry")
  
  if (!is.null( pixbuf )) checkPtrType(pixbuf, "GdkPixbuf")

  w <- .RGtkCall("S_gtk_entry_set_icon_from_pixbuf", object, icon.pos, pixbuf, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntrySetIconFromStock <-
function(object, icon.pos, stock.id = NULL)
{
  checkPtrType(object, "GtkEntry")
  
  if (!is.null( stock.id )) stock.id <- as.character(stock.id)

  w <- .RGtkCall("S_gtk_entry_set_icon_from_stock", object, icon.pos, stock.id, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntrySetIconFromIconName <-
function(object, icon.pos, icon.name = NULL)
{
  checkPtrType(object, "GtkEntry")
  
  if (!is.null( icon.name )) icon.name <- as.character(icon.name)

  w <- .RGtkCall("S_gtk_entry_set_icon_from_icon_name", object, icon.pos, icon.name, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntrySetIconFromGicon <-
function(object, icon.pos, icon = NULL)
{
  checkPtrType(object, "GtkEntry")
  
  if (!is.null( icon )) checkPtrType(icon, "GIcon")

  w <- .RGtkCall("S_gtk_entry_set_icon_from_gicon", object, icon.pos, icon, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryGetIconStorageType <-
function(object, icon.pos)
{
  checkPtrType(object, "GtkEntry")
  

  w <- .RGtkCall("S_gtk_entry_get_icon_storage_type", object, icon.pos, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryGetIconPixbuf <-
function(object, icon.pos)
{
  checkPtrType(object, "GtkEntry")
  

  w <- .RGtkCall("S_gtk_entry_get_icon_pixbuf", object, icon.pos, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryGetIconStock <-
function(object, icon.pos)
{
  checkPtrType(object, "GtkEntry")
  

  w <- .RGtkCall("S_gtk_entry_get_icon_stock", object, icon.pos, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryGetIconName <-
function(object, icon.pos)
{
  checkPtrType(object, "GtkEntry")
  

  w <- .RGtkCall("S_gtk_entry_get_icon_name", object, icon.pos, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryGetIconGicon <-
function(object, icon.pos)
{
  checkPtrType(object, "GtkEntry")
  

  w <- .RGtkCall("S_gtk_entry_get_icon_gicon", object, icon.pos, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntrySetIconActivatable <-
function(object, icon.pos, activatable)
{
  checkPtrType(object, "GtkEntry")
  
  activatable <- as.logical(activatable)

  w <- .RGtkCall("S_gtk_entry_set_icon_activatable", object, icon.pos, activatable, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryGetIconActivatable <-
function(object, icon.pos)
{
  checkPtrType(object, "GtkEntry")
  

  w <- .RGtkCall("S_gtk_entry_get_icon_activatable", object, icon.pos, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntrySetIconSensitive <-
function(object, icon.pos, sensitive)
{
  checkPtrType(object, "GtkEntry")
  
  sensitive <- as.logical(sensitive)

  w <- .RGtkCall("S_gtk_entry_set_icon_sensitive", object, icon.pos, sensitive, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryGetIconSensitive <-
function(object, icon.pos)
{
  checkPtrType(object, "GtkEntry")
  

  w <- .RGtkCall("S_gtk_entry_get_icon_sensitive", object, icon.pos, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryGetIconAtPos <-
function(object, x, y)
{
  checkPtrType(object, "GtkEntry")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_gtk_entry_get_icon_at_pos", object, x, y, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntrySetIconTooltipText <-
function(object, icon.pos, tooltip = NULL)
{
  checkPtrType(object, "GtkEntry")
  
  if (!is.null( tooltip )) tooltip <- as.character(tooltip)

  w <- .RGtkCall("S_gtk_entry_set_icon_tooltip_text", object, icon.pos, tooltip, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryGetIconTooltipText <-
function(object, icon.pos)
{
  checkPtrType(object, "GtkEntry")
  

  w <- .RGtkCall("S_gtk_entry_get_icon_tooltip_text", object, icon.pos, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntrySetIconTooltipMarkup <-
function(object, icon.pos, tooltip = NULL)
{
  checkPtrType(object, "GtkEntry")
  
  if (!is.null( tooltip )) tooltip <- as.character(tooltip)

  w <- .RGtkCall("S_gtk_entry_set_icon_tooltip_markup", object, icon.pos, tooltip, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryGetIconTooltipMarkup <-
function(object, icon.pos)
{
  checkPtrType(object, "GtkEntry")
  

  w <- .RGtkCall("S_gtk_entry_get_icon_tooltip_markup", object, icon.pos, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryUnsetInvisibleChar <-
function(object)
{
  checkPtrType(object, "GtkEntry")

  w <- .RGtkCall("S_gtk_entry_unset_invisible_char", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntrySetIconDragSource <-
function(object, icon.pos, target.list, actions)
{
  checkPtrType(object, "GtkEntry")
  
  checkPtrType(target.list, "GtkTargetList")
  

  w <- .RGtkCall("S_gtk_entry_set_icon_drag_source", object, icon.pos, target.list, actions, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryGetCurrentIconDragSource <-
function(object)
{
  checkPtrType(object, "GtkEntry")

  w <- .RGtkCall("S_gtk_entry_get_current_icon_drag_source", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkImageMenuItemSetAlwaysShowImage <-
function(object, always.show)
{
  checkPtrType(object, "GtkImageMenuItem")
  always.show <- as.logical(always.show)

  w <- .RGtkCall("S_gtk_image_menu_item_set_always_show_image", object, always.show, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkImageMenuItemGetAlwaysShowImage <-
function(object)
{
  checkPtrType(object, "GtkImageMenuItem")

  w <- .RGtkCall("S_gtk_image_menu_item_get_always_show_image", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkImageMenuItemSetUseStock <-
function(object, use.stock)
{
  checkPtrType(object, "GtkImageMenuItem")
  use.stock <- as.logical(use.stock)

  w <- .RGtkCall("S_gtk_image_menu_item_set_use_stock", object, use.stock, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkImageMenuItemGetUseStock <-
function(object)
{
  checkPtrType(object, "GtkImageMenuItem")

  w <- .RGtkCall("S_gtk_image_menu_item_get_use_stock", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkImageMenuItemSetAccelGroup <-
function(object, accel.group)
{
  checkPtrType(object, "GtkImageMenuItem")
  checkPtrType(accel.group, "GtkAccelGroup")

  w <- .RGtkCall("S_gtk_image_menu_item_set_accel_group", object, accel.group, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIMMulticontextGetContextId <-
function(object)
{
  checkPtrType(object, "GtkIMMulticontext")

  w <- .RGtkCall("S_gtk_im_multicontext_get_context_id", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIMMulticontextSetContextId <-
function(object, context.id)
{
  checkPtrType(object, "GtkIMMulticontext")
  context.id <- as.character(context.id)

  w <- .RGtkCall("S_gtk_im_multicontext_set_context_id", object, context.id, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuItemSetLabel <-
function(object, label)
{
  checkPtrType(object, "GtkMenuItem")
  label <- as.character(label)

  w <- .RGtkCall("S_gtk_menu_item_set_label", object, label, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuItemGetLabel <-
function(object)
{
  checkPtrType(object, "GtkMenuItem")

  w <- .RGtkCall("S_gtk_menu_item_get_label", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkMenuItemSetUseUnderline <-
function(object, setting)
{
  checkPtrType(object, "GtkMenuItem")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_menu_item_set_use_underline", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuItemGetUseUnderline <-
function(object)
{
  checkPtrType(object, "GtkMenuItem")

  w <- .RGtkCall("S_gtk_menu_item_get_use_underline", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkOrientableGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_orientable_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkOrientableSetOrientation <-
function(object, orientation)
{
  checkPtrType(object, "GtkOrientable")
  

  w <- .RGtkCall("S_gtk_orientable_set_orientation", object, orientation, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkOrientableGetOrientation <-
function(object)
{
  checkPtrType(object, "GtkOrientable")

  w <- .RGtkCall("S_gtk_orientable_get_orientation", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintOperationDrawPageFinish <-
function(object)
{
  checkPtrType(object, "GtkPrintOperation")

  w <- .RGtkCall("S_gtk_print_operation_draw_page_finish", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintOperationSetDeferDrawing <-
function(object)
{
  checkPtrType(object, "GtkPrintOperation")

  w <- .RGtkCall("S_gtk_print_operation_set_defer_drawing", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintSettingsGetResolutionX <-
function(object)
{
  checkPtrType(object, "GtkPrintSettings")

  w <- .RGtkCall("S_gtk_print_settings_get_resolution_x", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsGetResolutionY <-
function(object)
{
  checkPtrType(object, "GtkPrintSettings")

  w <- .RGtkCall("S_gtk_print_settings_get_resolution_y", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsSetResolutionXy <-
function(object, resolution.x, resolution.y)
{
  checkPtrType(object, "GtkPrintSettings")
  resolution.x <- as.integer(resolution.x)
  resolution.y <- as.integer(resolution.y)

  w <- .RGtkCall("S_gtk_print_settings_set_resolution_xy", object, resolution.x, resolution.y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintSettingsGetPrinterLpi <-
function(object)
{
  checkPtrType(object, "GtkPrintSettings")

  w <- .RGtkCall("S_gtk_print_settings_get_printer_lpi", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintSettingsSetPrinterLpi <-
function(object, lpi)
{
  checkPtrType(object, "GtkPrintSettings")
  lpi <- as.numeric(lpi)

  w <- .RGtkCall("S_gtk_print_settings_set_printer_lpi", object, lpi, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkScaleAddMark <-
function(object, value, position, markup = NULL)
{
  checkPtrType(object, "GtkScale")
  value <- as.numeric(value)
  
  if (!is.null( markup )) markup <- as.character(markup)

  w <- .RGtkCall("S_gtk_scale_add_mark", object, value, position, markup, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkScaleClearMarks <-
function(object)
{
  checkPtrType(object, "GtkScale")

  w <- .RGtkCall("S_gtk_scale_clear_marks", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkSelectionDataGetSelection <-
function(object)
{
  checkPtrType(object, "GtkSelectionData")

  w <- .RGtkCall("S_gtk_selection_data_get_selection", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkStatusIconGetHasTooltip <-
function(object)
{
  checkPtrType(object, "GtkStatusIcon")

  w <- .RGtkCall("S_gtk_status_icon_get_has_tooltip", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkStatusIconGetTooltipText <-
function(object)
{
  checkPtrType(object, "GtkStatusIcon")

  w <- .RGtkCall("S_gtk_status_icon_get_tooltip_text", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkStatusIconGetTooltipMarkup <-
function(object)
{
  checkPtrType(object, "GtkStatusIcon")

  w <- .RGtkCall("S_gtk_status_icon_get_tooltip_markup", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkStatusIconSetHasTooltip <-
function(object, has.tooltip)
{
  checkPtrType(object, "GtkStatusIcon")
  has.tooltip <- as.logical(has.tooltip)

  w <- .RGtkCall("S_gtk_status_icon_set_has_tooltip", object, has.tooltip, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkStatusIconSetTooltipText <-
function(object, text)
{
  checkPtrType(object, "GtkStatusIcon")
  text <- as.character(text)

  w <- .RGtkCall("S_gtk_status_icon_set_tooltip_text", object, text, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkStatusIconSetTooltipMarkup <-
function(object, markup = NULL)
{
  checkPtrType(object, "GtkStatusIcon")
  if (!is.null( markup )) markup <- as.character(markup)

  w <- .RGtkCall("S_gtk_status_icon_set_tooltip_markup", object, markup, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkStyleGetStyleProperty <-
function(object, widget.type, property.name)
{
  checkPtrType(object, "GtkStyle")
  widget.type <- as.GType(widget.type)
  property.name <- as.character(property.name)

  w <- .RGtkCall("S_gtk_style_get_style_property", object, widget.type, property.name, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowGetDefaultIconName <-
function()
{
  

  w <- .RGtkCall("S_gtk_window_get_default_icon_name", PACKAGE = "RGtk2")

  return(w)
} 


gtkCellRendererSetAlignment <-
function(object, xalign, yalign)
{
  checkPtrType(object, "GtkCellRenderer")
  xalign <- as.numeric(xalign)
  yalign <- as.numeric(yalign)

  w <- .RGtkCall("S_gtk_cell_renderer_set_alignment", object, xalign, yalign, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCellRendererGetAlignment <-
function(object, xalign, yalign)
{
  checkPtrType(object, "GtkCellRenderer")
  xalign <- as.list(as.numeric(xalign))
  yalign <- as.list(as.numeric(yalign))

  w <- .RGtkCall("S_gtk_cell_renderer_get_alignment", object, xalign, yalign, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCellRendererSetPadding <-
function(object, xpad, ypad)
{
  checkPtrType(object, "GtkCellRenderer")
  xpad <- as.integer(xpad)
  ypad <- as.integer(ypad)

  w <- .RGtkCall("S_gtk_cell_renderer_set_padding", object, xpad, ypad, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCellRendererGetPadding <-
function(object, xpad, ypad)
{
  checkPtrType(object, "GtkCellRenderer")
  xpad <- as.list(as.integer(xpad))
  ypad <- as.list(as.integer(ypad))

  w <- .RGtkCall("S_gtk_cell_renderer_get_padding", object, xpad, ypad, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCellRendererSetVisible <-
function(object, visible)
{
  checkPtrType(object, "GtkCellRenderer")
  visible <- as.logical(visible)

  w <- .RGtkCall("S_gtk_cell_renderer_set_visible", object, visible, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCellRendererGetVisible <-
function(object)
{
  checkPtrType(object, "GtkCellRenderer")

  w <- .RGtkCall("S_gtk_cell_renderer_get_visible", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkCellRendererSetSensitive <-
function(object, sensitive)
{
  checkPtrType(object, "GtkCellRenderer")
  sensitive <- as.logical(sensitive)

  w <- .RGtkCall("S_gtk_cell_renderer_set_sensitive", object, sensitive, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCellRendererGetSensitive <-
function(object)
{
  checkPtrType(object, "GtkCellRenderer")

  w <- .RGtkCall("S_gtk_cell_renderer_get_sensitive", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkCellRendererToggleGetActivatable <-
function(object)
{
  checkPtrType(object, "GtkCellRendererToggle")

  w <- .RGtkCall("S_gtk_cell_renderer_toggle_get_activatable", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkCellRendererToggleSetActivatable <-
function(object, setting)
{
  checkPtrType(object, "GtkCellRendererToggle")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_cell_renderer_toggle_set_activatable", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryNewWithBuffer <-
function(buffer, show = TRUE)
{
  checkPtrType(buffer, "GtkEntryBuffer")

  w <- .RGtkCall("S_gtk_entry_new_with_buffer", buffer, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkEntryGetBuffer <-
function(object)
{
  checkPtrType(object, "GtkEntry")

  w <- .RGtkCall("S_gtk_entry_get_buffer", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntrySetBuffer <-
function(object, buffer)
{
  checkPtrType(object, "GtkEntry")
  checkPtrType(buffer, "GtkEntryBuffer")

  w <- .RGtkCall("S_gtk_entry_set_buffer", object, buffer, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryBufferGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_entry_buffer_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryBufferNew <-
function(initial.chars = NULL, n.initial.chars = -1)
{
  if (!is.null( initial.chars )) initial.chars <- as.character(initial.chars)
  n.initial.chars <- as.integer(n.initial.chars)

  w <- .RGtkCall("S_gtk_entry_buffer_new", initial.chars, n.initial.chars, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryBufferGetBytes <-
function(object)
{
  checkPtrType(object, "GtkEntryBuffer")

  w <- .RGtkCall("S_gtk_entry_buffer_get_bytes", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryBufferGetLength <-
function(object)
{
  checkPtrType(object, "GtkEntryBuffer")

  w <- .RGtkCall("S_gtk_entry_buffer_get_length", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryBufferGetText <-
function(object)
{
  checkPtrType(object, "GtkEntryBuffer")

  w <- .RGtkCall("S_gtk_entry_buffer_get_text", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryBufferSetText <-
function(object, chars, n.chars)
{
  checkPtrType(object, "GtkEntryBuffer")
  chars <- as.character(chars)
  n.chars <- as.integer(n.chars)

  w <- .RGtkCall("S_gtk_entry_buffer_set_text", object, chars, n.chars, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryBufferSetMaxLength <-
function(object, max.length)
{
  checkPtrType(object, "GtkEntryBuffer")
  max.length <- as.integer(max.length)

  w <- .RGtkCall("S_gtk_entry_buffer_set_max_length", object, max.length, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryBufferGetMaxLength <-
function(object)
{
  checkPtrType(object, "GtkEntryBuffer")

  w <- .RGtkCall("S_gtk_entry_buffer_get_max_length", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryBufferInsertText <-
function(object, position, chars, n.chars)
{
  checkPtrType(object, "GtkEntryBuffer")
  position <- as.numeric(position)
  chars <- as.character(chars)
  n.chars <- as.integer(n.chars)

  w <- .RGtkCall("S_gtk_entry_buffer_insert_text", object, position, chars, n.chars, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryBufferDeleteText <-
function(object, position, n.chars)
{
  checkPtrType(object, "GtkEntryBuffer")
  position <- as.numeric(position)
  n.chars <- as.integer(n.chars)

  w <- .RGtkCall("S_gtk_entry_buffer_delete_text", object, position, n.chars, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryBufferEmitInsertedText <-
function(object, position, chars, n.chars)
{
  checkPtrType(object, "GtkEntryBuffer")
  position <- as.numeric(position)
  chars <- as.character(chars)
  n.chars <- as.numeric(n.chars)

  w <- .RGtkCall("S_gtk_entry_buffer_emit_inserted_text", object, position, chars, n.chars, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkEntryBufferEmitDeletedText <-
function(object, position, n.chars)
{
  checkPtrType(object, "GtkEntryBuffer")
  position <- as.numeric(position)
  n.chars <- as.numeric(n.chars)

  w <- .RGtkCall("S_gtk_entry_buffer_emit_deleted_text", object, position, n.chars, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFileChooserSetCreateFolders <-
function(object, create.folders)
{
  checkPtrType(object, "GtkFileChooser")
  create.folders <- as.logical(create.folders)

  w <- .RGtkCall("S_gtk_file_chooser_set_create_folders", object, create.folders, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkFileChooserGetCreateFolders <-
function(object)
{
  checkPtrType(object, "GtkFileChooser")

  w <- .RGtkCall("S_gtk_file_chooser_get_create_folders", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkIconViewSetItemPadding <-
function(object, item.padding)
{
  checkPtrType(object, "GtkIconView")
  item.padding <- as.integer(item.padding)

  w <- .RGtkCall("S_gtk_icon_view_set_item_padding", object, item.padding, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkIconViewGetItemPadding <-
function(object)
{
  checkPtrType(object, "GtkIconView")

  w <- .RGtkCall("S_gtk_icon_view_get_item_padding", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkInfoBarGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_info_bar_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkInfoBarNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_info_bar_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkInfoBarGetActionArea <-
function(object)
{
  checkPtrType(object, "GtkInfoBar")

  w <- .RGtkCall("S_gtk_info_bar_get_action_area", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkInfoBarGetContentArea <-
function(object)
{
  checkPtrType(object, "GtkInfoBar")

  w <- .RGtkCall("S_gtk_info_bar_get_content_area", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkInfoBarAddActionWidget <-
function(object, child, response.id)
{
  checkPtrType(object, "GtkInfoBar")
  checkPtrType(child, "GtkWidget")
  response.id <- as.integer(response.id)

  w <- .RGtkCall("S_gtk_info_bar_add_action_widget", object, child, response.id, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkInfoBarAddButton <-
function(object, button.text, response.id)
{
  checkPtrType(object, "GtkInfoBar")
  button.text <- as.character(button.text)
  response.id <- as.integer(response.id)

  w <- .RGtkCall("S_gtk_info_bar_add_button", object, button.text, response.id, PACKAGE = "RGtk2")

  return(w)
} 


gtkInfoBarSetResponseSensitive <-
function(object, response.id, setting)
{
  checkPtrType(object, "GtkInfoBar")
  response.id <- as.integer(response.id)
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_info_bar_set_response_sensitive", object, response.id, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkInfoBarSetDefaultResponse <-
function(object, response.id)
{
  checkPtrType(object, "GtkInfoBar")
  response.id <- as.integer(response.id)

  w <- .RGtkCall("S_gtk_info_bar_set_default_response", object, response.id, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkInfoBarResponse <-
function(object, response.id)
{
  checkPtrType(object, "GtkInfoBar")
  response.id <- as.integer(response.id)

  w <- .RGtkCall("S_gtk_info_bar_response", object, response.id, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkInfoBarSetMessageType <-
function(object, message.type)
{
  checkPtrType(object, "GtkInfoBar")
  

  w <- .RGtkCall("S_gtk_info_bar_set_message_type", object, message.type, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkInfoBarGetMessageType <-
function(object)
{
  checkPtrType(object, "GtkInfoBar")

  w <- .RGtkCall("S_gtk_info_bar_get_message_type", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkLabelGetCurrentUri <-
function(object)
{
  checkPtrType(object, "GtkLabel")

  w <- .RGtkCall("S_gtk_label_get_current_uri", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkLabelSetTrackVisitedLinks <-
function(object, track.links)
{
  checkPtrType(object, "GtkLabel")
  track.links <- as.logical(track.links)

  w <- .RGtkCall("S_gtk_label_set_track_visited_links", object, track.links, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkLabelGetTrackVisitedLinks <-
function(object)
{
  checkPtrType(object, "GtkLabel")

  w <- .RGtkCall("S_gtk_label_get_track_visited_links", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkMenuSetReserveToggleSize <-
function(object, reserve.toggle.size)
{
  checkPtrType(object, "GtkMenu")
  reserve.toggle.size <- as.logical(reserve.toggle.size)

  w <- .RGtkCall("S_gtk_menu_set_reserve_toggle_size", object, reserve.toggle.size, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkMenuGetReserveToggleSize <-
function(object)
{
  checkPtrType(object, "GtkMenu")

  w <- .RGtkCall("S_gtk_menu_get_reserve_toggle_size", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintOperationSetSupportSelection <-
function(object, support.selection)
{
  checkPtrType(object, "GtkPrintOperation")
  support.selection <- as.logical(support.selection)

  w <- .RGtkCall("S_gtk_print_operation_set_support_selection", object, support.selection, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintOperationGetSupportSelection <-
function(object)
{
  checkPtrType(object, "GtkPrintOperation")

  w <- .RGtkCall("S_gtk_print_operation_get_support_selection", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintOperationSetHasSelection <-
function(object, has.selection)
{
  checkPtrType(object, "GtkPrintOperation")
  has.selection <- as.logical(has.selection)

  w <- .RGtkCall("S_gtk_print_operation_set_has_selection", object, has.selection, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintOperationGetHasSelection <-
function(object)
{
  checkPtrType(object, "GtkPrintOperation")

  w <- .RGtkCall("S_gtk_print_operation_get_has_selection", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintOperationSetEmbedPageSetup <-
function(object, embed)
{
  checkPtrType(object, "GtkPrintOperation")
  embed <- as.logical(embed)

  w <- .RGtkCall("S_gtk_print_operation_set_embed_page_setup", object, embed, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintOperationGetEmbedPageSetup <-
function(object)
{
  checkPtrType(object, "GtkPrintOperation")

  w <- .RGtkCall("S_gtk_print_operation_get_embed_page_setup", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkPrintOperationGetNPagesToPrint <-
function(object)
{
  checkPtrType(object, "GtkPrintOperation")

  w <- .RGtkCall("S_gtk_print_operation_get_n_pages_to_print", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRangeSetFlippable <-
function(object, flippable)
{
  checkPtrType(object, "GtkRange")
  flippable <- as.logical(flippable)

  w <- .RGtkCall("S_gtk_range_set_flippable", object, flippable, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRangeGetFlippable <-
function(object)
{
  checkPtrType(object, "GtkRange")

  w <- .RGtkCall("S_gtk_range_get_flippable", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkStatusIconSetTitle <-
function(object, title)
{
  checkPtrType(object, "GtkStatusIcon")
  title <- as.character(title)

  w <- .RGtkCall("S_gtk_status_icon_set_title", object, title, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkStatusIconGetTitle <-
function(object)
{
  checkPtrType(object, "GtkStatusIcon")

  w <- .RGtkCall("S_gtk_status_icon_get_title", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetGetAllocation <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_allocation", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetSetAllocation <-
function(object, allocation)
{
  checkPtrType(object, "GtkWidget")
  allocation <- as.GtkAllocation(allocation)

  w <- .RGtkCall("S_gtk_widget_set_allocation", object, allocation, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetGetAppPaintable <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_app_paintable", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetSetCanDefault <-
function(object, can.default)
{
  checkPtrType(object, "GtkWidget")
  can.default <- as.logical(can.default)

  w <- .RGtkCall("S_gtk_widget_set_can_default", object, can.default, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetGetCanDefault <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_can_default", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetSetCanFocus <-
function(object, can.focus)
{
  checkPtrType(object, "GtkWidget")
  can.focus <- as.logical(can.focus)

  w <- .RGtkCall("S_gtk_widget_set_can_focus", object, can.focus, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetGetCanFocus <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_can_focus", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetGetDoubleBuffered <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_double_buffered", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetSetHasWindow <-
function(object, has.window)
{
  checkPtrType(object, "GtkWidget")
  has.window <- as.logical(has.window)

  w <- .RGtkCall("S_gtk_widget_set_has_window", object, has.window, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetGetHasWindow <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_has_window", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetSetReceivesDefault <-
function(object, receives.default)
{
  checkPtrType(object, "GtkWidget")
  receives.default <- as.logical(receives.default)

  w <- .RGtkCall("S_gtk_widget_set_receives_default", object, receives.default, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetGetReceivesDefault <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_receives_default", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetGetSensitive <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_sensitive", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetGetState <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_state", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetSetVisible <-
function(object, visible)
{
  checkPtrType(object, "GtkWidget")
  visible <- as.logical(visible)

  w <- .RGtkCall("S_gtk_widget_set_visible", object, visible, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetGetVisible <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_visible", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetSetWindow <-
function(object, window)
{
  checkPtrType(object, "GtkWidget")
  checkPtrType(window, "GdkWindow")

  w <- .RGtkCall("S_gtk_widget_set_window", object, window, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetHasDefault <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_has_default", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetHasFocus <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_has_focus", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetHasGrab <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_has_grab", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetIsSensitive <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_is_sensitive", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetIsToplevel <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_is_toplevel", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetIsDrawable <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_is_drawable", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkHSVGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_hsv_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkHSVNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_hsv_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkHSVSetColor <-
function(object, h, s, v)
{
  checkPtrType(object, "GtkHSV")
  h <- as.numeric(h)
  s <- as.numeric(s)
  v <- as.numeric(v)

  w <- .RGtkCall("S_gtk_hsv_set_color", object, h, s, v, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkHSVGetColor <-
function(object, h, s, v)
{
  checkPtrType(object, "GtkHSV")
  h <- as.list(as.numeric(h))
  s <- as.list(as.numeric(s))
  v <- as.list(as.numeric(v))

  w <- .RGtkCall("S_gtk_hsv_get_color", object, h, s, v, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkHSVSetMetrics <-
function(object, size, ring.width)
{
  checkPtrType(object, "GtkHSV")
  size <- as.integer(size)
  ring.width <- as.integer(ring.width)

  w <- .RGtkCall("S_gtk_hsv_set_metrics", object, size, ring.width, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkHSVGetMetrics <-
function(object, size, ring.width)
{
  checkPtrType(object, "GtkHSV")
  size <- as.list(as.integer(size))
  ring.width <- as.list(as.integer(ring.width))

  w <- .RGtkCall("S_gtk_hsv_get_metrics", object, size, ring.width, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkHSVIsAdjusting <-
function(object)
{
  checkPtrType(object, "GtkHSV")

  w <- .RGtkCall("S_gtk_hsv_is_adjusting", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkHSVToRgb <-
function(h, s, r, g, b)
{
  h <- as.numeric(h)
  s <- as.numeric(s)
  r <- as.list(as.numeric(r))
  g <- as.list(as.numeric(g))
  b <- as.list(as.numeric(b))

  w <- .RGtkCall("S_gtk_hsv_to_rgb", h, s, r, g, b, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRgbToHsv <-
function(r, g, h, s, v)
{
  r <- as.numeric(r)
  g <- as.numeric(g)
  h <- as.list(as.numeric(h))
  s <- as.list(as.numeric(s))
  v <- as.list(as.numeric(v))

  w <- .RGtkCall("S_gtk_rgb_to_hsv", r, g, h, s, v, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolPaletteGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_tool_palette_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkToolPaletteNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_tool_palette_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkToolPaletteSetGroupPosition <-
function(object, group, position)
{
  checkPtrType(object, "GtkToolPalette")
  checkPtrType(group, "GtkToolItemGroup")
  position <- as.integer(position)

  w <- .RGtkCall("S_gtk_tool_palette_set_group_position", object, group, position, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolPaletteSetExclusive <-
function(object, group, exclusive)
{
  checkPtrType(object, "GtkToolPalette")
  checkPtrType(group, "GtkToolItemGroup")
  exclusive <- as.logical(exclusive)

  w <- .RGtkCall("S_gtk_tool_palette_set_exclusive", object, group, exclusive, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolPaletteSetExpand <-
function(object, group, expand)
{
  checkPtrType(object, "GtkToolPalette")
  checkPtrType(group, "GtkToolItemGroup")
  expand <- as.logical(expand)

  w <- .RGtkCall("S_gtk_tool_palette_set_expand", object, group, expand, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolPaletteGetGroupPosition <-
function(object, group)
{
  checkPtrType(object, "GtkToolPalette")
  checkPtrType(group, "GtkToolItemGroup")

  w <- .RGtkCall("S_gtk_tool_palette_get_group_position", object, group, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolPaletteGetExclusive <-
function(object, group)
{
  checkPtrType(object, "GtkToolPalette")
  checkPtrType(group, "GtkToolItemGroup")

  w <- .RGtkCall("S_gtk_tool_palette_get_exclusive", object, group, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolPaletteGetExpand <-
function(object, group)
{
  checkPtrType(object, "GtkToolPalette")
  checkPtrType(group, "GtkToolItemGroup")

  w <- .RGtkCall("S_gtk_tool_palette_get_expand", object, group, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolPaletteSetIconSize <-
function(object, icon.size)
{
  checkPtrType(object, "GtkToolPalette")
  

  w <- .RGtkCall("S_gtk_tool_palette_set_icon_size", object, icon.size, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolPaletteUnsetIconSize <-
function(object)
{
  checkPtrType(object, "GtkToolPalette")

  w <- .RGtkCall("S_gtk_tool_palette_unset_icon_size", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolPaletteSetStyle <-
function(object, style)
{
  checkPtrType(object, "GtkToolPalette")
  

  w <- .RGtkCall("S_gtk_tool_palette_set_style", object, style, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolPaletteUnsetStyle <-
function(object)
{
  checkPtrType(object, "GtkToolPalette")

  w <- .RGtkCall("S_gtk_tool_palette_unset_style", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolPaletteGetIconSize <-
function(object)
{
  checkPtrType(object, "GtkToolPalette")

  w <- .RGtkCall("S_gtk_tool_palette_get_icon_size", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolPaletteGetStyle <-
function(object)
{
  checkPtrType(object, "GtkToolPalette")

  w <- .RGtkCall("S_gtk_tool_palette_get_style", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolPaletteGetDropItem <-
function(object, x, y)
{
  checkPtrType(object, "GtkToolPalette")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_gtk_tool_palette_get_drop_item", object, x, y, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolPaletteGetDropGroup <-
function(object, x, y)
{
  checkPtrType(object, "GtkToolPalette")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_gtk_tool_palette_get_drop_group", object, x, y, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolPaletteGetDragItem <-
function(object, selection)
{
  checkPtrType(object, "GtkToolPalette")
  checkPtrType(selection, "GtkSelectionData")

  w <- .RGtkCall("S_gtk_tool_palette_get_drag_item", object, selection, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolPaletteSetDragSource <-
function(object, targets)
{
  checkPtrType(object, "GtkToolPalette")
  

  w <- .RGtkCall("S_gtk_tool_palette_set_drag_source", object, targets, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolPaletteAddDragDest <-
function(object, widget, flags, targets, actions)
{
  checkPtrType(object, "GtkToolPalette")
  checkPtrType(widget, "GtkWidget")
  
  
  

  w <- .RGtkCall("S_gtk_tool_palette_add_drag_dest", object, widget, flags, targets, actions, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolPaletteGetHadjustment <-
function(object)
{
  checkPtrType(object, "GtkToolPalette")

  w <- .RGtkCall("S_gtk_tool_palette_get_hadjustment", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolPaletteGetVadjustment <-
function(object)
{
  checkPtrType(object, "GtkToolPalette")

  w <- .RGtkCall("S_gtk_tool_palette_get_vadjustment", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolPaletteGetDragTargetItem <-
function()
{
  

  w <- .RGtkCall("S_gtk_tool_palette_get_drag_target_item", PACKAGE = "RGtk2")

  return(w)
} 


gtkToolPaletteGetDragTargetGroup <-
function()
{
  

  w <- .RGtkCall("S_gtk_tool_palette_get_drag_target_group", PACKAGE = "RGtk2")

  return(w)
} 


gtkToolItemGetEllipsizeMode <-
function(object)
{
  checkPtrType(object, "GtkToolItem")

  w <- .RGtkCall("S_gtk_tool_item_get_ellipsize_mode", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolItemGetTextAlignment <-
function(object)
{
  checkPtrType(object, "GtkToolItem")

  w <- .RGtkCall("S_gtk_tool_item_get_text_alignment", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolItemGetTextOrientation <-
function(object)
{
  checkPtrType(object, "GtkToolItem")

  w <- .RGtkCall("S_gtk_tool_item_get_text_orientation", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolItemGetTextSizeGroup <-
function(object)
{
  checkPtrType(object, "GtkToolItem")

  w <- .RGtkCall("S_gtk_tool_item_get_text_size_group", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolItemGroupGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_tool_item_group_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkToolItemGroupNew <-
function(label, show = TRUE)
{
  label <- as.character(label)

  w <- .RGtkCall("S_gtk_tool_item_group_new", label, PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkToolItemGroupSetLabel <-
function(object, label)
{
  checkPtrType(object, "GtkToolItemGroup")
  label <- as.character(label)

  w <- .RGtkCall("S_gtk_tool_item_group_set_label", object, label, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolItemGroupSetLabelWidget <-
function(object, label.widget)
{
  checkPtrType(object, "GtkToolItemGroup")
  checkPtrType(label.widget, "GtkWidget")

  w <- .RGtkCall("S_gtk_tool_item_group_set_label_widget", object, label.widget, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolItemGroupSetCollapsed <-
function(object, collapsed)
{
  checkPtrType(object, "GtkToolItemGroup")
  collapsed <- as.logical(collapsed)

  w <- .RGtkCall("S_gtk_tool_item_group_set_collapsed", object, collapsed, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolItemGroupSetEllipsize <-
function(object, ellipsize)
{
  checkPtrType(object, "GtkToolItemGroup")
  

  w <- .RGtkCall("S_gtk_tool_item_group_set_ellipsize", object, ellipsize, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolItemGroupSetHeaderRelief <-
function(object, style)
{
  checkPtrType(object, "GtkToolItemGroup")
  

  w <- .RGtkCall("S_gtk_tool_item_group_set_header_relief", object, style, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolItemGroupGetLabel <-
function(object)
{
  checkPtrType(object, "GtkToolItemGroup")

  w <- .RGtkCall("S_gtk_tool_item_group_get_label", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolItemGroupGetLabelWidget <-
function(object)
{
  checkPtrType(object, "GtkToolItemGroup")

  w <- .RGtkCall("S_gtk_tool_item_group_get_label_widget", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolItemGroupGetCollapsed <-
function(object)
{
  checkPtrType(object, "GtkToolItemGroup")

  w <- .RGtkCall("S_gtk_tool_item_group_get_collapsed", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolItemGroupGetEllipsize <-
function(object)
{
  checkPtrType(object, "GtkToolItemGroup")

  w <- .RGtkCall("S_gtk_tool_item_group_get_ellipsize", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolItemGroupGetHeaderRelief <-
function(object)
{
  checkPtrType(object, "GtkToolItemGroup")

  w <- .RGtkCall("S_gtk_tool_item_group_get_header_relief", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolItemGroupInsert <-
function(object, item, position)
{
  checkPtrType(object, "GtkToolItemGroup")
  checkPtrType(item, "GtkToolItem")
  position <- as.integer(position)

  w <- .RGtkCall("S_gtk_tool_item_group_insert", object, item, position, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolItemGroupSetItemPosition <-
function(object, item, position)
{
  checkPtrType(object, "GtkToolItemGroup")
  checkPtrType(item, "GtkToolItem")
  position <- as.integer(position)

  w <- .RGtkCall("S_gtk_tool_item_group_set_item_position", object, item, position, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkToolItemGroupGetItemPosition <-
function(object, item)
{
  checkPtrType(object, "GtkToolItemGroup")
  checkPtrType(item, "GtkToolItem")

  w <- .RGtkCall("S_gtk_tool_item_group_get_item_position", object, item, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolItemGroupGetNItems <-
function(object)
{
  checkPtrType(object, "GtkToolItemGroup")

  w <- .RGtkCall("S_gtk_tool_item_group_get_n_items", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolItemGroupGetNthItem <-
function(object, index)
{
  checkPtrType(object, "GtkToolItemGroup")
  index <- as.numeric(index)

  w <- .RGtkCall("S_gtk_tool_item_group_get_nth_item", object, index, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolItemGroupGetDropItem <-
function(object, x, y)
{
  checkPtrType(object, "GtkToolItemGroup")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_gtk_tool_item_group_get_drop_item", object, x, y, PACKAGE = "RGtk2")

  return(w)
} 


gtkSpinnerGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_spinner_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkSpinnerNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_spinner_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkSpinnerStart <-
function(object)
{
  checkPtrType(object, "GtkSpinner")

  w <- .RGtkCall("S_gtk_spinner_start", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkSpinnerStop <-
function(object)
{
  checkPtrType(object, "GtkSpinner")

  w <- .RGtkCall("S_gtk_spinner_stop", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkCellRendererSpinnerGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_cell_renderer_spinner_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkCellRendererSpinnerNew <-
function()
{
  

  w <- .RGtkCall("S_gtk_cell_renderer_spinner_new", PACKAGE = "RGtk2")

  return(w)
} 


gtkActionGetAlwaysShowImage <-
function(object)
{
  checkPtrType(object, "GtkAction")

  w <- .RGtkCall("S_gtk_action_get_always_show_image", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkActionSetAlwaysShowImage <-
function(object, always.show)
{
  checkPtrType(object, "GtkAction")
  always.show <- as.logical(always.show)

  w <- .RGtkCall("S_gtk_action_set_always_show_image", object, always.show, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkDialogGetWidgetForResponse <-
function(object, response.id)
{
  checkPtrType(object, "GtkDialog")
  response.id <- as.integer(response.id)

  w <- .RGtkCall("S_gtk_dialog_get_widget_for_response", object, response.id, PACKAGE = "RGtk2")

  return(w)
} 


gtkOffscreenWindowGetType <-
function()
{
  

  w <- .RGtkCall("S_gtk_offscreen_window_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gtkOffscreenWindowNew <-
function(show = TRUE)
{
  

  w <- .RGtkCall("S_gtk_offscreen_window_new", PACKAGE = "RGtk2")

  if(show)
    gtkWidgetShowAll(w)

  return(w)
} 


gtkOffscreenWindowGetPixmap <-
function(object)
{
  checkPtrType(object, "GtkOffscreenWindow")

  w <- .RGtkCall("S_gtk_offscreen_window_get_pixmap", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkOffscreenWindowGetPixbuf <-
function(object)
{
  checkPtrType(object, "GtkOffscreenWindow")

  w <- .RGtkCall("S_gtk_offscreen_window_get_pixbuf", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryGetIconWindow <-
function(object, icon.pos)
{
  checkPtrType(object, "GtkEntry")
  

  w <- .RGtkCall("S_gtk_entry_get_icon_window", object, icon.pos, PACKAGE = "RGtk2")

  return(w)
} 


gtkEntryGetTextWindow <-
function(object)
{
  checkPtrType(object, "GtkEntry")

  w <- .RGtkCall("S_gtk_entry_get_text_window", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkNotebookGetActionWidget <-
function(object, pack.type)
{
  checkPtrType(object, "GtkNotebook")
  

  w <- .RGtkCall("S_gtk_notebook_get_action_widget", object, pack.type, PACKAGE = "RGtk2")

  return(w)
} 


gtkNotebookSetActionWidget <-
function(object, widget, pack.type)
{
  checkPtrType(object, "GtkNotebook")
  checkPtrType(widget, "GtkWidget")
  

  w <- .RGtkCall("S_gtk_notebook_set_action_widget", object, widget, pack.type, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPaintSpinner <-
function(object, window, state.type, area, widget, detail, step, x, y, width, height)
{
  checkPtrType(object, "GtkStyle")
  checkPtrType(window, "GdkWindow")
  
  area <- as.GdkRectangle(area)
  checkPtrType(widget, "GtkWidget")
  detail <- as.character(detail)
  step <- as.numeric(step)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gtk_paint_spinner", object, window, state.type, area, widget, detail, step, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPanedGetHandleWindow <-
function(object)
{
  checkPtrType(object, "GtkPaned")

  w <- .RGtkCall("S_gtk_paned_get_handle_window", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRangeGetMinSliderSize <-
function(object)
{
  checkPtrType(object, "GtkRange")

  w <- .RGtkCall("S_gtk_range_get_min_slider_size", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRangeGetRangeRect <-
function(object)
{
  checkPtrType(object, "GtkRange")

  w <- .RGtkCall("S_gtk_range_get_range_rect", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRangeGetSliderRange <-
function(object)
{
  checkPtrType(object, "GtkRange")

  w <- .RGtkCall("S_gtk_range_get_slider_range", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRangeGetSliderSizeFixed <-
function(object)
{
  checkPtrType(object, "GtkRange")

  w <- .RGtkCall("S_gtk_range_get_slider_size_fixed", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkRangeSetMinSliderSize <-
function(object, min.size)
{
  checkPtrType(object, "GtkRange")
  min.size <- as.logical(min.size)

  w <- .RGtkCall("S_gtk_range_set_min_slider_size", object, min.size, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkRangeSetSliderSizeFixed <-
function(object, size.fixed)
{
  checkPtrType(object, "GtkRange")
  size.fixed <- as.logical(size.fixed)

  w <- .RGtkCall("S_gtk_range_set_slider_size_fixed", object, size.fixed, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkStatusbarGetMessageArea <-
function(object)
{
  checkPtrType(object, "GtkStatusbar")

  w <- .RGtkCall("S_gtk_statusbar_get_message_area", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkStatusIconSetName <-
function(object, name)
{
  checkPtrType(object, "GtkStatusIcon")
  name <- as.character(name)

  w <- .RGtkCall("S_gtk_status_icon_set_name", object, name, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkViewportGetBinWindow <-
function(object)
{
  checkPtrType(object, "GtkViewport")

  w <- .RGtkCall("S_gtk_viewport_get_bin_window", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetSetRealized <-
function(object, realized)
{
  checkPtrType(object, "GtkWidget")
  realized <- as.logical(realized)

  w <- .RGtkCall("S_gtk_widget_set_realized", object, realized, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetGetRealized <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_realized", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetSetMapped <-
function(object, mapped)
{
  checkPtrType(object, "GtkWidget")
  mapped <- as.logical(mapped)

  w <- .RGtkCall("S_gtk_widget_set_mapped", object, mapped, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetGetMapped <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_mapped", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetGetRequisition <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_get_requisition", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWidgetStyleAttach <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_style_attach", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWidgetHasRcStyle <-
function(object)
{
  checkPtrType(object, "GtkWidget")

  w <- .RGtkCall("S_gtk_widget_has_rc_style", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowSetMnemonicsVisible <-
function(object, setting)
{
  checkPtrType(object, "GtkWindow")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_gtk_window_set_mnemonics_visible", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkWindowGetMnemonicsVisible <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_get_mnemonics_visible", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkWindowGetWindowType <-
function(object)
{
  checkPtrType(object, "GtkWindow")

  w <- .RGtkCall("S_gtk_window_get_window_type", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkTooltipSetIconFromGicon <-
function(object, gicon, size)
{
  checkPtrType(object, "GtkTooltip")
  checkPtrType(gicon, "GIcon")
  

  w <- .RGtkCall("S_gtk_tooltip_set_icon_from_gicon", object, gicon, size, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gtkPrintContextGetHardMargins <-
function(object)
{
  checkPtrType(object, "GtkPrintContext")

  w <- .RGtkCall("S_gtk_print_context_get_hard_margins", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolShellGetTextOrientation <-
function(object)
{
  checkPtrType(object, "GtkToolShell")

  w <- .RGtkCall("S_gtk_tool_shell_get_text_orientation", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolShellGetTextAlignment <-
function(object)
{
  checkPtrType(object, "GtkToolShell")

  w <- .RGtkCall("S_gtk_tool_shell_get_text_alignment", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolShellGetEllipsizeMode <-
function(object)
{
  checkPtrType(object, "GtkToolShell")

  w <- .RGtkCall("S_gtk_tool_shell_get_ellipsize_mode", object, PACKAGE = "RGtk2")

  return(w)
} 


gtkToolShellGetTextSizeGroup <-
function(object)
{
  checkPtrType(object, "GtkToolShell")

  w <- .RGtkCall("S_gtk_tool_shell_get_text_size_group", object, PACKAGE = "RGtk2")

  return(w)
} 

