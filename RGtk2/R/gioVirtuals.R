if(!exists('.virtuals')) .virtuals <- new.env()
assign("GAppInfo", c("dup", "equal", "get_id", "get_name", "get_description", "get_executable", "get_icon", "launch", "supports_uris", "supports_files", "launch_uris", "should_show", "set_as_default_for_type", "set_as_default_for_extension", "add_supports_type", "can_remove_supports_type", "remove_supports_type", "get_commandline"), .virtuals)
assign("GAppLaunchContext", c("get_display", "get_startup_notify_id", "launch_failed"), .virtuals)
assign("GAsyncResult", c("get_user_data", "get_source_object"), .virtuals)
assign("GBufferedInputStream", c("fill", "fill_async", "fill_finish"), .virtuals)
assign("GDrive", c("get_name", "get_icon", "has_volumes", "get_volumes", "is_media_removable", "has_media", "is_media_check_automatic", "can_poll_for_media", "can_eject", "eject", "eject_finish", "poll_for_media", "poll_for_media_finish", "get_identifier", "enumerate_identifiers", "get_start_stop_type", "start", "start_finish", "stop", "stop_finish", "can_start", "can_start_degraded", "can_stop", "eject_with_operation", "eject_with_operation_finish"), .virtuals)
assign("GFileEnumerator", c("next_file", "close_fn", "next_files_async", "next_files_finish", "close_async", "close_finish"), .virtuals)
assign("GFile", c("dup", "equal", "get_basename", "get_path", "get_uri", "get_parse_name", "get_parent", "get_child_for_display_name", "prefix_matches", "get_relative_path", "resolve_relative_path", "is_native", "has_uri_scheme", "get_uri_scheme", "read_fn", "read_async", "read_finish", "append_to", "create", "replace", "append_to_async", "append_to_finish", "create_async", "create_finish", "replace_async", "replace_finish", "query_info", "query_info_async", "query_info_finish", "query_filesystem_info", "query_filesystem_info_async", "query_filesystem_info_finish", "find_enclosing_mount", "find_enclosing_mount_async", "find_enclosing_mount_finish", "enumerate_children", "enumerate_children_async", "enumerate_children_finish", "set_display_name", "set_display_name_async", "set_display_name_finish", "delete_file", "trash", "copy", "copy_async", "copy_finish", "move", "make_directory", "make_symbolic_link", "query_settable_attributes", "query_writable_namespaces", "set_attribute", "set_attributes_from_info", "set_attributes_async", "set_attributes_finish", "mount_enclosing_volume", "mount_enclosing_volume_finish", "mount_mountable", "mount_mountable_finish", "unmount_mountable", "unmount_mountable_finish", "eject_mountable", "eject_mountable_finish", "monitor_dir", "monitor_file", "create_readwrite", "create_readwrite_async", "create_readwrite_finish", "eject_mountable_with_operation", "eject_mountable_with_operation_finish", "open_readwrite", "open_readwrite_async", "open_readwrite_finish", "poll_mountable", "poll_mountable_finish", "replace_readwrite", "replace_readwrite_async", "replace_readwrite_finish", "start_mountable", "start_mountable_finish", "stop_mountable", "stop_mountable_finish", "unmount_mountable_with_operation", "unmount_mountable_with_operation_finish"), .virtuals)
assign("GFileInputStream", c("query_info", "query_info_async", "query_info_finish"), .virtuals)
assign("GFileMonitor", c("cancel"), .virtuals)
assign("GFileOutputStream", c("query_info", "query_info_async", "query_info_finish", "get_etag"), .virtuals)
assign("GIcon", c("hash", "equal"), .virtuals)
assign("GInputStream", c("skip", "close_fn", "read_finish", "skip_async", "skip_finish", "close_async", "close_finish"), .virtuals)
assign("GLoadableIcon", c("load", "load_async", "load_finish"), .virtuals)
assign("GMount", c("get_root", "get_name", "get_icon", "get_uuid", "get_volume", "get_drive", "can_unmount", "can_eject", "unmount", "unmount_finish", "eject", "eject_finish", "remount", "remount_finish", "guess_content_type", "guess_content_type_finish", "guess_content_type_sync", "unmount_with_operation", "unmount_with_operation_finish", "eject_with_operation", "eject_with_operation_finish"), .virtuals)
assign("GOutputStream", c("write_fn", "splice", "flush", "close_fn", "write_async", "write_finish", "splice_async", "splice_finish", "flush_async", "flush_finish", "close_async", "close_finish"), .virtuals)
assign("GSeekable", c("tell", "can_seek", "seek", "can_truncate", "truncate_fn"), .virtuals)
assign("GVfs", c("is_active", "get_file_for_path", "get_file_for_uri", "parse_name", "get_supported_uri_schemes"), .virtuals)
assign("GVolume", c("get_name", "get_icon", "get_uuid", "get_drive", "get_mount", "can_mount", "can_eject", "should_automount", "mount_fn", "mount_finish", "eject", "eject_finish", "get_identifier", "enumerate_identifiers", "get_activation_root", "eject_with_operation", "eject_with_operation_finish"), .virtuals)
assign("GVolumeMonitor", c("get_connected_drives", "get_volumes", "get_mounts", "get_volume_for_uuid", "get_mount_for_uuid"), .virtuals)
assign("GAsyncInitable", c("init_async", "init_finish"), .virtuals)
assign("GFileIOStream", c("query_info", "query_info_async", "query_info_finish", "get_etag"), .virtuals)
assign("GInetAddress", c("to_string", "to_bytes"), .virtuals)
assign("GInitable", c("init"), .virtuals)
assign("GIOStream", c("get_input_stream", "get_output_stream", "close_fn", "close_async", "close_finish"), .virtuals)
assign("GResolver", c("lookup_by_name", "lookup_by_name_async", "lookup_by_name_finish", "lookup_by_address", "lookup_by_address_async", "lookup_by_address_finish", "lookup_service", "lookup_service_async", "lookup_service_finish"), .virtuals)
assign("GSocketAddressEnumerator", c("next", "next_async", "next_finish"), .virtuals)
assign("GSocketAddress", c("get_family", "to_native", "get_native_size"), .virtuals)
assign("GSocketConnectable", c("enumerate"), .virtuals)
assign("GSocketControlMessage", c("get_size", "get_level", "get_type", "serialize"), .virtuals)
assign("GSocketListener", c("changed"), .virtuals)


gAppInfoIfaceDup <-
function(object.class, object)
{
  checkPtrType(object.class, "GAppInfoIface")
  checkPtrType(object, "GAppInfo")

  w <- .RGtkCall("S_gapp_info_iface_dup", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gAppInfoIfaceEqual <-
function(object.class, object, appinfo2)
{
  checkPtrType(object.class, "GAppInfoIface")
  checkPtrType(object, "GAppInfo")
  checkPtrType(appinfo2, "GAppInfo")

  w <- .RGtkCall("S_gapp_info_iface_equal", object.class, object, appinfo2, PACKAGE = "RGtk2")

  return(w)
}

gAppInfoIfaceGetId <-
function(object.class, object)
{
  checkPtrType(object.class, "GAppInfoIface")
  checkPtrType(object, "GAppInfo")

  w <- .RGtkCall("S_gapp_info_iface_get_id", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gAppInfoIfaceGetName <-
function(object.class, object)
{
  checkPtrType(object.class, "GAppInfoIface")
  checkPtrType(object, "GAppInfo")

  w <- .RGtkCall("S_gapp_info_iface_get_name", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gAppInfoIfaceGetDescription <-
function(object.class, object)
{
  checkPtrType(object.class, "GAppInfoIface")
  checkPtrType(object, "GAppInfo")

  w <- .RGtkCall("S_gapp_info_iface_get_description", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gAppInfoIfaceGetExecutable <-
function(object.class, object)
{
  checkPtrType(object.class, "GAppInfoIface")
  checkPtrType(object, "GAppInfo")

  w <- .RGtkCall("S_gapp_info_iface_get_executable", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gAppInfoIfaceGetIcon <-
function(object.class, object)
{
  checkPtrType(object.class, "GAppInfoIface")
  checkPtrType(object, "GAppInfo")

  w <- .RGtkCall("S_gapp_info_iface_get_icon", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gAppInfoIfaceLaunch <-
function(object.class, object, files, launch.context, .errwarn = TRUE)
{
  checkPtrType(object.class, "GAppInfoIface")
  checkPtrType(object, "GAppInfo")
  files <- lapply(files, function(x) { x <- as.GList(x); x })
  checkPtrType(launch.context, "GAppLaunchContext")

  w <- .RGtkCall("S_gapp_info_iface_launch", object.class, object, files, launch.context, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gAppInfoIfaceSupportsUris <-
function(object.class, object)
{
  checkPtrType(object.class, "GAppInfoIface")
  checkPtrType(object, "GAppInfo")

  w <- .RGtkCall("S_gapp_info_iface_supports_uris", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gAppInfoIfaceSupportsFiles <-
function(object.class, object)
{
  checkPtrType(object.class, "GAppInfoIface")
  checkPtrType(object, "GAppInfo")

  w <- .RGtkCall("S_gapp_info_iface_supports_files", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gAppInfoIfaceLaunchUris <-
function(object.class, object, uris, launch.context, .errwarn = TRUE)
{
  checkPtrType(object.class, "GAppInfoIface")
  checkPtrType(object, "GAppInfo")
  uris <- lapply(uris, function(x) { x <- as.GList(x); x })
  checkPtrType(launch.context, "GAppLaunchContext")

  w <- .RGtkCall("S_gapp_info_iface_launch_uris", object.class, object, uris, launch.context, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gAppInfoIfaceShouldShow <-
function(object.class, object)
{
  checkPtrType(object.class, "GAppInfoIface")
  checkPtrType(object, "GAppInfo")

  w <- .RGtkCall("S_gapp_info_iface_should_show", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gAppInfoIfaceSetAsDefaultForType <-
function(object.class, object, content.type, .errwarn = TRUE)
{
  checkPtrType(object.class, "GAppInfoIface")
  checkPtrType(object, "GAppInfo")
  content.type <- as.character(content.type)

  w <- .RGtkCall("S_gapp_info_iface_set_as_default_for_type", object.class, object, content.type, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gAppInfoIfaceSetAsDefaultForExtension <-
function(object.class, object, extension, .errwarn = TRUE)
{
  checkPtrType(object.class, "GAppInfoIface")
  checkPtrType(object, "GAppInfo")
  extension <- as.character(extension)

  w <- .RGtkCall("S_gapp_info_iface_set_as_default_for_extension", object.class, object, extension, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gAppInfoIfaceAddSupportsType <-
function(object.class, object, content.type, .errwarn = TRUE)
{
  checkPtrType(object.class, "GAppInfoIface")
  checkPtrType(object, "GAppInfo")
  content.type <- as.character(content.type)

  w <- .RGtkCall("S_gapp_info_iface_add_supports_type", object.class, object, content.type, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gAppInfoIfaceCanRemoveSupportsType <-
function(object.class, object)
{
  checkPtrType(object.class, "GAppInfoIface")
  checkPtrType(object, "GAppInfo")

  w <- .RGtkCall("S_gapp_info_iface_can_remove_supports_type", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gAppInfoIfaceRemoveSupportsType <-
function(object.class, object, content.type, .errwarn = TRUE)
{
  checkPtrType(object.class, "GAppInfoIface")
  checkPtrType(object, "GAppInfo")
  content.type <- as.character(content.type)

  w <- .RGtkCall("S_gapp_info_iface_remove_supports_type", object.class, object, content.type, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gAppLaunchContextClassGetDisplay <-
function(object.class, object, info, files)
{
  checkPtrType(object.class, "GAppLaunchContextClass")
  checkPtrType(object, "GAppLaunchContext")
  checkPtrType(info, "GAppInfo")
  files <- lapply(files, function(x) { x <- as.GList(x); x })

  w <- .RGtkCall("S_gapp_launch_context_class_get_display", object.class, object, info, files, PACKAGE = "RGtk2")

  return(w)
}

gAppLaunchContextClassGetStartupNotifyId <-
function(object.class, object, info, files)
{
  checkPtrType(object.class, "GAppLaunchContextClass")
  checkPtrType(object, "GAppLaunchContext")
  checkPtrType(info, "GAppInfo")
  files <- lapply(files, function(x) { x <- as.GList(x); x })

  w <- .RGtkCall("S_gapp_launch_context_class_get_startup_notify_id", object.class, object, info, files, PACKAGE = "RGtk2")

  return(w)
}

gAppLaunchContextClassLaunchFailed <-
function(object.class, object, startup.notify.id)
{
  checkPtrType(object.class, "GAppLaunchContextClass")
  checkPtrType(object, "GAppLaunchContext")
  startup.notify.id <- as.character(startup.notify.id)

  w <- .RGtkCall("S_gapp_launch_context_class_launch_failed", object.class, object, startup.notify.id, PACKAGE = "RGtk2")

  return(invisible(w))
}

gAsyncResultIfaceGetUserData <-
function(object.class, object)
{
  checkPtrType(object.class, "GAsyncResultIface")
  checkPtrType(object, "GAsyncResult")

  w <- .RGtkCall("S_gasync_result_iface_get_user_data", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gAsyncResultIfaceGetSourceObject <-
function(object.class, object)
{
  checkPtrType(object.class, "GAsyncResultIface")
  checkPtrType(object, "GAsyncResult")

  w <- .RGtkCall("S_gasync_result_iface_get_source_object", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gBufferedInputStreamClassFill <-
function(object.class, object, count, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GBufferedInputStreamClass")
  checkPtrType(object, "GBufferedInputStream")
  count <- as.integer(count)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gbuffered_input_stream_class_fill", object.class, object, count, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gBufferedInputStreamClassFillAsync <-
function(object.class, object, count, io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GBufferedInputStreamClass")
  checkPtrType(object, "GBufferedInputStream")
  count <- as.integer(count)
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gbuffered_input_stream_class_fill_async", object.class, object, count, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gBufferedInputStreamClassFillFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GBufferedInputStreamClass")
  checkPtrType(object, "GBufferedInputStream")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_gbuffered_input_stream_class_fill_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gDriveIfaceGetName <-
function(object.class, object)
{
  checkPtrType(object.class, "GDriveIface")
  checkPtrType(object, "GDrive")

  w <- .RGtkCall("S_gdrive_iface_get_name", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gDriveIfaceGetIcon <-
function(object.class, object)
{
  checkPtrType(object.class, "GDriveIface")
  checkPtrType(object, "GDrive")

  w <- .RGtkCall("S_gdrive_iface_get_icon", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gDriveIfaceHasVolumes <-
function(object.class, object)
{
  checkPtrType(object.class, "GDriveIface")
  checkPtrType(object, "GDrive")

  w <- .RGtkCall("S_gdrive_iface_has_volumes", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gDriveIfaceGetVolumes <-
function(object.class, object)
{
  checkPtrType(object.class, "GDriveIface")
  checkPtrType(object, "GDrive")

  w <- .RGtkCall("S_gdrive_iface_get_volumes", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gDriveIfaceIsMediaRemovable <-
function(object.class, object)
{
  checkPtrType(object.class, "GDriveIface")
  checkPtrType(object, "GDrive")

  w <- .RGtkCall("S_gdrive_iface_is_media_removable", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gDriveIfaceHasMedia <-
function(object.class, object)
{
  checkPtrType(object.class, "GDriveIface")
  checkPtrType(object, "GDrive")

  w <- .RGtkCall("S_gdrive_iface_has_media", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gDriveIfaceIsMediaCheckAutomatic <-
function(object.class, object)
{
  checkPtrType(object.class, "GDriveIface")
  checkPtrType(object, "GDrive")

  w <- .RGtkCall("S_gdrive_iface_is_media_check_automatic", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gDriveIfaceCanPollForMedia <-
function(object.class, object)
{
  checkPtrType(object.class, "GDriveIface")
  checkPtrType(object, "GDrive")

  w <- .RGtkCall("S_gdrive_iface_can_poll_for_media", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gDriveIfaceCanEject <-
function(object.class, object)
{
  checkPtrType(object.class, "GDriveIface")
  checkPtrType(object, "GDrive")

  w <- .RGtkCall("S_gdrive_iface_can_eject", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gDriveIfaceEject <-
function(object.class, object, flags = "G_MOUNT_UNMOUNT_NONE", cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GDriveIface")
  checkPtrType(object, "GDrive")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gdrive_iface_eject", object.class, object, flags, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gDriveIfaceEjectFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GDriveIface")
  checkPtrType(object, "GDrive")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_gdrive_iface_eject_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gDriveIfacePollForMedia <-
function(object.class, object, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GDriveIface")
  checkPtrType(object, "GDrive")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gdrive_iface_poll_for_media", object.class, object, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gDriveIfacePollForMediaFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GDriveIface")
  checkPtrType(object, "GDrive")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_gdrive_iface_poll_for_media_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gDriveIfaceGetIdentifier <-
function(object.class, object, kind)
{
  checkPtrType(object.class, "GDriveIface")
  checkPtrType(object, "GDrive")
  kind <- as.character(kind)

  w <- .RGtkCall("S_gdrive_iface_get_identifier", object.class, object, kind, PACKAGE = "RGtk2")

  return(w)
}

gDriveIfaceEnumerateIdentifiers <-
function(object.class, object)
{
  checkPtrType(object.class, "GDriveIface")
  checkPtrType(object, "GDrive")

  w <- .RGtkCall("S_gdrive_iface_enumerate_identifiers", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gFileEnumeratorClassNextFile <-
function(object.class, object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileEnumeratorClass")
  checkPtrType(object, "GFileEnumerator")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gfile_enumerator_class_next_file", object.class, object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileEnumeratorClassCloseFn <-
function(object.class, object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileEnumeratorClass")
  checkPtrType(object, "GFileEnumerator")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gfile_enumerator_class_close_fn", object.class, object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileEnumeratorClassNextFilesAsync <-
function(object.class, object, num.files, io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GFileEnumeratorClass")
  checkPtrType(object, "GFileEnumerator")
  num.files <- as.integer(num.files)
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gfile_enumerator_class_next_files_async", object.class, object, num.files, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gFileEnumeratorClassNextFilesFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileEnumeratorClass")
  checkPtrType(object, "GFileEnumerator")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_gfile_enumerator_class_next_files_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileEnumeratorClassCloseAsync <-
function(object.class, object, io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GFileEnumeratorClass")
  checkPtrType(object, "GFileEnumerator")
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gfile_enumerator_class_close_async", object.class, object, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gFileEnumeratorClassCloseFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileEnumeratorClass")
  checkPtrType(object, "GFileEnumerator")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_gfile_enumerator_class_close_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceDup <-
function(object.class, object)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")

  w <- .RGtkCall("S_gfile_iface_dup", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gFileIfaceEqual <-
function(object.class, object, file2)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  checkPtrType(file2, "GFile")

  w <- .RGtkCall("S_gfile_iface_equal", object.class, object, file2, PACKAGE = "RGtk2")

  return(w)
}

gFileIfaceGetBasename <-
function(object.class, object)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")

  w <- .RGtkCall("S_gfile_iface_get_basename", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gFileIfaceGetPath <-
function(object.class, object)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")

  w <- .RGtkCall("S_gfile_iface_get_path", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gFileIfaceGetUri <-
function(object.class, object)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")

  w <- .RGtkCall("S_gfile_iface_get_uri", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gFileIfaceGetParseName <-
function(object.class, object)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")

  w <- .RGtkCall("S_gfile_iface_get_parse_name", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gFileIfaceGetParent <-
function(object.class, object)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")

  w <- .RGtkCall("S_gfile_iface_get_parent", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gFileIfaceGetChildForDisplayName <-
function(object.class, object, display.name, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  display.name <- as.character(display.name)

  w <- .RGtkCall("S_gfile_iface_get_child_for_display_name", object.class, object, display.name, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfacePrefixMatches <-
function(object.class, object, file)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  checkPtrType(file, "GFile")

  w <- .RGtkCall("S_gfile_iface_prefix_matches", object.class, object, file, PACKAGE = "RGtk2")

  return(w)
}

gFileIfaceGetRelativePath <-
function(object.class, object, descendant)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  checkPtrType(descendant, "GFile")

  w <- .RGtkCall("S_gfile_iface_get_relative_path", object.class, object, descendant, PACKAGE = "RGtk2")

  return(w)
}

gFileIfaceResolveRelativePath <-
function(object.class, object, relative.path)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  relative.path <- as.character(relative.path)

  w <- .RGtkCall("S_gfile_iface_resolve_relative_path", object.class, object, relative.path, PACKAGE = "RGtk2")

  return(w)
}

gFileIfaceIsNative <-
function(object.class, object)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")

  w <- .RGtkCall("S_gfile_iface_is_native", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gFileIfaceHasUriScheme <-
function(object.class, object, uri.scheme)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  uri.scheme <- as.character(uri.scheme)

  w <- .RGtkCall("S_gfile_iface_has_uri_scheme", object.class, object, uri.scheme, PACKAGE = "RGtk2")

  return(w)
}

gFileIfaceGetUriScheme <-
function(object.class, object)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")

  w <- .RGtkCall("S_gfile_iface_get_uri_scheme", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gFileIfaceReadFn <-
function(object.class, object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gfile_iface_read_fn", object.class, object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceReadAsync <-
function(object.class, object, io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gfile_iface_read_async", object.class, object, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gFileIfaceReadFinish <-
function(object.class, object, res, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  checkPtrType(res, "GAsyncResult")

  w <- .RGtkCall("S_gfile_iface_read_finish", object.class, object, res, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceAppendTo <-
function(object.class, object, flags = "G_FILE_CREATE_NONE", cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gfile_iface_append_to", object.class, object, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceCreate <-
function(object.class, object, flags = "G_FILE_CREATE_NONE", cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gfile_iface_create", object.class, object, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceReplace <-
function(object.class, object, etag, make.backup, flags = "G_FILE_CREATE_NONE", cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  etag <- as.character(etag)
  make.backup <- as.logical(make.backup)
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gfile_iface_replace", object.class, object, etag, make.backup, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceAppendToAsync <-
function(object.class, object, flags = "G_FILE_CREATE_NONE", io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gfile_iface_append_to_async", object.class, object, flags, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gFileIfaceAppendToFinish <-
function(object.class, object, res, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  checkPtrType(res, "GAsyncResult")

  w <- .RGtkCall("S_gfile_iface_append_to_finish", object.class, object, res, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceCreateAsync <-
function(object.class, object, flags = "G_FILE_CREATE_NONE", io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gfile_iface_create_async", object.class, object, flags, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gFileIfaceCreateFinish <-
function(object.class, object, res, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  checkPtrType(res, "GAsyncResult")

  w <- .RGtkCall("S_gfile_iface_create_finish", object.class, object, res, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceReplaceAsync <-
function(object.class, object, etag, make.backup, flags = "G_FILE_CREATE_NONE", io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  etag <- as.character(etag)
  make.backup <- as.logical(make.backup)
  
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gfile_iface_replace_async", object.class, object, etag, make.backup, flags, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gFileIfaceReplaceFinish <-
function(object.class, object, res, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  checkPtrType(res, "GAsyncResult")

  w <- .RGtkCall("S_gfile_iface_replace_finish", object.class, object, res, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceQueryInfo <-
function(object.class, object, attributes, flags = "G_FILE_QUERY_INFO_NONE", cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  attributes <- as.character(attributes)
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gfile_iface_query_info", object.class, object, attributes, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceQueryInfoAsync <-
function(object.class, object, attributes, flags = "G_FILE_QUERY_INFO_NONE", io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  attributes <- as.character(attributes)
  
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gfile_iface_query_info_async", object.class, object, attributes, flags, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gFileIfaceQueryInfoFinish <-
function(object.class, object, res, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  checkPtrType(res, "GAsyncResult")

  w <- .RGtkCall("S_gfile_iface_query_info_finish", object.class, object, res, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceQueryFilesystemInfo <-
function(object.class, object, attributes, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  attributes <- as.character(attributes)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gfile_iface_query_filesystem_info", object.class, object, attributes, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceQueryFilesystemInfoAsync <-
function(object.class, object, attributes, io.priority, cancellable, callback, user.data = NULL)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  attributes <- as.character(attributes)
  io.priority <- as.integer(io.priority)
  checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gfile_iface_query_filesystem_info_async", object.class, object, attributes, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gFileIfaceQueryFilesystemInfoFinish <-
function(object.class, object, res, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  checkPtrType(res, "GAsyncResult")

  w <- .RGtkCall("S_gfile_iface_query_filesystem_info_finish", object.class, object, res, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceFindEnclosingMount <-
function(object.class, object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gfile_iface_find_enclosing_mount", object.class, object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceFindEnclosingMountAsync <-
function(object.class, object, io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gfile_iface_find_enclosing_mount_async", object.class, object, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gFileIfaceFindEnclosingMountFinish <-
function(object.class, object, res, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  checkPtrType(res, "GAsyncResult")

  w <- .RGtkCall("S_gfile_iface_find_enclosing_mount_finish", object.class, object, res, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceEnumerateChildren <-
function(object.class, object, attributes, flags = "G_FILE_QUERY_INFO_NONE", cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  attributes <- as.character(attributes)
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gfile_iface_enumerate_children", object.class, object, attributes, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceEnumerateChildrenAsync <-
function(object.class, object, attributes, flags = "G_FILE_QUERY_INFO_NONE", io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  attributes <- as.character(attributes)
  
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gfile_iface_enumerate_children_async", object.class, object, attributes, flags, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gFileIfaceEnumerateChildrenFinish <-
function(object.class, object, res, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  checkPtrType(res, "GAsyncResult")

  w <- .RGtkCall("S_gfile_iface_enumerate_children_finish", object.class, object, res, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceSetDisplayName <-
function(object.class, object, display.name, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  display.name <- as.character(display.name)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gfile_iface_set_display_name", object.class, object, display.name, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceSetDisplayNameAsync <-
function(object.class, object, display.name, io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  display.name <- as.character(display.name)
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gfile_iface_set_display_name_async", object.class, object, display.name, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gFileIfaceSetDisplayNameFinish <-
function(object.class, object, res, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  checkPtrType(res, "GAsyncResult")

  w <- .RGtkCall("S_gfile_iface_set_display_name_finish", object.class, object, res, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceDeleteFile <-
function(object.class, object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gfile_iface_delete_file", object.class, object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceTrash <-
function(object.class, object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gfile_iface_trash", object.class, object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceCopy <-
function(object.class, object, destination, flags = "G_FILE_COPY_NONE", cancellable = NULL, progress.callback, progress.callback.data, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  checkPtrType(destination, "GFile")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  progress.callback <- as.function(progress.callback)
  

  w <- .RGtkCall("S_gfile_iface_copy", object.class, object, destination, flags, cancellable, progress.callback, progress.callback.data, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceCopyAsync <-
function(object.class, object, destination, flags = "G_FILE_COPY_NONE", io.priority = 0, cancellable = NULL, progress.callback, progress.callback.data, callback, user.data = NULL)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  checkPtrType(destination, "GFile")
  
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  progress.callback <- as.function(progress.callback)
  
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gfile_iface_copy_async", object.class, object, destination, flags, io.priority, cancellable, progress.callback, progress.callback.data, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gFileIfaceCopyFinish <-
function(object.class, object, res, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  checkPtrType(res, "GAsyncResult")

  w <- .RGtkCall("S_gfile_iface_copy_finish", object.class, object, res, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceMove <-
function(object.class, object, destination, flags = "G_FILE_COPY_NONE", cancellable = NULL, progress.callback, progress.callback.data, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  checkPtrType(destination, "GFile")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  progress.callback <- as.function(progress.callback)
  

  w <- .RGtkCall("S_gfile_iface_move", object.class, object, destination, flags, cancellable, progress.callback, progress.callback.data, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceMakeDirectory <-
function(object.class, object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gfile_iface_make_directory", object.class, object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceMakeSymbolicLink <-
function(object.class, object, symlink.value, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  symlink.value <- as.character(symlink.value)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gfile_iface_make_symbolic_link", object.class, object, symlink.value, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceQuerySettableAttributes <-
function(object.class, object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gfile_iface_query_settable_attributes", object.class, object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceQueryWritableNamespaces <-
function(object.class, object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gfile_iface_query_writable_namespaces", object.class, object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceSetAttribute <-
function(object.class, object, attribute, type, value.p, flags = "G_FILE_QUERY_INFO_NONE", cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  attribute <- as.character(attribute)
  
  
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gfile_iface_set_attribute", object.class, object, attribute, type, value.p, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceSetAttributesFromInfo <-
function(object.class, object, info, flags = "G_FILE_QUERY_INFO_NONE", cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  checkPtrType(info, "GFileInfo")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gfile_iface_set_attributes_from_info", object.class, object, info, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceSetAttributesAsync <-
function(object.class, object, info, flags = "G_FILE_QUERY_INFO_NONE", io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  checkPtrType(info, "GFileInfo")
  
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gfile_iface_set_attributes_async", object.class, object, info, flags, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gFileIfaceSetAttributesFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_gfile_iface_set_attributes_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceMountEnclosingVolume <-
function(object.class, object, flags, mount.operation, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  
  checkPtrType(mount.operation, "GMountOperation")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gfile_iface_mount_enclosing_volume", object.class, object, flags, mount.operation, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gFileIfaceMountEnclosingVolumeFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_gfile_iface_mount_enclosing_volume_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceMountMountable <-
function(object.class, object, flags, mount.operation, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  
  checkPtrType(mount.operation, "GMountOperation")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gfile_iface_mount_mountable", object.class, object, flags, mount.operation, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gFileIfaceMountMountableFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_gfile_iface_mount_mountable_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceUnmountMountable <-
function(object.class, object, flags = "G_MOUNT_UNMOUNT_NONE", cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gfile_iface_unmount_mountable", object.class, object, flags, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gFileIfaceUnmountMountableFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_gfile_iface_unmount_mountable_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceEjectMountable <-
function(object.class, object, flags = "G_MOUNT_UNMOUNT_NONE", cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gfile_iface_eject_mountable", object.class, object, flags, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gFileIfaceEjectMountableFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_gfile_iface_eject_mountable_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceMonitorDir <-
function(object.class, object, flags = "G_FILE_MONITOR_NONE", cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gfile_iface_monitor_dir", object.class, object, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceMonitorFile <-
function(object.class, object, flags = "G_FILE_MONITOR_NONE", cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gfile_iface_monitor_file", object.class, object, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileInputStreamClassQueryInfo <-
function(object.class, object, attributes, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileInputStreamClass")
  checkPtrType(object, "GFileInputStream")
  attributes <- as.character(attributes)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gfile_input_stream_class_query_info", object.class, object, attributes, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileInputStreamClassQueryInfoAsync <-
function(object.class, object, attributes, io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GFileInputStreamClass")
  checkPtrType(object, "GFileInputStream")
  attributes <- as.character(attributes)
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gfile_input_stream_class_query_info_async", object.class, object, attributes, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gFileInputStreamClassQueryInfoFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileInputStreamClass")
  checkPtrType(object, "GFileInputStream")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_gfile_input_stream_class_query_info_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileMonitorClassCancel <-
function(object.class, object)
{
  checkPtrType(object.class, "GFileMonitorClass")
  checkPtrType(object, "GFileMonitor")

  w <- .RGtkCall("S_gfile_monitor_class_cancel", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gFileOutputStreamClassQueryInfo <-
function(object.class, object, attributes, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileOutputStreamClass")
  checkPtrType(object, "GFileOutputStream")
  attributes <- as.character(attributes)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gfile_output_stream_class_query_info", object.class, object, attributes, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileOutputStreamClassQueryInfoAsync <-
function(object.class, object, attributes, io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GFileOutputStreamClass")
  checkPtrType(object, "GFileOutputStream")
  attributes <- as.character(attributes)
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gfile_output_stream_class_query_info_async", object.class, object, attributes, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gFileOutputStreamClassQueryInfoFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileOutputStreamClass")
  checkPtrType(object, "GFileOutputStream")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_gfile_output_stream_class_query_info_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileOutputStreamClassGetEtag <-
function(object.class, object)
{
  checkPtrType(object.class, "GFileOutputStreamClass")
  checkPtrType(object, "GFileOutputStream")

  w <- .RGtkCall("S_gfile_output_stream_class_get_etag", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gIconIfaceHash <-
function(object.class, object)
{
  checkPtrType(object.class, "GIconIface")
  checkPtrType(object, "GIcon")

  w <- .RGtkCall("S_gicon_iface_hash", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gIconIfaceEqual <-
function(object.class, object, icon2)
{
  checkPtrType(object.class, "GIconIface")
  checkPtrType(object, "GIcon")
  checkPtrType(icon2, "GIcon")

  w <- .RGtkCall("S_gicon_iface_equal", object.class, object, icon2, PACKAGE = "RGtk2")

  return(w)
}

gInputStreamClassSkip <-
function(object.class, object, count, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GInputStreamClass")
  checkPtrType(object, "GInputStream")
  count <- as.numeric(count)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_ginput_stream_class_skip", object.class, object, count, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gInputStreamClassCloseFn <-
function(object.class, object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GInputStreamClass")
  checkPtrType(object, "GInputStream")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_ginput_stream_class_close_fn", object.class, object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gInputStreamClassReadFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GInputStreamClass")
  checkPtrType(object, "GInputStream")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_ginput_stream_class_read_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gInputStreamClassSkipAsync <-
function(object.class, object, count, io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GInputStreamClass")
  checkPtrType(object, "GInputStream")
  count <- as.numeric(count)
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_ginput_stream_class_skip_async", object.class, object, count, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gInputStreamClassSkipFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GInputStreamClass")
  checkPtrType(object, "GInputStream")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_ginput_stream_class_skip_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gInputStreamClassCloseAsync <-
function(object.class, object, io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GInputStreamClass")
  checkPtrType(object, "GInputStream")
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_ginput_stream_class_close_async", object.class, object, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gInputStreamClassCloseFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GInputStreamClass")
  checkPtrType(object, "GInputStream")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_ginput_stream_class_close_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gLoadableIconIfaceLoad <-
function(object.class, object, size, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GLoadableIconIface")
  checkPtrType(object, "GLoadableIcon")
  size <- as.integer(size)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gloadable_icon_iface_load", object.class, object, size, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gLoadableIconIfaceLoadAsync <-
function(object.class, object, size, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GLoadableIconIface")
  checkPtrType(object, "GLoadableIcon")
  size <- as.integer(size)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gloadable_icon_iface_load_async", object.class, object, size, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gLoadableIconIfaceLoadFinish <-
function(object.class, object, res, type, .errwarn = TRUE)
{
  checkPtrType(object.class, "GLoadableIconIface")
  checkPtrType(object, "GLoadableIcon")
  checkPtrType(res, "GAsyncResult")
  type <- as.list(as.character(type))

  w <- .RGtkCall("S_gloadable_icon_iface_load_finish", object.class, object, res, type, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gMountIfaceGetRoot <-
function(object.class, object)
{
  checkPtrType(object.class, "GMountIface")
  checkPtrType(object, "GMount")

  w <- .RGtkCall("S_gmount_iface_get_root", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gMountIfaceGetName <-
function(object.class, object)
{
  checkPtrType(object.class, "GMountIface")
  checkPtrType(object, "GMount")

  w <- .RGtkCall("S_gmount_iface_get_name", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gMountIfaceGetIcon <-
function(object.class, object)
{
  checkPtrType(object.class, "GMountIface")
  checkPtrType(object, "GMount")

  w <- .RGtkCall("S_gmount_iface_get_icon", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gMountIfaceGetUuid <-
function(object.class, object)
{
  checkPtrType(object.class, "GMountIface")
  checkPtrType(object, "GMount")

  w <- .RGtkCall("S_gmount_iface_get_uuid", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gMountIfaceGetVolume <-
function(object.class, object)
{
  checkPtrType(object.class, "GMountIface")
  checkPtrType(object, "GMount")

  w <- .RGtkCall("S_gmount_iface_get_volume", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gMountIfaceGetDrive <-
function(object.class, object)
{
  checkPtrType(object.class, "GMountIface")
  checkPtrType(object, "GMount")

  w <- .RGtkCall("S_gmount_iface_get_drive", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gMountIfaceCanUnmount <-
function(object.class, object)
{
  checkPtrType(object.class, "GMountIface")
  checkPtrType(object, "GMount")

  w <- .RGtkCall("S_gmount_iface_can_unmount", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gMountIfaceCanEject <-
function(object.class, object)
{
  checkPtrType(object.class, "GMountIface")
  checkPtrType(object, "GMount")

  w <- .RGtkCall("S_gmount_iface_can_eject", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gMountIfaceUnmount <-
function(object.class, object, flags = "G_MOUNT_UNMOUNT_NONE", cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GMountIface")
  checkPtrType(object, "GMount")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gmount_iface_unmount", object.class, object, flags, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gMountIfaceUnmountFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GMountIface")
  checkPtrType(object, "GMount")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_gmount_iface_unmount_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gMountIfaceEject <-
function(object.class, object, flags = "G_MOUNT_UNMOUNT_NONE", cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GMountIface")
  checkPtrType(object, "GMount")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gmount_iface_eject", object.class, object, flags, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gMountIfaceEjectFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GMountIface")
  checkPtrType(object, "GMount")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_gmount_iface_eject_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gMountIfaceRemount <-
function(object.class, object, flags, mount.operation, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GMountIface")
  checkPtrType(object, "GMount")
  
  checkPtrType(mount.operation, "GMountOperation")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gmount_iface_remount", object.class, object, flags, mount.operation, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gMountIfaceRemountFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GMountIface")
  checkPtrType(object, "GMount")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_gmount_iface_remount_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gOutputStreamClassWriteFn <-
function(object.class, object, buffer, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GOutputStreamClass")
  checkPtrType(object, "GOutputStream")
  buffer <- as.list(as.raw(buffer))
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_goutput_stream_class_write_fn", object.class, object, buffer, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gOutputStreamClassSplice <-
function(object.class, object, source, flags = "G_OUTPUT_STREAM_SPLICE_NONE", cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GOutputStreamClass")
  checkPtrType(object, "GOutputStream")
  checkPtrType(source, "GInputStream")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_goutput_stream_class_splice", object.class, object, source, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gOutputStreamClassFlush <-
function(object.class, object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GOutputStreamClass")
  checkPtrType(object, "GOutputStream")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_goutput_stream_class_flush", object.class, object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gOutputStreamClassCloseFn <-
function(object.class, object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GOutputStreamClass")
  checkPtrType(object, "GOutputStream")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_goutput_stream_class_close_fn", object.class, object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gOutputStreamClassWriteAsync <-
function(object.class, object, buffer, io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GOutputStreamClass")
  checkPtrType(object, "GOutputStream")
  buffer <- as.list(as.raw(buffer))
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_goutput_stream_class_write_async", object.class, object, buffer, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gOutputStreamClassWriteFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GOutputStreamClass")
  checkPtrType(object, "GOutputStream")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_goutput_stream_class_write_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gOutputStreamClassSpliceAsync <-
function(object.class, object, source, flags = "G_OUTPUT_STREAM_SPLICE_NONE", io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GOutputStreamClass")
  checkPtrType(object, "GOutputStream")
  checkPtrType(source, "GInputStream")
  
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_goutput_stream_class_splice_async", object.class, object, source, flags, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gOutputStreamClassSpliceFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GOutputStreamClass")
  checkPtrType(object, "GOutputStream")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_goutput_stream_class_splice_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gOutputStreamClassFlushAsync <-
function(object.class, object, io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GOutputStreamClass")
  checkPtrType(object, "GOutputStream")
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_goutput_stream_class_flush_async", object.class, object, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gOutputStreamClassFlushFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GOutputStreamClass")
  checkPtrType(object, "GOutputStream")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_goutput_stream_class_flush_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gOutputStreamClassCloseAsync <-
function(object.class, object, io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GOutputStreamClass")
  checkPtrType(object, "GOutputStream")
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_goutput_stream_class_close_async", object.class, object, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gOutputStreamClassCloseFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GOutputStreamClass")
  checkPtrType(object, "GOutputStream")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_goutput_stream_class_close_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gSeekableIfaceTell <-
function(object.class, object)
{
  checkPtrType(object.class, "GSeekableIface")
  checkPtrType(object, "GSeekable")

  w <- .RGtkCall("S_gseekable_iface_tell", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gSeekableIfaceCanSeek <-
function(object.class, object)
{
  checkPtrType(object.class, "GSeekableIface")
  checkPtrType(object, "GSeekable")

  w <- .RGtkCall("S_gseekable_iface_can_seek", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gSeekableIfaceSeek <-
function(object.class, object, offset, type = "G_SEEK_SET", cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GSeekableIface")
  checkPtrType(object, "GSeekable")
  offset <- as.numeric(offset)
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gseekable_iface_seek", object.class, object, offset, type, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gSeekableIfaceCanTruncate <-
function(object.class, object)
{
  checkPtrType(object.class, "GSeekableIface")
  checkPtrType(object, "GSeekable")

  w <- .RGtkCall("S_gseekable_iface_can_truncate", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gSeekableIfaceTruncateFn <-
function(object.class, object, offset, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GSeekableIface")
  checkPtrType(object, "GSeekable")
  offset <- as.numeric(offset)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gseekable_iface_truncate_fn", object.class, object, offset, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gVfsClassIsActive <-
function(object.class, object)
{
  checkPtrType(object.class, "GVfsClass")
  checkPtrType(object, "GVfs")

  w <- .RGtkCall("S_gvfs_class_is_active", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gVfsClassGetFileForPath <-
function(object.class, object, path)
{
  checkPtrType(object.class, "GVfsClass")
  checkPtrType(object, "GVfs")
  path <- as.character(path)

  w <- .RGtkCall("S_gvfs_class_get_file_for_path", object.class, object, path, PACKAGE = "RGtk2")

  return(w)
}

gVfsClassGetFileForUri <-
function(object.class, object, uri)
{
  checkPtrType(object.class, "GVfsClass")
  checkPtrType(object, "GVfs")
  uri <- as.character(uri)

  w <- .RGtkCall("S_gvfs_class_get_file_for_uri", object.class, object, uri, PACKAGE = "RGtk2")

  return(w)
}

gVfsClassParseName <-
function(object.class, object, parse.name)
{
  checkPtrType(object.class, "GVfsClass")
  checkPtrType(object, "GVfs")
  parse.name <- as.character(parse.name)

  w <- .RGtkCall("S_gvfs_class_parse_name", object.class, object, parse.name, PACKAGE = "RGtk2")

  return(w)
}

gVfsClassGetSupportedUriSchemes <-
function(object.class, object)
{
  checkPtrType(object.class, "GVfsClass")
  checkPtrType(object, "GVfs")

  w <- .RGtkCall("S_gvfs_class_get_supported_uri_schemes", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gVolumeIfaceGetName <-
function(object.class, object)
{
  checkPtrType(object.class, "GVolumeIface")
  checkPtrType(object, "GVolume")

  w <- .RGtkCall("S_gvolume_iface_get_name", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gVolumeIfaceGetIcon <-
function(object.class, object)
{
  checkPtrType(object.class, "GVolumeIface")
  checkPtrType(object, "GVolume")

  w <- .RGtkCall("S_gvolume_iface_get_icon", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gVolumeIfaceGetUuid <-
function(object.class, object)
{
  checkPtrType(object.class, "GVolumeIface")
  checkPtrType(object, "GVolume")

  w <- .RGtkCall("S_gvolume_iface_get_uuid", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gVolumeIfaceGetDrive <-
function(object.class, object)
{
  checkPtrType(object.class, "GVolumeIface")
  checkPtrType(object, "GVolume")

  w <- .RGtkCall("S_gvolume_iface_get_drive", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gVolumeIfaceGetMount <-
function(object.class, object)
{
  checkPtrType(object.class, "GVolumeIface")
  checkPtrType(object, "GVolume")

  w <- .RGtkCall("S_gvolume_iface_get_mount", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gVolumeIfaceCanMount <-
function(object.class, object)
{
  checkPtrType(object.class, "GVolumeIface")
  checkPtrType(object, "GVolume")

  w <- .RGtkCall("S_gvolume_iface_can_mount", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gVolumeIfaceCanEject <-
function(object.class, object)
{
  checkPtrType(object.class, "GVolumeIface")
  checkPtrType(object, "GVolume")

  w <- .RGtkCall("S_gvolume_iface_can_eject", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gVolumeIfaceShouldAutomount <-
function(object.class, object)
{
  checkPtrType(object.class, "GVolumeIface")
  checkPtrType(object, "GVolume")

  w <- .RGtkCall("S_gvolume_iface_should_automount", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gVolumeIfaceMountFn <-
function(object.class, object, flags, mount.operation, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GVolumeIface")
  checkPtrType(object, "GVolume")
  
  checkPtrType(mount.operation, "GMountOperation")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gvolume_iface_mount_fn", object.class, object, flags, mount.operation, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gVolumeIfaceMountFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GVolumeIface")
  checkPtrType(object, "GVolume")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_gvolume_iface_mount_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gVolumeIfaceEject <-
function(object.class, object, flags = "G_MOUNT_UNMOUNT_NONE", cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GVolumeIface")
  checkPtrType(object, "GVolume")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gvolume_iface_eject", object.class, object, flags, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gVolumeIfaceEjectFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GVolumeIface")
  checkPtrType(object, "GVolume")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_gvolume_iface_eject_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gVolumeMonitorClassGetConnectedDrives <-
function(object.class, object)
{
  checkPtrType(object.class, "GVolumeMonitorClass")
  checkPtrType(object, "GVolumeMonitor")

  w <- .RGtkCall("S_gvolume_monitor_class_get_connected_drives", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gVolumeMonitorClassGetVolumes <-
function(object.class, object)
{
  checkPtrType(object.class, "GVolumeMonitorClass")
  checkPtrType(object, "GVolumeMonitor")

  w <- .RGtkCall("S_gvolume_monitor_class_get_volumes", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gVolumeMonitorClassGetMounts <-
function(object.class, object)
{
  checkPtrType(object.class, "GVolumeMonitorClass")
  checkPtrType(object, "GVolumeMonitor")

  w <- .RGtkCall("S_gvolume_monitor_class_get_mounts", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gVolumeMonitorClassGetVolumeForUuid <-
function(object.class, object, uuid)
{
  checkPtrType(object.class, "GVolumeMonitorClass")
  checkPtrType(object, "GVolumeMonitor")
  uuid <- as.character(uuid)

  w <- .RGtkCall("S_gvolume_monitor_class_get_volume_for_uuid", object.class, object, uuid, PACKAGE = "RGtk2")

  return(w)
}

gVolumeMonitorClassGetMountForUuid <-
function(object.class, object, uuid)
{
  checkPtrType(object.class, "GVolumeMonitorClass")
  checkPtrType(object, "GVolumeMonitor")
  uuid <- as.character(uuid)

  w <- .RGtkCall("S_gvolume_monitor_class_get_mount_for_uuid", object.class, object, uuid, PACKAGE = "RGtk2")

  return(w)
}

gMountIfaceGuessContentType <-
function(object.class, object, force.rescan, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GMountIface")
  checkPtrType(object, "GMount")
  force.rescan <- as.logical(force.rescan)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gmount_iface_guess_content_type", object.class, object, force.rescan, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gMountIfaceGuessContentTypeFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GMountIface")
  checkPtrType(object, "GMount")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_gmount_iface_guess_content_type_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gMountIfaceGuessContentTypeSync <-
function(object.class, object, force.rescan, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GMountIface")
  checkPtrType(object, "GMount")
  force.rescan <- as.logical(force.rescan)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gmount_iface_guess_content_type_sync", object.class, object, force.rescan, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gVolumeIfaceGetIdentifier <-
function(object.class, object, kind)
{
  checkPtrType(object.class, "GVolumeIface")
  checkPtrType(object, "GVolume")
  kind <- as.character(kind)

  w <- .RGtkCall("S_gvolume_iface_get_identifier", object.class, object, kind, PACKAGE = "RGtk2")

  return(w)
}

gVolumeIfaceEnumerateIdentifiers <-
function(object.class, object)
{
  checkPtrType(object.class, "GVolumeIface")
  checkPtrType(object, "GVolume")

  w <- .RGtkCall("S_gvolume_iface_enumerate_identifiers", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gVolumeIfaceGetActivationRoot <-
function(object.class, object)
{
  checkPtrType(object.class, "GVolumeIface")
  checkPtrType(object, "GVolume")

  w <- .RGtkCall("S_gvolume_iface_get_activation_root", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gAppInfoIfaceGetCommandline <-
function(object.class, object)
{
  checkPtrType(object.class, "GAppInfoIface")
  checkPtrType(object, "GAppInfo")

  w <- .RGtkCall("S_gapp_info_iface_get_commandline", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gAsyncInitableIfaceInitAsync <-
function(object.class, object, io.priority, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GAsyncInitableIface")
  checkPtrType(object, "GAsyncInitable")
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gasync_initable_iface_init_async", object.class, object, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gAsyncInitableIfaceInitFinish <-
function(object.class, object, res, .errwarn = TRUE)
{
  checkPtrType(object.class, "GAsyncInitableIface")
  checkPtrType(object, "GAsyncInitable")
  checkPtrType(res, "GAsyncResult")

  w <- .RGtkCall("S_gasync_initable_iface_init_finish", object.class, object, res, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gDriveIfaceGetStartStopType <-
function(object.class, object)
{
  checkPtrType(object.class, "GDriveIface")
  checkPtrType(object, "GDrive")

  w <- .RGtkCall("S_gdrive_iface_get_start_stop_type", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gDriveIfaceStart <-
function(object.class, object, flags, mount.operation, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GDriveIface")
  checkPtrType(object, "GDrive")
  
  checkPtrType(mount.operation, "GMountOperation")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gdrive_iface_start", object.class, object, flags, mount.operation, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gDriveIfaceStartFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GDriveIface")
  checkPtrType(object, "GDrive")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_gdrive_iface_start_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gDriveIfaceStop <-
function(object.class, object, flags, mount.operation, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GDriveIface")
  checkPtrType(object, "GDrive")
  
  checkPtrType(mount.operation, "GMountOperation")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gdrive_iface_stop", object.class, object, flags, mount.operation, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gDriveIfaceStopFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GDriveIface")
  checkPtrType(object, "GDrive")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_gdrive_iface_stop_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gDriveIfaceCanStart <-
function(object.class, object)
{
  checkPtrType(object.class, "GDriveIface")
  checkPtrType(object, "GDrive")

  w <- .RGtkCall("S_gdrive_iface_can_start", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gDriveIfaceCanStartDegraded <-
function(object.class, object)
{
  checkPtrType(object.class, "GDriveIface")
  checkPtrType(object, "GDrive")

  w <- .RGtkCall("S_gdrive_iface_can_start_degraded", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gDriveIfaceCanStop <-
function(object.class, object)
{
  checkPtrType(object.class, "GDriveIface")
  checkPtrType(object, "GDrive")

  w <- .RGtkCall("S_gdrive_iface_can_stop", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gDriveIfaceEjectWithOperation <-
function(object.class, object, flags, mount.operation, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GDriveIface")
  checkPtrType(object, "GDrive")
  
  checkPtrType(mount.operation, "GMountOperation")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gdrive_iface_eject_with_operation", object.class, object, flags, mount.operation, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gDriveIfaceEjectWithOperationFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GDriveIface")
  checkPtrType(object, "GDrive")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_gdrive_iface_eject_with_operation_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceCreateReadwrite <-
function(object.class, object, flags, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gfile_iface_create_readwrite", object.class, object, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceCreateReadwriteAsync <-
function(object.class, object, flags, io.priority, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gfile_iface_create_readwrite_async", object.class, object, flags, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gFileIfaceCreateReadwriteFinish <-
function(object.class, object, res, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  checkPtrType(res, "GAsyncResult")

  w <- .RGtkCall("S_gfile_iface_create_readwrite_finish", object.class, object, res, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceEjectMountableWithOperation <-
function(object.class, object, flags, mount.operation, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  
  checkPtrType(mount.operation, "GMountOperation")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gfile_iface_eject_mountable_with_operation", object.class, object, flags, mount.operation, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gFileIfaceEjectMountableWithOperationFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_gfile_iface_eject_mountable_with_operation_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceOpenReadwrite <-
function(object.class, object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gfile_iface_open_readwrite", object.class, object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceOpenReadwriteAsync <-
function(object.class, object, io.priority, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gfile_iface_open_readwrite_async", object.class, object, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gFileIfaceOpenReadwriteFinish <-
function(object.class, object, res, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  checkPtrType(res, "GAsyncResult")

  w <- .RGtkCall("S_gfile_iface_open_readwrite_finish", object.class, object, res, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfacePollMountable <-
function(object.class, object, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gfile_iface_poll_mountable", object.class, object, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gFileIfacePollMountableFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_gfile_iface_poll_mountable_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceReplaceReadwrite <-
function(object.class, object, etag, make.backup, flags, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  etag <- as.character(etag)
  make.backup <- as.logical(make.backup)
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gfile_iface_replace_readwrite", object.class, object, etag, make.backup, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceReplaceReadwriteAsync <-
function(object.class, object, etag, make.backup, flags, io.priority, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  etag <- as.character(etag)
  make.backup <- as.logical(make.backup)
  
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gfile_iface_replace_readwrite_async", object.class, object, etag, make.backup, flags, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gFileIfaceReplaceReadwriteFinish <-
function(object.class, object, res, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  checkPtrType(res, "GAsyncResult")

  w <- .RGtkCall("S_gfile_iface_replace_readwrite_finish", object.class, object, res, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceStartMountable <-
function(object.class, object, flags, start.operation, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  
  checkPtrType(start.operation, "GMountOperation")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gfile_iface_start_mountable", object.class, object, flags, start.operation, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gFileIfaceStartMountableFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_gfile_iface_start_mountable_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceStopMountable <-
function(object.class, object, flags, mount.operation, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  
  checkPtrType(mount.operation, "GMountOperation")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gfile_iface_stop_mountable", object.class, object, flags, mount.operation, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gFileIfaceStopMountableFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_gfile_iface_stop_mountable_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIfaceUnmountMountableWithOperation <-
function(object.class, object, flags, mount.operation, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  
  checkPtrType(mount.operation, "GMountOperation")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gfile_iface_unmount_mountable_with_operation", object.class, object, flags, mount.operation, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gFileIfaceUnmountMountableWithOperationFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIface")
  checkPtrType(object, "GFile")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_gfile_iface_unmount_mountable_with_operation_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIOStreamClassQueryInfo <-
function(object.class, object, attributes, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIOStreamClass")
  checkPtrType(object, "GFileIOStream")
  attributes <- as.character(attributes)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gfile_iostream_class_query_info", object.class, object, attributes, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIOStreamClassQueryInfoAsync <-
function(object.class, object, attributes, io.priority, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GFileIOStreamClass")
  checkPtrType(object, "GFileIOStream")
  attributes <- as.character(attributes)
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gfile_iostream_class_query_info_async", object.class, object, attributes, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gFileIOStreamClassQueryInfoFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GFileIOStreamClass")
  checkPtrType(object, "GFileIOStream")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_gfile_iostream_class_query_info_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gFileIOStreamClassGetEtag <-
function(object.class, object)
{
  checkPtrType(object.class, "GFileIOStreamClass")
  checkPtrType(object, "GFileIOStream")

  w <- .RGtkCall("S_gfile_iostream_class_get_etag", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gInetAddressClassToString <-
function(object.class, object)
{
  checkPtrType(object.class, "GInetAddressClass")
  checkPtrType(object, "GInetAddress")

  w <- .RGtkCall("S_ginet_address_class_to_string", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gInetAddressClassToBytes <-
function(object.class, object)
{
  checkPtrType(object.class, "GInetAddressClass")
  checkPtrType(object, "GInetAddress")

  w <- .RGtkCall("S_ginet_address_class_to_bytes", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gInitableIfaceInit <-
function(object.class, object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GInitableIface")
  checkPtrType(object, "GInitable")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_ginitable_iface_init", object.class, object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gIOStreamClassGetInputStream <-
function(object.class, object)
{
  checkPtrType(object.class, "GIOStreamClass")
  checkPtrType(object, "GIOStream")

  w <- .RGtkCall("S_giostream_class_get_input_stream", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gIOStreamClassGetOutputStream <-
function(object.class, object)
{
  checkPtrType(object.class, "GIOStreamClass")
  checkPtrType(object, "GIOStream")

  w <- .RGtkCall("S_giostream_class_get_output_stream", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gIOStreamClassCloseFn <-
function(object.class, object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GIOStreamClass")
  checkPtrType(object, "GIOStream")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_giostream_class_close_fn", object.class, object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gIOStreamClassCloseAsync <-
function(object.class, object, io.priority, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GIOStreamClass")
  checkPtrType(object, "GIOStream")
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_giostream_class_close_async", object.class, object, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gIOStreamClassCloseFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GIOStreamClass")
  checkPtrType(object, "GIOStream")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_giostream_class_close_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gMountIfaceUnmountWithOperation <-
function(object.class, object, flags, mount.operation, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GMountIface")
  checkPtrType(object, "GMount")
  
  checkPtrType(mount.operation, "GMountOperation")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gmount_iface_unmount_with_operation", object.class, object, flags, mount.operation, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gMountIfaceUnmountWithOperationFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GMountIface")
  checkPtrType(object, "GMount")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_gmount_iface_unmount_with_operation_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gMountIfaceEjectWithOperation <-
function(object.class, object, flags, mount.operation, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GMountIface")
  checkPtrType(object, "GMount")
  
  checkPtrType(mount.operation, "GMountOperation")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gmount_iface_eject_with_operation", object.class, object, flags, mount.operation, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gMountIfaceEjectWithOperationFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GMountIface")
  checkPtrType(object, "GMount")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_gmount_iface_eject_with_operation_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gResolverClassLookupByName <-
function(object.class, object, hostname, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GResolverClass")
  checkPtrType(object, "GResolver")
  hostname <- as.character(hostname)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gresolver_class_lookup_by_name", object.class, object, hostname, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gResolverClassLookupByNameAsync <-
function(object.class, object, hostname, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GResolverClass")
  checkPtrType(object, "GResolver")
  hostname <- as.character(hostname)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gresolver_class_lookup_by_name_async", object.class, object, hostname, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gResolverClassLookupByNameFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GResolverClass")
  checkPtrType(object, "GResolver")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_gresolver_class_lookup_by_name_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gResolverClassLookupByAddress <-
function(object.class, object, address, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GResolverClass")
  checkPtrType(object, "GResolver")
  checkPtrType(address, "GInetAddress")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gresolver_class_lookup_by_address", object.class, object, address, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gResolverClassLookupByAddressAsync <-
function(object.class, object, address, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GResolverClass")
  checkPtrType(object, "GResolver")
  checkPtrType(address, "GInetAddress")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gresolver_class_lookup_by_address_async", object.class, object, address, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gResolverClassLookupByAddressFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GResolverClass")
  checkPtrType(object, "GResolver")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_gresolver_class_lookup_by_address_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gResolverClassLookupService <-
function(object.class, object, rrname, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GResolverClass")
  checkPtrType(object, "GResolver")
  rrname <- as.character(rrname)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gresolver_class_lookup_service", object.class, object, rrname, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gResolverClassLookupServiceAsync <-
function(object.class, object, rrname, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GResolverClass")
  checkPtrType(object, "GResolver")
  rrname <- as.character(rrname)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gresolver_class_lookup_service_async", object.class, object, rrname, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gResolverClassLookupServiceFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GResolverClass")
  checkPtrType(object, "GResolver")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_gresolver_class_lookup_service_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gSocketAddressEnumeratorClassNext <-
function(object.class, object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object.class, "GSocketAddressEnumeratorClass")
  checkPtrType(object, "GSocketAddressEnumerator")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_gsocket_address_enumerator_class_next", object.class, object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gSocketAddressEnumeratorClassNextAsync <-
function(object.class, object, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GSocketAddressEnumeratorClass")
  checkPtrType(object, "GSocketAddressEnumerator")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gsocket_address_enumerator_class_next_async", object.class, object, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gSocketAddressEnumeratorClassNextFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GSocketAddressEnumeratorClass")
  checkPtrType(object, "GSocketAddressEnumerator")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_gsocket_address_enumerator_class_next_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gSocketAddressClassGetFamily <-
function(object.class, object)
{
  checkPtrType(object.class, "GSocketAddressClass")
  checkPtrType(object, "GSocketAddress")

  w <- .RGtkCall("S_gsocket_address_class_get_family", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gSocketAddressClassToNative <-
function(object.class, object, dest, destlen, .errwarn = TRUE)
{
  checkPtrType(object.class, "GSocketAddressClass")
  checkPtrType(object, "GSocketAddress")
  
  destlen <- as.numeric(destlen)

  w <- .RGtkCall("S_gsocket_address_class_to_native", object.class, object, dest, destlen, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}

gSocketAddressClassGetNativeSize <-
function(object.class, object)
{
  checkPtrType(object.class, "GSocketAddressClass")
  checkPtrType(object, "GSocketAddress")

  w <- .RGtkCall("S_gsocket_address_class_get_native_size", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gSocketConnectableIfaceEnumerate <-
function(object.class, object)
{
  checkPtrType(object.class, "GSocketConnectableIface")
  checkPtrType(object, "GSocketConnectable")

  w <- .RGtkCall("S_gsocket_connectable_iface_enumerate", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gSocketControlMessageClassGetSize <-
function(object.class, object)
{
  checkPtrType(object.class, "GSocketControlMessageClass")
  checkPtrType(object, "GSocketControlMessage")

  w <- .RGtkCall("S_gsocket_control_message_class_get_size", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gSocketControlMessageClassGetLevel <-
function(object.class, object)
{
  checkPtrType(object.class, "GSocketControlMessageClass")
  checkPtrType(object, "GSocketControlMessage")

  w <- .RGtkCall("S_gsocket_control_message_class_get_level", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gSocketControlMessageClassGetType <-
function(object.class, object)
{
  checkPtrType(object.class, "GSocketControlMessageClass")
  checkPtrType(object, "GSocketControlMessage")

  w <- .RGtkCall("S_gsocket_control_message_class_get_type", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gSocketControlMessageClassSerialize <-
function(object.class, object, data)
{
  checkPtrType(object.class, "GSocketControlMessageClass")
  checkPtrType(object, "GSocketControlMessage")
  

  w <- .RGtkCall("S_gsocket_control_message_class_serialize", object.class, object, data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gSocketListenerClassChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "GSocketListenerClass")
  checkPtrType(object, "GSocketListener")

  w <- .RGtkCall("S_gsocket_listener_class_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gVolumeIfaceEjectWithOperation <-
function(object.class, object, flags, mount.operation, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object.class, "GVolumeIface")
  checkPtrType(object, "GVolume")
  
  checkPtrType(mount.operation, "GMountOperation")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_gvolume_iface_eject_with_operation", object.class, object, flags, mount.operation, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
}

gVolumeIfaceEjectWithOperationFinish <-
function(object.class, object, result, .errwarn = TRUE)
{
  checkPtrType(object.class, "GVolumeIface")
  checkPtrType(object, "GVolume")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_gvolume_iface_eject_with_operation_finish", object.class, object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
}