
gAppInfoGetType <-
function()
{
  

  w <- .RGtkCall("S_g_app_info_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gAppInfoLaunchDefaultForUri <-
function(uri, launch.context, .errwarn = TRUE)
{
  uri <- as.character(uri)
  checkPtrType(launch.context, "GAppLaunchContext")

  w <- .RGtkCall("S_g_app_info_launch_default_for_uri", uri, launch.context, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gAppLaunchContextGetType <-
function()
{
  

  w <- .RGtkCall("S_g_app_launch_context_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gAppInfoCreateFromCommandline <-
function(commandline, application.name = NULL, flags = "G_APP_INFO_CREATE_NONE", .errwarn = TRUE)
{
  commandline <- as.character(commandline)
  if (!is.null( application.name )) application.name <- as.character(application.name)
  

  w <- .RGtkCall("S_g_app_info_create_from_commandline", commandline, application.name, flags, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gAppInfoDup <-
function(object)
{
  checkPtrType(object, "GAppInfo")

  w <- .RGtkCall("S_g_app_info_dup", object, PACKAGE = "RGtk2")

  return(w)
} 


gAppInfoEqual <-
function(object, appinfo2)
{
  checkPtrType(object, "GAppInfo")
  checkPtrType(appinfo2, "GAppInfo")

  w <- .RGtkCall("S_g_app_info_equal", object, appinfo2, PACKAGE = "RGtk2")

  return(w)
} 


gAppInfoGetId <-
function(object)
{
  checkPtrType(object, "GAppInfo")

  w <- .RGtkCall("S_g_app_info_get_id", object, PACKAGE = "RGtk2")

  return(w)
} 


gAppInfoGetName <-
function(object)
{
  checkPtrType(object, "GAppInfo")

  w <- .RGtkCall("S_g_app_info_get_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gAppInfoGetDescription <-
function(object)
{
  checkPtrType(object, "GAppInfo")

  w <- .RGtkCall("S_g_app_info_get_description", object, PACKAGE = "RGtk2")

  return(w)
} 


gAppInfoGetExecutable <-
function(object)
{
  checkPtrType(object, "GAppInfo")

  w <- .RGtkCall("S_g_app_info_get_executable", object, PACKAGE = "RGtk2")

  return(w)
} 


gAppInfoGetIcon <-
function(object)
{
  checkPtrType(object, "GAppInfo")

  w <- .RGtkCall("S_g_app_info_get_icon", object, PACKAGE = "RGtk2")

  return(w)
} 


gAppInfoLaunch <-
function(object, files, launch.context, .errwarn = TRUE)
{
  checkPtrType(object, "GAppInfo")
  files <- as.GList(files)
  checkPtrType(launch.context, "GAppLaunchContext")

  w <- .RGtkCall("S_g_app_info_launch", object, files, launch.context, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gAppInfoSupportsUris <-
function(object)
{
  checkPtrType(object, "GAppInfo")

  w <- .RGtkCall("S_g_app_info_supports_uris", object, PACKAGE = "RGtk2")

  return(w)
} 


gAppInfoSupportsFiles <-
function(object)
{
  checkPtrType(object, "GAppInfo")

  w <- .RGtkCall("S_g_app_info_supports_files", object, PACKAGE = "RGtk2")

  return(w)
} 


gAppInfoLaunchUris <-
function(object, uris, launch.context, .errwarn = TRUE)
{
  checkPtrType(object, "GAppInfo")
  uris <- as.GList(uris)
  checkPtrType(launch.context, "GAppLaunchContext")

  w <- .RGtkCall("S_g_app_info_launch_uris", object, uris, launch.context, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gAppInfoShouldShow <-
function(object)
{
  checkPtrType(object, "GAppInfo")

  w <- .RGtkCall("S_g_app_info_should_show", object, PACKAGE = "RGtk2")

  return(w)
} 


gAppInfoSetAsDefaultForType <-
function(object, content.type, .errwarn = TRUE)
{
  checkPtrType(object, "GAppInfo")
  content.type <- as.character(content.type)

  w <- .RGtkCall("S_g_app_info_set_as_default_for_type", object, content.type, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gAppInfoSetAsDefaultForExtension <-
function(object, extension, .errwarn = TRUE)
{
  checkPtrType(object, "GAppInfo")
  extension <- as.character(extension)

  w <- .RGtkCall("S_g_app_info_set_as_default_for_extension", object, extension, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gAppInfoAddSupportsType <-
function(object, content.type, .errwarn = TRUE)
{
  checkPtrType(object, "GAppInfo")
  content.type <- as.character(content.type)

  w <- .RGtkCall("S_g_app_info_add_supports_type", object, content.type, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gAppInfoCanRemoveSupportsType <-
function(object)
{
  checkPtrType(object, "GAppInfo")

  w <- .RGtkCall("S_g_app_info_can_remove_supports_type", object, PACKAGE = "RGtk2")

  return(w)
} 


gAppInfoRemoveSupportsType <-
function(object, content.type, .errwarn = TRUE)
{
  checkPtrType(object, "GAppInfo")
  content.type <- as.character(content.type)

  w <- .RGtkCall("S_g_app_info_remove_supports_type", object, content.type, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gAppInfoGetAll <-
function()
{
  

  w <- .RGtkCall("S_g_app_info_get_all", PACKAGE = "RGtk2")

  return(w)
} 


gAppInfoGetAllForType <-
function(content.type)
{
  content.type <- as.character(content.type)

  w <- .RGtkCall("S_g_app_info_get_all_for_type", content.type, PACKAGE = "RGtk2")

  return(w)
} 


gAppInfoGetDefaultForType <-
function(content.type, must.support.uris)
{
  content.type <- as.character(content.type)
  must.support.uris <- as.logical(must.support.uris)

  w <- .RGtkCall("S_g_app_info_get_default_for_type", content.type, must.support.uris, PACKAGE = "RGtk2")

  return(w)
} 


gAppInfoGetDefaultForUriScheme <-
function(uri.scheme)
{
  uri.scheme <- as.character(uri.scheme)

  w <- .RGtkCall("S_g_app_info_get_default_for_uri_scheme", uri.scheme, PACKAGE = "RGtk2")

  return(w)
} 


gAppLaunchContextNew <-
function()
{
  

  w <- .RGtkCall("S_g_app_launch_context_new", PACKAGE = "RGtk2")

  return(w)
} 


gAppLaunchContextGetDisplay <-
function(object, info, files)
{
  checkPtrType(object, "GAppLaunchContext")
  checkPtrType(info, "GAppInfo")
  files <- as.GList(files)

  w <- .RGtkCall("S_g_app_launch_context_get_display", object, info, files, PACKAGE = "RGtk2")

  return(w)
} 


gAppLaunchContextGetStartupNotifyId <-
function(object, info, files)
{
  checkPtrType(object, "GAppLaunchContext")
  checkPtrType(info, "GAppInfo")
  files <- as.GList(files)

  w <- .RGtkCall("S_g_app_launch_context_get_startup_notify_id", object, info, files, PACKAGE = "RGtk2")

  return(w)
} 


gAppLaunchContextLaunchFailed <-
function(object, startup.notify.id)
{
  checkPtrType(object, "GAppLaunchContext")
  startup.notify.id <- as.character(startup.notify.id)

  w <- .RGtkCall("S_g_app_launch_context_launch_failed", object, startup.notify.id, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gAsyncResultGetType <-
function()
{
  

  w <- .RGtkCall("S_g_async_result_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gAsyncResultGetUserData <-
function(object)
{
  checkPtrType(object, "GAsyncResult")

  w <- .RGtkCall("S_g_async_result_get_user_data", object, PACKAGE = "RGtk2")

  return(w)
} 


gAsyncResultGetSourceObject <-
function(object)
{
  checkPtrType(object, "GAsyncResult")

  w <- .RGtkCall("S_g_async_result_get_source_object", object, PACKAGE = "RGtk2")

  return(w)
} 


gBufferedInputStreamGetType <-
function()
{
  

  w <- .RGtkCall("S_g_buffered_input_stream_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gBufferedInputStreamNew <-
function(base.stream = NULL)
{
  

  w <- .RGtkCall("S_g_buffered_input_stream_new", base.stream, PACKAGE = "RGtk2")

  return(w)
} 


gBufferedInputStreamNewSized <-
function(base.stream, size)
{
  checkPtrType(base.stream, "GInputStream")
  size <- as.numeric(size)

  w <- .RGtkCall("S_g_buffered_input_stream_new_sized", base.stream, size, PACKAGE = "RGtk2")

  return(w)
} 


gBufferedInputStreamGetBufferSize <-
function(object)
{
  checkPtrType(object, "GBufferedInputStream")

  w <- .RGtkCall("S_g_buffered_input_stream_get_buffer_size", object, PACKAGE = "RGtk2")

  return(w)
} 


gBufferedInputStreamSetBufferSize <-
function(object, size)
{
  checkPtrType(object, "GBufferedInputStream")
  size <- as.numeric(size)

  w <- .RGtkCall("S_g_buffered_input_stream_set_buffer_size", object, size, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gBufferedInputStreamGetAvailable <-
function(object)
{
  checkPtrType(object, "GBufferedInputStream")

  w <- .RGtkCall("S_g_buffered_input_stream_get_available", object, PACKAGE = "RGtk2")

  return(w)
} 


gBufferedInputStreamPeekBuffer <-
function(object)
{
  checkPtrType(object, "GBufferedInputStream")

  w <- .RGtkCall("S_g_buffered_input_stream_peek_buffer", object, PACKAGE = "RGtk2")

  return(w)
} 


gBufferedInputStreamFill <-
function(object, count, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GBufferedInputStream")
  count <- as.integer(count)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_buffered_input_stream_fill", object, count, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gBufferedInputStreamFillAsync <-
function(object, count, io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GBufferedInputStream")
  count <- as.integer(count)
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_buffered_input_stream_fill_async", object, count, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gBufferedInputStreamFillFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GBufferedInputStream")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_buffered_input_stream_fill_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gBufferedInputStreamReadByte <-
function(object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GBufferedInputStream")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_buffered_input_stream_read_byte", object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gBufferedOutputStreamGetType <-
function()
{
  

  w <- .RGtkCall("S_g_buffered_output_stream_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gBufferedOutputStreamNew <-
function(base.stream = NULL)
{
  

  w <- .RGtkCall("S_g_buffered_output_stream_new", base.stream, PACKAGE = "RGtk2")

  return(w)
} 


gBufferedOutputStreamNewSized <-
function(base.stream, size)
{
  checkPtrType(base.stream, "GOutputStream")
  size <- as.numeric(size)

  w <- .RGtkCall("S_g_buffered_output_stream_new_sized", base.stream, size, PACKAGE = "RGtk2")

  return(w)
} 


gBufferedOutputStreamGetBufferSize <-
function(object)
{
  checkPtrType(object, "GBufferedOutputStream")

  w <- .RGtkCall("S_g_buffered_output_stream_get_buffer_size", object, PACKAGE = "RGtk2")

  return(w)
} 


gBufferedOutputStreamSetBufferSize <-
function(object, size)
{
  checkPtrType(object, "GBufferedOutputStream")
  size <- as.numeric(size)

  w <- .RGtkCall("S_g_buffered_output_stream_set_buffer_size", object, size, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gBufferedOutputStreamGetAutoGrow <-
function(object)
{
  checkPtrType(object, "GBufferedOutputStream")

  w <- .RGtkCall("S_g_buffered_output_stream_get_auto_grow", object, PACKAGE = "RGtk2")

  return(w)
} 


gBufferedOutputStreamSetAutoGrow <-
function(object, auto.grow)
{
  checkPtrType(object, "GBufferedOutputStream")
  auto.grow <- as.logical(auto.grow)

  w <- .RGtkCall("S_g_buffered_output_stream_set_auto_grow", object, auto.grow, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gCancellableGetType <-
function()
{
  

  w <- .RGtkCall("S_g_cancellable_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gCancellableNew <-
function()
{
  

  w <- .RGtkCall("S_g_cancellable_new", PACKAGE = "RGtk2")

  return(w)
} 


gCancellableIsCancelled <-
function(object)
{
  checkPtrType(object, "GCancellable")

  w <- .RGtkCall("S_g_cancellable_is_cancelled", object, PACKAGE = "RGtk2")

  return(w)
} 


gCancellableSetErrorIfCancelled <-
function(object, .errwarn = TRUE)
{
  checkPtrType(object, "GCancellable")

  w <- .RGtkCall("S_g_cancellable_set_error_if_cancelled", object, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gCancellableGetFd <-
function(object)
{
  checkPtrType(object, "GCancellable")

  w <- .RGtkCall("S_g_cancellable_get_fd", object, PACKAGE = "RGtk2")

  return(w)
} 


gCancellableGetCurrent <-
function()
{
  

  w <- .RGtkCall("S_g_cancellable_get_current", PACKAGE = "RGtk2")

  return(w)
} 


gCancellablePushCurrent <-
function(object)
{
  checkPtrType(object, "GCancellable")

  w <- .RGtkCall("S_g_cancellable_push_current", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gCancellablePopCurrent <-
function(object)
{
  checkPtrType(object, "GCancellable")

  w <- .RGtkCall("S_g_cancellable_pop_current", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gCancellableReset <-
function(object)
{
  checkPtrType(object, "GCancellable")

  w <- .RGtkCall("S_g_cancellable_reset", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gCancellableCancel <-
function(object)
{
  checkPtrType(object, "GCancellable")

  w <- .RGtkCall("S_g_cancellable_cancel", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gContentTypeEquals <-
function(type1, type2)
{
  type1 <- as.character(type1)
  type2 <- as.character(type2)

  w <- .RGtkCall("S_g_content_type_equals", type1, type2, PACKAGE = "RGtk2")

  return(w)
} 


gContentTypeIsA <-
function(type, supertype)
{
  type <- as.character(type)
  supertype <- as.character(supertype)

  w <- .RGtkCall("S_g_content_type_is_a", type, supertype, PACKAGE = "RGtk2")

  return(w)
} 


gContentTypeIsUnknown <-
function(type)
{
  type <- as.character(type)

  w <- .RGtkCall("S_g_content_type_is_unknown", type, PACKAGE = "RGtk2")

  return(w)
} 


gContentTypeGetDescription <-
function(type)
{
  type <- as.character(type)

  w <- .RGtkCall("S_g_content_type_get_description", type, PACKAGE = "RGtk2")

  return(w)
} 


gContentTypeGetMimeType <-
function(type)
{
  type <- as.character(type)

  w <- .RGtkCall("S_g_content_type_get_mime_type", type, PACKAGE = "RGtk2")

  return(w)
} 


gContentTypeGetIcon <-
function(type)
{
  type <- as.character(type)

  w <- .RGtkCall("S_g_content_type_get_icon", type, PACKAGE = "RGtk2")

  return(w)
} 


gContentTypeCanBeExecutable <-
function(type)
{
  type <- as.character(type)

  w <- .RGtkCall("S_g_content_type_can_be_executable", type, PACKAGE = "RGtk2")

  return(w)
} 


gContentTypeGuess <-
function(filename, data)
{
  filename <- as.character(filename)
  data <- as.list(as.raw(data))

  w <- .RGtkCall("S_g_content_type_guess", filename, data, PACKAGE = "RGtk2")

  return(w)
} 


gContentTypesGetRegistered <-
function()
{
  

  w <- .RGtkCall("S_g_content_types_get_registered", PACKAGE = "RGtk2")

  return(w)
} 


gDataInputStreamGetType <-
function()
{
  

  w <- .RGtkCall("S_g_data_input_stream_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gDataInputStreamNew <-
function(base.stream = NULL)
{
  

  w <- .RGtkCall("S_g_data_input_stream_new", base.stream, PACKAGE = "RGtk2")

  return(w)
} 


gDataInputStreamSetByteOrder <-
function(object, order)
{
  checkPtrType(object, "GDataInputStream")
  

  w <- .RGtkCall("S_g_data_input_stream_set_byte_order", object, order, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gDataInputStreamGetByteOrder <-
function(object)
{
  checkPtrType(object, "GDataInputStream")

  w <- .RGtkCall("S_g_data_input_stream_get_byte_order", object, PACKAGE = "RGtk2")

  return(w)
} 


gDataInputStreamSetNewlineType <-
function(object, type)
{
  checkPtrType(object, "GDataInputStream")
  

  w <- .RGtkCall("S_g_data_input_stream_set_newline_type", object, type, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gDataInputStreamGetNewlineType <-
function(object)
{
  checkPtrType(object, "GDataInputStream")

  w <- .RGtkCall("S_g_data_input_stream_get_newline_type", object, PACKAGE = "RGtk2")

  return(w)
} 


gDataInputStreamReadByte <-
function(object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GDataInputStream")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_data_input_stream_read_byte", object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gDataInputStreamReadInt16 <-
function(object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GDataInputStream")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_data_input_stream_read_int16", object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gDataInputStreamReadUint16 <-
function(object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GDataInputStream")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_data_input_stream_read_uint16", object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gDataInputStreamReadInt32 <-
function(object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GDataInputStream")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_data_input_stream_read_int32", object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gDataInputStreamReadUint32 <-
function(object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GDataInputStream")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_data_input_stream_read_uint32", object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gDataInputStreamReadInt64 <-
function(object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GDataInputStream")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_data_input_stream_read_int64", object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gDataInputStreamReadUint64 <-
function(object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GDataInputStream")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_data_input_stream_read_uint64", object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gDataInputStreamReadLine <-
function(object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GDataInputStream")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_data_input_stream_read_line", object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gDataInputStreamReadUntil <-
function(object, stop.chars, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GDataInputStream")
  stop.chars <- as.character(stop.chars)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_data_input_stream_read_until", object, stop.chars, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gDataOutputStreamGetType <-
function()
{
  

  w <- .RGtkCall("S_g_data_output_stream_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gDataOutputStreamNew <-
function(base.stream = NULL)
{
  

  w <- .RGtkCall("S_g_data_output_stream_new", base.stream, PACKAGE = "RGtk2")

  return(w)
} 


gDataOutputStreamSetByteOrder <-
function(object, order)
{
  checkPtrType(object, "GDataOutputStream")
  

  w <- .RGtkCall("S_g_data_output_stream_set_byte_order", object, order, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gDataOutputStreamGetByteOrder <-
function(object)
{
  checkPtrType(object, "GDataOutputStream")

  w <- .RGtkCall("S_g_data_output_stream_get_byte_order", object, PACKAGE = "RGtk2")

  return(w)
} 


gDataOutputStreamPutByte <-
function(object, data, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GDataOutputStream")
  data <- as.raw(data)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_data_output_stream_put_byte", object, data, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gDataOutputStreamPutInt16 <-
function(object, data, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GDataOutputStream")
  data <- as.integer(data)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_data_output_stream_put_int16", object, data, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gDataOutputStreamPutUint16 <-
function(object, data, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GDataOutputStream")
  data <- as.integer(data)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_data_output_stream_put_uint16", object, data, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gDataOutputStreamPutInt32 <-
function(object, data, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GDataOutputStream")
  data <- as.integer(data)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_data_output_stream_put_int32", object, data, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gDataOutputStreamPutUint32 <-
function(object, data, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GDataOutputStream")
  data <- as.numeric(data)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_data_output_stream_put_uint32", object, data, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gDataOutputStreamPutInt64 <-
function(object, data, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GDataOutputStream")
  data <- as.numeric(data)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_data_output_stream_put_int64", object, data, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gDataOutputStreamPutUint64 <-
function(object, data, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GDataOutputStream")
  data <- as.numeric(data)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_data_output_stream_put_uint64", object, data, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gDataOutputStreamPutString <-
function(object, str, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GDataOutputStream")
  str <- as.character(str)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_data_output_stream_put_string", object, str, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gDriveGetType <-
function()
{
  

  w <- .RGtkCall("S_g_drive_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gDriveGetName <-
function(object)
{
  checkPtrType(object, "GDrive")

  w <- .RGtkCall("S_g_drive_get_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gDriveGetIcon <-
function(object)
{
  checkPtrType(object, "GDrive")

  w <- .RGtkCall("S_g_drive_get_icon", object, PACKAGE = "RGtk2")

  return(w)
} 


gDriveHasVolumes <-
function(object)
{
  checkPtrType(object, "GDrive")

  w <- .RGtkCall("S_g_drive_has_volumes", object, PACKAGE = "RGtk2")

  return(w)
} 


gDriveGetVolumes <-
function(object)
{
  checkPtrType(object, "GDrive")

  w <- .RGtkCall("S_g_drive_get_volumes", object, PACKAGE = "RGtk2")

  return(w)
} 


gDriveIsMediaRemovable <-
function(object)
{
  checkPtrType(object, "GDrive")

  w <- .RGtkCall("S_g_drive_is_media_removable", object, PACKAGE = "RGtk2")

  return(w)
} 


gDriveHasMedia <-
function(object)
{
  checkPtrType(object, "GDrive")

  w <- .RGtkCall("S_g_drive_has_media", object, PACKAGE = "RGtk2")

  return(w)
} 


gDriveIsMediaCheckAutomatic <-
function(object)
{
  checkPtrType(object, "GDrive")

  w <- .RGtkCall("S_g_drive_is_media_check_automatic", object, PACKAGE = "RGtk2")

  return(w)
} 


gDriveCanPollForMedia <-
function(object)
{
  checkPtrType(object, "GDrive")

  w <- .RGtkCall("S_g_drive_can_poll_for_media", object, PACKAGE = "RGtk2")

  return(w)
} 


gDriveCanEject <-
function(object)
{
  checkPtrType(object, "GDrive")

  w <- .RGtkCall("S_g_drive_can_eject", object, PACKAGE = "RGtk2")

  return(w)
} 


gDriveEject <-
function(object, flags = "G_MOUNT_UNMOUNT_NONE", cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GDrive")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_drive_eject", object, flags, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gDriveEjectFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GDrive")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_drive_eject_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gDrivePollForMedia <-
function(object, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GDrive")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_drive_poll_for_media", object, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gDrivePollForMediaFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GDrive")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_drive_poll_for_media_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gDriveGetIdentifier <-
function(object, kind)
{
  checkPtrType(object, "GDrive")
  kind <- as.character(kind)

  w <- .RGtkCall("S_g_drive_get_identifier", object, kind, PACKAGE = "RGtk2")

  return(w)
} 


gDriveEnumerateIdentifiers <-
function(object)
{
  checkPtrType(object, "GDrive")

  w <- .RGtkCall("S_g_drive_enumerate_identifiers", object, PACKAGE = "RGtk2")

  return(w)
} 


gFileAttributeInfoListNew <-
function()
{
  

  w <- .RGtkCall("S_g_file_attribute_info_list_new", PACKAGE = "RGtk2")

  return(w)
} 


gFileAttributeInfoListLookup <-
function(object, name)
{
  checkPtrType(object, "GFileAttributeInfoList")
  name <- as.character(name)

  w <- .RGtkCall("S_g_file_attribute_info_list_lookup", object, name, PACKAGE = "RGtk2")

  return(w)
} 


gFileAttributeInfoListAdd <-
function(object, name, type, flags = "G_FILE_ATTRIBUTE_INFO_NONE")
{
  checkPtrType(object, "GFileAttributeInfoList")
  name <- as.character(name)
  
  

  w <- .RGtkCall("S_g_file_attribute_info_list_add", object, name, type, flags, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileEnumeratorGetType <-
function()
{
  

  w <- .RGtkCall("S_g_file_enumerator_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gFileEnumeratorNextFile <-
function(object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFileEnumerator")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_enumerator_next_file", object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileEnumeratorClose <-
function(object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFileEnumerator")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_enumerator_close", object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileEnumeratorNextFilesAsync <-
function(object, num.files, io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GFileEnumerator")
  num.files <- as.integer(num.files)
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_file_enumerator_next_files_async", object, num.files, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileEnumeratorNextFilesFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GFileEnumerator")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_file_enumerator_next_files_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileEnumeratorCloseAsync <-
function(object, io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GFileEnumerator")
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_file_enumerator_close_async", object, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileEnumeratorCloseFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GFileEnumerator")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_file_enumerator_close_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileEnumeratorIsClosed <-
function(object)
{
  checkPtrType(object, "GFileEnumerator")

  w <- .RGtkCall("S_g_file_enumerator_is_closed", object, PACKAGE = "RGtk2")

  return(w)
} 


gFileEnumeratorHasPending <-
function(object)
{
  checkPtrType(object, "GFileEnumerator")

  w <- .RGtkCall("S_g_file_enumerator_has_pending", object, PACKAGE = "RGtk2")

  return(w)
} 


gFileEnumeratorSetPending <-
function(object, pending)
{
  checkPtrType(object, "GFileEnumerator")
  pending <- as.logical(pending)

  w <- .RGtkCall("S_g_file_enumerator_set_pending", object, pending, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileGetType <-
function()
{
  

  w <- .RGtkCall("S_g_file_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gFileNewForPath <-
function(path)
{
  path <- as.character(path)

  w <- .RGtkCall("S_g_file_new_for_path", path, PACKAGE = "RGtk2")

  return(w)
} 


gFileNewForUri <-
function(uri)
{
  uri <- as.character(uri)

  w <- .RGtkCall("S_g_file_new_for_uri", uri, PACKAGE = "RGtk2")

  return(w)
} 


gFileNewForCommandlineArg <-
function(arg)
{
  arg <- as.character(arg)

  w <- .RGtkCall("S_g_file_new_for_commandline_arg", arg, PACKAGE = "RGtk2")

  return(w)
} 


gFileParseName <-
function(parse.name)
{
  parse.name <- as.character(parse.name)

  w <- .RGtkCall("S_g_file_parse_name", parse.name, PACKAGE = "RGtk2")

  return(w)
} 


gFileDup <-
function(object)
{
  checkPtrType(object, "GFile")

  w <- .RGtkCall("S_g_file_dup", object, PACKAGE = "RGtk2")

  return(w)
} 


gFileHash <-
function(file)
{
  

  w <- .RGtkCall("S_g_file_hash", file, PACKAGE = "RGtk2")

  return(w)
} 


gFileEqual <-
function(object, file2)
{
  checkPtrType(object, "GFile")
  checkPtrType(file2, "GFile")

  w <- .RGtkCall("S_g_file_equal", object, file2, PACKAGE = "RGtk2")

  return(w)
} 


gFileGetBasename <-
function(object)
{
  checkPtrType(object, "GFile")

  w <- .RGtkCall("S_g_file_get_basename", object, PACKAGE = "RGtk2")

  return(w)
} 


gFileGetPath <-
function(object)
{
  checkPtrType(object, "GFile")

  w <- .RGtkCall("S_g_file_get_path", object, PACKAGE = "RGtk2")

  return(w)
} 


gFileGetUri <-
function(object)
{
  checkPtrType(object, "GFile")

  w <- .RGtkCall("S_g_file_get_uri", object, PACKAGE = "RGtk2")

  return(w)
} 


gFileGetParseName <-
function(object)
{
  checkPtrType(object, "GFile")

  w <- .RGtkCall("S_g_file_get_parse_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gFileGetParent <-
function(object)
{
  checkPtrType(object, "GFile")

  w <- .RGtkCall("S_g_file_get_parent", object, PACKAGE = "RGtk2")

  return(w)
} 


gFileGetChild <-
function(object, name)
{
  checkPtrType(object, "GFile")
  name <- as.character(name)

  w <- .RGtkCall("S_g_file_get_child", object, name, PACKAGE = "RGtk2")

  return(w)
} 


gFileGetChildForDisplayName <-
function(object, display.name, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  display.name <- as.character(display.name)

  w <- .RGtkCall("S_g_file_get_child_for_display_name", object, display.name, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileHasPrefix <-
function(object, descendant)
{
  checkPtrType(object, "GFile")
  checkPtrType(descendant, "GFile")

  w <- .RGtkCall("S_g_file_has_prefix", object, descendant, PACKAGE = "RGtk2")

  return(w)
} 


gFileGetRelativePath <-
function(object, descendant)
{
  checkPtrType(object, "GFile")
  checkPtrType(descendant, "GFile")

  w <- .RGtkCall("S_g_file_get_relative_path", object, descendant, PACKAGE = "RGtk2")

  return(w)
} 


gFileResolveRelativePath <-
function(object, relative.path)
{
  checkPtrType(object, "GFile")
  relative.path <- as.character(relative.path)

  w <- .RGtkCall("S_g_file_resolve_relative_path", object, relative.path, PACKAGE = "RGtk2")

  return(w)
} 


gFileIsNative <-
function(object)
{
  checkPtrType(object, "GFile")

  w <- .RGtkCall("S_g_file_is_native", object, PACKAGE = "RGtk2")

  return(w)
} 


gFileHasUriScheme <-
function(object, uri.scheme)
{
  checkPtrType(object, "GFile")
  uri.scheme <- as.character(uri.scheme)

  w <- .RGtkCall("S_g_file_has_uri_scheme", object, uri.scheme, PACKAGE = "RGtk2")

  return(w)
} 


gFileGetUriScheme <-
function(object)
{
  checkPtrType(object, "GFile")

  w <- .RGtkCall("S_g_file_get_uri_scheme", object, PACKAGE = "RGtk2")

  return(w)
} 


gFileRead <-
function(object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_read", object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileReadAsync <-
function(object, io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GFile")
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_file_read_async", object, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileReadFinish <-
function(object, res, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  checkPtrType(res, "GAsyncResult")

  w <- .RGtkCall("S_g_file_read_finish", object, res, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileAppendTo <-
function(object, flags = "G_FILE_CREATE_NONE", cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_append_to", object, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileCreate <-
function(object, flags = "G_FILE_CREATE_NONE", cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_create", object, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileReplace <-
function(object, etag, make.backup, flags = "G_FILE_CREATE_NONE", cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  etag <- as.character(etag)
  make.backup <- as.logical(make.backup)
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_replace", object, etag, make.backup, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileAppendToAsync <-
function(object, flags = "G_FILE_CREATE_NONE", io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GFile")
  
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_file_append_to_async", object, flags, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileAppendToFinish <-
function(object, res, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  checkPtrType(res, "GAsyncResult")

  w <- .RGtkCall("S_g_file_append_to_finish", object, res, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileCreateAsync <-
function(object, flags = "G_FILE_CREATE_NONE", io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GFile")
  
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_file_create_async", object, flags, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileCreateFinish <-
function(object, res, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  checkPtrType(res, "GAsyncResult")

  w <- .RGtkCall("S_g_file_create_finish", object, res, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileReplaceAsync <-
function(object, etag, make.backup, flags = "G_FILE_CREATE_NONE", io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GFile")
  etag <- as.character(etag)
  make.backup <- as.logical(make.backup)
  
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_file_replace_async", object, etag, make.backup, flags, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileReplaceFinish <-
function(object, res, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  checkPtrType(res, "GAsyncResult")

  w <- .RGtkCall("S_g_file_replace_finish", object, res, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileQueryExists <-
function(object, cancellable = NULL)
{
  checkPtrType(object, "GFile")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_query_exists", object, cancellable, PACKAGE = "RGtk2")

  return(w)
} 


gFileQueryInfo <-
function(object, attributes, flags = "G_FILE_QUERY_INFO_NONE", cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  attributes <- as.character(attributes)
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_query_info", object, attributes, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileQueryInfoAsync <-
function(object, attributes, flags = "G_FILE_QUERY_INFO_NONE", io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GFile")
  attributes <- as.character(attributes)
  
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_file_query_info_async", object, attributes, flags, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileQueryInfoFinish <-
function(object, res, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  checkPtrType(res, "GAsyncResult")

  w <- .RGtkCall("S_g_file_query_info_finish", object, res, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileQueryFilesystemInfo <-
function(object, attributes, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  attributes <- as.character(attributes)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_query_filesystem_info", object, attributes, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileQueryFilesystemInfoAsync <-
function(object, attributes, io.priority, cancellable, callback, user.data = NULL)
{
  checkPtrType(object, "GFile")
  attributes <- as.character(attributes)
  io.priority <- as.integer(io.priority)
  checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_file_query_filesystem_info_async", object, attributes, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileQueryFilesystemInfoFinish <-
function(object, res, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  checkPtrType(res, "GAsyncResult")

  w <- .RGtkCall("S_g_file_query_filesystem_info_finish", object, res, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileFindEnclosingMount <-
function(object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_find_enclosing_mount", object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileFindEnclosingMountAsync <-
function(object, io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GFile")
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_file_find_enclosing_mount_async", object, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileFindEnclosingMountFinish <-
function(object, res, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  checkPtrType(res, "GAsyncResult")

  w <- .RGtkCall("S_g_file_find_enclosing_mount_finish", object, res, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileEnumerateChildren <-
function(object, attributes, flags = "G_FILE_QUERY_INFO_NONE", cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  attributes <- as.character(attributes)
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_enumerate_children", object, attributes, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileEnumerateChildrenAsync <-
function(object, attributes, flags = "G_FILE_QUERY_INFO_NONE", io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GFile")
  attributes <- as.character(attributes)
  
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_file_enumerate_children_async", object, attributes, flags, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileEnumerateChildrenFinish <-
function(object, res, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  checkPtrType(res, "GAsyncResult")

  w <- .RGtkCall("S_g_file_enumerate_children_finish", object, res, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileSetDisplayName <-
function(object, display.name, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  display.name <- as.character(display.name)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_set_display_name", object, display.name, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileSetDisplayNameAsync <-
function(object, display.name, io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GFile")
  display.name <- as.character(display.name)
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_file_set_display_name_async", object, display.name, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileSetDisplayNameFinish <-
function(object, res, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  checkPtrType(res, "GAsyncResult")

  w <- .RGtkCall("S_g_file_set_display_name_finish", object, res, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileDelete <-
function(object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_delete", object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileTrash <-
function(object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_trash", object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileCopy <-
function(object, destination, flags = "G_FILE_COPY_NONE", cancellable = NULL, progress.callback, progress.callback.data, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  checkPtrType(destination, "GFile")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  progress.callback <- as.function(progress.callback)
  

  w <- .RGtkCall("S_g_file_copy", object, destination, flags, cancellable, progress.callback, progress.callback.data, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileCopyAsync <-
function(object, destination, flags = "G_FILE_COPY_NONE", io.priority = 0, cancellable = NULL, progress.callback, progress.callback.data, callback, user.data = NULL)
{
  checkPtrType(object, "GFile")
  checkPtrType(destination, "GFile")
  
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  progress.callback <- as.function(progress.callback)
  
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_file_copy_async", object, destination, flags, io.priority, cancellable, progress.callback, progress.callback.data, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileCopyFinish <-
function(object, res, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  checkPtrType(res, "GAsyncResult")

  w <- .RGtkCall("S_g_file_copy_finish", object, res, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileMove <-
function(object, destination, flags = "G_FILE_COPY_NONE", cancellable = NULL, progress.callback, progress.callback.data, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  checkPtrType(destination, "GFile")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  progress.callback <- as.function(progress.callback)
  

  w <- .RGtkCall("S_g_file_move", object, destination, flags, cancellable, progress.callback, progress.callback.data, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileMakeDirectory <-
function(object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_make_directory", object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileMakeSymbolicLink <-
function(object, symlink.value, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  symlink.value <- as.character(symlink.value)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_make_symbolic_link", object, symlink.value, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileQuerySettableAttributes <-
function(object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_query_settable_attributes", object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileQueryWritableNamespaces <-
function(object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_query_writable_namespaces", object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileSetAttribute <-
function(object, attribute, type, value.p, flags = "G_FILE_QUERY_INFO_NONE", cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  attribute <- as.character(attribute)
  
  
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_set_attribute", object, attribute, type, value.p, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileSetAttributesFromInfo <-
function(object, info, flags = "G_FILE_QUERY_INFO_NONE", cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  checkPtrType(info, "GFileInfo")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_set_attributes_from_info", object, info, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileSetAttributesAsync <-
function(object, info, flags = "G_FILE_QUERY_INFO_NONE", io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GFile")
  checkPtrType(info, "GFileInfo")
  
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_file_set_attributes_async", object, info, flags, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileSetAttributesFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_file_set_attributes_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileSetAttributeString <-
function(object, attribute, value, flags = "G_FILE_QUERY_INFO_NONE", cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  attribute <- as.character(attribute)
  value <- as.character(value)
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_set_attribute_string", object, attribute, value, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileSetAttributeByteString <-
function(object, attribute, value, flags = "G_FILE_QUERY_INFO_NONE", cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  attribute <- as.character(attribute)
  value <- as.character(value)
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_set_attribute_byte_string", object, attribute, value, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileSetAttributeUint32 <-
function(object, attribute, value, flags = "G_FILE_QUERY_INFO_NONE", cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  attribute <- as.character(attribute)
  value <- as.numeric(value)
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_set_attribute_uint32", object, attribute, value, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileSetAttributeInt32 <-
function(object, attribute, value, flags = "G_FILE_QUERY_INFO_NONE", cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  attribute <- as.character(attribute)
  value <- as.integer(value)
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_set_attribute_int32", object, attribute, value, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileSetAttributeUint64 <-
function(object, attribute, value, flags = "G_FILE_QUERY_INFO_NONE", cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  attribute <- as.character(attribute)
  value <- as.numeric(value)
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_set_attribute_uint64", object, attribute, value, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileSetAttributeInt64 <-
function(object, attribute, value, flags = "G_FILE_QUERY_INFO_NONE", cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  attribute <- as.character(attribute)
  value <- as.numeric(value)
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_set_attribute_int64", object, attribute, value, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileMountEnclosingVolume <-
function(object, flags, mount.operation, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GFile")
  
  checkPtrType(mount.operation, "GMountOperation")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_file_mount_enclosing_volume", object, flags, mount.operation, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileMountEnclosingVolumeFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_file_mount_enclosing_volume_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileMountMountable <-
function(object, flags, mount.operation, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GFile")
  
  checkPtrType(mount.operation, "GMountOperation")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_file_mount_mountable", object, flags, mount.operation, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileMountMountableFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_file_mount_mountable_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileUnmountMountable <-
function(object, flags = "G_MOUNT_UNMOUNT_NONE", cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GFile")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_file_unmount_mountable", object, flags, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileUnmountMountableFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_file_unmount_mountable_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileEjectMountable <-
function(object, flags = "G_MOUNT_UNMOUNT_NONE", cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GFile")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_file_eject_mountable", object, flags, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileEjectMountableFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_file_eject_mountable_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileCopyAttributes <-
function(object, destination, flags = "G_FILE_COPY_NONE", cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  checkPtrType(destination, "GFile")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_copy_attributes", object, destination, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileMonitorDirectory <-
function(object, flags = "G_FILE_MONITOR_NONE", cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_monitor_directory", object, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileMonitorFile <-
function(object, flags = "G_FILE_MONITOR_NONE", cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_monitor_file", object, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileQueryDefaultHandler <-
function(object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_query_default_handler", object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileLoadContents <-
function(object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_load_contents", object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileLoadContentsAsync <-
function(object, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GFile")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_file_load_contents_async", object, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileLoadContentsFinish <-
function(object, res, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  checkPtrType(res, "GAsyncResult")

  w <- .RGtkCall("S_g_file_load_contents_finish", object, res, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileLoadPartialContentsFinish <-
function(object, res, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  checkPtrType(res, "GAsyncResult")

  w <- .RGtkCall("S_g_file_load_partial_contents_finish", object, res, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileReplaceContents <-
function(object, contents, length, etag, make.backup, flags = "G_FILE_CREATE_NONE", cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  contents <- as.character(contents)
  length <- as.numeric(length)
  etag <- as.character(etag)
  make.backup <- as.logical(make.backup)
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_replace_contents", object, contents, length, etag, make.backup, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileReplaceContentsAsync <-
function(object, contents, length, etag, make.backup, flags = "G_FILE_CREATE_NONE", cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GFile")
  contents <- as.character(contents)
  length <- as.numeric(length)
  etag <- as.character(etag)
  make.backup <- as.logical(make.backup)
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_file_replace_contents_async", object, contents, length, etag, make.backup, flags, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileReplaceContentsFinish <-
function(object, res, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  checkPtrType(res, "GAsyncResult")

  w <- .RGtkCall("S_g_file_replace_contents_finish", object, res, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileIconGetType <-
function()
{
  

  w <- .RGtkCall("S_g_file_icon_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gFileIconNew <-
function(file)
{
  checkPtrType(file, "GFile")

  w <- .RGtkCall("S_g_file_icon_new", file, PACKAGE = "RGtk2")

  return(w)
} 


gFileIconGetFile <-
function(object)
{
  checkPtrType(object, "GFileIcon")

  w <- .RGtkCall("S_g_file_icon_get_file", object, PACKAGE = "RGtk2")

  return(w)
} 


gFileInfoGetType <-
function()
{
  

  w <- .RGtkCall("S_g_file_info_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gFileInfoNew <-
function()
{
  

  w <- .RGtkCall("S_g_file_info_new", PACKAGE = "RGtk2")

  return(w)
} 


gFileInfoDup <-
function(object)
{
  checkPtrType(object, "GFileInfo")

  w <- .RGtkCall("S_g_file_info_dup", object, PACKAGE = "RGtk2")

  return(w)
} 


gFileInfoCopyInto <-
function(object, dest.info)
{
  checkPtrType(object, "GFileInfo")
  checkPtrType(dest.info, "GFileInfo")

  w <- .RGtkCall("S_g_file_info_copy_into", object, dest.info, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileInfoHasAttribute <-
function(object, attribute)
{
  checkPtrType(object, "GFileInfo")
  attribute <- as.character(attribute)

  w <- .RGtkCall("S_g_file_info_has_attribute", object, attribute, PACKAGE = "RGtk2")

  return(w)
} 


gFileInfoListAttributes <-
function(object, name.space)
{
  checkPtrType(object, "GFileInfo")
  name.space <- as.character(name.space)

  w <- .RGtkCall("S_g_file_info_list_attributes", object, name.space, PACKAGE = "RGtk2")

  return(w)
} 


gFileInfoGetAttributeData <-
function(object, attribute)
{
  checkPtrType(object, "GFileInfo")
  attribute <- as.character(attribute)

  w <- .RGtkCall("S_g_file_info_get_attribute_data", object, attribute, PACKAGE = "RGtk2")

  return(w)
} 


gFileInfoGetAttributeType <-
function(object, attribute)
{
  checkPtrType(object, "GFileInfo")
  attribute <- as.character(attribute)

  w <- .RGtkCall("S_g_file_info_get_attribute_type", object, attribute, PACKAGE = "RGtk2")

  return(w)
} 


gFileInfoRemoveAttribute <-
function(object, attribute)
{
  checkPtrType(object, "GFileInfo")
  attribute <- as.character(attribute)

  w <- .RGtkCall("S_g_file_info_remove_attribute", object, attribute, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileInfoGetAttributeStatus <-
function(object, attribute)
{
  checkPtrType(object, "GFileInfo")
  attribute <- as.character(attribute)

  w <- .RGtkCall("S_g_file_info_get_attribute_status", object, attribute, PACKAGE = "RGtk2")

  return(w)
} 


gFileInfoGetAttributeAsString <-
function(object, attribute)
{
  checkPtrType(object, "GFileInfo")
  attribute <- as.character(attribute)

  w <- .RGtkCall("S_g_file_info_get_attribute_as_string", object, attribute, PACKAGE = "RGtk2")

  return(w)
} 


gFileInfoGetAttributeString <-
function(object, attribute)
{
  checkPtrType(object, "GFileInfo")
  attribute <- as.character(attribute)

  w <- .RGtkCall("S_g_file_info_get_attribute_string", object, attribute, PACKAGE = "RGtk2")

  return(w)
} 


gFileInfoGetAttributeByteString <-
function(object, attribute)
{
  checkPtrType(object, "GFileInfo")
  attribute <- as.character(attribute)

  w <- .RGtkCall("S_g_file_info_get_attribute_byte_string", object, attribute, PACKAGE = "RGtk2")

  return(w)
} 


gFileInfoGetAttributeBoolean <-
function(object, attribute)
{
  checkPtrType(object, "GFileInfo")
  attribute <- as.character(attribute)

  w <- .RGtkCall("S_g_file_info_get_attribute_boolean", object, attribute, PACKAGE = "RGtk2")

  return(w)
} 


gFileInfoGetAttributeUint32 <-
function(object, attribute)
{
  checkPtrType(object, "GFileInfo")
  attribute <- as.character(attribute)

  w <- .RGtkCall("S_g_file_info_get_attribute_uint32", object, attribute, PACKAGE = "RGtk2")

  return(w)
} 


gFileInfoGetAttributeInt32 <-
function(object, attribute)
{
  checkPtrType(object, "GFileInfo")
  attribute <- as.character(attribute)

  w <- .RGtkCall("S_g_file_info_get_attribute_int32", object, attribute, PACKAGE = "RGtk2")

  return(w)
} 


gFileInfoGetAttributeUint64 <-
function(object, attribute)
{
  checkPtrType(object, "GFileInfo")
  attribute <- as.character(attribute)

  w <- .RGtkCall("S_g_file_info_get_attribute_uint64", object, attribute, PACKAGE = "RGtk2")

  return(w)
} 


gFileInfoGetAttributeInt64 <-
function(object, attribute)
{
  checkPtrType(object, "GFileInfo")
  attribute <- as.character(attribute)

  w <- .RGtkCall("S_g_file_info_get_attribute_int64", object, attribute, PACKAGE = "RGtk2")

  return(w)
} 


gFileInfoGetAttributeObject <-
function(object, attribute)
{
  checkPtrType(object, "GFileInfo")
  attribute <- as.character(attribute)

  w <- .RGtkCall("S_g_file_info_get_attribute_object", object, attribute, PACKAGE = "RGtk2")

  return(w)
} 


gFileInfoSetAttribute <-
function(object, attribute, type, value.p)
{
  checkPtrType(object, "GFileInfo")
  attribute <- as.character(attribute)
  
  

  w <- .RGtkCall("S_g_file_info_set_attribute", object, attribute, type, value.p, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileInfoSetAttributeString <-
function(object, attribute, attr.value)
{
  checkPtrType(object, "GFileInfo")
  attribute <- as.character(attribute)
  attr.value <- as.character(attr.value)

  w <- .RGtkCall("S_g_file_info_set_attribute_string", object, attribute, attr.value, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileInfoSetAttributeByteString <-
function(object, attribute, attr.value)
{
  checkPtrType(object, "GFileInfo")
  attribute <- as.character(attribute)
  attr.value <- as.character(attr.value)

  w <- .RGtkCall("S_g_file_info_set_attribute_byte_string", object, attribute, attr.value, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileInfoSetAttributeBoolean <-
function(object, attribute, attr.value)
{
  checkPtrType(object, "GFileInfo")
  attribute <- as.character(attribute)
  attr.value <- as.logical(attr.value)

  w <- .RGtkCall("S_g_file_info_set_attribute_boolean", object, attribute, attr.value, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileInfoSetAttributeUint32 <-
function(object, attribute, attr.value)
{
  checkPtrType(object, "GFileInfo")
  attribute <- as.character(attribute)
  attr.value <- as.numeric(attr.value)

  w <- .RGtkCall("S_g_file_info_set_attribute_uint32", object, attribute, attr.value, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileInfoSetAttributeInt32 <-
function(object, attribute, attr.value)
{
  checkPtrType(object, "GFileInfo")
  attribute <- as.character(attribute)
  attr.value <- as.integer(attr.value)

  w <- .RGtkCall("S_g_file_info_set_attribute_int32", object, attribute, attr.value, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileInfoSetAttributeUint64 <-
function(object, attribute, attr.value)
{
  checkPtrType(object, "GFileInfo")
  attribute <- as.character(attribute)
  attr.value <- as.numeric(attr.value)

  w <- .RGtkCall("S_g_file_info_set_attribute_uint64", object, attribute, attr.value, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileInfoSetAttributeInt64 <-
function(object, attribute, attr.value)
{
  checkPtrType(object, "GFileInfo")
  attribute <- as.character(attribute)
  attr.value <- as.numeric(attr.value)

  w <- .RGtkCall("S_g_file_info_set_attribute_int64", object, attribute, attr.value, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileInfoSetAttributeObject <-
function(object, attribute, attr.value)
{
  checkPtrType(object, "GFileInfo")
  attribute <- as.character(attribute)
  checkPtrType(attr.value, "GObject")

  w <- .RGtkCall("S_g_file_info_set_attribute_object", object, attribute, attr.value, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileInfoClearStatus <-
function(object)
{
  checkPtrType(object, "GFileInfo")

  w <- .RGtkCall("S_g_file_info_clear_status", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileInfoGetFileType <-
function(object)
{
  checkPtrType(object, "GFileInfo")

  w <- .RGtkCall("S_g_file_info_get_file_type", object, PACKAGE = "RGtk2")

  return(w)
} 


gFileInfoGetIsHidden <-
function(object)
{
  checkPtrType(object, "GFileInfo")

  w <- .RGtkCall("S_g_file_info_get_is_hidden", object, PACKAGE = "RGtk2")

  return(w)
} 


gFileInfoGetIsBackup <-
function(object)
{
  checkPtrType(object, "GFileInfo")

  w <- .RGtkCall("S_g_file_info_get_is_backup", object, PACKAGE = "RGtk2")

  return(w)
} 


gFileInfoGetIsSymlink <-
function(object)
{
  checkPtrType(object, "GFileInfo")

  w <- .RGtkCall("S_g_file_info_get_is_symlink", object, PACKAGE = "RGtk2")

  return(w)
} 


gFileInfoGetName <-
function(object)
{
  checkPtrType(object, "GFileInfo")

  w <- .RGtkCall("S_g_file_info_get_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gFileInfoGetDisplayName <-
function(object)
{
  checkPtrType(object, "GFileInfo")

  w <- .RGtkCall("S_g_file_info_get_display_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gFileInfoGetEditName <-
function(object)
{
  checkPtrType(object, "GFileInfo")

  w <- .RGtkCall("S_g_file_info_get_edit_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gFileInfoGetIcon <-
function(object)
{
  checkPtrType(object, "GFileInfo")

  w <- .RGtkCall("S_g_file_info_get_icon", object, PACKAGE = "RGtk2")

  return(w)
} 


gFileInfoGetContentType <-
function(object)
{
  checkPtrType(object, "GFileInfo")

  w <- .RGtkCall("S_g_file_info_get_content_type", object, PACKAGE = "RGtk2")

  return(w)
} 


gFileInfoGetSize <-
function(object)
{
  checkPtrType(object, "GFileInfo")

  w <- .RGtkCall("S_g_file_info_get_size", object, PACKAGE = "RGtk2")

  return(w)
} 


gFileInfoGetModificationTime <-
function(object, result)
{
  checkPtrType(object, "GFileInfo")
  result <- as.GTimeVal(result)

  w <- .RGtkCall("S_g_file_info_get_modification_time", object, result, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileInfoGetSymlinkTarget <-
function(object)
{
  checkPtrType(object, "GFileInfo")

  w <- .RGtkCall("S_g_file_info_get_symlink_target", object, PACKAGE = "RGtk2")

  return(w)
} 


gFileInfoGetEtag <-
function(object)
{
  checkPtrType(object, "GFileInfo")

  w <- .RGtkCall("S_g_file_info_get_etag", object, PACKAGE = "RGtk2")

  return(w)
} 


gFileInfoGetSortOrder <-
function(object)
{
  checkPtrType(object, "GFileInfo")

  w <- .RGtkCall("S_g_file_info_get_sort_order", object, PACKAGE = "RGtk2")

  return(w)
} 


gFileInfoSetAttributeMask <-
function(object, mask)
{
  checkPtrType(object, "GFileInfo")
  checkPtrType(mask, "GFileAttributeMatcher")

  w <- .RGtkCall("S_g_file_info_set_attribute_mask", object, mask, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileInfoUnsetAttributeMask <-
function(object)
{
  checkPtrType(object, "GFileInfo")

  w <- .RGtkCall("S_g_file_info_unset_attribute_mask", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileInfoSetFileType <-
function(object, type)
{
  checkPtrType(object, "GFileInfo")
  

  w <- .RGtkCall("S_g_file_info_set_file_type", object, type, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileInfoSetIsHidden <-
function(object, is.hidden)
{
  checkPtrType(object, "GFileInfo")
  is.hidden <- as.logical(is.hidden)

  w <- .RGtkCall("S_g_file_info_set_is_hidden", object, is.hidden, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileInfoSetIsSymlink <-
function(object, is.symlink)
{
  checkPtrType(object, "GFileInfo")
  is.symlink <- as.logical(is.symlink)

  w <- .RGtkCall("S_g_file_info_set_is_symlink", object, is.symlink, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileInfoSetName <-
function(object, name)
{
  checkPtrType(object, "GFileInfo")
  name <- as.character(name)

  w <- .RGtkCall("S_g_file_info_set_name", object, name, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileInfoSetDisplayName <-
function(object, display.name)
{
  checkPtrType(object, "GFileInfo")
  display.name <- as.character(display.name)

  w <- .RGtkCall("S_g_file_info_set_display_name", object, display.name, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileInfoSetEditName <-
function(object, edit.name)
{
  checkPtrType(object, "GFileInfo")
  edit.name <- as.character(edit.name)

  w <- .RGtkCall("S_g_file_info_set_edit_name", object, edit.name, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileInfoSetIcon <-
function(object, icon)
{
  checkPtrType(object, "GFileInfo")
  checkPtrType(icon, "GIcon")

  w <- .RGtkCall("S_g_file_info_set_icon", object, icon, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileInfoSetContentType <-
function(object, content.type)
{
  checkPtrType(object, "GFileInfo")
  content.type <- as.character(content.type)

  w <- .RGtkCall("S_g_file_info_set_content_type", object, content.type, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileInfoSetSize <-
function(object, size)
{
  checkPtrType(object, "GFileInfo")
  size <- as.numeric(size)

  w <- .RGtkCall("S_g_file_info_set_size", object, size, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileInfoSetModificationTime <-
function(object, mtime)
{
  checkPtrType(object, "GFileInfo")
  mtime <- as.GTimeVal(mtime)

  w <- .RGtkCall("S_g_file_info_set_modification_time", object, mtime, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileInfoSetSymlinkTarget <-
function(object, symlink.target)
{
  checkPtrType(object, "GFileInfo")
  symlink.target <- as.character(symlink.target)

  w <- .RGtkCall("S_g_file_info_set_symlink_target", object, symlink.target, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileInfoSetSortOrder <-
function(object, sort.order)
{
  checkPtrType(object, "GFileInfo")
  sort.order <- as.integer(sort.order)

  w <- .RGtkCall("S_g_file_info_set_sort_order", object, sort.order, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileAttributeMatcherNew <-
function(attributes)
{
  attributes <- as.character(attributes)

  w <- .RGtkCall("S_g_file_attribute_matcher_new", attributes, PACKAGE = "RGtk2")

  return(w)
} 


gFileAttributeMatcherMatches <-
function(object, attribute)
{
  checkPtrType(object, "GFileAttributeMatcher")
  attribute <- as.character(attribute)

  w <- .RGtkCall("S_g_file_attribute_matcher_matches", object, attribute, PACKAGE = "RGtk2")

  return(w)
} 


gFileAttributeMatcherMatchesOnly <-
function(object, attribute)
{
  checkPtrType(object, "GFileAttributeMatcher")
  attribute <- as.character(attribute)

  w <- .RGtkCall("S_g_file_attribute_matcher_matches_only", object, attribute, PACKAGE = "RGtk2")

  return(w)
} 


gFileAttributeMatcherEnumerateNamespace <-
function(object, ns)
{
  checkPtrType(object, "GFileAttributeMatcher")
  ns <- as.character(ns)

  w <- .RGtkCall("S_g_file_attribute_matcher_enumerate_namespace", object, ns, PACKAGE = "RGtk2")

  return(w)
} 


gFileAttributeMatcherEnumerateNext <-
function(object)
{
  checkPtrType(object, "GFileAttributeMatcher")

  w <- .RGtkCall("S_g_file_attribute_matcher_enumerate_next", object, PACKAGE = "RGtk2")

  return(w)
} 


gFileInputStreamGetType <-
function()
{
  

  w <- .RGtkCall("S_g_file_input_stream_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gFileInputStreamQueryInfo <-
function(object, attributes, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFileInputStream")
  attributes <- as.character(attributes)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_input_stream_query_info", object, attributes, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileInputStreamQueryInfoAsync <-
function(object, attributes, io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GFileInputStream")
  attributes <- as.character(attributes)
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_file_input_stream_query_info_async", object, attributes, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileInputStreamQueryInfoFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GFileInputStream")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_file_input_stream_query_info_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileMonitorGetType <-
function()
{
  

  w <- .RGtkCall("S_g_file_monitor_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gFileMonitorCancel <-
function(object)
{
  checkPtrType(object, "GFileMonitor")

  w <- .RGtkCall("S_g_file_monitor_cancel", object, PACKAGE = "RGtk2")

  return(w)
} 


gFileMonitorIsCancelled <-
function(object)
{
  checkPtrType(object, "GFileMonitor")

  w <- .RGtkCall("S_g_file_monitor_is_cancelled", object, PACKAGE = "RGtk2")

  return(w)
} 


gFileMonitorSetRateLimit <-
function(object, limit.msecs)
{
  checkPtrType(object, "GFileMonitor")
  limit.msecs <- as.integer(limit.msecs)

  w <- .RGtkCall("S_g_file_monitor_set_rate_limit", object, limit.msecs, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileMonitorEmitEvent <-
function(object, file, other.file, event.type)
{
  checkPtrType(object, "GFileMonitor")
  checkPtrType(file, "GFile")
  checkPtrType(other.file, "GFile")
  

  w <- .RGtkCall("S_g_file_monitor_emit_event", object, file, other.file, event.type, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFilenameCompleterGetType <-
function()
{
  

  w <- .RGtkCall("S_g_filename_completer_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gFilenameCompleterNew <-
function()
{
  

  w <- .RGtkCall("S_g_filename_completer_new", PACKAGE = "RGtk2")

  return(w)
} 


gFilenameCompleterGetCompletionSuffix <-
function(object, initial.text)
{
  checkPtrType(object, "GFilenameCompleter")
  initial.text <- as.character(initial.text)

  w <- .RGtkCall("S_g_filename_completer_get_completion_suffix", object, initial.text, PACKAGE = "RGtk2")

  return(w)
} 


gFilenameCompleterGetCompletions <-
function(object, initial.text)
{
  checkPtrType(object, "GFilenameCompleter")
  initial.text <- as.character(initial.text)

  w <- .RGtkCall("S_g_filename_completer_get_completions", object, initial.text, PACKAGE = "RGtk2")

  return(w)
} 


gFilenameCompleterSetDirsOnly <-
function(object, dirs.only)
{
  checkPtrType(object, "GFilenameCompleter")
  dirs.only <- as.logical(dirs.only)

  w <- .RGtkCall("S_g_filename_completer_set_dirs_only", object, dirs.only, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileOutputStreamGetType <-
function()
{
  

  w <- .RGtkCall("S_g_file_output_stream_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gFileOutputStreamQueryInfo <-
function(object, attributes, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFileOutputStream")
  attributes <- as.character(attributes)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_output_stream_query_info", object, attributes, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileOutputStreamQueryInfoAsync <-
function(object, attributes, io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GFileOutputStream")
  attributes <- as.character(attributes)
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_file_output_stream_query_info_async", object, attributes, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileOutputStreamQueryInfoFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GFileOutputStream")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_file_output_stream_query_info_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileOutputStreamGetEtag <-
function(object)
{
  checkPtrType(object, "GFileOutputStream")

  w <- .RGtkCall("S_g_file_output_stream_get_etag", object, PACKAGE = "RGtk2")

  return(w)
} 


gFilterInputStreamGetType <-
function()
{
  

  w <- .RGtkCall("S_g_filter_input_stream_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gFilterInputStreamGetBaseStream <-
function(object)
{
  checkPtrType(object, "GFilterInputStream")

  w <- .RGtkCall("S_g_filter_input_stream_get_base_stream", object, PACKAGE = "RGtk2")

  return(w)
} 


gFilterOutputStreamGetType <-
function()
{
  

  w <- .RGtkCall("S_g_filter_output_stream_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gFilterOutputStreamGetBaseStream <-
function(object)
{
  checkPtrType(object, "GFilterOutputStream")

  w <- .RGtkCall("S_g_filter_output_stream_get_base_stream", object, PACKAGE = "RGtk2")

  return(w)
} 


gIconGetType <-
function()
{
  

  w <- .RGtkCall("S_g_icon_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gIconHash <-
function(icon)
{
  

  w <- .RGtkCall("S_g_icon_hash", icon, PACKAGE = "RGtk2")

  return(w)
} 


gIconEqual <-
function(object, icon2)
{
  checkPtrType(object, "GIcon")
  checkPtrType(icon2, "GIcon")

  w <- .RGtkCall("S_g_icon_equal", object, icon2, PACKAGE = "RGtk2")

  return(w)
} 


gInputStreamGetType <-
function()
{
  

  w <- .RGtkCall("S_g_input_stream_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gInputStreamRead <-
function(object, count, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GInputStream")
  count <- as.numeric(count)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_input_stream_read", object, count, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gInputStreamReadAll <-
function(object, count, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GInputStream")
  count <- as.numeric(count)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_input_stream_read_all", object, count, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gInputStreamSkip <-
function(object, count, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GInputStream")
  count <- as.numeric(count)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_input_stream_skip", object, count, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gInputStreamClose <-
function(object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GInputStream")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_input_stream_close", object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gInputStreamReadAsync <-
function(object, count, io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GInputStream")
  count <- as.numeric(count)
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_input_stream_read_async", object, count, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(w)
} 


gInputStreamReadFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GInputStream")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_input_stream_read_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gInputStreamSkipAsync <-
function(object, count, io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GInputStream")
  count <- as.numeric(count)
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_input_stream_skip_async", object, count, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gInputStreamSkipFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GInputStream")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_input_stream_skip_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gInputStreamCloseAsync <-
function(object, io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GInputStream")
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_input_stream_close_async", object, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gInputStreamCloseFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GInputStream")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_input_stream_close_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gInputStreamIsClosed <-
function(object)
{
  checkPtrType(object, "GInputStream")

  w <- .RGtkCall("S_g_input_stream_is_closed", object, PACKAGE = "RGtk2")

  return(w)
} 


gInputStreamHasPending <-
function(object)
{
  checkPtrType(object, "GInputStream")

  w <- .RGtkCall("S_g_input_stream_has_pending", object, PACKAGE = "RGtk2")

  return(w)
} 


gInputStreamSetPending <-
function(object, .errwarn = TRUE)
{
  checkPtrType(object, "GInputStream")

  w <- .RGtkCall("S_g_input_stream_set_pending", object, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gInputStreamClearPending <-
function(object)
{
  checkPtrType(object, "GInputStream")

  w <- .RGtkCall("S_g_input_stream_clear_pending", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gAppInfoCreateFlagsGetType <-
function()
{
  

  w <- .RGtkCall("S_g_app_info_create_flags_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gDataStreamByteOrderGetType <-
function()
{
  

  w <- .RGtkCall("S_g_data_stream_byte_order_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gDataStreamNewlineTypeGetType <-
function()
{
  

  w <- .RGtkCall("S_g_data_stream_newline_type_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gFileQueryInfoFlagsGetType <-
function()
{
  

  w <- .RGtkCall("S_g_file_query_info_flags_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gFileCreateFlagsGetType <-
function()
{
  

  w <- .RGtkCall("S_g_file_create_flags_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gFileCopyFlagsGetType <-
function()
{
  

  w <- .RGtkCall("S_g_file_copy_flags_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gFileMonitorFlagsGetType <-
function()
{
  

  w <- .RGtkCall("S_g_file_monitor_flags_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gFileAttributeTypeGetType <-
function()
{
  

  w <- .RGtkCall("S_g_file_attribute_type_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gFileAttributeInfoFlagsGetType <-
function()
{
  

  w <- .RGtkCall("S_g_file_attribute_info_flags_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gFileAttributeStatusGetType <-
function()
{
  

  w <- .RGtkCall("S_g_file_attribute_status_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gFileTypeGetType <-
function()
{
  

  w <- .RGtkCall("S_g_file_type_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gFileMonitorEventGetType <-
function()
{
  

  w <- .RGtkCall("S_g_file_monitor_event_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gIoErrorEnumGetType <-
function()
{
  

  w <- .RGtkCall("S_g_io_error_enum_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gAskPasswordFlagsGetType <-
function()
{
  

  w <- .RGtkCall("S_g_ask_password_flags_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gPasswordSaveGetType <-
function()
{
  

  w <- .RGtkCall("S_g_password_save_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gOutputStreamSpliceFlagsGetType <-
function()
{
  

  w <- .RGtkCall("S_g_output_stream_splice_flags_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gIoErrorQuark <-
function()
{
  

  w <- .RGtkCall("S_g_io_error_quark", PACKAGE = "RGtk2")

  return(w)
} 


gIoErrorFromErrno <-
function(err.no)
{
  err.no <- as.integer(err.no)

  w <- .RGtkCall("S_g_io_error_from_errno", err.no, PACKAGE = "RGtk2")

  return(w)
} 


gIOModuleGetType <-
function()
{
  

  w <- .RGtkCall("S_g_io_module_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gIOModuleNew <-
function(filename)
{
  filename <- as.character(filename)

  w <- .RGtkCall("S_g_io_module_new", filename, PACKAGE = "RGtk2")

  return(w)
} 


gIOModulesLoadAllInDirectory <-
function(dirname)
{
  dirname <- as.character(dirname)

  w <- .RGtkCall("S_g_io_modules_load_all_in_directory", dirname, PACKAGE = "RGtk2")

  return(w)
} 


gIoSchedulerCancelAllJobs <-
function()
{
  

  w <- .RGtkCall("S_g_io_scheduler_cancel_all_jobs", PACKAGE = "RGtk2")

  return(w)
} 


gIoSchedulerJobSendToMainloop <-
function(object, func, user.data = NULL)
{
  checkPtrType(object, "GIOSchedulerJob")
  func <- as.function(func)
  

  w <- .RGtkCall("S_g_io_scheduler_job_send_to_mainloop", object, func, user.data, PACKAGE = "RGtk2")

  return(w)
} 


gIoSchedulerJobSendToMainloopAsync <-
function(object, func, user.data = NULL)
{
  checkPtrType(object, "GIOSchedulerJob")
  func <- as.function(func)
  

  w <- .RGtkCall("S_g_io_scheduler_job_send_to_mainloop_async", object, func, user.data, PACKAGE = "RGtk2")

  return(w)
} 


gLoadableIconGetType <-
function()
{
  

  w <- .RGtkCall("S_g_loadable_icon_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gLoadableIconLoad <-
function(object, size, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GLoadableIcon")
  size <- as.integer(size)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_loadable_icon_load", object, size, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gLoadableIconLoadAsync <-
function(object, size, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GLoadableIcon")
  size <- as.integer(size)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_loadable_icon_load_async", object, size, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gLoadableIconLoadFinish <-
function(object, res, type, .errwarn = TRUE)
{
  checkPtrType(object, "GLoadableIcon")
  checkPtrType(res, "GAsyncResult")
  type <- as.list(as.character(type))

  w <- .RGtkCall("S_g_loadable_icon_load_finish", object, res, type, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gMemoryInputStreamGetType <-
function()
{
  

  w <- .RGtkCall("S_g_memory_input_stream_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gMemoryInputStreamNew <-
function()
{
  

  w <- .RGtkCall("S_g_memory_input_stream_new", PACKAGE = "RGtk2")

  return(w)
} 


gMemoryInputStreamNewFromData <-
function(data)
{
  data <- as.list(as.raw(data))

  w <- .RGtkCall("S_g_memory_input_stream_new_from_data", data, PACKAGE = "RGtk2")

  return(w)
} 


gMemoryInputStreamAddData <-
function(object, data)
{
  checkPtrType(object, "GMemoryInputStream")
  data <- as.list(as.raw(data))

  w <- .RGtkCall("S_g_memory_input_stream_add_data", object, data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gMemoryOutputStreamGetType <-
function()
{
  

  w <- .RGtkCall("S_g_memory_output_stream_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gMemoryOutputStreamNew <-
function(len)
{
  len <- as.numeric(len)

  w <- .RGtkCall("S_g_memory_output_stream_new", len, PACKAGE = "RGtk2")

  return(w)
} 


gMemoryOutputStreamGetData <-
function(object)
{
  checkPtrType(object, "GMemoryOutputStream")

  w <- .RGtkCall("S_g_memory_output_stream_get_data", object, PACKAGE = "RGtk2")

  return(w)
} 


gMemoryOutputStreamGetSize <-
function(object)
{
  checkPtrType(object, "GMemoryOutputStream")

  w <- .RGtkCall("S_g_memory_output_stream_get_size", object, PACKAGE = "RGtk2")

  return(w)
} 


gMountGetType <-
function()
{
  

  w <- .RGtkCall("S_g_mount_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gMountGetRoot <-
function(object)
{
  checkPtrType(object, "GMount")

  w <- .RGtkCall("S_g_mount_get_root", object, PACKAGE = "RGtk2")

  return(w)
} 


gMountGetName <-
function(object)
{
  checkPtrType(object, "GMount")

  w <- .RGtkCall("S_g_mount_get_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gMountGetIcon <-
function(object)
{
  checkPtrType(object, "GMount")

  w <- .RGtkCall("S_g_mount_get_icon", object, PACKAGE = "RGtk2")

  return(w)
} 


gMountGetUuid <-
function(object)
{
  checkPtrType(object, "GMount")

  w <- .RGtkCall("S_g_mount_get_uuid", object, PACKAGE = "RGtk2")

  return(w)
} 


gMountGetVolume <-
function(object)
{
  checkPtrType(object, "GMount")

  w <- .RGtkCall("S_g_mount_get_volume", object, PACKAGE = "RGtk2")

  return(w)
} 


gMountGetDrive <-
function(object)
{
  checkPtrType(object, "GMount")

  w <- .RGtkCall("S_g_mount_get_drive", object, PACKAGE = "RGtk2")

  return(w)
} 


gMountCanUnmount <-
function(object)
{
  checkPtrType(object, "GMount")

  w <- .RGtkCall("S_g_mount_can_unmount", object, PACKAGE = "RGtk2")

  return(w)
} 


gMountCanEject <-
function(object)
{
  checkPtrType(object, "GMount")

  w <- .RGtkCall("S_g_mount_can_eject", object, PACKAGE = "RGtk2")

  return(w)
} 


gMountUnmount <-
function(object, flags = "G_MOUNT_UNMOUNT_NONE", cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GMount")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_mount_unmount", object, flags, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gMountUnmountFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GMount")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_mount_unmount_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gMountEject <-
function(object, flags = "G_MOUNT_UNMOUNT_NONE", cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GMount")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_mount_eject", object, flags, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gMountEjectFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GMount")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_mount_eject_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gMountRemount <-
function(object, flags, mount.operation, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GMount")
  
  checkPtrType(mount.operation, "GMountOperation")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_mount_remount", object, flags, mount.operation, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gMountRemountFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GMount")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_mount_remount_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gMountOperationGetType <-
function()
{
  

  w <- .RGtkCall("S_g_mount_operation_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gMountOperationNew <-
function()
{
  

  w <- .RGtkCall("S_g_mount_operation_new", PACKAGE = "RGtk2")

  return(w)
} 


gMountOperationGetUsername <-
function(object)
{
  checkPtrType(object, "GMountOperation")

  w <- .RGtkCall("S_g_mount_operation_get_username", object, PACKAGE = "RGtk2")

  return(w)
} 


gMountOperationSetUsername <-
function(object, username)
{
  checkPtrType(object, "GMountOperation")
  username <- as.character(username)

  w <- .RGtkCall("S_g_mount_operation_set_username", object, username, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gMountOperationGetPassword <-
function(object)
{
  checkPtrType(object, "GMountOperation")

  w <- .RGtkCall("S_g_mount_operation_get_password", object, PACKAGE = "RGtk2")

  return(w)
} 


gMountOperationSetPassword <-
function(object, password)
{
  checkPtrType(object, "GMountOperation")
  password <- as.character(password)

  w <- .RGtkCall("S_g_mount_operation_set_password", object, password, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gMountOperationGetAnonymous <-
function(object)
{
  checkPtrType(object, "GMountOperation")

  w <- .RGtkCall("S_g_mount_operation_get_anonymous", object, PACKAGE = "RGtk2")

  return(w)
} 


gMountOperationSetAnonymous <-
function(object, anonymous)
{
  checkPtrType(object, "GMountOperation")
  anonymous <- as.logical(anonymous)

  w <- .RGtkCall("S_g_mount_operation_set_anonymous", object, anonymous, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gMountOperationGetDomain <-
function(object)
{
  checkPtrType(object, "GMountOperation")

  w <- .RGtkCall("S_g_mount_operation_get_domain", object, PACKAGE = "RGtk2")

  return(w)
} 


gMountOperationSetDomain <-
function(object, domain)
{
  checkPtrType(object, "GMountOperation")
  domain <- as.character(domain)

  w <- .RGtkCall("S_g_mount_operation_set_domain", object, domain, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gMountOperationGetPasswordSave <-
function(object)
{
  checkPtrType(object, "GMountOperation")

  w <- .RGtkCall("S_g_mount_operation_get_password_save", object, PACKAGE = "RGtk2")

  return(w)
} 


gMountOperationSetPasswordSave <-
function(object, save)
{
  checkPtrType(object, "GMountOperation")
  

  w <- .RGtkCall("S_g_mount_operation_set_password_save", object, save, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gMountOperationGetChoice <-
function(object)
{
  checkPtrType(object, "GMountOperation")

  w <- .RGtkCall("S_g_mount_operation_get_choice", object, PACKAGE = "RGtk2")

  return(w)
} 


gMountOperationSetChoice <-
function(object, choice)
{
  checkPtrType(object, "GMountOperation")
  choice <- as.integer(choice)

  w <- .RGtkCall("S_g_mount_operation_set_choice", object, choice, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gMountOperationReply <-
function(object, result)
{
  checkPtrType(object, "GMountOperation")
  

  w <- .RGtkCall("S_g_mount_operation_reply", object, result, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gNativeVolumeMonitorGetType <-
function()
{
  

  w <- .RGtkCall("S_g_native_volume_monitor_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gOutputStreamGetType <-
function()
{
  

  w <- .RGtkCall("S_g_output_stream_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gOutputStreamWrite <-
function(object, buffer, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GOutputStream")
  buffer <- as.list(as.raw(buffer))
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_output_stream_write", object, buffer, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gOutputStreamWriteAll <-
function(object, buffer, bytes.written, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GOutputStream")
  buffer <- as.list(as.raw(buffer))
  bytes.written <- as.list(as.numeric(bytes.written))
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_output_stream_write_all", object, buffer, bytes.written, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gOutputStreamSplice <-
function(object, source, flags = "G_OUTPUT_STREAM_SPLICE_NONE", cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GOutputStream")
  checkPtrType(source, "GInputStream")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_output_stream_splice", object, source, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gOutputStreamFlush <-
function(object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GOutputStream")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_output_stream_flush", object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gOutputStreamClose <-
function(object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GOutputStream")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_output_stream_close", object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gOutputStreamWriteAsync <-
function(object, buffer, io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GOutputStream")
  buffer <- as.list(as.raw(buffer))
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_output_stream_write_async", object, buffer, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(w)
} 


gOutputStreamWriteFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GOutputStream")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_output_stream_write_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gOutputStreamSpliceAsync <-
function(object, source, flags = "G_OUTPUT_STREAM_SPLICE_NONE", io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GOutputStream")
  checkPtrType(source, "GInputStream")
  
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_output_stream_splice_async", object, source, flags, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gOutputStreamSpliceFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GOutputStream")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_output_stream_splice_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gOutputStreamFlushAsync <-
function(object, io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GOutputStream")
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_output_stream_flush_async", object, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gOutputStreamFlushFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GOutputStream")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_output_stream_flush_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gOutputStreamCloseAsync <-
function(object, io.priority = 0, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GOutputStream")
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_output_stream_close_async", object, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gOutputStreamCloseFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GOutputStream")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_output_stream_close_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gOutputStreamIsClosed <-
function(object)
{
  checkPtrType(object, "GOutputStream")

  w <- .RGtkCall("S_g_output_stream_is_closed", object, PACKAGE = "RGtk2")

  return(w)
} 


gOutputStreamHasPending <-
function(object)
{
  checkPtrType(object, "GOutputStream")

  w <- .RGtkCall("S_g_output_stream_has_pending", object, PACKAGE = "RGtk2")

  return(w)
} 


gOutputStreamSetPending <-
function(object, .errwarn = TRUE)
{
  checkPtrType(object, "GOutputStream")

  w <- .RGtkCall("S_g_output_stream_set_pending", object, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gOutputStreamClearPending <-
function(object)
{
  checkPtrType(object, "GOutputStream")

  w <- .RGtkCall("S_g_output_stream_clear_pending", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gSeekableGetType <-
function()
{
  

  w <- .RGtkCall("S_g_seekable_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gSeekableTell <-
function(object)
{
  checkPtrType(object, "GSeekable")

  w <- .RGtkCall("S_g_seekable_tell", object, PACKAGE = "RGtk2")

  return(w)
} 


gSeekableCanSeek <-
function(object)
{
  checkPtrType(object, "GSeekable")

  w <- .RGtkCall("S_g_seekable_can_seek", object, PACKAGE = "RGtk2")

  return(w)
} 


gSeekableSeek <-
function(object, offset, type = "G_SEEK_SET", cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GSeekable")
  offset <- as.numeric(offset)
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_seekable_seek", object, offset, type, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSeekableCanTruncate <-
function(object)
{
  checkPtrType(object, "GSeekable")

  w <- .RGtkCall("S_g_seekable_can_truncate", object, PACKAGE = "RGtk2")

  return(w)
} 


gSeekableTruncate <-
function(object, offset, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GSeekable")
  offset <- as.numeric(offset)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_seekable_truncate", object, offset, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSimpleAsyncResultGetType <-
function()
{
  

  w <- .RGtkCall("S_g_simple_async_result_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gSimpleAsyncResultNew <-
function(source.object, callback, user.data = NULL, source.tag)
{
  checkPtrType(source.object, "GObject")
  callback <- as.function(callback)
  
  

  w <- .RGtkCall("S_g_simple_async_result_new", source.object, callback, user.data, source.tag, PACKAGE = "RGtk2")

  return(w)
} 


gSimpleAsyncResultNewFromError <-
function(source.object, callback, user.data = NULL)
{
  checkPtrType(source.object, "GObject")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_simple_async_result_new_from_error", source.object, callback, user.data, PACKAGE = "RGtk2")

  return(w)
} 


gSimpleAsyncResultSetOpResGpointer <-
function(object, op.res)
{
  checkPtrType(object, "GSimpleAsyncResult")
  

  w <- .RGtkCall("S_g_simple_async_result_set_op_res_gpointer", object, op.res, PACKAGE = "RGtk2")

  return(w)
} 


gSimpleAsyncResultGetOpResGpointer <-
function(object)
{
  checkPtrType(object, "GSimpleAsyncResult")

  w <- .RGtkCall("S_g_simple_async_result_get_op_res_gpointer", object, PACKAGE = "RGtk2")

  return(w)
} 


gSimpleAsyncResultSetOpResGssize <-
function(object, op.res)
{
  checkPtrType(object, "GSimpleAsyncResult")
  op.res <- as.integer(op.res)

  w <- .RGtkCall("S_g_simple_async_result_set_op_res_gssize", object, op.res, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gSimpleAsyncResultGetOpResGssize <-
function(object)
{
  checkPtrType(object, "GSimpleAsyncResult")

  w <- .RGtkCall("S_g_simple_async_result_get_op_res_gssize", object, PACKAGE = "RGtk2")

  return(w)
} 


gSimpleAsyncResultSetOpResGboolean <-
function(object, op.res)
{
  checkPtrType(object, "GSimpleAsyncResult")
  op.res <- as.logical(op.res)

  w <- .RGtkCall("S_g_simple_async_result_set_op_res_gboolean", object, op.res, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gSimpleAsyncResultGetOpResGboolean <-
function(object)
{
  checkPtrType(object, "GSimpleAsyncResult")

  w <- .RGtkCall("S_g_simple_async_result_get_op_res_gboolean", object, PACKAGE = "RGtk2")

  return(w)
} 


gSimpleAsyncResultGetSourceTag <-
function(object)
{
  checkPtrType(object, "GSimpleAsyncResult")

  w <- .RGtkCall("S_g_simple_async_result_get_source_tag", object, PACKAGE = "RGtk2")

  return(w)
} 


gSimpleAsyncResultSetHandleCancellation <-
function(object, handle.cancellation)
{
  checkPtrType(object, "GSimpleAsyncResult")
  handle.cancellation <- as.logical(handle.cancellation)

  w <- .RGtkCall("S_g_simple_async_result_set_handle_cancellation", object, handle.cancellation, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gSimpleAsyncResultComplete <-
function(object)
{
  checkPtrType(object, "GSimpleAsyncResult")

  w <- .RGtkCall("S_g_simple_async_result_complete", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gSimpleAsyncResultCompleteInIdle <-
function(object)
{
  checkPtrType(object, "GSimpleAsyncResult")

  w <- .RGtkCall("S_g_simple_async_result_complete_in_idle", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gSimpleAsyncResultSetFromError <-
function(object)
{
  checkPtrType(object, "GSimpleAsyncResult")

  w <- .RGtkCall("S_g_simple_async_result_set_from_error", object, PACKAGE = "RGtk2")

  return(w)
} 


gSimpleAsyncResultPropagateError <-
function(object, .errwarn = TRUE)
{
  checkPtrType(object, "GSimpleAsyncResult")

  w <- .RGtkCall("S_g_simple_async_result_propagate_error", object, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSimpleAsyncReportGerrorInIdle <-
function(object, callback, user.data = NULL)
{
  checkPtrType(object, "GObject")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_simple_async_report_gerror_in_idle", object, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gThemedIconGetType <-
function()
{
  

  w <- .RGtkCall("S_g_themed_icon_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gThemedIconNew <-
function(iconname = NULL)
{
  

  w <- .RGtkCall("S_g_themed_icon_new", iconname, PACKAGE = "RGtk2")

  return(w)
} 


gThemedIconNewWithDefaultFallbacks <-
function(iconname)
{
  iconname <- as.character(iconname)

  w <- .RGtkCall("S_g_themed_icon_new_with_default_fallbacks", iconname, PACKAGE = "RGtk2")

  return(w)
} 


gThemedIconNewFromNames <-
function(iconnames, len)
{
  iconnames <- as.list(as.character(iconnames))
  len <- as.integer(len)

  w <- .RGtkCall("S_g_themed_icon_new_from_names", iconnames, len, PACKAGE = "RGtk2")

  return(w)
} 


gThemedIconGetNames <-
function(object)
{
  checkPtrType(object, "GThemedIcon")

  w <- .RGtkCall("S_g_themed_icon_get_names", object, PACKAGE = "RGtk2")

  return(w)
} 


gThemedIconAppendName <-
function(object, iconname)
{
  checkPtrType(object, "GThemedIcon")
  iconname <- as.character(iconname)

  w <- .RGtkCall("S_g_themed_icon_append_name", object, iconname, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gVfsGetType <-
function()
{
  

  w <- .RGtkCall("S_g_vfs_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gVfsIsActive <-
function(object)
{
  checkPtrType(object, "GVfs")

  w <- .RGtkCall("S_g_vfs_is_active", object, PACKAGE = "RGtk2")

  return(w)
} 


gVfsGetFileForPath <-
function(object, path)
{
  checkPtrType(object, "GVfs")
  path <- as.character(path)

  w <- .RGtkCall("S_g_vfs_get_file_for_path", object, path, PACKAGE = "RGtk2")

  return(w)
} 


gVfsGetFileForUri <-
function(object, uri)
{
  checkPtrType(object, "GVfs")
  uri <- as.character(uri)

  w <- .RGtkCall("S_g_vfs_get_file_for_uri", object, uri, PACKAGE = "RGtk2")

  return(w)
} 


gVfsParseName <-
function(object, parse.name)
{
  checkPtrType(object, "GVfs")
  parse.name <- as.character(parse.name)

  w <- .RGtkCall("S_g_vfs_parse_name", object, parse.name, PACKAGE = "RGtk2")

  return(w)
} 


gVfsGetDefault <-
function()
{
  

  w <- .RGtkCall("S_g_vfs_get_default", PACKAGE = "RGtk2")

  return(w)
} 


gVfsGetLocal <-
function()
{
  

  w <- .RGtkCall("S_g_vfs_get_local", PACKAGE = "RGtk2")

  return(w)
} 


gVfsGetSupportedUriSchemes <-
function(object)
{
  checkPtrType(object, "GVfs")

  w <- .RGtkCall("S_g_vfs_get_supported_uri_schemes", object, PACKAGE = "RGtk2")

  return(w)
} 


gVolumeGetType <-
function()
{
  

  w <- .RGtkCall("S_g_volume_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gVolumeGetName <-
function(object)
{
  checkPtrType(object, "GVolume")

  w <- .RGtkCall("S_g_volume_get_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gVolumeGetIcon <-
function(object)
{
  checkPtrType(object, "GVolume")

  w <- .RGtkCall("S_g_volume_get_icon", object, PACKAGE = "RGtk2")

  return(w)
} 


gVolumeGetUuid <-
function(object)
{
  checkPtrType(object, "GVolume")

  w <- .RGtkCall("S_g_volume_get_uuid", object, PACKAGE = "RGtk2")

  return(w)
} 


gVolumeGetDrive <-
function(object)
{
  checkPtrType(object, "GVolume")

  w <- .RGtkCall("S_g_volume_get_drive", object, PACKAGE = "RGtk2")

  return(w)
} 


gVolumeGetMount <-
function(object)
{
  checkPtrType(object, "GVolume")

  w <- .RGtkCall("S_g_volume_get_mount", object, PACKAGE = "RGtk2")

  return(w)
} 


gVolumeCanMount <-
function(object)
{
  checkPtrType(object, "GVolume")

  w <- .RGtkCall("S_g_volume_can_mount", object, PACKAGE = "RGtk2")

  return(w)
} 


gVolumeCanEject <-
function(object)
{
  checkPtrType(object, "GVolume")

  w <- .RGtkCall("S_g_volume_can_eject", object, PACKAGE = "RGtk2")

  return(w)
} 


gVolumeShouldAutomount <-
function(object)
{
  checkPtrType(object, "GVolume")

  w <- .RGtkCall("S_g_volume_should_automount", object, PACKAGE = "RGtk2")

  return(w)
} 


gVolumeMount <-
function(object, flags, mount.operation, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GVolume")
  
  checkPtrType(mount.operation, "GMountOperation")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_volume_mount", object, flags, mount.operation, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gVolumeMountFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GVolume")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_volume_mount_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gVolumeEject <-
function(object, flags = "G_MOUNT_UNMOUNT_NONE", cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GVolume")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_volume_eject", object, flags, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gVolumeEjectFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GVolume")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_volume_eject_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gVolumeMonitorGetType <-
function()
{
  

  w <- .RGtkCall("S_g_volume_monitor_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gVolumeMonitorGet <-
function()
{
  

  w <- .RGtkCall("S_g_volume_monitor_get", PACKAGE = "RGtk2")

  return(w)
} 


gVolumeMonitorGetConnectedDrives <-
function(object)
{
  checkPtrType(object, "GVolumeMonitor")

  w <- .RGtkCall("S_g_volume_monitor_get_connected_drives", object, PACKAGE = "RGtk2")

  return(w)
} 


gVolumeMonitorGetVolumes <-
function(object)
{
  checkPtrType(object, "GVolumeMonitor")

  w <- .RGtkCall("S_g_volume_monitor_get_volumes", object, PACKAGE = "RGtk2")

  return(w)
} 


gVolumeMonitorGetMounts <-
function(object)
{
  checkPtrType(object, "GVolumeMonitor")

  w <- .RGtkCall("S_g_volume_monitor_get_mounts", object, PACKAGE = "RGtk2")

  return(w)
} 


gVolumeMonitorGetVolumeForUuid <-
function(object, uuid)
{
  checkPtrType(object, "GVolumeMonitor")
  uuid <- as.character(uuid)

  w <- .RGtkCall("S_g_volume_monitor_get_volume_for_uuid", object, uuid, PACKAGE = "RGtk2")

  return(w)
} 


gVolumeMonitorGetMountForUuid <-
function(object, uuid)
{
  checkPtrType(object, "GVolumeMonitor")
  uuid <- as.character(uuid)

  w <- .RGtkCall("S_g_volume_monitor_get_mount_for_uuid", object, uuid, PACKAGE = "RGtk2")

  return(w)
} 


gVolumeMonitorAdoptOrphanMount <-
function(mount)
{
  checkPtrType(mount, "GMount")

  w <- .RGtkCall("S_g_volume_monitor_adopt_orphan_mount", mount, PACKAGE = "RGtk2")

  return(w)
} 


gIoExtensionPointRegister <-
function(name)
{
  name <- as.character(name)

  w <- .RGtkCall("S_g_io_extension_point_register", name, PACKAGE = "RGtk2")

  return(w)
} 


gIoExtensionPointLookup <-
function(name)
{
  name <- as.character(name)

  w <- .RGtkCall("S_g_io_extension_point_lookup", name, PACKAGE = "RGtk2")

  return(w)
} 


gIoExtensionPointSetRequiredType <-
function(object, type)
{
  checkPtrType(object, "GIOExtensionPoint")
  type <- as.GType(type)

  w <- .RGtkCall("S_g_io_extension_point_set_required_type", object, type, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gIoExtensionPointGetRequiredType <-
function(object)
{
  checkPtrType(object, "GIOExtensionPoint")

  w <- .RGtkCall("S_g_io_extension_point_get_required_type", object, PACKAGE = "RGtk2")

  return(w)
} 


gIoExtensionPointGetExtensions <-
function(object)
{
  checkPtrType(object, "GIOExtensionPoint")

  w <- .RGtkCall("S_g_io_extension_point_get_extensions", object, PACKAGE = "RGtk2")

  return(w)
} 


gIoExtensionPointGetExtensionByName <-
function(object, name)
{
  checkPtrType(object, "GIOExtensionPoint")
  name <- as.character(name)

  w <- .RGtkCall("S_g_io_extension_point_get_extension_by_name", object, name, PACKAGE = "RGtk2")

  return(w)
} 


gIoExtensionPointImplement <-
function(extension.point.name, type, extension.name, priority)
{
  extension.point.name <- as.character(extension.point.name)
  type <- as.GType(type)
  extension.name <- as.character(extension.name)
  priority <- as.integer(priority)

  w <- .RGtkCall("S_g_io_extension_point_implement", extension.point.name, type, extension.name, priority, PACKAGE = "RGtk2")

  return(w)
} 


gIoExtensionGetType <-
function(object)
{
  checkPtrType(object, "GIOExtension")

  w <- .RGtkCall("S_g_io_extension_get_type", object, PACKAGE = "RGtk2")

  return(w)
} 


gIoExtensionGetName <-
function(object)
{
  checkPtrType(object, "GIOExtension")

  w <- .RGtkCall("S_g_io_extension_get_name", object, PACKAGE = "RGtk2")

  return(w)
} 


gIoExtensionGetPriority <-
function(object)
{
  checkPtrType(object, "GIOExtension")

  w <- .RGtkCall("S_g_io_extension_get_priority", object, PACKAGE = "RGtk2")

  return(w)
} 


gIoExtensionRefClass <-
function(object)
{
  checkPtrType(object, "GIOExtension")

  w <- .RGtkCall("S_g_io_extension_ref_class", object, PACKAGE = "RGtk2")

  return(w)
} 


gContentTypeFromMimeType <-
function(mime.type)
{
  mime.type <- as.character(mime.type)

  w <- .RGtkCall("S_g_content_type_from_mime_type", mime.type, PACKAGE = "RGtk2")

  return(w)
} 


gContentTypeGuessForTree <-
function(root)
{
  checkPtrType(root, "GFile")

  w <- .RGtkCall("S_g_content_type_guess_for_tree", root, PACKAGE = "RGtk2")

  return(w)
} 


gEmblemedIconGetType <-
function()
{
  

  w <- .RGtkCall("S_g_emblemed_icon_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gEmblemedIconNew <-
function(icon, emblem)
{
  checkPtrType(icon, "GIcon")
  checkPtrType(emblem, "GEmblem")

  w <- .RGtkCall("S_g_emblemed_icon_new", icon, emblem, PACKAGE = "RGtk2")

  return(w)
} 


gEmblemedIconGetIcon <-
function(object)
{
  checkPtrType(object, "GEmblemedIcon")

  w <- .RGtkCall("S_g_emblemed_icon_get_icon", object, PACKAGE = "RGtk2")

  return(w)
} 


gEmblemedIconGetEmblems <-
function(object)
{
  checkPtrType(object, "GEmblemedIcon")

  w <- .RGtkCall("S_g_emblemed_icon_get_emblems", object, PACKAGE = "RGtk2")

  return(w)
} 


gEmblemedIconAddEmblem <-
function(object, emblem)
{
  checkPtrType(object, "GEmblemedIcon")
  checkPtrType(emblem, "GEmblem")

  w <- .RGtkCall("S_g_emblemed_icon_add_emblem", object, emblem, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gEmblemGetType <-
function()
{
  

  w <- .RGtkCall("S_g_emblem_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gEmblemNew <-
function(icon = NULL, origin = NULL)
{
  

  w <- .RGtkCall("S_g_emblem_new", icon, origin, PACKAGE = "RGtk2")

  return(w)
} 


gEmblemNewWithOrigin <-
function(icon, origin)
{
  checkPtrType(icon, "GIcon")
  

  w <- .RGtkCall("S_g_emblem_new_with_origin", icon, origin, PACKAGE = "RGtk2")

  return(w)
} 


gEmblemGetIcon <-
function(object)
{
  checkPtrType(object, "GEmblem")

  w <- .RGtkCall("S_g_emblem_get_icon", object, PACKAGE = "RGtk2")

  return(w)
} 


gEmblemGetOrigin <-
function(object)
{
  checkPtrType(object, "GEmblem")

  w <- .RGtkCall("S_g_emblem_get_origin", object, PACKAGE = "RGtk2")

  return(w)
} 


gFileQueryFileType <-
function(object, flags, cancellable = NULL)
{
  checkPtrType(object, "GFile")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_query_file_type", object, flags, cancellable, PACKAGE = "RGtk2")

  return(w)
} 


gFileMakeDirectoryWithParents <-
function(object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_make_directory_with_parents", object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileMonitor <-
function(object, flags = "G_FILE_MONITOR_NONE", cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_monitor", object, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gMemoryOutputStreamGetDataSize <-
function(object)
{
  checkPtrType(object, "GMemoryOutputStream")

  w <- .RGtkCall("S_g_memory_output_stream_get_data_size", object, PACKAGE = "RGtk2")

  return(w)
} 


gMountGuessContentType <-
function(object, force.rescan, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GMount")
  force.rescan <- as.logical(force.rescan)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_mount_guess_content_type", object, force.rescan, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gMountGuessContentTypeFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GMount")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_mount_guess_content_type_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gMountGuessContentTypeSync <-
function(object, force.rescan, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GMount")
  force.rescan <- as.logical(force.rescan)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_mount_guess_content_type_sync", object, force.rescan, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gThemedIconPrependName <-
function(object, iconname)
{
  checkPtrType(object, "GThemedIcon")
  iconname <- as.character(iconname)

  w <- .RGtkCall("S_g_themed_icon_prepend_name", object, iconname, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gVolumeGetIdentifier <-
function(object, kind)
{
  checkPtrType(object, "GVolume")
  kind <- as.character(kind)

  w <- .RGtkCall("S_g_volume_get_identifier", object, kind, PACKAGE = "RGtk2")

  return(w)
} 


gVolumeEnumerateIdentifiers <-
function(object)
{
  checkPtrType(object, "GVolume")

  w <- .RGtkCall("S_g_volume_enumerate_identifiers", object, PACKAGE = "RGtk2")

  return(w)
} 


gVolumeGetActivationRoot <-
function(object)
{
  checkPtrType(object, "GVolume")

  w <- .RGtkCall("S_g_volume_get_activation_root", object, PACKAGE = "RGtk2")

  return(w)
} 


gFileEnumeratorGetContainer <-
function(object)
{
  checkPtrType(object, "GFileEnumerator")

  w <- .RGtkCall("S_g_file_enumerator_get_container", object, PACKAGE = "RGtk2")

  return(w)
} 


gAppInfoResetTypeAssociations <-
function(content.type)
{
  content.type <- as.character(content.type)

  w <- .RGtkCall("S_g_app_info_reset_type_associations", content.type, PACKAGE = "RGtk2")

  return(w)
} 


gAppInfoCanDelete <-
function(object)
{
  checkPtrType(object, "GAppInfo")

  w <- .RGtkCall("S_g_app_info_can_delete", object, PACKAGE = "RGtk2")

  return(w)
} 


gAppInfoDelete <-
function(object)
{
  checkPtrType(object, "GAppInfo")

  w <- .RGtkCall("S_g_app_info_delete", object, PACKAGE = "RGtk2")

  return(w)
} 


gAppInfoGetCommandline <-
function(object)
{
  checkPtrType(object, "GAppInfo")

  w <- .RGtkCall("S_g_app_info_get_commandline", object, PACKAGE = "RGtk2")

  return(w)
} 


gDataInputStreamReadUntilAsync <-
function(object, stop.chars, io.priority, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GDataInputStream")
  stop.chars <- as.character(stop.chars)
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_data_input_stream_read_until_async", object, stop.chars, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gDataInputStreamReadUntilFinish <-
function(object, result, length, .errwarn = TRUE)
{
  checkPtrType(object, "GDataInputStream")
  checkPtrType(result, "GAsyncResult")
  length <- as.list(as.numeric(length))

  w <- .RGtkCall("S_g_data_input_stream_read_until_finish", object, result, length, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gDataInputStreamReadLineAsync <-
function(object, io.priority, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GDataInputStream")
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_data_input_stream_read_line_async", object, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gDataInputStreamReadLineFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GDataInputStream")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_data_input_stream_read_line_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gIconToString <-
function(object)
{
  checkPtrType(object, "GIcon")

  w <- .RGtkCall("S_g_icon_to_string", object, PACKAGE = "RGtk2")

  return(w)
} 


gIconNewForString <-
function(str, .errwarn = TRUE)
{
  str <- as.character(str)

  w <- .RGtkCall("S_g_icon_new_for_string", str, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gMountIsShadowed <-
function(object)
{
  checkPtrType(object, "GMount")

  w <- .RGtkCall("S_g_mount_is_shadowed", object, PACKAGE = "RGtk2")

  return(w)
} 


gMountShadow <-
function(object)
{
  checkPtrType(object, "GMount")

  w <- .RGtkCall("S_g_mount_shadow", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gMountUnshadow <-
function(object)
{
  checkPtrType(object, "GMount")

  w <- .RGtkCall("S_g_mount_unshadow", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFilterInputStreamGetCloseBaseStream <-
function(object)
{
  checkPtrType(object, "GFilterInputStream")

  w <- .RGtkCall("S_g_filter_input_stream_get_close_base_stream", object, PACKAGE = "RGtk2")

  return(w)
} 


gFilterInputStreamSetCloseBaseStream <-
function(object, close.base)
{
  checkPtrType(object, "GFilterInputStream")
  close.base <- as.logical(close.base)

  w <- .RGtkCall("S_g_filter_input_stream_set_close_base_stream", object, close.base, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFilterOutputStreamGetCloseBaseStream <-
function(object)
{
  checkPtrType(object, "GFilterOutputStream")

  w <- .RGtkCall("S_g_filter_output_stream_get_close_base_stream", object, PACKAGE = "RGtk2")

  return(w)
} 


gFilterOutputStreamSetCloseBaseStream <-
function(object, close.base)
{
  checkPtrType(object, "GFilterOutputStream")
  close.base <- as.logical(close.base)

  w <- .RGtkCall("S_g_filter_output_stream_set_close_base_stream", object, close.base, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gAsyncInitableGetType <-
function()
{
  

  w <- .RGtkCall("S_g_async_initable_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gAsyncInitableInitAsync <-
function(object, io.priority, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GAsyncInitable")
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_async_initable_init_async", object, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gAsyncInitableInitFinish <-
function(object, res, .errwarn = TRUE)
{
  checkPtrType(object, "GAsyncInitable")
  checkPtrType(res, "GAsyncResult")

  w <- .RGtkCall("S_g_async_initable_init_finish", object, res, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gAsyncInitableNewFinish <-
function(object, res, .errwarn = TRUE)
{
  checkPtrType(object, "GAsyncInitable")
  checkPtrType(res, "GAsyncResult")

  w <- .RGtkCall("S_g_async_initable_new_finish", object, res, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gCancellableDisconnect <-
function(object, handler.id)
{
  checkPtrType(object, "GCancellable")
  handler.id <- as.numeric(handler.id)

  w <- .RGtkCall("S_g_cancellable_disconnect", object, handler.id, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gCancellableReleaseFd <-
function(object)
{
  checkPtrType(object, "GCancellable")

  w <- .RGtkCall("S_g_cancellable_release_fd", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gDriveCanStart <-
function(object)
{
  checkPtrType(object, "GDrive")

  w <- .RGtkCall("S_g_drive_can_start", object, PACKAGE = "RGtk2")

  return(w)
} 


gDriveCanStartDegraded <-
function(object)
{
  checkPtrType(object, "GDrive")

  w <- .RGtkCall("S_g_drive_can_start_degraded", object, PACKAGE = "RGtk2")

  return(w)
} 


gDriveCanStop <-
function(object)
{
  checkPtrType(object, "GDrive")

  w <- .RGtkCall("S_g_drive_can_stop", object, PACKAGE = "RGtk2")

  return(w)
} 


gDriveEjectWithOperation <-
function(object, flags, mount.operation, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GDrive")
  
  checkPtrType(mount.operation, "GMountOperation")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_drive_eject_with_operation", object, flags, mount.operation, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gDriveEjectWithOperationFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GDrive")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_drive_eject_with_operation_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gDriveGetStartStopType <-
function(object)
{
  checkPtrType(object, "GDrive")

  w <- .RGtkCall("S_g_drive_get_start_stop_type", object, PACKAGE = "RGtk2")

  return(w)
} 


gDriveStart <-
function(object, flags, mount.operation, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GDrive")
  
  checkPtrType(mount.operation, "GMountOperation")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_drive_start", object, flags, mount.operation, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gDriveStartFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GDrive")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_drive_start_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gDriveStop <-
function(object, flags, mount.operation, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GDrive")
  
  checkPtrType(mount.operation, "GMountOperation")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_drive_stop", object, flags, mount.operation, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gDriveStopFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GDrive")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_drive_stop_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileCreateReadwrite <-
function(object, flags, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_create_readwrite", object, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileCreateReadwriteAsync <-
function(object, flags, io.priority, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GFile")
  
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_file_create_readwrite_async", object, flags, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileCreateReadwriteFinish <-
function(object, res, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  checkPtrType(res, "GAsyncResult")

  w <- .RGtkCall("S_g_file_create_readwrite_finish", object, res, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileEjectMountableWithOperation <-
function(object, flags, mount.operation, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GFile")
  
  checkPtrType(mount.operation, "GMountOperation")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_file_eject_mountable_with_operation", object, flags, mount.operation, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileEjectMountableWithOperationFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_file_eject_mountable_with_operation_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileOpenReadwrite <-
function(object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_open_readwrite", object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileOpenReadwriteAsync <-
function(object, io.priority, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GFile")
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_file_open_readwrite_async", object, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileOpenReadwriteFinish <-
function(object, res, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  checkPtrType(res, "GAsyncResult")

  w <- .RGtkCall("S_g_file_open_readwrite_finish", object, res, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFilePollMountable <-
function(object, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GFile")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_file_poll_mountable", object, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFilePollMountableFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_file_poll_mountable_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileReplaceReadwrite <-
function(object, etag, make.backup, flags, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  etag <- as.character(etag)
  make.backup <- as.logical(make.backup)
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_replace_readwrite", object, etag, make.backup, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileReplaceReadwriteAsync <-
function(object, etag, make.backup, flags, io.priority, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GFile")
  etag <- as.character(etag)
  make.backup <- as.logical(make.backup)
  
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_file_replace_readwrite_async", object, etag, make.backup, flags, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileReplaceReadwriteFinish <-
function(object, res, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  checkPtrType(res, "GAsyncResult")

  w <- .RGtkCall("S_g_file_replace_readwrite_finish", object, res, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileStartMountable <-
function(object, flags, start.operation, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GFile")
  
  checkPtrType(start.operation, "GMountOperation")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_file_start_mountable", object, flags, start.operation, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileStartMountableFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_file_start_mountable_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileStopMountable <-
function(object, flags, mount.operation, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GFile")
  
  checkPtrType(mount.operation, "GMountOperation")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_file_stop_mountable", object, flags, mount.operation, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileStopMountableFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_file_stop_mountable_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileSupportsThreadContexts <-
function(object)
{
  checkPtrType(object, "GFile")

  w <- .RGtkCall("S_g_file_supports_thread_contexts", object, PACKAGE = "RGtk2")

  return(w)
} 


gFileUnmountMountableWithOperation <-
function(object, flags, mount.operation, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GFile")
  
  checkPtrType(mount.operation, "GMountOperation")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_file_unmount_mountable_with_operation", object, flags, mount.operation, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileUnmountMountableWithOperationFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GFile")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_file_unmount_mountable_with_operation_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileInfoHasNamespace <-
function(object, name.space)
{
  checkPtrType(object, "GFileInfo")
  name.space <- as.character(name.space)

  w <- .RGtkCall("S_g_file_info_has_namespace", object, name.space, PACKAGE = "RGtk2")

  return(w)
} 


gFileInfoSetAttributeStatus <-
function(object, attribute, status)
{
  checkPtrType(object, "GFileInfo")
  attribute <- as.character(attribute)
  

  w <- .RGtkCall("S_g_file_info_set_attribute_status", object, attribute, status, PACKAGE = "RGtk2")

  return(w)
} 


gFileInfoGetAttributeStringv <-
function(object, attribute)
{
  checkPtrType(object, "GFileInfo")
  attribute <- as.character(attribute)

  w <- .RGtkCall("S_g_file_info_get_attribute_stringv", object, attribute, PACKAGE = "RGtk2")

  return(w)
} 


gFileInfoSetAttributeStringv <-
function(object, attribute, attr.value)
{
  checkPtrType(object, "GFileInfo")
  attribute <- as.character(attribute)
  attr.value <- as.list(as.character(attr.value))

  w <- .RGtkCall("S_g_file_info_set_attribute_stringv", object, attribute, attr.value, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileIOStreamGetType <-
function()
{
  

  w <- .RGtkCall("S_g_file_io_stream_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gFileIOStreamQueryInfo <-
function(object, attributes, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GFileIOStream")
  attributes <- as.character(attributes)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_file_io_stream_query_info", object, attributes, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileIOStreamQueryInfoAsync <-
function(object, attributes, io.priority, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GFileIOStream")
  attributes <- as.character(attributes)
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_file_io_stream_query_info_async", object, attributes, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gFileIOStreamQueryInfoFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GFileIOStream")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_file_io_stream_query_info_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gFileIOStreamGetEtag <-
function(object)
{
  checkPtrType(object, "GFileIOStream")

  w <- .RGtkCall("S_g_file_io_stream_get_etag", object, PACKAGE = "RGtk2")

  return(w)
} 


gInetAddressGetType <-
function()
{
  

  w <- .RGtkCall("S_g_inet_address_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gInetAddressNewFromString <-
function(string)
{
  string <- as.character(string)

  w <- .RGtkCall("S_g_inet_address_new_from_string", string, PACKAGE = "RGtk2")

  return(w)
} 


gInetAddressNewFromBytes <-
function(bytes, family)
{
  bytes <- as.list(as.raw(bytes))
  

  w <- .RGtkCall("S_g_inet_address_new_from_bytes", bytes, family, PACKAGE = "RGtk2")

  return(w)
} 


gInetAddressNewLoopback <-
function(family)
{
  

  w <- .RGtkCall("S_g_inet_address_new_loopback", family, PACKAGE = "RGtk2")

  return(w)
} 


gInetAddressNewAny <-
function(family)
{
  

  w <- .RGtkCall("S_g_inet_address_new_any", family, PACKAGE = "RGtk2")

  return(w)
} 


gInetAddressToString <-
function(object)
{
  checkPtrType(object, "GInetAddress")

  w <- .RGtkCall("S_g_inet_address_to_string", object, PACKAGE = "RGtk2")

  return(w)
} 


gInetAddressToBytes <-
function(object)
{
  checkPtrType(object, "GInetAddress")

  w <- .RGtkCall("S_g_inet_address_to_bytes", object, PACKAGE = "RGtk2")

  return(w)
} 


gInetAddressGetNativeSize <-
function(object)
{
  checkPtrType(object, "GInetAddress")

  w <- .RGtkCall("S_g_inet_address_get_native_size", object, PACKAGE = "RGtk2")

  return(w)
} 


gInetAddressGetFamily <-
function(object)
{
  checkPtrType(object, "GInetAddress")

  w <- .RGtkCall("S_g_inet_address_get_family", object, PACKAGE = "RGtk2")

  return(w)
} 


gInetAddressGetIsAny <-
function(object)
{
  checkPtrType(object, "GInetAddress")

  w <- .RGtkCall("S_g_inet_address_get_is_any", object, PACKAGE = "RGtk2")

  return(w)
} 


gInetAddressGetIsLoopback <-
function(object)
{
  checkPtrType(object, "GInetAddress")

  w <- .RGtkCall("S_g_inet_address_get_is_loopback", object, PACKAGE = "RGtk2")

  return(w)
} 


gInetAddressGetIsLinkLocal <-
function(object)
{
  checkPtrType(object, "GInetAddress")

  w <- .RGtkCall("S_g_inet_address_get_is_link_local", object, PACKAGE = "RGtk2")

  return(w)
} 


gInetAddressGetIsSiteLocal <-
function(object)
{
  checkPtrType(object, "GInetAddress")

  w <- .RGtkCall("S_g_inet_address_get_is_site_local", object, PACKAGE = "RGtk2")

  return(w)
} 


gInetAddressGetIsMulticast <-
function(object)
{
  checkPtrType(object, "GInetAddress")

  w <- .RGtkCall("S_g_inet_address_get_is_multicast", object, PACKAGE = "RGtk2")

  return(w)
} 


gInetAddressGetIsMcGlobal <-
function(object)
{
  checkPtrType(object, "GInetAddress")

  w <- .RGtkCall("S_g_inet_address_get_is_mc_global", object, PACKAGE = "RGtk2")

  return(w)
} 


gInetAddressGetIsMcLinkLocal <-
function(object)
{
  checkPtrType(object, "GInetAddress")

  w <- .RGtkCall("S_g_inet_address_get_is_mc_link_local", object, PACKAGE = "RGtk2")

  return(w)
} 


gInetAddressGetIsMcNodeLocal <-
function(object)
{
  checkPtrType(object, "GInetAddress")

  w <- .RGtkCall("S_g_inet_address_get_is_mc_node_local", object, PACKAGE = "RGtk2")

  return(w)
} 


gInetAddressGetIsMcOrgLocal <-
function(object)
{
  checkPtrType(object, "GInetAddress")

  w <- .RGtkCall("S_g_inet_address_get_is_mc_org_local", object, PACKAGE = "RGtk2")

  return(w)
} 


gInetAddressGetIsMcSiteLocal <-
function(object)
{
  checkPtrType(object, "GInetAddress")

  w <- .RGtkCall("S_g_inet_address_get_is_mc_site_local", object, PACKAGE = "RGtk2")

  return(w)
} 


gInitableGetType <-
function()
{
  

  w <- .RGtkCall("S_g_initable_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gInitableInit <-
function(object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GInitable")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_initable_init", object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gIOStreamGetType <-
function()
{
  

  w <- .RGtkCall("S_g_io_stream_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gIOStreamGetInputStream <-
function(object)
{
  checkPtrType(object, "GIOStream")

  w <- .RGtkCall("S_g_io_stream_get_input_stream", object, PACKAGE = "RGtk2")

  return(w)
} 


gIOStreamGetOutputStream <-
function(object)
{
  checkPtrType(object, "GIOStream")

  w <- .RGtkCall("S_g_io_stream_get_output_stream", object, PACKAGE = "RGtk2")

  return(w)
} 


gIOStreamClose <-
function(object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GIOStream")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_io_stream_close", object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gIOStreamCloseAsync <-
function(object, io.priority, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GIOStream")
  io.priority <- as.integer(io.priority)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_io_stream_close_async", object, io.priority, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gIOStreamCloseFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GIOStream")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_io_stream_close_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gIOStreamIsClosed <-
function(object)
{
  checkPtrType(object, "GIOStream")

  w <- .RGtkCall("S_g_io_stream_is_closed", object, PACKAGE = "RGtk2")

  return(w)
} 


gIOStreamHasPending <-
function(object)
{
  checkPtrType(object, "GIOStream")

  w <- .RGtkCall("S_g_io_stream_has_pending", object, PACKAGE = "RGtk2")

  return(w)
} 


gIOStreamSetPending <-
function(object, .errwarn = TRUE)
{
  checkPtrType(object, "GIOStream")

  w <- .RGtkCall("S_g_io_stream_set_pending", object, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gIOStreamClearPending <-
function(object)
{
  checkPtrType(object, "GIOStream")

  w <- .RGtkCall("S_g_io_stream_clear_pending", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gMountUnmountWithOperation <-
function(object, flags, mount.operation, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GMount")
  
  checkPtrType(mount.operation, "GMountOperation")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_mount_unmount_with_operation", object, flags, mount.operation, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gMountUnmountWithOperationFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GMount")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_mount_unmount_with_operation_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gMountEjectWithOperation <-
function(object, flags, mount.operation, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GMount")
  
  checkPtrType(mount.operation, "GMountOperation")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_mount_eject_with_operation", object, flags, mount.operation, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gMountEjectWithOperationFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GMount")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_mount_eject_with_operation_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gNetworkAddressGetType <-
function()
{
  

  w <- .RGtkCall("S_g_network_address_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gNetworkAddressNew <-
function(hostname, port)
{
  hostname <- as.character(hostname)
  port <- as.integer(port)

  w <- .RGtkCall("S_g_network_address_new", hostname, port, PACKAGE = "RGtk2")

  return(w)
} 


gNetworkAddressParse <-
function(host.and.port, default.port, .errwarn = TRUE)
{
  host.and.port <- as.character(host.and.port)
  default.port <- as.integer(default.port)

  w <- .RGtkCall("S_g_network_address_parse", host.and.port, default.port, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gNetworkAddressGetHostname <-
function(object)
{
  checkPtrType(object, "GNetworkAddress")

  w <- .RGtkCall("S_g_network_address_get_hostname", object, PACKAGE = "RGtk2")

  return(w)
} 


gNetworkAddressGetPort <-
function(object)
{
  checkPtrType(object, "GNetworkAddress")

  w <- .RGtkCall("S_g_network_address_get_port", object, PACKAGE = "RGtk2")

  return(w)
} 


gNetworkServiceGetType <-
function()
{
  

  w <- .RGtkCall("S_g_network_service_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gNetworkServiceNew <-
function(service, protocol, domain)
{
  service <- as.character(service)
  protocol <- as.character(protocol)
  domain <- as.character(domain)

  w <- .RGtkCall("S_g_network_service_new", service, protocol, domain, PACKAGE = "RGtk2")

  return(w)
} 


gNetworkServiceGetService <-
function(object)
{
  checkPtrType(object, "GNetworkService")

  w <- .RGtkCall("S_g_network_service_get_service", object, PACKAGE = "RGtk2")

  return(w)
} 


gNetworkServiceGetProtocol <-
function(object)
{
  checkPtrType(object, "GNetworkService")

  w <- .RGtkCall("S_g_network_service_get_protocol", object, PACKAGE = "RGtk2")

  return(w)
} 


gNetworkServiceGetDomain <-
function(object)
{
  checkPtrType(object, "GNetworkService")

  w <- .RGtkCall("S_g_network_service_get_domain", object, PACKAGE = "RGtk2")

  return(w)
} 


gResolverGetType <-
function()
{
  

  w <- .RGtkCall("S_g_resolver_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gResolverGetDefault <-
function()
{
  

  w <- .RGtkCall("S_g_resolver_get_default", PACKAGE = "RGtk2")

  return(w)
} 


gResolverSetDefault <-
function(object)
{
  checkPtrType(object, "GResolver")

  w <- .RGtkCall("S_g_resolver_set_default", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gResolverLookupByName <-
function(object, hostname, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GResolver")
  hostname <- as.character(hostname)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_resolver_lookup_by_name", object, hostname, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gResolverLookupByNameAsync <-
function(object, hostname, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GResolver")
  hostname <- as.character(hostname)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_resolver_lookup_by_name_async", object, hostname, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gResolverLookupByNameFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GResolver")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_resolver_lookup_by_name_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gResolverFreeAddresses <-
function(addresses)
{
  addresses <- as.GList(addresses)

  w <- .RGtkCall("S_g_resolver_free_addresses", addresses, PACKAGE = "RGtk2")

  return(w)
} 


gResolverLookupByAddress <-
function(object, address, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GResolver")
  checkPtrType(address, "GInetAddress")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_resolver_lookup_by_address", object, address, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gResolverLookupByAddressAsync <-
function(object, address, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GResolver")
  checkPtrType(address, "GInetAddress")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_resolver_lookup_by_address_async", object, address, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gResolverLookupByAddressFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GResolver")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_resolver_lookup_by_address_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gResolverLookupService <-
function(object, service, protocol, domain, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GResolver")
  service <- as.character(service)
  protocol <- as.character(protocol)
  domain <- as.character(domain)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_resolver_lookup_service", object, service, protocol, domain, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gResolverLookupServiceAsync <-
function(object, service, protocol, domain, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GResolver")
  service <- as.character(service)
  protocol <- as.character(protocol)
  domain <- as.character(domain)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_resolver_lookup_service_async", object, service, protocol, domain, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gResolverLookupServiceFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GResolver")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_resolver_lookup_service_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gResolverFreeTargets <-
function(targets)
{
  targets <- as.GList(targets)

  w <- .RGtkCall("S_g_resolver_free_targets", targets, PACKAGE = "RGtk2")

  return(w)
} 


gResolverErrorQuark <-
function()
{
  

  w <- .RGtkCall("S_g_resolver_error_quark", PACKAGE = "RGtk2")

  return(w)
} 


gSocketAddressEnumeratorGetType <-
function()
{
  

  w <- .RGtkCall("S_g_socket_address_enumerator_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gSocketAddressEnumeratorNext <-
function(object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GSocketAddressEnumerator")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_socket_address_enumerator_next", object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketAddressEnumeratorNextAsync <-
function(object, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GSocketAddressEnumerator")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_socket_address_enumerator_next_async", object, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gSocketAddressEnumeratorNextFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GSocketAddressEnumerator")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_socket_address_enumerator_next_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketAddressGetType <-
function()
{
  

  w <- .RGtkCall("S_g_socket_address_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gSocketAddressGetFamily <-
function(object)
{
  checkPtrType(object, "GSocketAddress")

  w <- .RGtkCall("S_g_socket_address_get_family", object, PACKAGE = "RGtk2")

  return(w)
} 


gSocketAddressNewFromNative <-
function(native, len)
{
  
  len <- as.numeric(len)

  w <- .RGtkCall("S_g_socket_address_new_from_native", native, len, PACKAGE = "RGtk2")

  return(w)
} 


gSocketAddressToNative <-
function(object, dest, destlen, .errwarn = TRUE)
{
  checkPtrType(object, "GSocketAddress")
  
  destlen <- as.numeric(destlen)

  w <- .RGtkCall("S_g_socket_address_to_native", object, dest, destlen, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketAddressGetNativeSize <-
function(object)
{
  checkPtrType(object, "GSocketAddress")

  w <- .RGtkCall("S_g_socket_address_get_native_size", object, PACKAGE = "RGtk2")

  return(w)
} 


gSocketClientGetType <-
function()
{
  

  w <- .RGtkCall("S_g_socket_client_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gSocketClientNew <-
function()
{
  

  w <- .RGtkCall("S_g_socket_client_new", PACKAGE = "RGtk2")

  return(w)
} 


gSocketClientGetFamily <-
function(object)
{
  checkPtrType(object, "GSocketClient")

  w <- .RGtkCall("S_g_socket_client_get_family", object, PACKAGE = "RGtk2")

  return(w)
} 


gSocketClientSetFamily <-
function(object, family)
{
  checkPtrType(object, "GSocketClient")
  

  w <- .RGtkCall("S_g_socket_client_set_family", object, family, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gSocketClientGetSocketType <-
function(object)
{
  checkPtrType(object, "GSocketClient")

  w <- .RGtkCall("S_g_socket_client_get_socket_type", object, PACKAGE = "RGtk2")

  return(w)
} 


gSocketClientSetSocketType <-
function(object, type)
{
  checkPtrType(object, "GSocketClient")
  

  w <- .RGtkCall("S_g_socket_client_set_socket_type", object, type, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gSocketClientGetProtocol <-
function(object)
{
  checkPtrType(object, "GSocketClient")

  w <- .RGtkCall("S_g_socket_client_get_protocol", object, PACKAGE = "RGtk2")

  return(w)
} 


gSocketClientSetProtocol <-
function(object, protocol)
{
  checkPtrType(object, "GSocketClient")
  

  w <- .RGtkCall("S_g_socket_client_set_protocol", object, protocol, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gSocketClientGetLocalAddress <-
function(object)
{
  checkPtrType(object, "GSocketClient")

  w <- .RGtkCall("S_g_socket_client_get_local_address", object, PACKAGE = "RGtk2")

  return(w)
} 


gSocketClientSetLocalAddress <-
function(object, address)
{
  checkPtrType(object, "GSocketClient")
  checkPtrType(address, "GSocketAddress")

  w <- .RGtkCall("S_g_socket_client_set_local_address", object, address, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gSocketClientConnect <-
function(object, connectable, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GSocketClient")
  checkPtrType(connectable, "GSocketConnectable")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_socket_client_connect", object, connectable, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketClientConnectToHost <-
function(object, host.and.port, default.port, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GSocketClient")
  host.and.port <- as.character(host.and.port)
  default.port <- as.integer(default.port)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_socket_client_connect_to_host", object, host.and.port, default.port, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketClientConnectToService <-
function(object, domain, service, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GSocketClient")
  domain <- as.character(domain)
  service <- as.character(service)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_socket_client_connect_to_service", object, domain, service, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketClientConnectAsync <-
function(object, connectable, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GSocketClient")
  checkPtrType(connectable, "GSocketConnectable")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_socket_client_connect_async", object, connectable, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gSocketClientConnectFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GSocketClient")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_socket_client_connect_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketClientConnectToHostAsync <-
function(object, host.and.port, default.port, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GSocketClient")
  host.and.port <- as.character(host.and.port)
  default.port <- as.integer(default.port)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_socket_client_connect_to_host_async", object, host.and.port, default.port, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gSocketClientConnectToHostFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GSocketClient")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_socket_client_connect_to_host_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketClientConnectToServiceAsync <-
function(object, domain, service, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GSocketClient")
  domain <- as.character(domain)
  service <- as.character(service)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_socket_client_connect_to_service_async", object, domain, service, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gSocketClientConnectToServiceFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GSocketClient")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_socket_client_connect_to_service_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketConnectableGetType <-
function()
{
  

  w <- .RGtkCall("S_g_socket_connectable_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gSocketConnectableEnumerate <-
function(object)
{
  checkPtrType(object, "GSocketConnectable")

  w <- .RGtkCall("S_g_socket_connectable_enumerate", object, PACKAGE = "RGtk2")

  return(w)
} 


gSocketConnectionGetType <-
function()
{
  

  w <- .RGtkCall("S_g_socket_connection_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gSocketConnectionGetSocket <-
function(object)
{
  checkPtrType(object, "GSocketConnection")

  w <- .RGtkCall("S_g_socket_connection_get_socket", object, PACKAGE = "RGtk2")

  return(w)
} 


gSocketConnectionGetLocalAddress <-
function(object, .errwarn = TRUE)
{
  checkPtrType(object, "GSocketConnection")

  w <- .RGtkCall("S_g_socket_connection_get_local_address", object, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketConnectionGetRemoteAddress <-
function(object, .errwarn = TRUE)
{
  checkPtrType(object, "GSocketConnection")

  w <- .RGtkCall("S_g_socket_connection_get_remote_address", object, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketConnectionFactoryRegisterType <-
function(g.type, family, type, protocol)
{
  g.type <- as.GType(g.type)
  
  
  protocol <- as.integer(protocol)

  w <- .RGtkCall("S_g_socket_connection_factory_register_type", g.type, family, type, protocol, PACKAGE = "RGtk2")

  return(w)
} 


gSocketConnectionFactoryLookupType <-
function(family, type, protocol.id)
{
  
  
  protocol.id <- as.integer(protocol.id)

  w <- .RGtkCall("S_g_socket_connection_factory_lookup_type", family, type, protocol.id, PACKAGE = "RGtk2")

  return(w)
} 


gSocketConnectionFactoryCreateConnection <-
function(object)
{
  checkPtrType(object, "GSocket")

  w <- .RGtkCall("S_g_socket_connection_factory_create_connection", object, PACKAGE = "RGtk2")

  return(w)
} 


gSocketControlMessageGetType <-
function()
{
  

  w <- .RGtkCall("S_g_socket_control_message_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gSocketControlMessageGetSize <-
function(object)
{
  checkPtrType(object, "GSocketControlMessage")

  w <- .RGtkCall("S_g_socket_control_message_get_size", object, PACKAGE = "RGtk2")

  return(w)
} 


gSocketControlMessageGetLevel <-
function(object)
{
  checkPtrType(object, "GSocketControlMessage")

  w <- .RGtkCall("S_g_socket_control_message_get_level", object, PACKAGE = "RGtk2")

  return(w)
} 


gSocketControlMessageGetMsgType <-
function(object)
{
  checkPtrType(object, "GSocketControlMessage")

  w <- .RGtkCall("S_g_socket_control_message_get_msg_type", object, PACKAGE = "RGtk2")

  return(w)
} 


gSocketControlMessageSerialize <-
function(object, data)
{
  checkPtrType(object, "GSocketControlMessage")
  

  w <- .RGtkCall("S_g_socket_control_message_serialize", object, data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gSocketControlMessageDeserialize <-
function(level, type, size, data)
{
  level <- as.integer(level)
  type <- as.integer(type)
  size <- as.numeric(size)
  

  w <- .RGtkCall("S_g_socket_control_message_deserialize", level, type, size, data, PACKAGE = "RGtk2")

  return(w)
} 


gSocketGetType <-
function()
{
  

  w <- .RGtkCall("S_g_socket_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gSocketNew <-
function(family, type, protocol, .errwarn = TRUE)
{
  
  
  

  w <- .RGtkCall("S_g_socket_new", family, type, protocol, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketNewFromFd <-
function(fd, .errwarn = TRUE)
{
  fd <- as.integer(fd)

  w <- .RGtkCall("S_g_socket_new_from_fd", fd, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketGetFd <-
function(object)
{
  checkPtrType(object, "GSocket")

  w <- .RGtkCall("S_g_socket_get_fd", object, PACKAGE = "RGtk2")

  return(w)
} 


gSocketGetFamily <-
function(object)
{
  checkPtrType(object, "GSocket")

  w <- .RGtkCall("S_g_socket_get_family", object, PACKAGE = "RGtk2")

  return(w)
} 


gSocketGetSocketType <-
function(object)
{
  checkPtrType(object, "GSocket")

  w <- .RGtkCall("S_g_socket_get_socket_type", object, PACKAGE = "RGtk2")

  return(w)
} 


gSocketGetProtocol <-
function(object)
{
  checkPtrType(object, "GSocket")

  w <- .RGtkCall("S_g_socket_get_protocol", object, PACKAGE = "RGtk2")

  return(w)
} 


gSocketGetLocalAddress <-
function(object, .errwarn = TRUE)
{
  checkPtrType(object, "GSocket")

  w <- .RGtkCall("S_g_socket_get_local_address", object, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketGetRemoteAddress <-
function(object, .errwarn = TRUE)
{
  checkPtrType(object, "GSocket")

  w <- .RGtkCall("S_g_socket_get_remote_address", object, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketSetBlocking <-
function(object, blocking)
{
  checkPtrType(object, "GSocket")
  blocking <- as.logical(blocking)

  w <- .RGtkCall("S_g_socket_set_blocking", object, blocking, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gSocketGetBlocking <-
function(object)
{
  checkPtrType(object, "GSocket")

  w <- .RGtkCall("S_g_socket_get_blocking", object, PACKAGE = "RGtk2")

  return(w)
} 


gSocketSetKeepalive <-
function(object, keepalive)
{
  checkPtrType(object, "GSocket")
  keepalive <- as.logical(keepalive)

  w <- .RGtkCall("S_g_socket_set_keepalive", object, keepalive, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gSocketGetKeepalive <-
function(object)
{
  checkPtrType(object, "GSocket")

  w <- .RGtkCall("S_g_socket_get_keepalive", object, PACKAGE = "RGtk2")

  return(w)
} 


gSocketGetListenBacklog <-
function(object)
{
  checkPtrType(object, "GSocket")

  w <- .RGtkCall("S_g_socket_get_listen_backlog", object, PACKAGE = "RGtk2")

  return(w)
} 


gSocketSetListenBacklog <-
function(object, backlog)
{
  checkPtrType(object, "GSocket")
  backlog <- as.integer(backlog)

  w <- .RGtkCall("S_g_socket_set_listen_backlog", object, backlog, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gSocketIsConnected <-
function(object)
{
  checkPtrType(object, "GSocket")

  w <- .RGtkCall("S_g_socket_is_connected", object, PACKAGE = "RGtk2")

  return(w)
} 


gSocketBind <-
function(object, address, allow.reuse, .errwarn = TRUE)
{
  checkPtrType(object, "GSocket")
  checkPtrType(address, "GSocketAddress")
  allow.reuse <- as.logical(allow.reuse)

  w <- .RGtkCall("S_g_socket_bind", object, address, allow.reuse, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketConnect <-
function(object, address, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GSocket")
  checkPtrType(address, "GSocketAddress")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_socket_connect", object, address, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketCheckConnectResult <-
function(object, .errwarn = TRUE)
{
  checkPtrType(object, "GSocket")

  w <- .RGtkCall("S_g_socket_check_connect_result", object, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketConditionCheck <-
function(object, condition)
{
  checkPtrType(object, "GSocket")
  

  w <- .RGtkCall("S_g_socket_condition_check", object, condition, PACKAGE = "RGtk2")

  return(w)
} 


gSocketConditionWait <-
function(object, condition, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GSocket")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_socket_condition_wait", object, condition, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketAccept <-
function(object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GSocket")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_socket_accept", object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketListen <-
function(object, .errwarn = TRUE)
{
  checkPtrType(object, "GSocket")

  w <- .RGtkCall("S_g_socket_listen", object, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketReceive <-
function(object, size, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GSocket")
  size <- as.numeric(size)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_socket_receive", object, size, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketReceiveFrom <-
function(object, size, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GSocket")
  size <- as.numeric(size)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_socket_receive_from", object, size, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketSend <-
function(object, buffer, size, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GSocket")
  buffer <- as.character(buffer)
  size <- as.numeric(size)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_socket_send", object, buffer, size, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketSendTo <-
function(object, address, buffer, size, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GSocket")
  checkPtrType(address, "GSocketAddress")
  buffer <- as.character(buffer)
  size <- as.numeric(size)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_socket_send_to", object, address, buffer, size, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketReceiveMessage <-
function(object, flags = 0, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GSocket")
  flags <- as.list(as.integer(flags))
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_socket_receive_message", object, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketSendMessage <-
function(object, address, vectors, messages = NULL, flags = 0, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GSocket")
  checkPtrType(address, "GSocketAddress")
  vectors <- lapply(vectors, function(x) { checkPtrType(x, "GOutputVector"); x })
  messages <- lapply(messages, function(x) { checkPtrType(x, "GSocketControlMessage"); x })
  flags <- as.integer(flags)
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_socket_send_message", object, address, vectors, messages, flags, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketClose <-
function(object, .errwarn = TRUE)
{
  checkPtrType(object, "GSocket")

  w <- .RGtkCall("S_g_socket_close", object, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketShutdown <-
function(object, shutdown.read, shutdown.write, .errwarn = TRUE)
{
  checkPtrType(object, "GSocket")
  shutdown.read <- as.logical(shutdown.read)
  shutdown.write <- as.logical(shutdown.write)

  w <- .RGtkCall("S_g_socket_shutdown", object, shutdown.read, shutdown.write, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketIsClosed <-
function(object)
{
  checkPtrType(object, "GSocket")

  w <- .RGtkCall("S_g_socket_is_closed", object, PACKAGE = "RGtk2")

  return(w)
} 


gSocketCreateSource <-
function(object, condition, cancellable = NULL)
{
  checkPtrType(object, "GSocket")
  
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_socket_create_source", object, condition, cancellable, PACKAGE = "RGtk2")

  return(w)
} 


gSocketSpeaksIpv4 <-
function(object)
{
  checkPtrType(object, "GSocket")

  w <- .RGtkCall("S_g_socket_speaks_ipv4", object, PACKAGE = "RGtk2")

  return(w)
} 


gSocketListenerGetType <-
function()
{
  

  w <- .RGtkCall("S_g_socket_listener_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gSocketListenerNew <-
function()
{
  

  w <- .RGtkCall("S_g_socket_listener_new", PACKAGE = "RGtk2")

  return(w)
} 


gSocketListenerSetBacklog <-
function(object, listen.backlog)
{
  checkPtrType(object, "GSocketListener")
  listen.backlog <- as.integer(listen.backlog)

  w <- .RGtkCall("S_g_socket_listener_set_backlog", object, listen.backlog, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gSocketListenerAddSocket <-
function(object, socket, source.object, .errwarn = TRUE)
{
  checkPtrType(object, "GSocketListener")
  checkPtrType(socket, "GSocket")
  checkPtrType(source.object, "GObject")

  w <- .RGtkCall("S_g_socket_listener_add_socket", object, socket, source.object, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketListenerAddAddress <-
function(object, address, type, protocol, source.object, .errwarn = TRUE)
{
  checkPtrType(object, "GSocketListener")
  checkPtrType(address, "GSocketAddress")
  
  
  checkPtrType(source.object, "GObject")

  w <- .RGtkCall("S_g_socket_listener_add_address", object, address, type, protocol, source.object, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketListenerAddInetPort <-
function(object, port, source.object, .errwarn = TRUE)
{
  checkPtrType(object, "GSocketListener")
  port <- as.integer(port)
  checkPtrType(source.object, "GObject")

  w <- .RGtkCall("S_g_socket_listener_add_inet_port", object, port, source.object, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketListenerAcceptSocket <-
function(object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GSocketListener")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_socket_listener_accept_socket", object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketListenerAcceptSocketAsync <-
function(object, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GSocketListener")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_socket_listener_accept_socket_async", object, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gSocketListenerAcceptSocketFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GSocketListener")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_socket_listener_accept_socket_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketListenerAccept <-
function(object, cancellable = NULL, .errwarn = TRUE)
{
  checkPtrType(object, "GSocketListener")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")

  w <- .RGtkCall("S_g_socket_listener_accept", object, cancellable, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketListenerAcceptAsync <-
function(object, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GSocketListener")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_socket_listener_accept_async", object, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gSocketListenerAcceptFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GSocketListener")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_socket_listener_accept_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gSocketListenerClose <-
function(object)
{
  checkPtrType(object, "GSocketListener")

  w <- .RGtkCall("S_g_socket_listener_close", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gSocketServiceGetType <-
function()
{
  

  w <- .RGtkCall("S_g_socket_service_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gSocketServiceNew <-
function()
{
  

  w <- .RGtkCall("S_g_socket_service_new", PACKAGE = "RGtk2")

  return(w)
} 


gSocketServiceStart <-
function(object)
{
  checkPtrType(object, "GSocketService")

  w <- .RGtkCall("S_g_socket_service_start", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gSocketServiceStop <-
function(object)
{
  checkPtrType(object, "GSocketService")

  w <- .RGtkCall("S_g_socket_service_stop", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gSocketServiceIsActive <-
function(object)
{
  checkPtrType(object, "GSocketService")

  w <- .RGtkCall("S_g_socket_service_is_active", object, PACKAGE = "RGtk2")

  return(w)
} 


gSrvTargetGetType <-
function()
{
  

  w <- .RGtkCall("S_g_srv_target_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gSrvTargetNew <-
function(hostname, port, priority, weight)
{
  hostname <- as.character(hostname)
  port <- as.integer(port)
  priority <- as.integer(priority)
  weight <- as.integer(weight)

  w <- .RGtkCall("S_g_srv_target_new", hostname, port, priority, weight, PACKAGE = "RGtk2")

  return(w)
} 


gSrvTargetCopy <-
function(object)
{
  checkPtrType(object, "GSrvTarget")

  w <- .RGtkCall("S_g_srv_target_copy", object, PACKAGE = "RGtk2")

  return(w)
} 


gSrvTargetFree <-
function(object)
{
  checkPtrType(object, "GSrvTarget")

  w <- .RGtkCall("S_g_srv_target_free", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gSrvTargetGetHostname <-
function(object)
{
  checkPtrType(object, "GSrvTarget")

  w <- .RGtkCall("S_g_srv_target_get_hostname", object, PACKAGE = "RGtk2")

  return(w)
} 


gSrvTargetGetPort <-
function(object)
{
  checkPtrType(object, "GSrvTarget")

  w <- .RGtkCall("S_g_srv_target_get_port", object, PACKAGE = "RGtk2")

  return(w)
} 


gSrvTargetGetPriority <-
function(object)
{
  checkPtrType(object, "GSrvTarget")

  w <- .RGtkCall("S_g_srv_target_get_priority", object, PACKAGE = "RGtk2")

  return(w)
} 


gSrvTargetGetWeight <-
function(object)
{
  checkPtrType(object, "GSrvTarget")

  w <- .RGtkCall("S_g_srv_target_get_weight", object, PACKAGE = "RGtk2")

  return(w)
} 


gSrvTargetListSort <-
function(targets)
{
  targets <- as.GList(targets)

  w <- .RGtkCall("S_g_srv_target_list_sort", targets, PACKAGE = "RGtk2")

  return(w)
} 


gThreadedSocketServiceGetType <-
function()
{
  

  w <- .RGtkCall("S_g_threaded_socket_service_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gThreadedSocketServiceNew <-
function(max.threads)
{
  max.threads <- as.integer(max.threads)

  w <- .RGtkCall("S_g_threaded_socket_service_new", max.threads, PACKAGE = "RGtk2")

  return(w)
} 


gVolumeEjectWithOperation <-
function(object, flags, mount.operation, cancellable = NULL, callback, user.data = NULL)
{
  checkPtrType(object, "GVolume")
  
  checkPtrType(mount.operation, "GMountOperation")
  if (!is.null( cancellable )) checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  

  w <- .RGtkCall("S_g_volume_eject_with_operation", object, flags, mount.operation, cancellable, callback, user.data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gVolumeEjectWithOperationFinish <-
function(object, result, .errwarn = TRUE)
{
  checkPtrType(object, "GVolume")
  checkPtrType(result, "GAsyncResult")

  w <- .RGtkCall("S_g_volume_eject_with_operation_finish", object, result, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


gInetSocketAddressNew <-
function(address, port)
{
  checkPtrType(address, "GInetAddress")
  port <- as.integer(port)

  w <- .RGtkCall("S_g_inet_socket_address_new", address, port, PACKAGE = "RGtk2")

  return(w)
} 


gInetSocketAddressGetAddress <-
function(object)
{
  checkPtrType(object, "GInetSocketAddress")

  w <- .RGtkCall("S_g_inet_socket_address_get_address", object, PACKAGE = "RGtk2")

  return(w)
} 


gInetSocketAddressGetPort <-
function(object)
{
  checkPtrType(object, "GInetSocketAddress")

  w <- .RGtkCall("S_g_inet_socket_address_get_port", object, PACKAGE = "RGtk2")

  return(w)
} 


gTcpConnectionGetType <-
function()
{
  

  w <- .RGtkCall("S_g_tcp_connection_get_type", PACKAGE = "RGtk2")

  return(w)
} 


gTcpConnectionSetGracefulDisconnect <-
function(object, graceful.disconnect)
{
  checkPtrType(object, "GTcpConnection")
  graceful.disconnect <- as.logical(graceful.disconnect)

  w <- .RGtkCall("S_g_tcp_connection_set_graceful_disconnect", object, graceful.disconnect, PACKAGE = "RGtk2")

  return(invisible(w))
} 


gTcpConnectionGetGracefulDisconnect <-
function(object)
{
  checkPtrType(object, "GTcpConnection")

  w <- .RGtkCall("S_g_tcp_connection_get_graceful_disconnect", object, PACKAGE = "RGtk2")

  return(w)
} 

