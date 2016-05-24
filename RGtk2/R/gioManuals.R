
## use sprintf() to handle var args
gSimpleAsyncResultNewError <-
  function(source.object, callback, user.data, domain, code, format, ...)
{
  checkPtrType(source.object, "GObject")
  callback <- as.function(callback)
  
  domain <- as.GQuark(domain)
  code <- as.integer(code)
  format <- sprintf(format, ...)

  w <- .RGtkCall("S_g_simple_async_result_new_error", source.object, callback,
                 user.data, domain, code, format, PACKAGE = "RGtk2")

  return(w)
}
gSimpleAsyncResultSetError <-
  function(object, domain, code, format, ...)
{
  checkPtrType(object, "GSimpleAsyncResult")
  domain <- as.GQuark(domain)
  code <- as.integer(code)
  format <- sprintf(format, ...)

  w <- .RGtkCall("S_g_simple_async_result_set_error", object, domain, code,
                 format, PACKAGE = "RGtk2")

  return(invisible(w))
}
gSimpleAsyncReportErrorInIdle <-
  function(object, callback, user.data, domain, code, format, ...)
{
  checkPtrType(object, "GObject")
  callback <- as.function(callback)
  
  domain <- as.GQuark(domain)
  code <- as.integer(code)
  format <- sprintf(format, ...)

  w <- .RGtkCall("S_g_simple_async_report_error_in_idle", object, callback,
                 user.data, domain, code, format, PACKAGE = "RGtk2")

  return(w)
}

## varargs, property settings
gAsyncInitableNewAsync <-
  function(object.type, io.priority, cancellable, callback, user.data, ...)
{
  object.type <- as.GType(object.type)
  io.priority <- as.integer(io.priority)
  checkPtrType(cancellable, "GCancellable")
  callback <- as.function(callback)
  
  props <- list(...)
  if (is.null(names(props)))
    stop("Property values in '...' must be named")

  w <- .RGtkCall("S_g_async_initable_new_async", object.type, io.priority,
                 cancellable, callback, user.data, props, PACKAGE = "RGtk2")

  return(w)
}

gInitableNew <-
  function(object.type, cancellable, ..., .errwarn = TRUE)
{
  object.type <- as.GType(object.type)
  checkPtrType(cancellable, "GCancellable")
  
  props <- list(...)
  if (is.null(names(props)))
    stop("Property values in '...' must be named")

  w <- .RGtkCall("S_g_initable_new", object.type, cancellable, props,
                 PACKAGE = "RGtk2")
  
  w <- handleError(w, .errwarn)

  return(w)
}


## Vector of file attribute constants

GFileAttributeStandard <- c(type = "standard::type",
                            isHidden = "standard::is-hidden",
                            isBackup = "standard::is-backup",
                            isSymlink = "standard::is-symlink",
                            isVirtual = "standard::is-virtual",
                            name = "standard::name",
                            displayName = "standard::display-name",
                            editName = "standard::edit-name",
                            copyName = "standard::copy-name",
                            icon = "standard::icon",
                            contentType = "standard::content-type",
                            fastContentType = "standard::fast-content-type",
                            size = "standard::size",
                            allocatedSize = "standard::allocated-size",
                            symlinkTarget = "standard::symlink-target",
                            targetUri = "standard::target-uri",
                            sortOrder = "standard::sort-order",
                            description = "standard::description")
GFileAttributeEtag <- c(value = "etag::value")
GFileAttributeId <- c(file = "id::file",
                      filesystem = "id::filesystem")
GFileAttributeAccess <- c(canRead = "access:can-read",
                          canWrite = "access:can-write",
                          canExecute = "access:can-execute",
                          canDelete = "access:can-delete",
                          canTrash = "access:can-trash",
                          canRename = "access:can-rename")
GFileAttributeMountable <- c(canMount = "mountable::can-mount",
                             canUnmount = "mountable::can-unmount",
                             canEject = "mountable::can-eject",
                             unixDevice = "mountable::unix-device",
                             unixDeviceFile = "mountable::unix-device-file",
                             halUdi = "mountable::hal-udi",
                             canPoll = "mountable::can-poll",
                             isMediaCheckAutomatic =
                             "mountable::is-media-check-automatic",
                             canStart = "mountable::can-start",
                             canStartDegraded = "mountable::can-start-degraded",
                             canStop = "mountable::can-stop",
                             startStopType = "mountable::start-stop-type")
GFileAttributeTime <- c(modified = "time::modified",
                        modifiedUsec = "time::modified-usec",
                        access = "time::access",
                        accessUsec = "time::access-usec",
                        changed = "time::changed",
                        changedUsec = "time::changed-usec",
                        created = "time::created",
                        createdUsec = "time::created-usec")
GFileAttributeUnix <- c(device = "unix::device",
                        inode = "unix::inode",
                        mode = "unix::mode",
                        nlink = "unix::nlink",
                        uid = "unix::uid",
                        gid = "unix::gid",
                        rdev = "unix::rdev",
                        blockSize = "unix::block-size",
                        blocks = "unix::blocks",
                        isMountpoint = "unix::is-mountpoint")
GFileAttributeDos <- c(isArchive = "dos::is-archive",
                       isSystem = "dos::is-system")
GFileAttributeOwner <- c(user = "owner::user",
                         userReal = "owner::user-real",
                         group = "owner::group")
GFileAttributeThumbnail <- c(path = "thumbnail::path",
                             failed = "thumbnail::failed")
GFileAttributePreview <- c(icon = "preview::icon")
GFileAttributeFilesystem <- c(size = "filesystem::size",
                              free = "filesystem::free",
                              type = "filesystem::type",
                              readonly = "filesystem::readonly",
                              usePreview = "filesystem::use-preview")
GFileAttributeGvfs <- c(backend = "gvfs::backend")
GFileAttributeTrash <- c(itemCount = "trash::item-count")

## Some convenient constants for IO priority
GPriority <- c(default = 0L, high = -100L, highIdle = 100L, defaultIdle = 200L,
               low = 300L)
                        
                             
