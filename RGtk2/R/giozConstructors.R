gAppLaunchContext <- gAppLaunchContextNew

gCancellable <- gCancellableNew

gFilenameCompleter <- gFilenameCompleterNew

gFileInfo <- gFileInfoNew

gBufferedInputStream <-
function(base.stream, size)
{
  if (!missing(size)) {
    gBufferedInputStreamNewSized(base.stream, size)
  }
  else {
    gBufferedInputStreamNew(base.stream)
  }
}

gDataInputStream <- gDataInputStreamNew

gMemoryInputStream <-
function(data)
{
  if (!missing(data)) {
    gMemoryInputStreamNewFromData(data)
  }
  else {
    gMemoryInputStreamNew()
  }
}

gMountOperation <- gMountOperationNew

gMemoryOutputStream <- gMemoryOutputStreamNew

gBufferedOutputStream <-
function(base.stream, size)
{
  if (!missing(size)) {
    gBufferedOutputStreamNewSized(base.stream, size)
  }
  else {
    gBufferedOutputStreamNew(base.stream)
  }
}

gDataOutputStream <- gDataOutputStreamNew

gSimpleAsyncResult <-
function(source.object, callback, user.data = NULL, source.tag, domain, code, format, ...)
{
  if (!missing(source.tag)) {
    gSimpleAsyncResultNew(source.object, callback, user.data, source.tag)
  }
  else {
    if (!missing(domain)) {
      gSimpleAsyncResultNewError(source.object, callback, user.data, domain, code, format, ...)
    }
    else {
      gSimpleAsyncResultNewFromError(source.object, callback, user.data)
    }
  }
}

gFileIcon <- gFileIconNew

gThemedIcon <-
function(iconname, iconnames, len)
{
  if (!missing(iconname)) {
    gThemedIconNew(iconname)
  }
  else {
    gThemedIconNewFromNames(iconnames, len)
  }
}

gEmblem <-
function(icon, origin)
{
  gEmblemNew(icon, origin)
}

gEmblemedIcon <- gEmblemedIconNew

gNetworkAddress <- gNetworkAddressNew

gNetworkService <- gNetworkServiceNew

gSocket <- gSocketNew

gSocketClient <- gSocketClientNew

gSocketListener <- gSocketListenerNew

gSocketService <- gSocketServiceNew

gThreadedSocketService <- gThreadedSocketServiceNew

gInetSocketAddress <- gInetSocketAddressNew

