gFileAttributeInfoListGetInfos <-
function(obj)
{
  checkPtrType(obj, 'GFileAttributeInfoList')
  v <- .Call('S_GFileAttributeInfoListGetInfos', obj, PACKAGE = "RGtk2")
  v
} 
gFileAttributeInfoListGetNInfos <-
function(obj)
{
  checkPtrType(obj, 'GFileAttributeInfoList')
  v <- .Call('S_GFileAttributeInfoListGetNInfos', obj, PACKAGE = "RGtk2")
  v
} 
