vcov.heckit5rob <-
function(object, ...)
{
  ret=list()
  ret$regime1=object$vcov1
  ret$regime2=object$vcov2
  return(ret)
}
