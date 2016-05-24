fitted.heckit5rob <-
function(object, ...)
{
  res=list()
  res$regime1=fitted(object$stage21)
  res$regime2=fitted(object$stage22)
  return(res)
}
