residuals.heckit5rob <-
function(object, ...)
{
  res=list()
  res$regime1=resid(object$stage21)
  res$regime2=resid(object$stage22)
  return(res)
}
