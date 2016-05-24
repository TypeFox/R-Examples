coef.heckit5rob <-
function(object, ...)
{
  coeff=list()
  coeff$S=coef(object$stage1)
  coeff$O1=coef(object$stage21)
  coeff$O2=coef(object$stage22)
  class(coeff)<- c("coef.heckit5rob", class(coeff), class(coeff$S))
  return(coeff)
}
