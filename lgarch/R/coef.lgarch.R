coef.lgarch <-
function(object, arma=FALSE, ...)
{
  if(arma){
    result <- object$par.arma
  }else{
    result <- object$par
  }
  return(result)
}
