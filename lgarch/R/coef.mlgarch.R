coef.mlgarch <-
function(object, varma=FALSE, ...)
{
  if(varma){
    result <- object$par.varma
  }else{
    result <- object$par
  }
  return(result)
}
