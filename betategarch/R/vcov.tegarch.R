vcov.tegarch <-
function(object, ...)
{
if(is.null(object$hessian)){
  aux <- list()
  aux$asym <- object$model[2]
  aux$skew <- object$model[3]
  aux$iN <- length(object$y)
  aux$signnegy <- sign(-object$y)
  aux$u <- rep(0,aux$iN)

  #construct objective.f:
  if(object$model[1]==1){
    logl.penalty <- tegarchLogl(object$y, object$initial.values,
      lower=object$lower, upper=object$upper,
      lambda.initial=object$lambda.initial,
      logl.penalty=-1e+100, c.code=TRUE, aux=aux)
    objective.f <- function(pars, x=object$y){f <- -tegarchLogl(x,
      pars, lower=object$lower, upper=object$upper,
      lambda.initial=object$lambda.initial,
      logl.penalty=logl.penalty, c.code=TRUE, aux=aux); f}
  }else{
    logl.penalty <- tegarchLogl2(object$y, object$initial.values,
      lower=object$lower, upper=object$upper,
      lambda.initial=object$lambda.initial,
      logl.penalty=-1e+100, c.code=TRUE, aux=aux)
    objective.f <- function(pars, x=object$y){f <- -tegarchLogl2(x,
      pars, lower=object$lower, upper=object$upper,
      lambda.initial=object$lambda.initial,
      logl.penalty=logl.penalty, c.code=TRUE, aux=aux); f}
  }
  #compute Hessian vcovmat:
  hessian.numeric <- -optimHess(object$par, objective.f)
  vcovmat <- solve(-hessian.numeric)

}else{
   vcovmat <- solve(-object$hessian)
}
return(vcovmat)
}
