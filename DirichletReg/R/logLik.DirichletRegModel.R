logLik.DirichletRegModel <- function(object, ...){
  structure(object$logLik, df = sum(object$npar), class = "logLik")
}
