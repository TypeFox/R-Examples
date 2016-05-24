deviance.biglm<-function(object,...) object$qr$ss

formula.biglm <- function(x, ...)
{
   formula(x$terms)
}

deviance.bigglm<-function(object,...) object$deviance

family.bigglm<-function(object,...) object$family

AIC.biglm<-function(object,...,k=2) deviance(object) + k*(object$n-object$df.resid)
AIC.bigglm<-function(object,...,k=2) deviance(object) + k*(object$n-object$df.resid)


