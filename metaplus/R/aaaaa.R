setClass("mymle", contains="mle2")

setClass("summary.mymle", contains="summary.mle")

setClass("profile.mymle", contains="profile.mle")

setMethod("summary", "mymle", function(object, waldtest=TRUE, ...){
  cmat <- cbind(Estimate = object@coef,
                `Std. Error` = sqrt(diag(object@vcov)))
  zval <- cmat[,"Estimate"]/cmat[,"Std. Error"]
  pval <- 2*pnorm(-abs(zval))
  coefmat <- cbind(cmat,"z value"=zval,"Pr(z)"=pval)
  m2logL <- 2*object@min
  new("summary.mymle", call=object@call.orig, coef=coefmat, m2logL= m2logL)
})
