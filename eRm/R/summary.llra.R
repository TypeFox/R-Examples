summary.llra <- function(object, ...) UseMethod("summary.llra")

summary.llra <- function(object, level=0.95, ...)
  #summary for class llra
  {
    modi <- object$model
    calli <- deparse(object$call)
    logli <- object$loglik
    iti <- object$iter
    pari <- object$npar
    cii <- confint(object, "eta", level=level)
    se.eta <- object$se.eta
    names(se.eta) <- names(object$etapar)
    res <- list(etapar=object$etapar,se.eta=se.eta,ci=cii,iter=iti,model=modi,call=calli,npar=pari,loglik=logli,refGroup=object$refGroup)
    class(res) <- "summary.llra"
    res
 }
