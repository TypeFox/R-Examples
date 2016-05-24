Auc <- function(object,...){
  UseMethod("Auc",object=object)
}

Auc.default <- function(object,Spec){
  Sens <- object
  N <- length(Spec)
  if (length(unique(Spec)) < 2) return(NA)
  if (Spec[1]<Spec[N]) {Spec <- rev(Spec); Sens <- rev(Sens)}
  0.5 * sum(diff(c(0,1-Spec,1)) * (c(Sens,1) + c(0,Sens)))
}

Auc.Roc <- function(object,
                    digits=3,
                    print=TRUE){
  x <- object
  res <- c("Roc","AppRoc","BootcvRoc","NoInfRoc")
  names(res) <- c(x$method$name,c("App","BCV","NoInf"))
  found <- match(names(x),res,nomatch=FALSE)
  res <- res[found]
  out <- lapply(res,function(r){
    out <- sapply(1:length(x$models),function(w){
      Auc.default(object=x[[r]][[w]]$Sensitivity,
                  Spec=x[[r]][[w]]$Specificity)
    })
    names(out) <- names(x$models)
    out
  })
  if (print==TRUE){
    outMat <- do.call("cbind",out)
    outMat <- outMat[order(-out[[x$method$name]]),]
    print(outMat)
  }
  invisible(out)
}

