do.twinlm.strata <- function(x,fun,...) {
  res <- lapply(x$model,function(m) do.call(fun,c(list(m),list(...))))
  names(res) <- names(x$model)
  class(res) <- "do.twinlm.strata"
  newattr <- setdiff(names(attributes(x)),names(attributes(res)))
  attributes(res) <- c(attributes(res),attributes(x)[newattr])
  return(res)
}

##' @export
print.do.twinlm.strata <- function(x,...) {
  for (i in seq_len(length(x))) {    
    message(rep("-",60),sep="")
    message("Strata '",names(x)[i],"'",sep="")
    print(x[[i]])
  }
  if (!is.null(attributes(x)$time)) {
      message(rep("-",60),sep="")
      cat("\n")
      cat("Event of interest before time ", attributes(x)$time, "\n", sep="")      
  }
  return(invisible(x))
}


##' @export
plot.twinlm.strata <- function(x,...)
  suppressMessages(do.twinlm.strata(x,"plot",...))

##' @export
print.twinlm.strata <- function(x,...)
  print.do.twinlm.strata(x$model,...)

##' @export
summary.twinlm.strata <- function(object,...) 
  do.twinlm.strata(object,"summary",...)

##' @export
coef.twinlm.strata <- function(object,...) object$coef

##' @export
logLik.twinlm.strata <- function(object,indiv=FALSE,list=FALSE,...) {
  ll <- lapply(object$model,function(x) logLik(x,indiv=indiv,...))
  if (list) return(ll)
  if (!indiv) {
    res <- structure(sum(unlist(ll)),df=0,nall=0)
    for (i in seq(length(ll))) {
      attributes(res)$nall <- attributes(res)$nall+attributes(ll[[i]])$nall
      attributes(res)$df <- attributes(res)$df+attributes(ll[[i]])$df
    }
    ##  attributes(res)$nobs <- attributes(res)$nall-attributes(res)$df
    attributes(res)$nobs <- attributes(res)$nall
    class(res) <- "logLik"
    return(res)
  }
  return(unlist(ll))
}

##' @export
score.twinlm.strata <- function(x,...) {
  ss <- lapply(x$model,function(m) score(m,indiv=FALSE,...))
  return(unlist(ss))
}
