
logLik.iModel <- function(object,...)
  structure(object$fitinfo$logL, df=object$fitinfo$dimension["df"], class="logLik")

## Returns (df, AIC=-2logL + k df), so the objective is to mimimize this quantity
##
extractAIC.iModel <- function(fit, scale, k = 2, ...){
  c(fit$fitinfo$df, fit$fitinfo$dev - 2*fit$fitinfo$dimension["df"])
}

summary.iModel <- function(object, ...){
  glist <- object$glist

  isg   <- object$isGraphical
  isd   <- object$isDecomposable
  #cq    <- maxCliques(ugList(glist))$maxCliques
  cq    <- getCliques(ugList(glist))# $maxCliques
  ans   <- structure(list(glist=glist, isGraphical=isg, isDecomposable=isd, cliques=cq),
                     class="iModelsummary")
  ans
}

print.iModelsummary <- function(x,...){
  cat(sprintf("is graphical=%s; is decomposable=%s\n", x$isGraphical, x$isDecomposable))
  cat("generators (glist):\n")
  str(x$glist, give.head=FALSE, comp.str=" ", no.list=TRUE)
  #cat("EXPERIMENTAL: components: ", names(x),"\n")
  invisible(x)
}

.extractFIT <- function(object,...){
  c(object[[1]], object$df)
}

.glist2formula <- function (f) {
  if (inherits(f, "formula")) 
    return(f)
  ans <- try(as.formula(paste("~", paste(unlist(lapply(f, paste, collapse = "*")), 
                                         collapse = "+")), .GlobalEnv),silent=TRUE)
  if (inherits(ans, "try-error"))
    stop("Unable to create formula from list. \nCould be due to white space, strange characters etc. in variable names\n")
  ans
}

formula.iModel <- function(x,...){
	#list2rhsFormula(x$glist)
  .glist2formula(x$glist)
}

terms.iModel <- function(x,...){
	x$glist
}


