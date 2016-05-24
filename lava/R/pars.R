##' @export
`pars` <-
  function(x,...) UseMethod("pars")

##' @export
pars.default <- function(x,...) {
  if (!is.null(x$opt$estimate))
    return(x$opt$estimate)
  if (!is.null(x$coef))
    return(x$coef)
  return(coef(x))
}

##' @export
pars.lvm.missing <- function(x,reorder=FALSE,...) {
    res <- pars.default(x)
    if (reorder) {
        idx <- match(coef(Model(x)),names(coef(x)))
        return(res[idx])
    }
    return(res)
}



###{{{ pars.multigroupfit
## pars.multigroupfit <- function(x,...) {
##   res <- pars.default(x)
##   lapply(ee$model$lvm,coef))
##   coef()
##}
###}}}

###{{{ pars.lvm

##' @export
pars.lvm <- function(x, A, P, v, e, ...) {
  parres <- A[index(x)$M1==1]
  diagcorfree <- diag(P)[diag(index(x)$P1)==1]
  parres <- c(parres, diagcorfree)

  if (ncol(A)>1)
  for (i in seq_len(ncol(index(x)$P1)-1))
    for (j in seq(i+1,nrow(index(x)$P1))) {
      if (index(x)$P1[j,i]!=0) {
        parres <- c(parres, P[j,i])
      }
    }
  if (length(parres)>0)
  names(parres) <- paste0("p",seq_len(length(parres)))
  if (!missing(v)) {
    parres <- c( v[which(index(x)$v1==1)], parres)
  }
  if (!missing(e)) {
    parres <- c( parres, e[which(index(x)$e1==1)] )
  }
  return(parres)
}

###}}} pars.lvm
