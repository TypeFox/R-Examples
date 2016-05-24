##' @export
`exogenous` <-
function(x,...) UseMethod("exogenous")

##' @export
"exogenous<-" <- function(x,...,value) UseMethod("exogenous<-")

##' @export
`exogenous<-.lvm` <- function(x,silent=FALSE,
                              xfree=TRUE,
                              ...,value) {
  if (inherits(value,"formula")) {
    exogenous(x,...) <- all.vars(value)
    return(x)
  }
  not.in <- !(value%in%vars(x))
  if (any(not.in)) {
    addvar(x,reindex=FALSE) <- value[not.in]
  }
  xorg <- exogenous(x)
  x$exogenous <- value
  if (!is.null(value) & xfree) {
    notexo.idx <- xorg[which(!(xorg%in%value))]
    if (length(notexo.idx)>0) { ##  & mom) {
      if (length(notexo.idx)>1) {
        covariance(x,notexo.idx,pairwise=TRUE,exo=TRUE) <- NA
      }
      covariance(x,notexo.idx,vars(x),exo=TRUE) <- NA
      intercept(x,notexo.idx) <- x$mean[notexo.idx]
    }
  }
##  x$exogenous <- value
  index(x) <- reindex(x)
  return(x)
}

##' @export
`exogenous.lvm` <-
function(x,latent=FALSE,index=TRUE,...) {
  if (!index) {
    if (latent) {
      allvars <- vars(x)
    } else {
      allvars <- manifest(x)
    }
    M <- x$M
    res <- c()
    for (i in allvars)
      if (!any(M[,i]==1) & any(M[i,]==1))
        res <- c(res, i)
    return(res)
  }
  if (is.null(x$exogenous)) return(x$exogenous)
  if (all(!is.na(x$exogenous)) & !latent) {
    return(x$exogenous[x$exogenous%in%index(x)$manifest])
  }
  if (!latent)
    return(index(x)$exogenous)
  return(exogenous(x,latent=latent,index=FALSE,...))
}

##' @export
`exogenous.lvmfit` <-
function(x,...) {
  exogenous(Model(x),...)
}

##' @export
exogenous.list <- function(x,...) {
  exolist <- c()
  endolist <- c()
  for (i in seq_along(x)) {
    exolist <- c(exolist, exogenous(x[[i]]))
    endolist <- c(endolist, endogenous(x[[i]]))
  }
  endolist <- unique(endolist)
  exolist <- unique(exolist)
  return(exolist[!(exolist%in%endolist)])
}

##' @export
`exogenous.multigroup` <-
function(x,...) {
  exogenous(Model(x))
}
