Moments <- function(x,p,data,conditional=TRUE,...) {

}

##' @export
`moments` <-
  function(x,...) UseMethod("moments")

##' @export
moments.lvmfit <- function(x, p=pars(x),...) moments(Model(x),p=p,...)

##' @export
moments.lvm.missing <- function(x, p=pars(x), ...) {
    idx <- match(coef(Model(x)),names(coef(x)))
    moments.lvmfit(x,p=p[idx],...)
}


##' @export
moments.lvm <- function(x, p, debug=FALSE, conditional=FALSE, data=NULL, ...) {
##  moments.lvm <- function(x, p, meanpar=NULL, conditional=FALSE, debug=FALSE,...) {
### p: model-parameters as obtained from e.g. 'startvalues'.
###       (vector of regression parameters and variance parameters)
### meanpar: mean-parameters (optional)

  ii <- index(x)
  pp <- modelPar(x,p)
  AP <- with(pp, matrices(x,p,meanpar=meanpar,epars=p2,data=data,...))
  P <- AP$P
  v <- AP$v
  if (!is.null(v)) {
    names(v) <- ii$vars
  }

  J <- ii$J
  Jidx <- ii$obs.idx
  if (conditional) {
    J <- ii$Jy
    px <- ii$px
    exo <- exogenous(x)
    if (!is.null(v))
      v[exo] <- as.numeric(data[1,exo])
    P <-  px%*% tcrossprod(P, px)
    Jidx <- ii$endo.idx
  }

  Im <- diag(nrow=nrow(AP$A))
  if (ii$sparse) {
    IAi <- with(AP, as(Inverse(Im-t(A)),"sparseMatrix"))
    ##IAi <- as(solve(Matrix::Diagonal(nrow(A))-t(A)),"sparseMatrix")
    G <- as(J%*%IAi,"sparseMatrix")
  } else {
    IAi <- Inverse(Im-t(AP$A))
    G <- J%*%IAi
    ##G <- IAi[Jidx,,drop=FALSE]
  }

  xi <- NULL
  if (!is.null(v)) {
    xi <- G%*%v ## Model-specific mean vector
  }
  Cfull <- as.matrix(IAi %*% tcrossprod(P,IAi))
  C <- as.matrix(J %*% tcrossprod(Cfull,J))

  return(list(Cfull=Cfull, C=C, v=v, e=AP$e, xi=xi, A=AP$A, P=P, IAi=IAi, J=J, G=G, npar=ii$npar, npar.reg=ii$npar.reg, npar.mean=ii$npar.mean, npar.ex=ii$npar.ex, parval=AP$parval, constrain.idx=AP$constrain.idx, constrainpar=AP$constrainpar))
}
