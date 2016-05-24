CheckPars <- function(Y,J=1,B="pc",lam1=1,lam2=1,thresh=10^(-4),
                      maxiter=100,maxiter.B=1,maxiter.T=1) {
  if (!all(is.matrix(Y),is.double(Y))) {
    stop("'Y' must be a numeric matrix")
  }
  if (!all(length(J)==1,J<=ncol(Y),J>=1)) {
    stop("'J' must be an integer <= 'ncol(Y)' and > 0")
  }
  if (is.matrix(B)) {
    if (!is.double(B)) {
      stop("'B' must be a numeric matrix")
    }
    if (!all(dim(B)==c(nrow(Y),J))) {
      stop("'B' is of the wrong dimension")
    }
  } else {
    B <- match.arg(B,c("pc","rand"))
  }
  if (!all(length(thresh)==1,thresh>0)) {
    stop("'thresh' must be a real number > 0")
  }
  if (!all(length(maxiter)==1,maxiter>0)) {
    stop("'maxiter' must be an integer > 0")
  }
  if (!all(length(maxiter.B)==1,maxiter.B>0)) {
    stop("'maxiter.B' must be an integer > 0")
  }
  if (!all(length(maxiter.T)==1,maxiter.T>0)) {
    stop("'maxiter.T' must be an integer > 0")
  }
  if (!all(length(lam1)==1,is.double(lam1),length(lam2)==1,
           is.double(lam2))) {
    stop("'lam1' and 'lam2' must be real numbers")
  }
}
