orlm <-
function(formula, data, constr, rhs, nec, control=orlmcontrol()){
UseMethod("orlm")
}

orlm.formula <-
function(formula, data, constr, rhs, nec, control=orlmcontrol()){
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- cbind(model.response(mf, "numeric"))
  x <- model.matrix(mt, mf, contrasts)

  if (is.numeric(constr)) constr <- rbind(constr)
  if (!is.matrix(constr)) stop("constr needs to be a matrix.")
  if (ncol(y)*ncol(x) != ncol(constr)) stop(paste("constr has not correct dimensions.\nNumber of columns (",ncol(constr),") should equal the number of responses times the number of parameters: ",ncol(y), "*", ncol(x), "=",ncol(y)*ncol(x), sep=""))
  if (length(rhs) != nrow(constr)) stop(paste("rhs has a different number of elements than there are numbers of rows in constr (",length(rhs), " != ", nrow(constr), ")", sep=""))
  if (is.numeric(nec) & length(nec) != 1) stop("nec needs to be single a numeric value or a logical vector with the same length as the number of constraints.")
  if (is.logical(nec) & length(nec) != length(rhs)) stop("nec needs to be single a numeric value or a logical vector with the same length as the number of constraints.")
  if (is.logical(nec)){
    ord <- order(nec, decreasing=TRUE)
    constr <- constr[ord,,drop=FALSE]
    rhs <- rhs[ord]
    nec <- sum(nec)
  }
  if (nec < 0) stop("nec needs to be positive")
  if (nec > length(rhs)) stop(paste("nec is larger than the number of constraints. (",nec," > ",length(rhs),")", sep=""))    

  ########################
  ## unconstrained linear model
  unc <- lm(y ~ x-1)

  tBeta <- as.vector(coefficients(unc))
  Sigma <- (t(residuals(unc)) %*% residuals(unc))/nrow(x)
  detU <- 1

  ############################
  ## lin model with order restrictions

  orsolve <- function(tBeta, x, y, Constr, RHS, NEC){
    Sigma <- (t(y - x %*% matrix(tBeta, ncol=ncol(y))) %*% (y - x %*% matrix(tBeta, ncol=ncol(y))))/nrow(x)
    yVx <- kronecker(solve(Sigma), t(x)) %*% as.vector(y)
    dvec <- 2*yVx
    Dmat <- 2*kronecker(solve(Sigma), t(x) %*% x)
    Amat <- t(Constr)
    bvec <- RHS
    solve.QP(Dmat,dvec,Amat,bvec=bvec, meq=NEC)    
  }

  orBeta <- tBeta
  val <- 0
  for (i in 1:control$maxiter){
    sqp <- orsolve(orBeta, x, y, constr, rhs, nec)
    orBeta <- sqp$solution
    if (abs(sqp$value - val) <= control$absval) break else val <- sqp$value
  }
  if (i == control$maxiter & abs(sqp$value - val) > control$absval) warning("Maximum number of iterations reached without convergence!")
  orSigma <- (t(y - x %*% matrix(orBeta, ncol=ncol(y))) %*% (y - x %*% matrix(orBeta, ncol=ncol(y))))/nrow(x)

  N <- length(as.vector(y))
  loglik <- (-N/2.0)*log(2*pi) + (-1/2.0)*(nrow(x)*log(det(orSigma)) + ncol(y)*log(detU)) - (1/2.0)*N

  if (ncol(y) > 1){
    orBeta <- matrix(orBeta, ncol=ncol(y), dimnames=list(colnames(x), colnames(y)))
    tBeta <- matrix(tBeta, ncol=ncol(y), dimnames=list(colnames(x), colnames(y)))
  }
  names(orBeta) <- colnames(x)
  out <- list(call=cl, X=x, y=y, unccoefficients=tBeta, coefficients=orBeta, fitted=x %*% orBeta, residuals=y - x %*% orBeta, sigma=Sigma, orSigma=orSigma, logLik=loglik, constr=constr, rhs=rhs, nec=nec, Niter=i, iact=sqp$iact)
  class(out) <- "orlm"
  return(out)
}

