orgls <-
function(formula, data, constr, rhs, nec, weights=NULL, correlation=NULL, control=orlmcontrol()){
UseMethod("orgls")
}


orgls.formula <-
function(formula, data, constr, rhs, nec, weights=NULL, correlation=NULL, control=orlmcontrol()){
  cl <- match.call()
  if (!is.null(correlation)){
    groups <- getGroupsFormula(correlation)
  } else {
    groups <- NULL
  }
  glsSt <- glsStruct(corStruct = correlation, varStruct = varFunc(weights))
  model <- terms(formula, data = data)
  mfArgs <- list(formula = asOneFormula(formula(glsSt), formula, groups), data = data, na.action = na.fail)
  dataMod <- do.call("model.frame", mfArgs)
  origOrder <- row.names(dataMod)
  if (!is.null(groups)) {
    groups <- eval(parse(text = paste("~1", deparse(groups[[2]]), sep = "|")))
    grps <- getGroups(dataMod, groups, level = length(getGroupsFormula(groups, asList = TRUE)))
    ord <- order(grps)
    grps <- grps[ord]
    dataMod <- dataMod[ord, , drop = FALSE]
    revOrder <- match(origOrder, row.names(dataMod))
  } else {
    grps <- NULL
  }
  X <- model.frame(model, dataMod)
  contr <- lapply(X, function(el) if (inherits(el, "factor")) contrasts(el))
  contr <- contr[!unlist(lapply(contr, is.null))]
  x <- model.matrix(model, X)
  y <- eval(model[[2]], dataMod)

  if (is.numeric(constr)) constr <- rbind(constr)
  if (!is.matrix(constr)) stop("constr needs to be a matrix.")
  if (ncol(x) != ncol(constr)) stop(paste("constr has not correct dimensions.\nNumber of columns (",ncol(constr),") should equal the number of parameters: ", ncol(x), sep=""))
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
  unc <- gls(formula, data=dataMod, weights=weights, correlation=correlation, method="ML")

  ## extracting the variance-covariance structure
  if (is.null(unc$modelStruct$varStruct)){
    V <- diag(nrow(x))
  } else {
    V <- diag(attr(unc$modelStruct$varStruct, "weights"))
  }
  if (is.null(unc$modelStruct$corStruct)){
    crr <- diag(nrow(x))
  } else {
    cr <- corMatrix(unc$modelStruct$corStruct)
    crr <- if (is.matrix(cr)) cr else as.matrix(bdiag(cr))
  }
  W <- V %*% crr %*% V

  tBeta <- lm.gls(formula, data = dataMod, W=W)$coefficients
  
  # taken from lm.gls in package MASS
  # transforming X and y into a classical linear model framework
  eW <- eigen(W, TRUE)
  d <- eW$values
  if (any(d <= 0)) stop("'W' is not positive definite")
  eWv <- eW$vector
  A <- diag(sqrt(d)) %*% t(eWv)
  Ainv <- eWv %*% diag(1/sqrt(d))
  X <- A %*% x
  Y <- as.vector(A %*% y)  
  res <- Y - X %*% tBeta
  Sigma <- as.vector(t(res) %*% (res))/(nrow(x))  
  ############################
  ## lin model with order restrictions
  orsolve <- function(tBeta, X, Y, Constr, RHS, NEC){
    yVx <- t(X) %*% Y
    dvec <- 2*yVx
    Dmat <- 2*(t(X) %*% X)
    Amat <- t(Constr)
    bvec <- RHS
    solve.QP(Dmat,dvec,Amat,bvec=bvec, meq=NEC)    
  }
  orBeta <- tBeta
  val <- 0
  for (i in 1:control$maxiter){
    sqp <- orsolve(orBeta, X, Y, constr, rhs, nec)
    orBeta <- sqp$solution
    if (abs(sqp$value - val) <= control$absval) break else val <- sqp$value
  }
  if (i == control$maxiter & abs(sqp$value - val) > control$absval) warning("Maximum number of iterations reached without convergence!")

  ores <- (Y - X %*% orBeta)
  orSigma <- as.vector(t(ores) %*% (ores))/(nrow(x))
  Aores <- Ainv %*% (Y - X %*% orBeta)
  AorSigma <- as.vector(t(Aores) %*% diag(diag(W)) %*% (Aores))/(nrow(x))  
  p <- unc$dims$p
  N <- unc$dims$N
  Np <- N - p
  loglik <- (-N/2)*log(2*pi) + (-1/2)*(nrow(x)*log(AorSigma) + determinant(W, logarithm=TRUE)$modulus) - (1/2)*N + sum(log(diag(W)))

  names(orBeta) <- colnames(x)
  out <- list(call=cl, X=X, XW=x, y=Y, unccoefficients=tBeta, coefficients=orBeta, fitted=Ainv %*% (X %*% orBeta), residuals=Ainv %*% (Y - X %*% orBeta), sigma=Sigma, orSigma=orSigma, logLik=loglik, constr=constr, rhs=rhs, nec=nec, Niter=i, iact=sqp$iact, extrap=length(coef(unc[["modelStruct"]])), modelStruct=unc$modelStruct, W=W)
  class(out) <- "orgls"
  return(out)
}
