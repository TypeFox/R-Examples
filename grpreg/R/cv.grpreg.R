cv.grpreg <- function(X, y, group=1:ncol(X), ..., nfolds=10, seed, trace=FALSE) {
  if (!missing(seed)) set.seed(seed)
  fit <- grpreg(X=X, y=y, group=group, ...)
  multi <- FALSE
  if (is.matrix(y) && ncol(y) > 1) {
    multi <- TRUE
    m <- ncol(y)
    response.names <- if (is.null(colnames(y))) paste("Y",1:m,sep="") else colnames(y)
    y <- multiY(y)
    X <- multiX(X, m)
  }
  g <- if (multi) c(rep(0,m-1), fit$group) else fit$group
  E <- matrix(NA, nrow=length(y), ncol=length(fit$lambda))
  if (fit$family=="binomial") PE <- E

  n <- length(y)
  if (multi) {
    nn <- n/m
    cv.ind <- rep(ceiling(sample(1:nn)/nn*nfolds), each=m)
  } else if (fit$family=="binomial" & (min(table(y)) > nfolds)) {
    ind1 <- which(y==1)
    ind0 <- which(y==0)
    n1 <- length(ind1)
    n0 <- length(ind0)
    cv.ind1 <- ceiling(sample(1:n1)/n1*nfolds)
    cv.ind0 <- ceiling(sample(1:n0)/n0*nfolds)
    cv.ind <- numeric(n)
    cv.ind[y==1] <- cv.ind1
    cv.ind[y==0] <- cv.ind0
  } else {
    cv.ind <- ceiling(sample(1:n)/n*nfolds)
  }

  for (i in 1:nfolds) {
    if (trace) cat("Starting CV fold #",i,sep="","\n")
    X1 <- X[cv.ind!=i, , drop=FALSE]
    y1 <- y[cv.ind!=i]
    X2 <- X[cv.ind==i, , drop=FALSE]
    y2 <- y[cv.ind==i]

    args <- list(..., X=X1, y=y1, group=g)
    args$lambda <- fit$lambda
    args$warn <- FALSE
    fit.i <- do.call('grpreg', args)

    yhat <- matrix(predict(fit.i, X2, type="response"), length(y2))
    E[cv.ind==i, 1:length(fit.i$lambda)] <- loss.grpreg(y2, yhat, fit$family)
    if (fit$family=="binomial") PE[cv.ind==i, 1:length(fit.i$lambda)] <- (yhat < 0.5) == y2
  }

  ## Eliminate saturated lambda values, if any
  ind <- which(apply(is.finite(E), 2, all))
  E <- E[, ind, drop=FALSE]
  lambda <- fit$lambda[ind]

  ## Return
  cve <- apply(E, 2, mean)
  cvse <- apply(E, 2, sd) / sqrt(n)
  min <- which.min(cve)
  null.dev <- calcNullDev(X, y, group=g, family=fit$family)

  val <- list(cve=cve, cvse=cvse, lambda=lambda, fit=fit, min=min, lambda.min=lambda[min], null.dev=null.dev)
  if (fit$family=="binomial") val$pe <- apply(PE[,ind], 2, mean)
  structure(val, class="cv.grpreg")
}
