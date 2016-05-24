grpreg <- function(X, y, group=1:ncol(X), penalty=c("grLasso", "grMCP", "grSCAD", "gel", "cMCP", "gBridge", "gLasso", "gMCP"), family=c("gaussian","binomial", "poisson"), nlambda=100, lambda, lambda.min={if (nrow(X) > ncol(X)) 1e-4 else .05}, alpha=1, eps=.001, max.iter=1000, dfmax=p, gmax=J, gamma=3, tau=1/3, group.multiplier={if (strtrim(penalty,2)=="gr") sqrt(table(group[group!=0])) else rep(1,J)}, warn=TRUE, ...) {
  ## Error checking
  if (class(X) != "matrix") {
    tmp <- try(X <- as.matrix(X), silent=TRUE)
    if (class(tmp)[1] == "try-error") stop("X must be a matrix or able to be coerced to a matrix")
  }
  if (storage.mode(X)=="integer") X <- 1.0*X
  family <- match.arg(family)
  penalty <- match.arg(penalty)
  if (penalty=="gLasso") penalty <- "grLasso"
  if (penalty=="gMCP") penalty <- "cMCP"
  if (penalty=="gBridge") stop("gBridge has been divorced from the grpreg function; use the gBridge() function instead")
  if (alpha > 1 | alpha <= 0) stop("alpha must be in (0,1]")
  if (length(group)!=ncol(X)) stop("group does not match X")
  
  ## Reorder groups, if necessary
  xnames <- if (is.null(colnames(X))) paste("V",1:ncol(X),sep="") else colnames(X)
  if (any(order(group) != 1:length(group)) | !is.numeric(group)) {
    reorder.groups <- TRUE
    gf <- as.factor(group)
    if (any(levels(gf)=="0")) {
      gf <- relevel(gf, "0")
      g <- as.numeric(gf) - 1
      J <- max(g)
      tryCatch(names(group.multiplier) <- setdiff(levels(gf), "0"), finally="Length of group.multiplier must equal number of penalized groups")
    } else {
      g <- as.numeric(gf)
      J <- max(g)
      tryCatch(names(group.multiplier) <- levels(gf), finally="Length of group.multiplier must equal number of penalized groups")
    }
    g.ord <- order(g)
    g.ord.inv <- match(1:length(g), g.ord)
    g <- g[g.ord]
    X <- X[,g.ord]
  } else {
    reorder.groups <- FALSE
    g <- group
    J <- max(g)
    if (length(group.multiplier)!=max(g)) stop("Length of group.multiplier must equal number of penalized groups")
  }
  
  ## Set up XX, yy, lambda
  multi <- FALSE
  if (is.matrix(y) && ncol(y) > 1) {
    multi <- TRUE
    m <- ncol(y)
    response.names <- if (is.null(colnames(y))) paste("Y",1:m,sep="") else colnames(y)
    y <- multiY(y)
    X <- multiX(X, m)
    group <- g <- c(rep(0, m-1), rep(g, each=m))
    group.multiplier <- rep(1,J) ## Think about changing in the future
  }
  std <- .Call("standardize", X)
  XX <- std[[1]]
  center <- std[[2]]
  scale <- std[[3]]
  nz <- which(scale > 1e-6)
  zg <- setdiff(unique(g), unique(g[nz]))
  if (length(zg)) {
    J  <- J - length(zg)
    group.multiplier <- group.multiplier[-zg]
  }
  XX <- XX[ ,nz, drop=FALSE]
  g <- g[nz]
  if (strtrim(penalty,2)=="gr") {
    XX <- orthogonalize(XX, g)
    g <- attr(XX, "group")
  }
  K <- as.numeric(table(g))
  yy <- as.numeric(if (family=="gaussian") y - mean(y) else y)
  if (nrow(XX) != length(yy)) stop("X and y do not have the same number of observations")
  if (missing(lambda)) {
    lambda <- setupLambda(XX, yy, g, family, penalty, alpha, lambda.min, nlambda, group.multiplier)
    user.lambda <- FALSE
  } else {
    nlambda <- length(lambda)
    user.lambda <- TRUE
  }

  ## Fit
  n <- length(yy)
  p <- ncol(XX)
  K0 <- as.integer(if (min(g)==0) K[1] else 0)
  K1 <- as.integer(if (min(g)==0) cumsum(K) else c(0, cumsum(K)))
  if (K0) {
    lambda[1] <- lambda[1] + 1e-5
    user.lambda <- TRUE
  }
  if (family=="gaussian") {
    if (strtrim(penalty,2)=="gr") fit <- .Call("gdfit_gaussian", XX, yy, penalty, K1, K0, lambda, alpha, eps, as.integer(max.iter), gamma, as.double(group.multiplier), as.integer(dfmax), as.integer(gmax), as.integer(user.lambda))
    else fit <- .Call("lcdfit_gaussian", XX, yy, penalty, K1, K0, lambda, alpha, eps, 0, gamma, tau, as.integer(max.iter), as.double(group.multiplier), as.integer(dfmax), as.integer(gmax), as.integer(user.lambda))
    b <- rbind(mean(y), matrix(fit[[1]], nrow=p))
    iter <- fit[[2]]
    df <- fit[[3]] + 1 ## Intercept
    loss <- fit[[4]]
  }
  if (family=="binomial") {
    if (strtrim(penalty,2)=="gr") fit <- .Call("gdfit_binomial", XX, yy, penalty, K1, K0, lambda, alpha, eps, as.integer(max.iter), gamma, as.double(group.multiplier), as.integer(dfmax), as.integer(gmax), as.integer(warn), as.integer(user.lambda))
    else fit <- .Call("lcdfit_binomial", XX, yy, penalty, K1, K0, lambda, alpha, eps, 0, gamma, tau, as.integer(max.iter), as.double(group.multiplier), as.integer(dfmax), as.integer(gmax), as.integer(warn), as.integer(user.lambda))
    b <- rbind(fit[[1]], matrix(fit[[2]], nrow=p))
    iter <- fit[[3]]
    df <- fit[[4]]
    loss <- fit[[5]]
  }
  if (family=="poisson") {
    if (strtrim(penalty,2)=="gr") fit <- .Call("gdfit_poisson", XX, yy, penalty, K1, K0, lambda, alpha, eps, as.integer(max.iter), gamma, as.double(group.multiplier), as.integer(dfmax), as.integer(gmax), as.integer(warn), as.integer(user.lambda))
    else fit <- .Call("lcdfit_poisson", XX, yy, penalty, K1, K0, lambda, alpha, eps, 0, gamma, tau, as.integer(max.iter), as.double(group.multiplier), as.integer(dfmax), as.integer(gmax), as.integer(warn), as.integer(user.lambda))
    b <- rbind(fit[[1]], matrix(fit[[2]], nrow=p))
    iter <- fit[[3]]
    df <- fit[[4]]
    loss <- fit[[5]]
  }
  
  ## Eliminate saturated lambda values, if any
  ind <- !is.na(iter)
  b <- b[, ind, drop=FALSE]
  iter <- iter[ind]
  lambda <- lambda[ind]
  df <- df[ind]
  loss <- loss[ind]
  if (warn & any(iter==max.iter)) warning("Algorithm failed to converge for all values of lambda")

  ## Unstandardize
  if (strtrim(penalty,2)=="gr") b <- unorthogonalize(b, XX, g)
  b <- unstandardize(b, center[nz], scale[nz])
  beta <- matrix(0, nrow=(ncol(X)+1), ncol=length(lambda))
  beta[1,] <- b[1,]
  if (reorder.groups) {
    beta[nz+1,] <- b[-1,]
    beta[-1,] <- beta[1+g.ord.inv,]
  } else {
    beta[nz+1,] <- b[-1,]
  }
  
  ## Names
  varnames <- c("(Intercept)", xnames)
  if (multi) {
    beta[2:m,] <- sweep(beta[2:m,,drop=FALSE], 2, beta[1,], FUN="+")
    beta <- array(beta, dim=c(m, nrow(beta)/m, ncol(beta)))
    group <- group[-(1:(m-1))]
    dimnames(beta) <- list(response.names, varnames, round(lambda,digits=4))
  } else {
    dimnames(beta) <- list(varnames, round(lambda,digits=4))
  }
  
  val <- structure(list(beta=beta,
                        family=family,
                        group=group,
                        lambda=lambda,
                        alpha=alpha,
                        loss = loss,
                        n = n,
                        penalty=penalty,
                        df=df,
                        iter=iter,
                        group.multiplier=group.multiplier),
                   class = "grpreg")
  if (family=="poisson") val$y <- y
  val
}
