BC <-
function(Y, X = NULL, model = 'null', covlist = NULL, condition = NULL, its = 100, ...) {
  if (!is.matrix(Y)) {
    stop("Y must be a matrix")
  }
  if (!(is.null(X) | is.matrix(X))) {
    stop("X must be a matrix or null")
  }
  
  n <- nrow(Y)
  nsp <- ncol(Y)
  Xname <- colnames(X)
  
  if (model %in% c('null', 'environment', 'community', 'full')) {
    if (model == 'null') {
      X <- matrix(1, n, 1)
      updateR <- F
      covlist <- NULL
    } else if (model == 'environment') {
      X <- cbind(rep(1, n), X)
      updateR <- F
    } else if (model == 'community') {
      X <- matrix(1, n, 1)
      updateR <- T
      covlist <- NULL
    } else if (model == 'full') {
      X <- cbind(rep(1, n), X)
      updateR <- T
    }
  } else {
    stop ("model must be one of: 'null', 'environment', 'community' or  'full'")
  }
  
  k <- ncol(X)
  # matrix column names
  if (is.matrix(X) & is.null(Xname)) {
    colnames(X) <- paste("cov", 0:(k-1), sep = "")
  }
  colnames(X)[1] <- "intercept"
  if(is.null(colnames(Y))) {
    colnames(Y) <- paste("sp", 1:ncol(Y), sep = "")
  }
  
  # covariate list
  if (is.null(covlist)) {  # if null, add them all
    for (i in 1:nsp) {
      covlist[i] <- list(1:k)
    }
  } else {
    if (!is.list(covlist)) {
      stop("covlist must be a list or null")
      }
    # otherwise add the intercept to the list
    covlist <- lapply(covlist, function(x) c(1, x+1) )
  }
  if (max(unlist(lapply(covlist, length))) > k) {
    stop ("too many elements in covlist")
  }
  if (max(unlist(lapply(covlist, max))) > k) {
    stop ("covlist index out of range of covariates")
  }
  
  if (!is.null(condition)) {
    if(!is.matrix(condition)){
      stop("condition must be a matrix or null")
    }
    if (is.null(colnames(condition))){
      colnames(condition) <- paste("cond", 1:ncol(condition), sep = ".")
    }
    # add them to X and adjust the covariate list
    ind <-  (ncol(X)+1):(ncol(X) + ncol(condition))
    X <- cbind(X, condition)
    covlist <- lapply(covlist, function(x) c(x, ind))
  }

  # initialise mu at species prevalence
  mu <- qnorm(rowMeans(t(Y)))
  mu <- matrix(rep(mu, n), n, nsp, byrow = TRUE) 
  
  # initialise z - just on the right side of 0
  trunc <- find_trunc(mu, t(Y))
  e <- matrix(0, n, nsp)
  e <- ifelse(trunc[, , 1] == -Inf, trunc[, , 2] - 1, trunc[, , 1] + 1)
  e <- sample_e(e, trunc, diag(nsp))
  z <- mu + e
  
  # initialise R
  W <- cov(e)
  if (!updateR) W <- diag(nsp)
  R <- cov2cor(W)
  
  reslis <- BCfit(Y, X, covlist, R, z, mu, updateR, its, ...)

  res <- NULL
  res$trace <- list(R = reslis$R, B = reslis$B, z = reslis$z)
  res$call <- list(model = model, Y = Y, X = X, covlist = covlist, its = its,
                   start = reslis$burn + 1, thin = reslis$thin)
  res$other <- list(mu = reslis$mu)
  class(res) <- "bayescomm"
  res
}
