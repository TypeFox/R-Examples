fanc <- function(
  x, factors, n.obs, cor.factor=FALSE, normalize=TRUE, rho.max, covmat,
  control=list())
{
  ## Organize control parameters
  con <- list(
    tol.em=1e-5, tol.cd=1e-5, tol.bfgs=1e-5, min.uniquevar=0.005,
    eta=0, zita=0, Delta=1e-5, init.coef=seq(0.3, 2, length=10),
    rho.max=ifelse(missing(rho.max), NA, rho.max),
    max.gamma=100, min.gamma=1.01,
    length.rho=30, length.gamma=9,
    maxit.em=10000, maxit.cd=500, maxit.bfgs=500,
    cor.factor=cor.factor, min.rhozero=FALSE,
    start="warm", ncand.initial=10, pmax_for_S=500, 
    maxit.initial=500,  
    progress=FALSE, openmp=FALSE, num.threads=0, gamma.ebic=1)
  con[names(control)] <- control
    
  ## Check error
  if (!missing(x)) {
    if (!is.matrix(x))
      stop('"x" must be a matrix.')
    if (!is.numeric(x))
      stop('"x" must be a matrix.')
    if (factors < 1 || factors >= ncol(x))
      stop('"factors" must be a positive integer less than the number of variables.')
  }
  if (!missing(covmat)) {
    if (!is.matrix(covmat))
      stop('"covmat" must be a covariance matrix.')
    if (factors < 1 || factors >= ncol(covmat))
      stop('"factors" must be a positive integer less than the number of variables.')
  }
  if (con$length.gamma < 1)
    stop('"length.gamma" in control must be a positive integer.')
  if (!missing(n.obs)) {
    if (n.obs < 2)
      stop('"n.obs" must be an integer greater than 2.')
  }
  if (!is.logical(cor.factor))
    stop('"cor.factor"  must be logical.')
  if (con$length.rho < 1)
    stop('"length.rho" in control must be a positive integer.')
  if (!missing(rho.max)) {
    if (rho.max <= 0)
      stop('"rho.max"  must be a positive real value.')
    if (rho.max > 8.25)
      warning('"rho.max" is greater than 8.25. In such cases, the reparametrization of the penalty funcion may be failed')
  }
  if (con$min.gamma <= 1)
    stop('"min.gamma" in control must be a real value greater than 1.')
  if (con$max.gamma <= con$min.gamma)
    stop('"max.gamma"  must be a real value greater than min.gamma.')
  if (con$eta < 0)
    stop('"eta"  must be a non-negative real vaule.')
  if (con$ncand.initial < 1)
    stop('"ncand.initial" in control must be a positive integer.')
  if (con$maxit.em < 1)
    stop('"maxit.em" in control must be a positive integer.')
  if (con$maxit.cd < 1)
    stop('"maxit.cd" in control must be a positive integer.')
  if (con$maxit.bfgs < 1)
    stop('"maxit.bfgs" in control must be a positive integer.')
  if (is.na(match(con$start, c("warm", "cold"))))
    stop('"start" in control must be "warm" or "cold".')
  if (con$Delta <= 0)
    stop('"Delta" in control must be a positive real value.')
  if (con$min.uniquevar <= 0)
    stop('"min.uniquevar" in control must be a positive real value.')
  if (con$tol.em <= 0)
    stop('"tol.em" in control must be a positive real value.')
  if (con$tol.cd <= 0)
    stop('"tol.cd" in control must be a positive real value.')
  if (con$tol.bfgs <= 0)
    stop('"tol.bfgs" in control must be a positive real value.')
  if (con$zita < 0)
    stop('"zita" in control must be a non-negative real vaule.')
  if (!is.logical(normalize))
    stop('"normalize" must be logical.')
  if (!is.logical(con$min.rhozero))
    stop('"min.rhozero" in control must be logical.')
  if (!is.logical(con$progress))
    stop('"progress" in control must be logical.')
  if (!is.logical(con$openmp))
    stop('"openmp" in control must be logical.')
  if (con$openmp && con$num.threads == 0) {
    con$num.threads <- parallel::detectCores()
  }

  con$tol.em <- as.double(con$tol.em)
  con$tol.cd <- as.double(con$tol.cd)
  con$tol.bfgs <- as.double(con$tol.bfgs)
  con$min.uniquevar <- as.double(con$min.uniquevar)
  con$eta <- as.double(con$eta)
  con$zita <- as.double(con$zita)
  con$Delta <- as.double(con$Delta)
  con$init.coef <- as.double(con$init.coef)
  con$rho.max <- as.double(con$rho.max)
  con$max.gamma <- as.double(con$max.gamma)
  con$min.gamma <- as.double(con$min.gamma)
  con$length.rho <- as.integer(con$length.rho)
  con$length.gamma <- as.integer(con$length.gamma)
  con$maxit.em <- as.integer(con$maxit.em)
  con$maxit.cd <- as.integer(con$maxit.cd)
  con$maxit.bfgs <- as.integer(con$maxit.bfgs)
  con$cor.factor <- as.integer(con$cor.factor)
  con$min.rhozero <- as.integer(con$min.rhozero)
  con$start <- as.character(con$start)
  con$ncand.initial <- as.integer(con$ncand.initial)
  con$pmax_for_S <- as.integer(con$pmax_for_S)
  #con$trace <- as.logical(con$trace)
  con$maxit.initial <- as.integer(con$maxit.initial)
  con$progress <- as.integer(con$progress)
  con$openmp <- as.integer(con$openmp)
  con$num.threads <- as.integer(con$num.threads)

  ## Prepare data
  missing.x <- missing(x)
  if (missing.x) {
    if (missing(covmat))
      stop("input data matrix or covariance matrix is needed.")
    p <- ncol(covmat)
    N <- ifelse(missing(n.obs), p + 1, n.obs)
    x <- 1
  } else {
    x.orig <- x
    p <- ncol(x)
    N <- nrow(x)
    x <- scale(x, scale=normalize)[1:N, 1:p]
    if (normalize) x <- x / sqrt(N - 1) * sqrt(N)
    if (p <= N || p <= con$pmax_for_S) {
      covmat <- cov(x) * (N - 1) / N
    } else {
      covmat <- diag(apply(x, 2, function(x) var(x) * (N - 1) / N))
    }
    x <- x[1 : N, 1 : p] / sqrt(N)
  }

  ## Run
  rslt <- .Call("RextCall_fanc",
                as.integer(p), as.integer(factors), as.integer(N),
                as.double(covmat), as.double(x), con)

  ## Convert Lambda to dgCMatrix
  rslt$loadings <- vector("list", con$length.gamma)
  for (gi in 1 : con$length.gamma) {
    rslt$loadings[[gi]] <- vector("list", con$length.rho)
    for (ri in 1 : con$length.rho) {
      arridx <- arrayInd(rslt$spi.loadings[[ri, gi]] + 1, .dim=c(p, factors))
      rslt$loadings[[gi]][[ri]] <-
        sparseMatrix(i=arridx[, 1], j=arridx[, 2],
                     x=rslt$spv.loadings[[ri, gi]], dims=c(p, factors))
    }
  }
  rslt$spi.loadings <- NULL
  rslt$spv.loadings <- NULL
  
  ## Convert fanc format (pivot diag.Psi and logF table, aggregate conv)
  rslt$uniquenesses <- array(apply(rslt$uniquenesses, 3, t),
                             dim=dim(rslt$uniquenesses)[c(2, 1, 3)])
  rslt$likelihood <- array(apply(rslt$likelihood, 3, t),
                           dim=dim(rslt$likelihood)[c(2, 1, 3)])
  rslt$convergence <- rowSums(matrix(rslt$convergence, 3))

  ## Attach rownames and colnames
  names.gamma <- paste("gamma", 1 : con$length.gamma, sep="")
  names.rho <- paste("rho", 1 : con$length.rho, sep="")
  names.var <-
    if (missing.x) {
      if (!is.null(colnames(covmat))) colnames(covmat)
      else paste("V", 1 : p, sep="")
    } else {
      if (!is.null(colnames(x))) colnames(x)
      else paste("V", 1 : p, sep="")
    }
  names.factor <- paste("Factor", 1 : factors, sep="")
  names(rslt$loadings) <- names.gamma
  for (gi in 1 : con$length.gamma) {
    names(rslt$loadings[[gi]]) <- names.rho
    for (ri in 1 : con$length.rho) {
      dimnames(rslt$loadings[[gi]][[ri]]) <- list(names.var, names.factor)
    }
  }
  dimnames(rslt$uniquenesses) <- list(names.rho, names.var, names.gamma)
  dimnames(rslt$Phi) <- list(names.factor, names.factor, names.rho, names.gamma)
  dimnames(rslt$likelihood) <-
    list(names.rho, c("logF", "penalty", "logF+penalty"), names.gamma)
  names(rslt$convergence) <- c("EM", "Coordinate descent", "BFGS")

  ## Attach other properties
  if (missing.x==FALSE) rslt$x <- x.orig
  rslt$call <- match.call()

  class(rslt) <- "fanc"
  rslt
}
