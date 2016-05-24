#### parse turns out to be quite time consuming
#### so better parse once and assign the expression to the body of a function
#### this is what fitdistr does in package MASS


#### get MLE
findMLE <- function(x, densfun, start, ...) {
  if (missing(densfun) || !(is.function(densfun)) )
    stop("'densfun' must be supplied as a function")
  log.arg <- "log" %in% names(formals(densfun))
  n <- NROW(x)
  
  myfn <- function(parm, ...) -sum(log(dens(parm, ...)))
  mylogfn <- function(parm, ...) -sum(dens(parm, ..., log = TRUE))

  ## correct the order of start to match the argument order in densfun
  nm <- names(start)
  f <- formals(densfun)
  args <- names(f)
  m <- match(nm, args)
  if (any(is.na(m))) 
    stop("'start' specifies names which are not arguments to 'densfun'")
  start <- start[order(m)]
  if (length(nm) == 1L) start <- c(unlist(start)) ## added for multivariate


  ## assuming for multivariate distribution, the parameters are stored
  ## in one parameter vector (e.g., param in dmymvnorm)
  dens <- function(parm, x, ...) densfun(x, parm, ...)
  if ((l <- length(nm)) > 1L) 
    body(dens) <- parse(text = paste("densfun(x,", paste("parm[", 1L:l, "]", collapse = ", "), ", ...)"))
  
  Call <- match.call(expand.dots = TRUE)
  Call[[1L]] <- as.name("optim")
  Call$densfun <- Call$start <- Call$x <- NULL
  Call$fn <- if ("log" %in% args) mylogfn else myfn
  Call$par <- start
  Call$x <- x
  ## Call$control <- list(maxit=20000) ## should not be hard coded as this.

  ## Call$hessian <- TRUE
  if (is.null(Call$method)) {
    if (any(c("lower", "upper") %in% names(Call))) 
      Call$method <- "L-BFGS-B"
##    else if (length(start) > 1L) 
##      Call$method <- "BFGS"
##    else if (length(start) == 1) Call$method <- "Brent"
    else Call$method <- "Nelder-Mead" ## for dim > 1, this works better
  }

  res <- eval.parent(Call)
  if (res$convergence > 0L) 
    stop("optimization failed")
  ## sds <- sqrt(diag(solve(res$hessian)))
  list(estimate = res$par,
       ## sd = sds,
       loglik = -res$value, 
       n = n)
}

gradControl <- function(eps=1e-4, d=0.0001,
                        zero.tol=sqrt(.Machine$double.eps/7e-7),
                        r=6, v=2, show.details=FALSE) {
  list(eps=eps, d = d, zero.tol = zero.tol, r = r, v = v, show.details = show.details)
}

gof.test <- function(x, distr, nsim, start,
                     simulation = c("mult", "pb"),
                     ng = 1000,
                     qgrid = c("empirical", "fitted"), ## only matters for dim = 1
                     gridStat = NULL,  ## for CvM statistics, dim >= 2
                     method.args = gradControl(), ## for jacobian,
                     derCdfWrtPar = NULL,         ## function dF.dTheta(theta, x)
                     derLogPdfWrtPar = NULL,      ## function dlogf.dTheta(theta, x)
                     ...) {
#### x: the data vector to be tested
#### distr: character string for the name of the distribution
#### nsim: bootstrap sample size
#### start: a named list to be passed to fitdistr
  n <- NROW(x)
  dim <- NCOL(x)
  if (is.null(gridStat)) gridStat <- if (dim == 1) TRUE else FALSE

  rdistr <- get(paste("r", distr, sep=""), mode = "function")
  pdistr <- get(paste("p", distr, sep=""), mode = "function")
  ddistr <- get(paste("d", distr, sep=""), mode = "function")
  if (dim == 1 && qgrid == "fitted") ## only for dim = 1
    qdistr <- get(paste("q", distr, sep=""), mode = "function") 
  
  if (!is.matrix(x)) x <- as.matrix(x)
  ## estimation
  fit <- findMLE(x, ddistr, start, ...)
  
  thetahat <- fit$estimate
  p <- length(thetahat)
  log.arg <- "log" %in% names(formals(ddistr))
  
  ## creat some calls to be evaluated later;
  ## ddistr.call, rdistr.call and pdistr.call are for dim == 1 only
  rdistr.call <- parse(text = paste("rdistr(n, ",  paste("theta[", 1L:p, "]", collapse = ", "), ")"))
  pdistr.call <- parse(text = paste("pdistr(x, ",  paste("theta[", 1L:p, "]", collapse = ", "), ")"))
  ddistr.call <- if (log.arg)  parse(text = paste("ddistr(x, ", paste("theta[", 1L:p, "]", collapse = ", "), ", log = TRUE)")) else parse(text = paste("log( ddistr(x, ", paste("theta[", 1L:p, "]", collapse = ", "), "))"))
  qdistr.call <- parse(text = paste("qdistr(p, ",  paste("theta[", 1L:p, "]", collapse = ", "), ")"))
  
  findMLE.call <- if (dim == 1) parse(text = paste("findMLE(x, ddistr, start = as.list(theta), ...)")) else parse(text = paste("findMLE(x, ddistr, start = list(param=theta), ...)"))

  ## for dim > 1, set up pdistr and rdistr
  dots <- names(list(...))  ## e.g., dispstr = "ex"
  pf <- formals(pdistr)
  pf.args <- names(pf)
  pf.m <- match(pf.args, dots)
  pf.m <- match(pf.args, dots)
  pf.m <- pf.m[!is.na(pf.m)]
  pdistr.my <- function(x, parm) {
    if (dim == 1) {
      ret <- eval(pdistr.call, list(theta=parm, x = x))
    }
    else {
      formals(pdistr)[names(list(...)[pf.m])] <- list(...)[pf.m]
      ret <- pdistr(x, parm)
    }
    ret 
  }
  ddistr.my <- function(x, parm) {
    if (dim == 1) {
      ret <- eval(ddistr.call, list(theta=parm, x = x))
    }
    else {
      formals(ddistr)[names(list(...)[pf.m])] <- list(...)[pf.m]
      if (log.arg) ret <- ddistr(x, parm, log=TRUE)
      else ret <- log(ddistr(x, parm))
    }
    ret 
  }
  rdistr.my <- function(n, parm) {
    if (dim == 1) {
      ret <- eval(rdistr.call, list(theta=parm, n = n))
    }
    else {
      formals(rdistr)[names(list(...)[pf.m])] <- list(...)[pf.m]
      formals(rdistr)$dim <- dim
      ret <- rdistr(n, parm)
    }
    ret 
  }
  
  ## generate a grid for numerigal integration if gridStat == TRUE
  xg <- NULL
  if (gridStat) {
    ng1 <- ceiling(ng^(1 / dim))
    pgrid <- (1:ng1) / (ng1 + 1)
    if (dim == 1 && qgrid == "fitted") {
      xglist <- eval(qdistr.call, list(theta = thetahat, p = pgrid))
      xg <- as.matrix(expand.grid(xglist))
    }
    else {
      xglist <- list()
      for (i in 1:NCOL(x)) xglist[[i]] <- quantile(x[,i], prob = pgrid)
      xg <- as.matrix(expand.grid(xglist))
    }
  }
  
  ## fxg <- ddistr.my(xg, thetahat) ## not needed
  if (dim > 1) ecdf <- function(x) mecdf(x, expand = NA, continuous = FALSE)
  ## need to be careful with mecdf which addes two fake points
  Fn <- ecdf(x)
  Ftheta <- function (x) pdistr.my(x, thetahat)
  
  ## test statistics
  stat <- getStat(x, n, Fn, Ftheta)
  names(stat) <- c("ksa", "cvma")

  if (gridStat) {
    stat <- c(getStat(xg, n, Fn, Ftheta), stat)
    names(stat) <- c("ks", "cvm", "ksa", "cvma")
  }

  pb <- function() {
    x.sim <- rdistr.my(n, thetahat)
    theta.sim <- eval(findMLE.call, list(theta = thetahat, x = x.sim))$estimate
    Fn.sim <- ecdf(x.sim)
    Ftheta.sim <- function(x) pdistr.my(x, theta.sim)
    stat <- getStat(x.sim, n, Fn.sim, Ftheta.sim)
    if (gridStat) {
      stat <- c(getStat(xg, n, Fn.sim, Ftheta.sim), stat)
    }
    stat
  }
  logdens <- function(theta) ddistr.my(x, theta)  
  llk <- function(theta) sum(logdens(theta))

  if (simulation == "pb") stat.sim <- replicate(nsim, pb())
  else {
    if (is.null(derLogPdfWrtPar))
      grad <- jacobian(logdens, thetahat, method.args = method.args)
    else grad <- derLogPdfWrtPar(thetahat, x)
    ## a matrix of n by p
    ## hinv <- solve( - hessian(llk, theta) / n ) ## a matrix p by p
    ## if (fisher == "G") hinv <- solve(var(t(grad))) else hinv <- solve( - hessian(llk, theta) / n )
    ## hinv <- solve(var(grad))
    ## psi <- hinv %*% t(grad)
    psi <- solve(var(grad), t(grad))
    stat.sim <- mult.sim(xg, x, pdistr.my, ddistr.my, thetahat, psi, nsim, gridStat = gridStat, method.args = method.args, derCdfWrtPar = derCdfWrtPar)
  }

  p.value <- rowSums(stat.sim >= stat) / (nsim + 1)
  list(statistic = stat, p.value = p.value, estimate = thetahat, stat.sim = stat.sim)
}

influFun <- function(xg, x, psi, F, theta, method.args, derCdfWrtPar) {
  ## xg: the grid point at which the function is to be evaluated
  ## x: observed data, n by d
  ## psi: influ of mle, a matrix of p by n
  ## F: cdf of X
  xg <- as.matrix(xg)
  x <- as.matrix(x)
  n <- NROW(x)
  d <- NCOL(x)
  ng <- NROW(xg)

  Ftheta <- function(theta) F(xg, theta)
  if (is.null(derCdfWrtPar))
    dF.dtheta <- jacobian(Ftheta, theta, method.args = method.args)
  else dF.dtheta <- derCdfWrtPar(theta, x)
  
  idx1 <- rep(1:n, ng)
  idx2 <- rep(1:ng, each = n)
  mat1 <- x[idx1,,drop = FALSE]
  mat2 <- xg[idx2,,drop = FALSE]
  bins <- matrix(ifelse(rowSums(mat1 <= mat2) == d, 1, 0), n, ng)
  ## influ.1 <- bins - matrix(Ftheta(theta), n, ng, byrow = TRUE) ## matrix of n by ng
  ## the second part matrix(...) is not needed when zcenter == TRUE
  influ.1 <- bins
  influ.2 <- t(dF.dtheta %*% psi)
  return(influ.1 - influ.2) ## a matrix of n by ng
}


influStat <- function(infl, z) {
  n <- ncol(z)
  simobj <- z %*% infl / n
  stat <- apply(simobj, 1, function(w) c(max(abs(w)), n * mean(w * w)))
  stat
}

mult.sim <- function(xg, x, pdistr.my, ddistr.my, thetahat, psi, nsim,
                     ## zcenter = TRUE, ## should always be TRUE
                     gridStat = TRUE,
                     method.args, derCdfWrtPar) {
  p <- length(thetahat)
  n <- NROW(x)
  d <- NCOL(x)
  ## xg is null unless gridStat == TRUE
  
  z <- matrix(rnorm(n * nsim), nrow = n)
  ## if (zcenter) z <- t(scale(z, center = TRUE, scale = FALSE)) ## nsim by n
  ## else z <- t(z)  ## nsim by n
  z <- t(scale(z, center = TRUE, scale = FALSE)) ## nsim by n

  infla <- influFun(x, x, psi, pdistr.my, thetahat, method.args, derCdfWrtPar)
  stat.sim <- influStat(infla, z)
  rownames(stat.sim) <- c("ksa", "cvma")
  if (gridStat) {
    infl <- influFun(xg, x, psi, pdistr.my, thetahat, method.args, derCdfWrtPar)
    stat.sim <- rbind(influStat(infl, z), stat.sim)
    rownames(stat.sim) <- c("ks", "cvm", "ksa", "cvma")
  }
  stat.sim
}

getStat <- function(x, n, Fn, Ftheta) {
  emp <- Fn(x) - Ftheta(x)
  ## n <- NROW(x)
  ##  c(ks = max(abs(emp)), cvm = n * weighted.mean(emp * emp, fxg),
  ##  wrong! because xg is not uniform grid
  c(ks = max(abs(emp)), cvm = n * mean(emp * emp))
}
