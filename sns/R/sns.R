sns <- function(x, fghEval, rnd=TRUE, gfit=NULL, mh.diag = FALSE, part = NULL
  , numderiv = 0, numderiv.method = c("Richardson", "simple"), numderiv.args = list()
  , ...)
{ 

  fitGaussian <- function(x, f, ...) 
  {
    ret <- f(x,...)                # Evaluate the function at 'x'
    Sigma <- solve(-ret$h)           
    mu <- x + Sigma %*% ret$g

    return (list(mu=as.vector(mu), # Newton method solution
               Sigma=Sigma,        # Inverse Hessian or Covariance matrix
               iSigma=-ret$h,      # Inverse covariance or Hessian
               f=ret$f,            # function value 
               g=ret$g))           # gradient 
  }
  
  numderiv <- as.integer(numderiv)
  if (numderiv > 0) {
    numderiv.method <- match.arg(numderiv.method)
    if (numderiv == 1) { # we need numeric hessian
      fghEval.int <- function(x, ...) {
        fg <- fghEval(x, ...)
        h <- hessian(func = function(x, ...) fghEval(x, ...)$f, x = x, ..., method = numderiv.method, method.args = numderiv.args)
        return (list(f = fg$f, g = fg$g, h = h))
      }
    } else { # we need numeric gradient and hessian
      fghEval.int <- function(x, ...) {
        f <- fghEval(x, ...)
        g <- grad(func = fghEval, x = x, ..., method = numderiv.method, method.args = numderiv.args)
        h <- hessian(func = fghEval, x = x, ..., method = numderiv.method, method.args = numderiv.args)
        return (list(f = f, g = g, h = h))
      }
    }
  } else {
    fghEval.int <- fghEval
  }
  
  if (!is.null(part)) { # code for 'state space partitioning', recursive call
    fghEval.part <- function(xsub, xfull, subset, ...) {
      xfull[subset] <- xsub
      ret <- fghEval.int(xfull, ...)
      return (list(f = ret$f, g = ret$g[subset], h = ret$h[subset, subset]))
    }
    npart <- length(part)
    accept <- logical(npart)
    if (mh.diag) {
      diag <- array(NA, dim = c(4, npart))
      dimnames(diag) <- list(c("log.p", "log.p.prop", "log.q", "log.q.prop"), 1:npart)
    }
    for (n in 1:npart) {
      subset <- part[[n]]
      rettmp <- sns(x = x[subset], fghEval = fghEval.part, rnd = rnd, gfit = NULL, part = NULL, mh.diag = mh.diag, xfull = x, subset = subset, ...)
      x[subset] <- rettmp
      accept[n] <- attr(rettmp, "accept")
      if (mh.diag) diag[, n] <- attr(rettmp, "mh.diag")
    }
    attr(x, "lp") <- attr(rettmp, "lp")
    attr(x, "accept") <- accept
    attr(x, "gfit") <- fitGaussian(x, fghEval.int, ...)
    if (mh.diag) attr(x, "mh.diag") <- mh.diag
    return (x)
  }

  # rnd: if FALSE, perform Newton's optimization (non-stochastic)
  # Fit Gaussian at x
  if (is.null(gfit)) gfit <- fitGaussian(x = x, f = fghEval.int, ...)
  mu     <- gfit$mu 
  Sigma  <- gfit$Sigma   # Covariance
  iSigma <- gfit$iSigma  # Inverse covariance

  K <- length(x)
      
  if (rnd) {
    # Draw sample from proposal distribution (Gaussian fit at x)
    x.prop <- as.vector(rmvnorm(n=1, mean=mu, sigma=Sigma))
  } else {
    # Run (non-stochastic) Newton optimization 
    rho <- 0.5; c <- 0.5;
    alphak <- 1; 
    d <- mu - x; # use newton's direction as step
    search_x <- as.vector(mu);
        
    fk <- gfit$f; # Values at the current point
    gk <- gfit$g; 
    fk1 <- fghEval.int(search_x, ...)$f; # Function value at searching point
    ls_iter <- 1;
    # Linesearch by backtracking from full Newton step
    while (fk1 < fk + c*alphak*(t(gk)%*%d) && ls_iter < 20) { 
      alphak <- alphak*rho; # if so, then go half way
      search_x <- x + alphak*d;
      fk1 <- fghEval.int(search_x, ...)$f;
      ls_iter <- ls_iter + 1;
    }
    x.prop <- as.vector(search_x);
  }

  log.q.prop <- dmvnorm(as.vector(x.prop), mu, Sigma, log=TRUE)
  
  # fit Gaussian at x.prop
  gfit.prop <- fitGaussian(x=x.prop,f=fghEval.int,...)
  mu.prop <- gfit.prop$mu
  Sigma.prop <- gfit.prop$Sigma
  iSigma.prop <- gfit.prop$iSigma

  # create MH acceptance ratio
  log.q <- dmvnorm(as.vector(x), mu.prop, Sigma.prop, log=TRUE)
  
  log.p <- gfit$f
  log.p.prop <- gfit.prop$f
  
  log.ratio <- (log.p.prop-log.p) + (log.q-log.q.prop)
  ratio <- min(1, exp(log.ratio))
    
  # perform acceptance test
  if (!rnd || ratio==1 || runif(1)<ratio) {
   gfit <- gfit.prop
	 x <- x.prop;
	 attr(x,"accept") <- TRUE
	 attr(x,"lp") <- log.p.prop
  } else {
	 attr(x,"accept") <- FALSE
	 attr(x,"lp") <- log.p
  }
  attr(x,"gfit") <- gfit

  if (mh.diag) {
    diag <- c(log.p, log.p.prop, log.q, log.q.prop)
    names(diag) <- c("log.p", "log.p.prop", "log.q", "log.q.prop")
    attr(x, "mh.diag") <- diag
  }

  return (x)
}

sns.run <- function(init, fghEval, niter = 100, nnr = min(10, round(niter/4))
  , mh.diag = FALSE, part = NULL, print.level = 0
  , report.progress = ceiling(niter/10)
  , numderiv = 0, numderiv.method = c("Richardson", "simple"), numderiv.args = list()
  , ...)
{
  fghEval.int <- sns.fghEval.numaug(fghEval, numderiv, numderiv.method, numderiv.args)

  # checking arguments
  stopifnot(niter >= 1)
  if (report.progress <= 0) {
      warning("invalid value specific for 'report.progress', using default.")
      report.progress <- ceiling(niter/10)
  }
  if (missing(init)) stop("starting point for MCMC chain must be specified")
  npart <- max(length(part),1)

  # initialization
  K <- length(init)
  x <- init
  accept <- matrix(logical(niter * npart), ncol = npart)
  lp <- double(niter)
  chain <- matrix( , nrow = niter, ncol = K)
  if (mh.diag) {
    diagnostic <- array(NA, dim = c(niter, 4, npart))
    dimnames(diagnostic) <- list(1:niter, c("log.p", "log.p.prop", "log.q", "log.q.prop"), 1:npart)
    if (niter > nnr && npart == 1) f.reldev <- rep(NA, niter - nnr)
  }

  # performing nr/mcmc iterations
  for (i in 1:niter) {
      x <- sns(x, fghEval.int, rnd = i > nnr, gfit = attr(x, "gfit"), mh.diag = mh.diag, part = part
               , numderiv = 0#, numderiv.method = numderiv.method, numderiv.args = numderiv.args
               , ...)
      accept[i, ] <- attr(x, "accept")
      lp[i] <- attr(x, "lp")
      chain[i, ] <- x
      if (mh.diag) {
        diagnostic[i, , ] <- attr(x, "mh.diag")
        if (npart == 1) {
          if (i == nnr) {
            f.ref <- attr(x, "gfit")$f
            x.ref <- attr(x, "gfit")$mu
          } else if (i > nnr) {
            if (attr(x, "gfit")$f < f.ref) f.reldev[i - nnr] <- (attr(x, "gfit")$f - f.ref + 0.5 * t(x - x.ref) %*% attr(x, "gfit")$iSigma %*% (x - x.ref)) / (f.ref - attr(x, "gfit")$f)
          }
        }
      }
      if (print.level && (i %% report.progress == 0))
          cat(paste0("finished iter ", i, " of ", niter, "\n"))
  }

  # assembling output
  attr(chain, "init") <- init
  attr(chain, "lp.init") <- fghEval.int(init, ...)$f
  attr(chain, "accept") <- accept
  attr(chain, "lp") <- lp
  attr(chain, "nnr") <- nnr
  attr(chain, "part") <- part
  if (mh.diag) {
    attr(chain, "mh.diag") <- drop(diagnostic)
    if (niter > nnr && npart == 1) attr(chain, "reldev") <- f.reldev
  }
  class(chain) <- "sns"
  return (chain)
}

