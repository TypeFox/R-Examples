mnp <- function(formula, data = parent.frame(), choiceX = NULL,
                cXnames = NULL, base = NULL, latent = FALSE,
                invcdf = FALSE, trace = TRUE, n.draws = 5000, p.var = "Inf", 
                p.df = n.dim+1, p.scale = 1, coef.start = 0,
                cov.start = 1, burnin = 0, thin = 0, verbose = FALSE) {   
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  mf$choiceX <- mf$cXnames <- mf$base <- mf$n.draws <- mf$latent <-
    mf$p.var <- mf$p.df <- mf$p.scale <- mf$coef.start <- mf$invcdf <-
      mf$trace <- mf$cov.start <- mf$verbose <- mf$burnin <- mf$thin <- NULL   
  mf[[1]] <- as.name("model.frame")
  mf$na.action <- 'na.pass'
  mf <- eval.parent(mf)

  ## fix this parameter
  p.alpha0 <- 1
  
  ## obtaining Y
  tmp <- ymatrix.mnp(mf, base=base, extra=TRUE, verbose=verbose)
  Y <- tmp$Y
  MoP <- tmp$MoP
  lev <- tmp$lev
  base <- tmp$base
  p <- tmp$p
  n.dim <- p - 1
  if(verbose)
    cat("\nThe base category is `", base, "'.\n\n", sep="") 
  if (p < 3)
    stop("The number of alternatives should be at least 3.")
  if(verbose) 
    cat("The total number of alternatives is ", p, ".\n\n", sep="") 
  if(verbose) {
    if (trace)
      cat("The trace restriction is used instead of the diagonal restriction.\n\n")
    else
      cat("The diagonal restriction is used instead of the trace restriction.\n\n")
  }
  
  ### obtaining X
  tmp <- xmatrix.mnp(formula, data=eval.parent(data),
                     choiceX=call$choiceX, cXnames=cXnames, 
                     base=base, n.dim=n.dim, lev=lev, MoP=MoP,
                     verbose=verbose, extra=TRUE)
  X <- tmp$X
  coefnames <- tmp$coefnames
  n.cov <- ncol(X) / n.dim
  
  ## listwise deletion for X
  na.ind <- apply(is.na(X), 1, sum)
  if (ncol(Y) == 1)
    na.ind <- na.ind + is.na(Y)
  Y <- Y[na.ind==0,]
  X <- X[na.ind==0,]
  n.obs <- nrow(X)
  
  if (verbose) {
    cat("The dimension of beta is ", n.cov, ".\n\n", sep="")
    cat("The number of observations is ", n.obs, ".\n\n", sep="")
    if (sum(na.ind>0)>0) {
      if (sum(na.ind>0)==1)
        cat("The observation ", (1:length(na.ind))[na.ind>0], " is dropped due to missing values.\n\n", sep="")
      else {
        cat("The following ", sum(na.ind>0), " observations are dropped due to missing values:\n", sep="")
        cat((1:length(na.ind))[na.ind>0], "\n\n")
      }
    }
  } 
  
  ## checking the prior for beta
  p.imp <- FALSE 
  if (p.var == Inf) {
    p.imp <- TRUE
    p.prec <- diag(0, n.cov)
    if (verbose)
      cat("Improper prior will be used for beta.\n\n")
  }
  else if (is.matrix(p.var)) {
    if (ncol(p.var) != n.cov || nrow(p.var) != n.cov)
      stop("The dimension of `p.var' should be ", n.cov, " x ", n.cov, sep="")
    if (sum(sign(eigen(p.var)$values) < 1) > 0)
      stop("`p.var' must be positive definite.")
    p.prec <- solve(p.var)
  }
  else {
    p.var <- diag(p.var, n.cov)
    p.prec <- solve(p.var)
  }
  p.mean <- rep(0, n.cov)

  ## checking prior for Sigma
  p.df <- eval(p.df)
  if (length(p.df) > 1)
    stop("`p.df' must be a positive integer.")
  if (p.df < n.dim)
    stop(paste("`p.df' must be at least ", n.dim, ".", sep=""))
  if (abs(as.integer(p.df) - p.df) > 0)
    stop("`p.df' must be a positive integer.")
  if (!is.matrix(p.scale))  
    p.scale <- diag(p.scale, n.dim)
  if (ncol(p.scale) != n.dim || nrow(p.scale) != n.dim)
    stop("`p.scale' must be ", n.dim, " x ", n.dim, sep="")
  if (sum(sign(eigen(p.scale)$values) < 1) > 0)
    stop("`p.scale' must be positive definite.")
  else if ((trace == FALSE) & (p.scale[1,1] != 1)) {
    p.scale[1,1] <- 1
    warning("p.scale[1,1] will be set to 1.")
  }
  Signames <- NULL
  for(j in 1:n.dim)
    for(k in 1:n.dim)
      if (j<=k)
        Signames <- c(Signames, paste(if(MoP) lev[j] else lev[j+1],
                                      ":", if(MoP) lev[k] else lev[k+1], sep="")) 

  ## checking starting values
  if (length(coef.start) == 1)
    coef.start <- rep(coef.start, n.cov)
  else if (length(coef.start) != n.cov)
    stop(paste("The dimenstion of `coef.start' must be  ",
               n.cov, ".", sep=""))
  if (!is.matrix(cov.start)) {
    cov.start <- diag(n.dim)*cov.start
    if (!trace)
      cov.start[1,1] <- 1
  }
  else if (ncol(cov.start) != n.dim || nrow(cov.start) != n.dim)
    stop("The dimension of `cov.start' must be ", n.dim, " x ", n.dim, sep="")
  else if (sum(sign(eigen(cov.start)$values) < 1) > 0)
    stop("`cov.start' must be a positive definite matrix.")
  else if ((trace == FALSE) & (cov.start[1,1] != 1)) {
    cov.start[1,1] <- 1
    warning("cov.start[1,1] will be set to 1.")
  }
  
  ## checking thinnig and burnin intervals
  if (burnin < 0)
    stop("`burnin' should be a non-negative integer.") 
  if (thin < 0)
    stop("`thin' should be a non-negative integer.")
  keep <- thin + 1
  
  ## running the algorithm
  if (latent)
    n.par <- n.cov + n.dim*(n.dim+1)/2 + n.dim*n.obs
  else
    n.par <- n.cov + n.dim*(n.dim+1)/2
  if(verbose)
    cat("Starting Gibbs sampler...\n")
  # recoding NA into -1
  Y[is.na(Y)] <- -1

  param <- .C("cMNPgibbs", as.integer(n.dim),
              as.integer(n.cov), as.integer(n.obs), as.integer(n.draws),
              as.double(p.mean), as.double(p.prec), as.integer(p.df),
              as.double(p.scale*p.alpha0), as.double(X), as.integer(Y), 
              as.double(coef.start), as.double(cov.start), 
              as.integer(p.imp), as.integer(invcdf),
              as.integer(burnin), as.integer(keep), as.integer(trace),
              as.integer(verbose), as.integer(MoP), as.integer(latent),
              pdStore = double(n.par*floor((n.draws-burnin)/keep)),
              PACKAGE="MNP")$pdStore 
  param <- matrix(param, ncol = n.par,
                  nrow = floor((n.draws-burnin)/keep), byrow=TRUE)
  if (latent) {
    W <- array(as.vector(t(param[,(n.par-n.dim*n.obs+1):n.par])),
               dim = c(n.dim, n.obs, floor((n.draws-burnin)/keep)),
               dimnames = list(lev[!(lev %in% base)], rownames(Y), NULL))
    param <- param[,1:(n.par-n.dim*n.obs)]
    }
  else
    W <- NULL
  colnames(param) <- c(coefnames, Signames)
    
  ##recoding -1 back into NA
  Y[Y==-1] <- NA

  ## returning the object
  res <- list(param = param, x = X, y = Y, w = W, call = call, alt = lev,
              n.alt = p, base = base, invcdf = invcdf, trace = trace,
              p.mean = if(p.imp) NULL else p.mean, p.var = p.var, 
              p.df = p.df, p.scale = p.scale, burnin = burnin, thin = thin) 
  class(res) <- "mnp"
  return(res)
}
  


