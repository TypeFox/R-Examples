"gelman.diag" <- function (x, confidence = 0.95, transform = FALSE,
                           autoburnin=TRUE, multivariate=TRUE) 
  ## Gelman and Rubin's diagnostic
  ## Gelman, A. and Rubin, D (1992). Inference from iterative simulation
  ## using multiple sequences.  Statistical Science, 7, 457-551.
  ##
  ## Correction and Multivariate generalization:
  ## Brooks, S.P. and Gelman, A. (1997) General methods for monitoring
  ## convergence of iterative simulations. Journal of Computational and
  ## Graphical Statistics, 7, 434-455.

{
  x <- as.mcmc.list(x)
  if (nchain(x) < 2) 
    stop("You need at least two chains")
  ## RGA added an autoburnin parameter here, because if I have already
  ## trimmed burn in, I don't want to do it again.
  if (autoburnin && start(x) < end(x)/2 ) 
    x <- window(x, start = end(x)/2 + 1)
  Niter <- niter(x)
  Nchain <- nchain(x)
  Nvar <- nvar(x)
  xnames <- varnames(x)
  
  if(transform)
    x <- gelman.transform(x)
  ##
  ## Estimate mean within-chain variance (W) and between-chain variance
  ## (B/Niter), and calculate sampling variances and covariance of the
  ## estimates (varW, varB, covWB)
  ##
  ## Multivariate (upper case)
  x <- lapply(x, as.matrix)
  S2 <- array(sapply(x, var, simplify=TRUE), dim=c(Nvar,Nvar,Nchain))
  W <- apply(S2, c(1,2), mean)
  xbar <- matrix(sapply(x, apply, 2, mean, simplify=TRUE), nrow=Nvar,
                 ncol=Nchain)
  B <- Niter * var(t(xbar))

  if(Nvar > 1 && multivariate) {
      ## We want the maximal eigenvalue of the square matrix X that
      ## solves WX = B. It is numerically easier to work with a
      ## symmetric matrix that has the same eigenvalues as X.
      if (is.R()) {
          CW <- chol(W)
          emax <- eigen(backsolve(CW, t(backsolve(CW, B, transpose=TRUE)),
                                  transpose=TRUE),
                        symmetric=TRUE, only.values=TRUE)$values[1]
      }
      else {
          emax <- eigen(qr.solve(W,B), symmetric=FALSE, only.values=TRUE)$values
      }
      mpsrf <- sqrt( (1 - 1/Niter) + (1 + 1/Nvar) * emax/Niter )
  }
  else
    mpsrf <- NULL
  ## Univariate (lower case)
  w <- diag(W)
  b <- diag(B)


  s2 <- matrix(apply(S2, 3, diag), nrow=Nvar, ncol=Nchain)
  muhat <- apply(xbar,1,mean)
  var.w <- apply(s2, 1, var)/Nchain              
  var.b <- (2 * b^2)/(Nchain - 1)      
  cov.wb <- (Niter/Nchain) * diag(var(t(s2), t(xbar^2)) -
                              2 * muhat * var(t(s2), t(xbar)))
  ##
  ## Posterior interval combines all uncertainties in a t interval with
  ## center muhat, scale sqrt(V), and df.V degrees of freedom.
  ##
  V <- (Niter - 1) * w / Niter  + (1 + 1/Nchain) * b/ Niter
  var.V <- ((Niter - 1)^2 * var.w + (1 + 1/Nchain)^2 * 
            var.b + 2 * (Niter - 1) * (1 + 1/Nchain) * cov.wb)/Niter^2
  df.V <- (2 * V^2)/var.V
  ##
  ## Potential scale reduction factor (that would be achieved by
  ## continuing simulations forever) is estimated by 
  ##   R = sqrt(V/W) * df.adj
  ## where df.adj is a degrees of freedom adjustment for the width
  ## of the t-interval.
  ##
  ## To calculate upper confidence interval we divide R2 = R^2 into two
  ## parts, fixed and random.  The upper limit of the random part is
  ## calculated assuming that B/W has an F distribution.
  ##
  df.adj <- (df.V + 3)/(df.V + 1)
  B.df <- Nchain - 1
  W.df <- (2 * w^2)/var.w
  R2.fixed <- (Niter - 1)/Niter
  R2.random <- (1 + 1/Nchain) * (1/Niter) * (b/w)
  R2.estimate <- R2.fixed + R2.random
  R2.upper <- R2.fixed + qf((1 + confidence)/2, B.df, W.df) * R2.random
  psrf <- cbind(sqrt(df.adj * R2.estimate), sqrt(df.adj * R2.upper))
  dimnames(psrf) <- list(xnames, c("Point est.", "Upper C.I."))
  
  out <- list(psrf = psrf, mpsrf=mpsrf)
  class(out) <- "gelman.diag"
  out
}

"gelman.transform" <- function(x)
  ## Gelman and Rubin diagnostic assumes a normal distribution. To
  ## improve the normal approximation, variables on [0, Inf) are log
  ## transformed, and variables on [0,1] are logit-transformed.
{
  if (!is.R())  {
    # in S-PLUS this function generates a superfluous warning,
    # so turn off all warnings during the function.
    oldWarn <- getOption("warn")
    options(warn=-1)
    on.exit(options (warn=oldWarn))
  }
  if (nvar(x) == 1) {
    z <- data.frame(lapply(x, unclass))
    if (min(z) > 0) {
      y <- if(max(z) < 1)
        log(z/(1-z))
      else log(z)
      for (j in 1:nchain(x)) x[[j]] <- y[,j]
    }
  }
  else for (i in 1:nvar(x)) {
    z <- data.frame(lapply(x[, i], unclass))
    if (min(z) > 0) {
      y <- if (max(z) < 1) 
        log(z/(1 - z))
      else log(z)
      for (j in 1:nchain(x)) x[[j]][, i] <- y[, j]
    }
  }
  return(x)
}

"gelman.mv.diag" <- function (x, confidence = 0.95, transform = FALSE)
{
  s2 <- sapply(x, var, simplify=TRUE)
  W <- matrix(apply(s2, 1, mean), nvar(x), nvar(x))
  xbar <- sapply(x, apply, 2, mean, simplify=TRUE)
  B <- niter(x) * var(t(xbar))
  emax <- eigen(qr.solve(W,B), symmetric=FALSE, only.values=TRUE)$values[1]
  mpsrf <- sqrt( (1 - 1/niter(x)) + (1 + 1/nvar(x)) * emax )
  return(mpsrf)
}

  
"print.gelman.diag" <-
  function (x, digits = 3, ...) 
{
  cat("Potential scale reduction factors:\n\n")
  print.default(x$psrf, digits = digits, ...)
  if(!is.null(x$mpsrf)) {
    cat("\nMultivariate psrf\n\n")
    cat(format(x$mpsrf,digits = digits))
  }
  cat("\n")
}

"gelman.plot" <-
  function (x, bin.width = 10, max.bins = 50, confidence = 0.95,
            transform = FALSE, autoburnin = TRUE, auto.layout = TRUE, ask,
            col = 1:2, lty = 1:2, xlab = "last iteration in chain",
            ylab = "shrink factor", type = "l", ...) 
{
  if (missing(ask)) {
    ask <- if (is.R()) {
      dev.interactive()
    }
    else {
      interactive()
    }
  }
  x <- as.mcmc.list(x)
  oldpar <- NULL
  on.exit(par(oldpar))
  if (auto.layout) 
    oldpar <- par(mfrow = set.mfrow(Nchains = nchain(x), Nparms = nvar(x)))
  y <- gelman.preplot(x, bin.width = bin.width, max.bins = max.bins, 
                      confidence = confidence, transform = transform,
                      autoburnin = autoburnin)
  all.na <- apply(is.na(y$shrink[, , 1, drop = FALSE]), 2, all)
  if (!any(all.na)) 
    for (j in 1:nvar(x)) {
      matplot(y$last.iter, y$shrink[, j, ], col = col, 
              lty = lty, xlab = xlab, ylab = ylab, type = type, 
              ...)
      abline(h = 1)
      ymax <- max(c(1, y$shrink[, j, ]), na.rm = TRUE)
      leg <- dimnames(y$shrink)[[3]]
      xmax <- max(y$last.iter)
      legend(xmax, ymax, legend = leg, lty = lty, bty = "n", 
             col = col, xjust = 1, yjust = 1)
      title(main = varnames(x)[j])
      if (j==1)
         oldpar <- c(oldpar, par(ask = ask))
    }
  return(invisible(y))
}

"gelman.preplot" <-
  function (x, bin.width = bin.width, max.bins = max.bins,
            confidence = confidence, transform = transform,
            autoburnin = autoburnin) 
{
  x <- as.mcmc.list(x)
  if (niter(x) <= 50) 
    stop("Less than 50 iterations in chain")
  nbin <- min(floor((niter(x) - 50)/thin(x)), max.bins)
  binw <- floor((niter(x) - 50)/nbin)
  last.iter <- c(seq(from = start(x) + 50 * thin(x), by = binw * 
                     thin(x), length = nbin), end(x))
  shrink <- array(dim = c(nbin + 1, nvar(x), 2))
  dimnames(shrink) <- list(last.iter, varnames(x),
                           c("median", paste(50 * (confidence + 1), "%",
                                             sep = ""))
                           )
  for (i in 1:(nbin + 1)) {
    shrink[i, , ] <- gelman.diag(window(x, end = last.iter[i]), 
                                 confidence = confidence,
                                 transform = transform,
                                 autoburnin = autoburnin)$psrf
  }
  all.na <- apply(is.na(shrink[, , 1, drop = FALSE]), 2, all)
  if (any(all.na)) {
    cat("\n******* Error: *******\n")
    cat("Cannot compute Gelman & Rubin's diagnostic for any chain \n")
    cat("segments for variables", varnames(x)[all.na], "\n")
    cat("This indicates convergence failure\n")
  }
  return(list(shrink = shrink, last.iter = last.iter))
}

if (!is.R()){

qr.solve <- function (a, b, tol = 1e-07) {
    if (!is.qr(a))
        a <- qr(a, tol = tol)
    nc <- ncol(a$qr)
    if (a$rank != nc)
        stop("singular matrix 'a' in solve")
    if (missing(b)) {
        if (nc != nrow(a$qr))
            stop("only square matrices can be inverted")
        b <- diag(1, nc)
    }
    return(qr.coef(a, b))
}

}
































