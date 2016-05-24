heavyFit <-
function(x, data, family = Student(df = 4), subset, na.action, control = heavy.control())
{
  Call <- match.call()
  if (missing(x))
    stop("'x' is not supplied")
  if (inherits(x, "formula")) {
    mt <- terms(x, data = data)
    if (attr(mt, "response") > 0)
      stop("response not allowed in formula")
    attr(mt, "intercept") <- 0
    mf <- match.call(expand.dots = FALSE)
    names(mf)[names(mf) == "x"] <- "formula"
    mf$family <- mf$control <- NULL
    mf[[1L]] <- as.name("model.frame")
    mf <- eval.parent(mf)
    na.act <- attr(mf, "na.action")
    z <- model.matrix(mt, mf)
  }
  else {
    z <- as.matrix(x)
    if (!missing(subset))
      z <- z[subset, , drop = FALSE]
    if (!missing(na.action))
      z <- na.omit(z)
    else
      z <- na.fail(z)
  }
  if (!is.numeric(z))
    stop("heavyFit applies only to numerical variables")
  znames <- dimnames(z)[[2]]
  dz <- dim(z)
  n <- dz[1]
  p <- dz[2]

  ## initial estimates
  center <- apply(z, 2, mean)
  f <- (n - 1) / n
  Scatter <- f * var(z)
  
  ## extract family info
  if (!inherits(family, "heavy.family"))
    stop("Use only with 'heavy.family' objects")
  if (is.null(family$family))
    stop("'family' not recognized")
  kind <- family$which
  if ((kind < 0) || (kind > 4))
    stop("not valid 'family' object")
  settings <- c(kind, family$npars, unlist(family$pars))

  ## set control values
  if (missing(control))
    control <- heavy.control()
  if (!control$algorithm)
    control$ncycles <- 1
  ctrl <- unlist(control)[1:4]
  ctrl <- c(ctrl, 0)

  ## Call fitter
  now <- proc.time()
  fit <- .C("mv_fit",
            z = as.double(t(z)),
            dims = as.integer(dz),
            settings = as.double(settings),
            center = as.double(center),
            Scatter = as.double(Scatter),
            distances = double(n),
            weights = as.double(rep(1, n)),
            logLik = double(1),
            acov = double(p^2),
            control = as.double(ctrl))
  speed <- proc.time() - now

  ## creating the output object
  out <- list(call = Call,
              x = z,
              dims = dz,
              family = family,
              settings = fit$settings,
              center = fit$center,
              Scatter = matrix(fit$Scatter, ncol = p),
              logLik = fit$logLik,
              numIter = fit$control[5],
              control = control,
              weights = fit$weights,
              distances = fit$distances,
              acov = matrix(fit$acov, ncol = p),
              speed = speed,
              converged = FALSE)
  names(out$center) <- znames
  dimnames(out$Scatter) <- list(znames, znames)
  if (!control$fix.shape) {
    if ((kind > 1) && (kind < 4)) {
      df <- signif(out$settings[3], 6)
      out$family$call <- call(out$family$family, df = df)
    }
  }
  if (out$numIter < control$maxIter)
    out$converged <- TRUE
  class(out) <- "heavyFit"
  out
}

print.heavyFit <-
function(x, digits = 4, ...)
{
  ## local functions
  print.symmetric <-
  function(z, digits = digits, ...)
  {
    ll <- lower.tri(z, diag = TRUE)
    z[ll] <- format(z[ll], ...)
    z[!ll] <- ""
    print(z, ..., quote = F)
  }
  cat("Call:\n")
  x$call$family <- x$family$call
  dput(x$call, control = NULL)
  if (x$converged)
    cat("Converged in", x$numIter, "iterations\n")
  else
    cat("Maximum number of iterations exceeded")
  cat("\nCenter:\n ")
  print(format(round(x$center, digits = digits)), quote = F, ...)
  cat("\nScatter matrix estimate:\n")
  print.symmetric(x$Scatter, digits = digits)
  nobs <- x$dims[1]
  cat("\nNumber of Observations:", nobs, "\n")
  invisible(x)
}

summary.heavyFit <-
function (object, ...)
{
  z <- object
  se <- sqrt(diag(z$acov))
  est <- z$center
  zval <- est / se
  ans <- z[c("call", "dims")]
  ans$family <- z$family
  ans$logLik <- z$logLik
  ans$Scatter <- z$Scatter
  ans$center <- cbind(est, se, zval, 2 * pnorm(abs(zval), lower.tail = FALSE))
  dimnames(ans$center) <- list(names(z$center),
        c("Estimate", "Std.Error", "Z value", "p-value"))
  ans$correlation <- z$acov / outer(se, se)
  dimnames(ans$correlation) <- dimnames(ans$center)[c(1,1)]
  class(ans) <- "summary.heavyFit"
  ans
}

print.summary.heavyFit <-
function(x, digits = 4, ...)
{
  ## local functions
  print.symmetric <-
  function(z, digits = digits, ...)
  {
    ll <- lower.tri(z, diag = TRUE)
    z[ll] <- format(z[ll], ...)
    z[!ll] <- ""
    print(z, ..., quote = F)
  }
  cat("Estimation under multivariate heavy-tailed distributions\n")
  cat(" Data:", paste(as.name(x$call$data), ";", sep = ""))
  print(x$family)
  cat("\nCenter:\n ")
  print(format(round(x$center, digits = digits)), quote = F, ...)
  cat("\nScatter matrix estimate:\n")
  print.symmetric(x$Scatter, digits = digits)
  nobs <- x$dims[1]
  cat("\nNumber of Observations:", nobs)
  p <- x$dims[2]
  cat("\nLog-likelihood:", format(x$logLik), "on", p * (p + 3) / 2, "degrees of freedom\n")
  invisible(x)
}
