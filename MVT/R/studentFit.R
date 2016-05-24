studentFit <-
function(x, data, family = Student(eta = .25), subset, na.action, control = mvt.control())
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
    stop("Student.fit applies only to numerical variables")
  znames <- dimnames(z)[[2]]
  dz <- dim(z)
  n <- dz[1]
  p <- dz[2]

  ## initial estimates
  center <- apply(z, 2, mean)
  f <- (n - 1) / n
  Scatter <- f * var(z)
  distances <- mahalanobis(z, center, Scatter)
  
  ## extract family info
  if (!inherits(family, "Student.family"))
    stop("Use only with 'Student.family' objects")
  if (is.null(family$family))
    stop("'family' not recognized")
  kind <- family$kind
  if ((kind < 0) || (kind > 1))
    stop("not valid 'family' object")
  settings <- c(kind, unlist(family$pars))

  ## set control values
  if (missing(control))
    control <- mvt.control()
  ctrl <- unlist(control)
  ctrl <- c(ctrl, 0)

  ## Call fitter
  now <- proc.time()
  fit <- .C("EM_fit",
            z = as.double(t(z)),
            dims = as.integer(dz),
            settings = as.double(settings),
            center = as.double(center),
            Scatter = as.double(Scatter),
            distances = as.double(distances),
            weights = as.double(rep(1, n)),
            logLik = double(1),
            control = as.double(ctrl))
  speed <- proc.time() - now

  ## creating the output object
  out <- list(call = Call,
              x = z,
              dims = dz,
              family = family,
              settings = fit$settings,
              start = list(center = center, Scatter = Scatter),
              center = fit$center,
              Scatter = matrix(fit$Scatter, ncol = p),
              logLik = fit$logLik,
              numIter = fit$control[4],
              control = control,
              weights = fit$weights,
              distances = fit$distances,
              speed = speed,
              converged = FALSE)
  names(out$center) <- znames
  dimnames(out$Scatter) <- list(znames, znames)
  eta <- out$settings[2]
  if (!control$fix.shape && kind != 0) {
    eta <- signif(eta, 4)
    out$family$call <- call(out$family$family, eta = eta)
  }
  out$eta <- out$settings[2]
  if (out$numIter < control$maxIter)
    out$converged <- TRUE
  class(out) <- "studentFit"
  out
}

print.studentFit <-
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
  if (x$dims[2] <= 5)
    print.symmetric(x$Scatter, digits = digits)
  else {
    print.symmetric(x$Scatter[1:5,1:5], digits = digits)
    cat("...")
  }
  nobs <- x$dims[1]
  cat("\nNumber of Observations:", nobs, "\n")
  invisible(x)
}
