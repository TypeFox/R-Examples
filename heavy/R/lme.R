heavyLme <-
function(fixed, random, groups, data = sys.frame(sys.parent()),
    family = Student(df = 4), subset, na.action = na.fail,
    control = heavy.control())
  UseMethod("heavyLme")

heavyLme.formula <-
function(fixed, random, groups, data = sys.frame(sys.parent()),
    family = Student(df = 4), subset, na.action = na.fail,
    control = heavy.control())
{
  ## local functions
  as.OneFormula <-
  function(..., omit = c(".", "pi"))
  { # combine formulas of a set of objects
    func <- function(x) { # should make all.vars generic
      if (is.list(x)) {
        return(unlist(lapply(x, all.vars)))
      } 
      all.vars(x)
    }
    names <- unique(unlist(lapply(list(...), func)))
    names <- names[is.na(match(names, omit))]
    if(length(names)) {
      eval(parse(text = paste("~", paste(names, collapse = "+")))[[1]])
    } else NULL
  }
  ##

  Call <- match.call()
  ## checking arguments
  if(!inherits(fixed, "formula") || length(fixed) != 3) {
    stop("\nFixed-effects model must be a formula of the form \"resp ~ pred\"")
  }
  if(missing(random)) {
    Call[["random"]] <- fixed[-2]
    random <- fixed[[3]]
  }
  groups <- asOneSidedFormula(groups)

  ## extract a data frame with enough information to evaluate 
  ## fixed, random and groups
  mfArgs <- list(formula = as.OneFormula(fixed, random, groups),
    data = data, na.action = na.action)
  if (!missing(subset)) {
    mfArgs[["subset"]] <- asOneSidedFormula(subset)[[2]]
  }
  dataMix <- do.call("model.frame", mfArgs)
  origOrder <- row.names(dataMix) # preserve the original order
  y <- eval(fixed[[2]], dataMix)
  ## sort the model.frame by groups and get the matrices and parameters
  ## used in the estimation procedures
  grp <- as.factor(eval(groups[[2]], dataMix))
  ord <- sort.list(grp) # order in which to sort the groups
  ## sorting the model.frame by groups and getting the matrices and parameters
  ## used in the estimation procedures
  grp <- grp[ord] # sorted groups
  glen <- tabulate(grp) # groups lengths
  ugrp <- as.character(unique(grp)) # unique groups
  n <- length(ugrp) # number of groups
  dataMix <- dataMix[ord,]  # sorted model.frame
  y <- y[ord] # sorted response vector
  N <- sum(glen)
  ZXrows <- nrow(dataMix)
  if(length(all.vars(fixed)) > 1) {
    X <- model.matrix(fixed, model.frame(fixed, dataMix))
  } else {
    X <- matrix(1, N, 1)
    dimnames(X) <- list(row.names(dataMix), "(Intercept)")
  }
  Xnames <- dimnames(X)[[2]]
  p <- ncol(X)
  if(length(all.vars(random))) {
    Z <- model.matrix(random, model.frame(random, dataMix))
  } else {
    Z <- X[,1,drop=F]
    dimnames(Z) <- list(dimnames(X)[[1]],"(Intercept)")
    Z[,1] <- 1
  }
  Znames <- dimnames(Z)[[2]];
  q <- ncol(Z)
  ZXcols <- p + q
  qsq <- q^2
  ## generating parameters used throughout the calculations
  ZX <- array(c(Z,X), c(ZXrows, ZXcols), list(rep("",N), c(Znames, Xnames)))
  qraux <- rep(0, n * ZXcols)
  DcLen <- pmin(glen, ZXcols)
  DcRows <- sum(DcLen)
  dims <- c(n, q, p, N, ZXrows, ZXcols, DcRows)

  ## initial estimates
  fit <- lsfit(X, y, intercept = FALSE)[1:2]
  scale <- sum(fit$residuals^2) / N
  theta <- apply(Z, 2, function(x) sum(x^2))
  theta <- 0.375 * sqrt(theta / ZXrows)
  theta <- diag(theta, nrow = length(theta))
  start <- c(fit$coefficients, theta, scale)
  names(start) <- NULL

  ## constructing the lmeData list
  lmeData <- list(ZX = ZX, y = y, qraux = qraux, dims = dims, glen = glen, 
                  DcLen = DcLen, grp = grp, ugrp = ugrp, origOrder = origOrder,
                  start = start)
  storage.mode(lmeData$ZX) <- "double"
  storage.mode(lmeData$y) <- "double"

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
    control <- heavy.control(algorithm = "NEM")
  if (!control$algorithm)
    control$ncycles <- 1
  ctrl <- unlist(control)
  ctrl <- c(ctrl, 0)

  ## call fitter
  now <- proc.time()
  fit <- .C("lme_fit",
            ZX = lmeData$ZX,
            y = lmeData$y,
            qraux = as.double(lmeData$qraux),
            dims = as.integer(lmeData$dims),
            glen = as.integer(lmeData$glen),
            DcLen = as.integer(lmeData$DcLen),
            settings = as.double(settings),
            coefficients = as.double(fit$coefficients),
            theta = as.double(theta),
            scale = as.double(scale),
            ranef = double(n * q),
            Root = double(n * qsq),
            distances = double(n),
            weights = as.double(rep(1, n)),
            logLik = double(1),
            control = as.double(ctrl))
  speed <- proc.time() - now
  
  ## compute fitted values and residuals
  Fitted <- .C("lme_fitted",
	       ZX = lmeData$ZX,
	       dims = as.integer(lmeData$dims),
	       glen = as.integer(lmeData$glen),
	       DcLen = as.integer(lmeData$DcLen),
	       coefficients = as.double(fit$coefficients),
	       ranef = as.double(fit$ranef),
	       conditional = as.double(rep(0, ZXrows)),
	       marginal = as.double(rep(0, ZXrows)))[c("marginal", "conditional")]
  Fitted <- as.data.frame(Fitted)
  Resid <- y - Fitted
  ## putting back in the original order of the data
  Fitted <- Fitted[origOrder,]
  Resid <- Resid[origOrder,]

  ## creating the output object
  lmeData$ZX <- fit$ZX
  lmeData$y <- fit$y
  lmeData$qraux <- fit$qraux
  out <- list(lmeData = lmeData,
              call = Call,
              family = family,
              coefficients = fit$coefficients,
              theta = matrix(fit$theta, ncol = q),
              scale = fit$scale,
              logLik = fit$logLik,
              numIter = fit$control[7],
              control = control,
              settings = fit$settings,
              ranef = matrix(fit$ranef, ncol = q, byrow = TRUE),
              weights = fit$weights,
              distances = fit$distances,
              Root.unscaled = fit$Root,
              Fitted = Fitted,
              Resid = Resid,
              speed = speed,
              converged = FALSE)
  if (!control$fix.shape) {
    if ((kind > 1) && (kind < 4)) {
      df <- signif(out$settings[3], 6)
      out$family$call <- call(out$family$family, df = df)
    }
  }
  if (out$numIter < control$maxIter)
    out$converged <- TRUE
  names(out$coefficients) <- Xnames
  dimnames(out$theta) <- list(Znames, Znames)
  dimnames(out$ranef) <- list(ugrp, Znames)
  names(out$weights) <- ugrp
  names(out$distances) <- ugrp
  class(out) <- "heavyLme"
  out
}

print.heavyLme <-
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
  cat("Linear mixed-effects model under heavy-tailed distributions\n")
  if (x$converged)
    cat(" Converged in", x$numIter, "iterations\n")
  else
    cat(" Maximum number of iterations exceeded\n")
  cat(" Data:", paste(as.name(x$call$data), ";", sep = ""))
  print(x$family)
  cat(" Log-likelihood:", format(x$logLik), "\n")
  cat("\nFixed:", deparse(x$call$fixed), "\n")
  print(x$coefficients)
  cat("\nRandom effects:\n")
  cat(" Formula:", paste(deparse(x$call$random), ";", sep = ""))
  cat(" Groups:", deparse(x$call$groups), "\n")
  cat(" Scale matrix estimate:\n")
  print.symmetric(x$theta)
  cat("Within-Group scale parameter:", format(x$scale), "\n")
  cat("\nNumber of Observations:", x$lmeData$dims[4], "\n")
  cat("Number of Groups:", x$lmeData$dims[1], "\n")
}

summary.heavyLme <-
function (object, ...)
{
  z <- object
  dims <- z$lmeData$dims
  p <- dims[3]
  ZX <- z$lmeData$ZX
  storage.mode(ZX) <- "double"
  control <- unlist(z$control)
  ## coefficients covariance matrix
  o <- .C("lme_acov",
          ZX = ZX,
          dims = as.integer(dims),
          glen = as.integer(z$lmeData$glen),
          DcLen = as.integer(z$lmeData$DcLen),
          settings = as.double(z$settings),
          Root = as.double(z$Root.unscaled),
          scale = as.double(z$scale),
          control = as.double(z$control),
          acov = double(p^2))
  est <- z$coefficients
  cov.unscaled <- solve(matrix(o$acov, ncol = p))
  acov <- z$scale * cov.unscaled
  se <- sqrt(diag(acov))
  zval <- est / se
  ans <- z[c("call", "family")]
  ans$dims <- z$lmeData$dims
  ans$logLik <- z$logLik
  ans$scale <- z$scale
  ans$theta <- z$theta
  ans$cov.unscaled <- cov.unscaled
  ans$coefficients <- cbind(est, se, zval, 2 * pnorm(abs(zval), lower.tail = FALSE))
  dimnames(ans$coefficients) <- list(names(z$coefficients),
        c("Estimate", "Std.Error", "Z-value", "p-value"))
  ans$correlation <- acov / outer(se, se)
  dimnames(ans$correlation) <- dimnames(ans$coefficients)[c(1,1)]
  class(ans) <- "summary.heavyLme"
  ans
}

print.summary.heavyLme <-
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
  cat("Linear mixed-effects model under heavy-tailed distributions\n")
  cat(" Data:", paste(as.name(x$call$data), ";", sep = ""))
  print(x$family)
  cat(" Log-likelihood:", format(x$logLik), "\n")
  cat("\nRandom effects:\n")
  cat(" Formula:", paste(deparse(x$call$random), ";", sep = ""))
  cat(" Groups:", deparse(x$call$groups), "\n")
  cat(" Scale matrix estimate:\n")
  print.symmetric(x$theta)
  cat("Within-Group scale parameter:", format(x$scale), "\n")
  cat("\nFixed:", deparse(x$call$fixed), "\n")
  print(format(round(x$coef, digits = digits)), quote = F, ...)
  cat("\nNumber of Observations:", x$dims[4], "\n")
  cat("Number of Groups:", x$dims[1], "\n")
  invisible(x)
}
