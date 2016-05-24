convInfo <- function(x, digits, show. = getOption("show.nls.convergence", TRUE)) 
{
  if (!is.null(x$convInfo)) 
    with(x$convInfo, {
      if (identical(x$call$algorithm, "port")) 
        cat("\nAlgorithm \"port\", convergence message: ", 
            stopMessage, "\n", sep = "")
      else {
        if (!isConv || show.) {
          cat("\nNumber of iterations", if (isConv) 
            "to convergence:"
            else "till stop:", finIter, "\nAchieved convergence tolerance:", 
            format(finTol, digits = digits))
          cat("\n")
        }
        if (!isConv) {
          cat("Reason stopped:", stopMessage)
          cat("\n")
        }
      }
    })
  invisible()
}

coef.onls <- function(object, ...)
{
  unlist(object$parONLS)
}

formula.onls <- function(x, ...)
{
  x$formula
}

weights.onls <- function(object, ...)
{
  object$weights
}

vcov.onls <- function(object, ...)
{
  sm <- summary(object)
  sm$cov.unscaled * sm$sigmaONLS^2  
}

df.residual.onls <- function (object, ...) 
{
  w <- object$weights
  n <- if (!is.null(w)) sum(w != 0) else length(residuals(object)) 
  n - length(coef(object))
}

residuals.onls <- function(object, ...)
{
  val <- as.vector(object$residONLS)
  if (!is.null(object$na.action)) 
    val <- napredict(object$na.action, val)
  lab <- "Vertical residuals from orthogonal model"
  attr(val, "label") <- lab
  val  
}

fitted.onls <- function(object, ...)
{
  val <- as.vector(object$fittedONLS)
  if (!is.null(object$na.action)) 
    val <- napredict(object$na.action, val)
  lab <- "Fitted values from orthogonal model"
  attr(val, "label") <- lab
  val  
}

deviance.onls <- function(object, ...)
{
  val <- sum(object$residONLS^2)  
  lab <- "Deviance (RSS) of vertical residuals from orthogonal model"
  attr(val, "label") <- lab
  val  
}

logLik.onls <- function(object, REML = FALSE, ...) 
{
  if (REML) 
    stop("cannot calculate REML log-likelihood for \"onls\" objects")
  res <- object$residONLS
  N <- length(res)
  if (is.null(w <- object$weights)) w <- rep_len(1, N)
  zw <- w == 0
  val <- -N * (log(2 * pi) + 1 - log(N) - sum(log(w + zw)) + 
                 log(sum(w * res^2)))/2
  attr(val, "df") <- 1L + length(coef(object))
  attr(val, "nobs") <- attr(val, "nall") <- sum(!zw)
  attr(val, "label") <- "Log-likelihood using vertical residuals from orthogonal model"
  class(val) <- "logLik.onls"
  val
}

print.logLik.onls <- function(x, digits = getOption("digits"), ...) 
{
  cat("'log Lik.' ", paste(format(c(x), digits = digits), collapse = ", "), 
      " (df=", format(attr(x, "df")), ")\n", 
      "\"", attr(x, "label"), "\"", sep = "")
  invisible(x)
}

predict.onls <- function(object, newdata, se.fit = FALSE, scale = NULL, df = Inf, 
                          interval = c("none", "confidence", "prediction"), level = 0.95, ...) 
{
  if (missing(newdata)) 
    return(as.vector(object$fittedONLS))
  if (!is.null(cl <- object$dataClasses)) 
    .checkMFClasses(cl, newdata)
  m <- object$model
  formula <- object$formula
  pred.name <- attr(object$pred, "name")
  m[pred.name] <- newdata
  eval(formula[[3]], envir = m)
}

print.onls <- function(x, ...) 
{
  digits <- max(3L, getOption("digits") - 3L)
  cat("Nonlinear orthogonal regression model\n")
  cat("  model: ", deparse(formula(x)), "\n", sep = "")
  cat("   data: ", deparse(x$data), "\n", sep = "")
  print(coef(x), digits = digits, ...)
  cat(" ", if (!is.null(x$weights) && diff(range(x$weights))) 
    "weighted ", "vertical residual sum-of-squares: ", format(deviance(x), 
                                                     digits = digits), "\n", sep = "")
  cat(" ", if (!is.null(x$weights) && diff(range(x$weights))) 
    "weighted ", "orthogonal residual sum-of-squares: ", format(deviance_o(x), 
                                                     digits = digits), "\n", sep = "")
  if (sum(x$ortho == TRUE) == length(x$ortho))
    cat(" PASSED: ", sum(x$ortho == TRUE), " out of ", length(x$ortho), " fitted points are orthogonal.\n", sep = "")
  else
    cat(" FAILED: Only ", sum(x$ortho == TRUE), " out of ", length(x$ortho), " fitted points are orthogonal.\n", sep = "")
  convInfo(x, digits = digits)
  invisible(x)
}

summary.onls <- function(object, correlation = FALSE, symbolic.cor = FALSE, ...) 
{
  residONLS <- as.vector(object$residONLS)
  w <- object$weights
  n <- if (!is.null(w)) sum(w > 0) else length(residONLS)
  param <- coef(object)
  pnames <- names(param)
  p <- length(param)
  rdf <- n - p  
  
  resvarONLS <- if (rdf <= 0) NaN else deviance(object)/rdf
  resvar_o <- if (rdf <= 0) NaN else deviance_o(object)/rdf 
  
  XtXinv <- chol2inv(qr.R(object$QR))
  dimnames(XtXinv) <- list(pnames, pnames)
  
  se <- sqrt(diag(XtXinv) * resvarONLS)
  tval <- param/se
  
  param <- cbind(param, se, tval, 2 * pt(abs(tval), rdf, lower.tail = FALSE))
  dimnames(param) <- list(pnames, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
  
  ans <- list(formula = formula(object), residONLS = residONLS, sigmaONLS = sqrt(resvarONLS),
              sigma_o = sqrt(resvar_o), df = c(p, rdf), cov.unscaled = XtXinv, 
              call = object$call, convInfo = object$convInfo, control = object$control, 
              na.action = object$na.action, coefficients = param, parameters = param)
  if (correlation && rdf > 0) {
    ans$correlation <- (XtXinv * resvarONLS)/outer(se, se)
    ans$symbolic.cor <- symbolic.cor
  }
  class(ans) <- "summary.onls"
  ans
}

print.summary.onls <- function(x, digits = max(3L, getOption("digits") - 3L), 
                               symbolic.cor = x$symbolic.cor, 
                                signif.stars = getOption("show.signif.stars"), ...) 
{
  cat("\nFormula: ", paste(deparse(x$formula), sep = "\n", collapse = "\n"), "\n", sep = "")
  df <- x$df
  rdf <- df[2L]
  cat("\nParameters:\n")
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars, ...)
  cat("\nResidual standard error of vertical distances:", format(signif(x$sigmaONLS, digits)), 
      "on", rdf, "degrees of freedom")
  cat("\nResidual standard error of orthogonal distances:", format(signif(x$sigma_o, digits)), 
      "on", rdf, "degrees of freedom")
  cat("\n")
  correl <- x$correlation
  if (!is.null(correl)) {
    p <- NCOL(correl)
    if (p > 1) {
      cat("\nCorrelation of Parameter Estimates:\n")
      if (is.logical(symbolic.cor) && symbolic.cor) {
        print(symnum(correl, abbr.colnames = NULL))
      }
      else {
        correl <- format(round(correl, 2), nsmall = 2L, 
                         digits = digits)
        correl[!lower.tri(correl)] <- ""
        print(correl[-1, -p, drop = FALSE], quote = FALSE)
      }
    }
  }
  convInfo(x, digits = digits)
  if (nzchar(mess <- naprint(x$na.action))) 
    cat("  (", mess, ")\n", sep = "")
  cat("\n")
  invisible(x)
}

plot.onls <- function(
  x, 
  fitted.nls = TRUE,
  fitted.onls = TRUE, 
  segments = TRUE,
  ...)
{
  object <- x
  X <- object$pred
  Y <- object$resp
  MODEL <- object$model
  
  ## display ONLS curve  
  if (fitted.onls) {
    MODEL <- object$model
    nameX <- attr(X, "name")    
    tempX <- MODEL[[nameX]]
    onlsX <- seq(min(tempX, na.rm = TRUE), max(tempX, na.rm = TRUE), length.out = 1000)
    MODEL[[nameX]] <- onlsX
    FORMULA <- object$call$formula[[3]]
    onlsY <- eval(FORMULA, envir = MODEL)     
  }  
  
  ## display NLS curve
  if (fitted.nls) {
    MODEL <- object$model
    nameX <- attr(X, "name")    
    tempX <- MODEL[[nameX]]
    nlsX <- seq(min(tempX, na.rm = TRUE), max(tempX, na.rm = TRUE), length.out = 1000)
    MODEL[[nameX]] <- nlsX
    FORMULA <- object$call$formula[[3]]
    m <- match(names(object$parNLS), names(MODEL))
    MODEL[m] <- object$parNLS
    nlsY <- eval(FORMULA, envir = MODEL)     
  }  
  
  plot(X, Y, xlab = attr(X, "name"), ylab = attr(Y, "name"), type = "p", ...)
  if (fitted.onls) lines(onlsX, onlsY, col = "red", ...)
  if (fitted.nls) lines(nlsX, nlsY, col = "blue", ...)
  
  ## display (x, y) -> (x0, y0) segments
  if (segments) {
    points(x0(object), y0(object), col = "red", pch = 16, cex = 0.5, ...)
    segments(x0(object), y0(object), X, Y, col = "red", ...)      
  }  
}

profile.onls <- function(fitted, which = 1L:npar, maxpts = 100, 
                         alphamax = 0.01, delta.t = cutoff/5, ...) 
{
  f.summary <- summary(fitted)
  std.err <- f.summary$coefficients[, "Std. Error"]
  nobs <- length(resid(fitted))
  origPars <- unlist(fitted$parONLS)
  npar <- length(origPars) 
  lower <- fitted$call$lower 
  lower <- rep_len(if (!is.null(lower)) as.double(lower) else -Inf, npar)
  upper <- fitted$call$upper
  upper <- rep_len(if (!is.null(upper)) as.double(upper) else Inf, npar)
  if (is.character(which))  which <- match(which, names(origPars), 0)
  which <- which[which >= 1 & which <= npar]
  cutoff <- sqrt(qf(1 - alphamax, 1L, nobs - npar))
  out <- vector("list", npar)
  S.hat <- deviance_o(fitted)  
  s2.hat <- summary(fitted)$sigma_o^2
    
  for (par in which) {
    pars <- origPars  
    sgn <- -1
    count <- 1
    varying <- rep.int(TRUE, npar)
    varying[par] <- FALSE
    FIXED <- !varying
    tau <- double(2 * maxpts)
    par.vals <- array(0, c(2L * maxpts, npar), list(NULL, names(pars)))
    tau[1L] <- 0
    par.vals[1, ] <- pars
    base <- pars[par]     
    profile.par.inc <- delta.t * std.err[par]
    pars[par] <- base - profile.par.inc    
    pars[par] <- pmin(upper[par], pmax(lower[par], pars[par]))
        
    while (count <= maxpts) {
      if (is.na(pars[par]) || isTRUE(all.equal(pars, par.vals[1, ])) 
          || pars[par] < lower[par] || pars[par] > upper[par] 
          || abs(pars[par] - base)/std.err[par] > 10 * cutoff) break    
      
      FIT1 <- update(fitted, start = pars, fixed = FIXED, verbose = FALSE)  
      dev1 <- deviance_o(FIT1)
      fstat <- (dev1 - S.hat)/s2.hat
      
      if (is.na(fstat) || fstat < 0) break
      newtau <- sgn * sqrt(fstat)
      
      if (abs(newtau - tau[count]) < 0.1) break
      count <- count + 1
      tau[count] <- newtau
      par.vals[count, ] <- pars <- coef(FIT1)[1L:npar]
      
      if (abs(tau[count]) > cutoff) break
      pars <- pars + ((pars - par.vals[count - 1, ]) * 
                        delta.t)/abs(tau[count] - tau[count - 1])
      pars[-par] <- pmin(upper[-par], pmax(lower[-par], pars[-par]))
    }
    
    ind <- seq_len(count)
    tau[ind] <- tau[rev(ind)]
    par.vals[ind, ] <- par.vals[rev(ind), ]
    sgn <- 1
    newmax <- count + maxpts
    pars <- par.vals[count, ]
    pars[par] <- base + profile.par.inc
    pars[par] <- pmin(upper[par], pmax(lower[par], pars[par]))
        
    while (count <= newmax) {      
      if (is.na(pars[par]) || isTRUE(all.equal(pars, par.vals[1, ])) 
          || pars[par] < lower[par] || pars[par] > 
            upper[par] || abs(pars[par] - base)/std.err[par] > 
            10 * cutoff) break
      
      FIT2 <- update(fitted, start = pars, fixed = FIXED, verbose = FALSE) 
      dev2 <- deviance_o(FIT2)
      fstat2 <- (dev2 - S.hat)/s2.hat
            
      if (is.na(fstat2) || fstat2 < 0) break
      newtau <- sgn * sqrt(fstat2)
      if (abs(newtau - tau[count]) < 0.1) break
      count <- count + 1
      tau[count] <- newtau
      par.vals[count, ] <- pars <- coef(FIT2)[1L:npar]
      if (abs(tau[count]) > cutoff) break
      pars <- pars + ((pars - par.vals[count - 1, ]) * 
                        delta.t)/abs(tau[count] - tau[count - 1])
      pars[-par] <- pmin(upper[-par], pmax(lower[-par], pars[-par]))
    }
    
    ind <- seq_len(count)    
    out[[par]] <- structure(list(tau = tau[ind], par.vals = par.vals[ind, , drop = FALSE]), 
                            class = "data.frame", row.names = as.character(ind), 
                            parameters = list(par = par, std.err = std.err[par]))
    
  }
  
  names(out)[which] <- names(coef(fitted))[which]
  out <- out[which]
  attr(out, "original.fit") <- fitted
  attr(out, "summary") <- f.summary
  class(out) <- c("profile.onls", "profile")    
  out
}

confint.onls <- function(object, parm, level = 0.95, ...) 
{
  pnames <- names(coef(object))
  if (missing(parm)) parm <- seq_along(pnames)
  if (is.numeric(parm)) parm <- pnames[parm]
  message("Waiting for profiling to be done...")
  utils::flush.console()
  object <- profile(object, which = parm, alphamax = (1 - level)/4)
  confint(object, parm = parm, level = level, ...)
}

confint.profile.onls <- function(object, parm = seq_along(pnames), level = 0.95, ...) 
{
  pnames <- names(object)
  ncoefs <- length(coef(attr(object, "original.fit")))
  of <- attr(object, "original.fit")
  if (is.numeric(parm)) parm <- pnames[parm]
  parm <- parm[parm %in% pnames]
  n <- length(fitted(of)) - length(coef(object))
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  pct <- paste(round(100 * a, 1), "%", sep = "")
  ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, pct))
  cutoff <- qt(a, n)
  for (pm in parm) {
    pro <- object[[pm]]
    sp <- if (ncoefs > 1) 
      spline(x = pro[, "par.vals"][, pm], y = pro$tau)
    else spline(x = pro[, "par.vals"], y = pro$tau)
    ci[pm, ] <- approx(sp$y, sp$x, xout = cutoff)$y
  }
  drop(ci)
}
