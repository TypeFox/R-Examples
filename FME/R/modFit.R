## -----------------------------------------------------------------------------
## Fitting the Model to Data
## -----------------------------------------------------------------------------

## a small utility function
overrule <- function(ini, new) if(is.null(new)) ini else new

modFit <- function(f, p, ..., lower = -Inf, upper = Inf,
                   method = c("Marq", "Port", "Newton", "Nelder-Mead", "BFGS", "CG",
                   "L-BFGS-B", "SANN", "Pseudo"), jac = NULL,
                   control = list(), hessian = TRUE) {

  ## check if valid input...
  np <- length(p)
  if (length(lower) != np & length(lower) != 1)
    stop("length of 'lower' should be either 1 or equal to number of parameters")
  if (length(upper) != np & length(upper) != 1)
    stop("length of 'upper' should be either 1 or equal to number of parameters")

  method <- match.arg(method)
  pnames <- names(p)

  Lower <- rep(lower, len = length(p))
  Upper <- rep(upper, len = length(p))
  if (any(p < lower) | any(p > upper))
    stop("values of 'p' should be inbetween 'lower' and 'upper'")


  ## are boundaries allowed in the method?
  bounds <- method %in% c("L-BFGS-B", "Port", "Pseudo")

  FF      <- NULL
  useCost <- method != "Marq" # marquardt uses residuals, others model cost

  Func <- function(p, ...) {
    ## The 'global' assignment ensures that, upon returning from the fitting
    ## the function need NOT be called once more, to obtain the residuals,
    ## sum of squares etc...
    FF <<- f(p, ...)
    cM <- class(FF) == "modCost"

    if (cM && useCost ) return(FF$model)
    if (cM)             return(FF$residuals$res)
    if (useCost) return (sum(FF^2)) else return(FF)
  }

  ## The methods cannot work with the "gradient function" - too confusing
  ## if jacobian is provided and the method is one of "", then the gradient
  ## function can be created from func and jac

  dots <- list(...)
  nmdots <- names(dots)
  if (method %in% c("BFGS", "CG", "L-BFGS") & "gr" %in% nmdots)
    stop ("cannot use gradient function 'gr' here; supply jacobian 'jac' instead")

  if (method == "Port" & "gradient" %in% nmdots)
    stop ("cannot use gradient function 'grad' here; supply jacobian 'jac' instead")

  grad <- NULL

  if (! is.null(jac))
    if (method %in% c("BFGS", "CG", "L-BFGS", "Port") & ! is.null(jac)){
      if (np ==1 )
        grad <- function (x) return(sum(2*f(x) * jac(x)))
      else
        grad <- function (x) return(colSums(2*f(x, ...) * jac(x, ...)))
    }

  Pars <- p
  estHess <- FALSE
  # Adapt function call if necessary
  if (bounds)                            # no need to change parameters
    Fun <- function(p, ...) Func(p, ...)
  else {
    # 1. Adapt initial parameter values
    lower <- -Inf
    upper <- Inf

    # lower and upper bounds...
    lu   <- which(is.finite(Lower) & is.finite(Upper))
    Pars[lu] <- tan(pi*((p[lu] - Lower[lu])/(Upper[lu] - Lower[lu]) - 0.5))

    # just lower bounds...
    l <- which(is.finite(Lower) & !is.finite(Upper))
    Pars[l] <-log(p[l]-Lower[l])

    # just upper bounds...
    u <- which(!is.finite(Lower) & is.finite(Upper))
    Pars[u] <- log(-p[u]+Upper[u])

    # 2. parameter transformation in function call
    Fun <- function(p, ...) {
      PP     <- p
      PP[lu] <- Lower[lu]+(Upper[lu] - Lower[lu]) * (atan(p[lu])/pi + 0.5)
      PP[l]  <- Lower[l] + exp(p[l])
      PP[u]  <- Upper[u] - exp(p[u])
      names(PP) <- pnames          # nlminb 'forgets' parameter names
      Func(PP, ...)
    }
  }
  limits <- (! bounds && length(c(lu, l, u)) > 0)

  ## optimise...
  if (method %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN")) {
    res <- optim(Pars, Fun, method = method, lower = lower, upper = upper,
                 control = control, hessian = hessian, gr = grad, ...)
    names(res)[2] <- "ssr"         # called "value" here
  }

  else if (method == "Port") {
    res <- nlminb(start = Pars, objective = Fun, lower = lower, upper = upper,
                  control = control, gradient = grad, ...)
    names(res)[2] <- "ssr"         # called "objective" here
    names(res)[6] <- "counts"
    if(hessian) estHess <- TRUE    # hessian not estimated by Port...
  }

  else if (method == "Marq") {
    Contr <- nls.lm.control()
    nmsC <- names(Contr)
    if (! "maxiter" %in% names(control))
      control$maxiter <- 100       # override too low default  (50)
    Contr[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC]) > 0)
      warning("unknown names in control: ", paste(noNms, collapse = ", "))

    res <- as.list(nls.lm(par = Pars, fn = Fun, control = Contr, ...)[])
    # renaming results for compatibility with other methods
    names(res)[7] <- "iterations"  # called "niter" here
    names(res)[9] <- "ssr"         # called "deviance" here
    names(res)[3] <- "residuals"
    res$hessian   <- 2*res$hessian # returns 0.5*hessian!
    Diag          <- unlist(res[6])
    res$diag      <- NULL
    res$diag      <- Diag
  }

  else if (method == "Newton") {
    nmsC <- c("typsize", "fscale", "print.level", "ndigit",
              "gradtol", "stepmax", "steptol", "iterlim", "check.analyticals")
    namc <- names(control)
    if (length(noNms <- namc[!namc %in% nmsC]) > 0)
       warning("unknown names in control: ", paste(noNms, collapse = ", "))

    typsize <- overrule( rep(1, length(p)), control$typsize)
    res <- nlm(p = Pars, f = Fun, ..., hessian = hessian,
               typsize = typsize,
               fscale = overrule(1, control$fscale),
               print.level = overrule(0, control$print.level),
               ndigit = overrule(12, control$ndigit),
               gradtol = overrule(1e-6, control$gradtol),
               stepmax = overrule(max(1000 * sqrt(sum((p/typsize)^2)), 1000),
                                  control$stepmax),
               steptol = overrule(1e-6, control$steptol),
               iterlim = overrule(100, control$iterlim),
               check.analyticals = overrule(TRUE, control$check.analyticals))
    # renaming results for compatibility with other methods
    names(res)[1] <- "ssr"  # called "minimum" here
    names(res)[2] <- "par"  # called "estimate" here
  }

  else if (method == "Pseudo") {
    res <- pseudoOptim(p = Pars, f = Fun, lower = lower, upper = upper,
                       control = control, ...)
    names(res)[2] <- "ssr"  #is called "cost" here
    if(hessian)
      estHess<-TRUE
  }


  if (limits) {
    respar      <- res$par
    res$par[lu] <- Lower[lu] + (Upper[lu]-Lower[lu])*(atan(respar[lu])/pi + 0.5)
    res$par[l]  <- Lower[l] + exp(respar[l])
    res$par[u]  <- Upper[u] - exp(respar[u])
    if (! bounds)
      estHess<-TRUE
  }

  names(res$par)<-names(p)
  class(res)    <- "modFit"
  # Karline: moved till here
  if (estHess)  {  # hessian (re)-estimated using backtransformed values
    useCost <- FALSE

    Fun <- function(p, ...) Func(p, ...)
    if (! is.null(jac))
      Jac <- jac(res$par)
    else Jac <- gradient(Fun, res$par, centered = TRUE, ...)
    res$hessian <- 2 * t(Jac) %*% Jac
  }

  if (!method == "Marq")
    if (class(FF) == "modCost")
      res$residuals <- FF$residuals$res
    else res$residuals  <- FF

  # mean square
  res$ms <- res$ssr/length(res$residuals)

  # mean square per varaible
  if (class(FF) == "modCost") {
    names(res$residuals)  <- FF$residuals$name
    res$var_ms            <- FF$var$SSR/FF$var$N
    res$var_ms_unscaled   <- FF$var$SSR.unscaled/FF$var$N
    res$var_ms_unweighted <- FF$var$SSR.unweighted/FF$var$N

    names(res$var_ms_unweighted) <- names(res$var_ms_unscaled) <-
      names(res$var_ms) <- FF$var$name
  } else res$var_ms <- res$var_ms_unweighted <- res$var_ms_unscaled <- NA

  res$rank <- np
  res$df.residual <- length(res$residuals) - res$rank
  if(!useCost & length(res$residuals) <= 1)
    stop ("Levenberg-Marquardt requires residuals - only one value is returned")

  res

}

## -----------------------------------------------------------------------------
## S3 methods of modFit
## -----------------------------------------------------------------------------

deviance.modFit <- function(object, ...)
  object$ssr

coef.modFit <- function(object, ...)
  unlist(object$par)

residuals.modFit <- function(object, ...)
  object$residuals

df.residual.modFit <- function(object, ...)
  object$df.residual

## -----------------------------------------------------------------------------

## inspired by summary.nls.lm
summary.modFit <- function (object, cov = TRUE,...) {
  param  <- object$par
  pnames <- names(param)
  p      <- length(param)
  covar  <- try(solve(0.5*object$hessian), silent = TRUE)   # unscaled covariance
  if (!is.numeric(covar)) {
    message <- "Cannot estimate covariance; system is singular"
    warning(message)
    covar <- matrix(data = NA, nrow = p, ncol = p)
  } else message <- "ok"

  rownames(covar) <- colnames(covar) <-pnames
  rdf    <- object$df.residual
  resvar <- object$ssr / rdf
  se     <- sqrt(diag(covar) * resvar)
  names(se) <- pnames
  tval      <- param / se
  modVariance <- object$ssr / length(object$residuals)

  param <- cbind(param, se, tval, 2 * pt(abs(tval), rdf, lower.tail = FALSE))
  dimnames(param) <- list(pnames, c("Estimate", "Std. Error",
                                    "t value", "Pr(>|t|)"))
  if(cov)
    ans <- list(residuals = object$residuals,
                residualVariance = resvar,
                sigma = sqrt(resvar),
                modVariance = modVariance,
                df = c(p, rdf), cov.unscaled = covar,
                cov.scaled = covar * resvar,
                info = object$info, niter = object$iterations,
                stopmess = message,
                par = param)
  else
    ans <- list(residuals = object$residuals,
                residualVariance = resvar,
                sigma = sqrt(resvar),
                modVariance = modVariance,
                df = c(p, rdf),
                info = object$info, niter = object$iterations,
                stopmess = message,
                par = param)
  class(ans) <- "summary.modFit"
  ans
}

## -----------------------------------------------------------------------------

print.summary.modFit <- function (x, digits = max(3, getOption("digits") - 3), ...) {
  df  <- x$df
  rdf <- df[2]
  cat("\nParameters:\n")
  printCoefmat(x$par, digits = digits, ...)
  cat("\nResidual standard error:",
      format(signif(x$sigma, digits)), "on", rdf, "degrees of freedom\n")

  printcor <- !is.null(x$cov.unscaled)
  if (printcor){
    Corr <- cov2cor(x$cov.unscaled)
    rownames(Corr) <- colnames(Corr) <- rownames(x$par)
    cat("\nParameter correlation:\n")
    print(Corr, digits = digits, ...)
  }

  invisible(x)
}

## -----------------------------------------------------------------------------

plot.modFit <- function(x, ask = NULL, ...) {
  vn  <- unique(names(x$residuals))
  lvn <- max(1, length(vn))

  np <- NP <- lvn
  if(! is.null(x$rsstrace)) np <- np + 1

  dots   <- list(...)
  nmdots <- names(dots)

  ## Set par mfrow and ask.
  ask <- setplotpar (nmdots, dots, np, ask)

  ## interactively wait if there are remaining figures
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }

  if(! is.null(x$rsstrace))
    plot(x$rsstrace, main = "residual sum of squares", xlab = "iteration",
         ylab = "-", log = "y")
  if(is.null(vn))
    plot(x$residuals, ylab = "-", main = "residuals", pch = 18, ...)
  else
    for(i in vn) {
      ii <- which(names(x$residuals) == i)
      plot(x$residuals[ii], main = i, ylab = "residuals", pch = 18, ...)
    }
}
