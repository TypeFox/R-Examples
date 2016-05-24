## some convenience functions generally applicable to survreg
## hopefully in survival soon

deviance.survreg <- function(object, ...)
  sum(residuals(object, type = "deviance")^2)

fitted.survreg <- function(object, ...)
  predict(object, type = "response", se.fit = FALSE)

nobs.survreg <- function(object, ...)
  length(object$linear.predictors)

weights.survreg <- function(object, ...)
  model.weights(model.frame(object))


## convenience tobit() interface to survreg()

tobit <- function(formula, left = 0, right = Inf, dist = "gaussian", subset = NULL, data = list(), ...)
{
  ## remember original environment
  oenv <- environment(formula)
  oformula <- eval(formula)
  
  ## process censoring
  stopifnot(all(left < right))
  lfin <- any(is.finite(left))
  rfin <- any(is.finite(right))
  
  ## formula processing: replace dependent variable
  ## original
  y <- formula[[2]]  
  if(lfin & rfin) { 
  ## interval censoring
    formula[[2]] <- call("Surv", call("ifelse", call(">=", y, substitute(right)), substitute(right), 
      call("ifelse", call("<=", y, substitute(left)), substitute(left), y)), time2 = substitute(right),
      call("ifelse", call(">=", y, substitute(right)), 0, call("ifelse", call("<=", y, substitute(left)), 2, 1)),
      type = "interval")
  } else if(!rfin) {
  ## left censoring
    formula[[2]] <- call("Surv", call("ifelse", call("<=", y, substitute(left)), substitute(left), y),
      call(">", y, substitute(left)) , type = "left")
  } else {
  ## right censoring
    formula[[2]] <- call("Surv", call("ifelse", call(">=", y, substitute(right)), substitute(right), y),
      call("<", y, substitute(right)) , type = "right")
  }
  
  ## call survreg
  cl <- ocl <- match.call()
  cl$formula <- formula
  cl$left <- NULL
  cl$right <- NULL
  cl$dist <- dist
  cl[[1]] <- as.name("survreg")
  rval <- eval(cl, oenv)

  ## slightly modify result
  class(rval) <- c("tobit", class(rval))
  ocl$formula <- oformula
  rval$call <- ocl
  rval$formula <- formula
  return(rval)
}


## add printing and summary methods that are more similar to 
## the corresponding methods for lm objects

print.tobit <- function(x, digits = max(3, getOption("digits") - 3), ...) 
{
  ## failure
  if(!is.null(x$fail)) {
    cat("tobit/survreg failed.", x$fail, "\n")
    return(invisible(x))
  }
  
  ## call
  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  ## coefficients
  coef <- x$coefficients
  if(any(nas <- is.na(coef))) {
    if (is.null(names(coef))) names(coef) <- paste("b", 1:length(coef), sep = "")
    cat("Coefficients: (", sum(nas), " not defined because of singularities)\n", sep = "")
  } else cat("Coefficients:\n")
  print.default(format(coef, digits = digits), print.gap = 2, quote = FALSE)

  ## scale
  if(nrow(x$var) == length(coef)) 
    cat("\nScale fixed at", format(x$scale, digits = digits), "\n")
  else if (length(x$scale) == 1) 
    cat("\nScale:", format(x$scale, digits = digits), "\n")
  else {
    cat("\nScale:\n")
    print(format(x$scale, digits = digits), ...)
  }
  
  ## return    
  cat("\n")
  invisible(x)
}

summary.tobit <- function(object, correlation = FALSE, symbolic.cor = FALSE, vcov. = NULL, ...) 
{
  ## failure
  if(!is.null(object$fail)) {
    warning("tobit/survreg failed.", object$fail, "   No summary provided\n")
    return(invisible(object))
  }
  
  ## rank
  if(all(is.na(object$coefficients))) {
    warning("This model has zero rank --- no summary is provided")
    return(invisible(object))
  }

  ## vcov
  if(is.null(vcov.)) vcov. <- vcov(object)
  else {
    if(is.function(vcov.)) vcov. <- vcov.(object)
  }
  
  ## coefmat
  coef <- coeftest(object, vcov. = vcov., ...)
  attr(coef, "method") <- NULL
  
  ## Wald test
  nc <- length(coef(object))
  has_intercept <- attr(terms(object), "intercept") > 0.5
  wald <- if(nc <= has_intercept) NULL else linearHypothesis(object,
    if(has_intercept) cbind(0, diag(nc-1)) else diag(nc),
    vcov. = vcov.)[2,3]
  ## instead of: waldtest(object, vcov = vcov.)

  ## correlation
  correlation <- if(correlation) cov2cor(vcov.) else NULL
    
  ## distribution
  dist <- object$dist
  if(is.character(dist)) sd <- survreg.distributions[[dist]]
    else sd <- dist
  if(length(object$parms)) pprint <- paste(sd$name, "distribution: parmameters =", object$parms)
    else pprint <- paste(sd$name, "distribution")

  ## number of observations
  nobs <- as.vector(table(factor(object$y[,2], levels = 2:0)))
  nobs <- c(sum(nobs), nobs)
  names(nobs) <- c("Total", "Left-censored", "Uncensored", "Right-censored")

  rval <- object[match(c("call", "df", "loglik", "iter", "na.action", "idf", "scale"),
    names(object), nomatch = 0)]
  rval <- c(rval, list(coefficients = coef, correlation = correlation,
    symbolic.cor = symbolic.cor, parms = pprint, n = nobs, wald = wald))

  class(rval) <- "summary.tobit"
  return(rval)
}

print.summary.tobit <- function(x, digits = max(3, getOption("digits") - 3), ...) 
{
  ## call
  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  ## observations and censoring
  if(length(x$na.action)) cat("Observations: (", naprint(x$na.action), ")\n", sep = "")
    else cat("Observations:\n")
  print(x$n)

  ## coefficients
  if(any(nas <- is.na(x$coefficients[,1])))
    cat("\nCoefficients: (", sum(nas), " not defined because of singularities)\n", sep = "")
  else cat("\nCoefficients:\n")
  printCoefmat(x$coefficients, digits = digits, ...)

  ## scale
  if("Log(scale)" %in% rownames(x$coefficients))
    cat("\nScale:", format(x$scale, digits = digits), "\n")
  else
    cat("\nScale fixed at", format(x$scale, digits = digits), "\n")

  ## logLik and Chi-squared test
  cat(paste("\n", x$parms, "\n", sep = ""))
  cat("Number of Newton-Raphson Iterations:", format(trunc(x$iter)), "\n")
  cat("Log-likelihood:", formatC(x$loglik[2], digits = digits), "on", x$df, "Df\n")
  if(!is.null(x$wald)) cat("Wald-statistic:", formatC(x$wald, digits = digits),
    "on", sum(x$df) - x$idf, "Df, p-value:",
    format.pval(pchisq(x$wald, sum(x$df) - x$idf, lower.tail = FALSE)), "\n")

  ## correlation
  correl <- x$correlation
  if (!is.null(correl)) {
    p <- NCOL(correl)
    if (p > 1) {
      cat("\nCorrelation of Coefficients:\n")
      if (is.logical(x$symbolic.cor) && x$symbolic.cor) {
        print(symnum(correl, abbr.colnames = NULL))
      }
      else {
        correl <- format(round(correl, 2), nsmall = 2, digits = digits)
    	correl[!lower.tri(correl)] <- ""
    	print(correl[-1, -p, drop = FALSE], quote = FALSE)
      }
    }
  }

  ## return
  cat("\n")
  invisible(x)
}


## as the apparent y ~ ... and actual Surv(y) ~ ... formula
## differ, some standard functionality has to be done by work-arounds

formula.tobit <- function(x, ...) x$formula

model.frame.tobit <- function(formula, ...)
{
  Call <- formula$call
  Call[[1]] <- as.name("model.frame")
  Call <- Call[match(c("", "formula", "data", "weights", "subset", "na.action"), names(Call), 0)]
  dots <- list(...)
  nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0)]
  Call[names(nargs)] <- nargs
  Call$formula <- formula$formula
  env <- environment(formula$terms)
  if(is.null(env)) env <- parent.frame()
  eval(Call, env)
}

update.tobit <- function(object, formula., ..., evaluate = TRUE)
{
  call <- object$call
  extras <- match.call(expand.dots = FALSE)$...
  if(!missing(formula.)) {
    ff <- formula(object)
    ff[[2]] <- call$formula[[2]]
    call$formula <- update.formula(ff, formula.)
  }
  if (length(extras) > 0) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  if(evaluate) eval(call, parent.frame())
  else call  
}

waldtest.tobit <- function(object, ..., test = c("Chisq", "F"), name = NULL)
{
  if(is.null(name)) name <- function(x) paste(deparse(x$call$formula), collapse="\n")
  waldtest.default(object, ..., test = match.arg(test), name = name)
}

lrtest.tobit <- function(object, ..., name = NULL)
{
  if(is.null(name)) name <- function(x) paste(deparse(x$call$formula), collapse="\n")
  lrtest.default(object, ..., name = name)
}

linearHypothesis.tobit <- function(model, hypothesis.matrix,
  rhs = NULL, vcov. = NULL, ...)
{
  if(compareVersion(packageDescription("car")$Version, "2.0-19") < 0L) {
    stop("The linearHypothesis() method requires package car version 2.0-19 or later.")
  }
  if(is.null(vcov.)) {
    vcov. <- vcov(model)
  } else {
    if(is.function(vcov.)) vcov. <- vcov.(model)
  }
  vcov. <- vcov.[-nrow(vcov.), -ncol(vcov.)]
  model$formula <- model$call$formula
  car::linearHypothesis.default(model,
    hypothesis.matrix = hypothesis.matrix, rhs = rhs, vcov. = vcov., ...)
}
