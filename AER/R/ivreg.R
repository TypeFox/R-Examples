ivreg <- function(formula, instruments, data, subset, na.action, weights, offset,
  contrasts = NULL, model = TRUE, y = TRUE, x = FALSE, ...)
{
  ## set up model.frame() call  
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  
  ## handle instruments for backward compatibility
  if(!missing(instruments)) {
    formula <- as.Formula(formula, instruments)
    cl$instruments <- NULL
    cl$formula <- formula(formula)
  } else {
    formula <- as.Formula(formula)
  }
  stopifnot(length(formula)[1] == 1L, length(formula)[2] %in% 1:2)
  
  ## try to handle dots in formula
  has_dot <- function(formula) inherits(try(terms(formula), silent = TRUE), "try-error")
  if(has_dot(formula)) {
    f1 <- formula(formula, rhs = 1)
    f2 <- formula(formula, lhs = 0, rhs = 2)
    if(!has_dot(f1) & has_dot(f2)) formula <- as.Formula(f1,
      update(formula(formula, lhs = 0, rhs = 1), f2))
  }
  
  ## call model.frame()
  mf$formula <- formula
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  
  ## extract response, terms, model matrices
  Y <- model.response(mf, "numeric")
  mt <- terms(formula, data = data)
  mtX <- terms(formula, data = data, rhs = 1)
  X <- model.matrix(mtX, mf, contrasts)
  if(length(formula)[2] < 2L) {
    mtZ <- NULL
    Z <- NULL
  } else {
    mtZ <- delete.response(terms(formula, data = data, rhs = 2))
    Z <- model.matrix(mtZ, mf, contrasts)
  }

  ## weights and offset
  weights <- model.weights(mf)
  offset <- model.offset(mf)
  if(is.null(offset)) offset <- 0
  if(length(offset) == 1) offset <- rep(offset, NROW(Y))
  offset <- as.vector(offset)

  ## call default interface
  rval <- ivreg.fit(X, Y, Z, weights, offset, ...)

  ## enhance information stored in fitted model object
  rval$call <- cl
  rval$formula <- formula(formula)
  rval$terms <- list(regressors = mtX, instruments = mtZ, full = mt)
  rval$na.action <- attr(mf, "na.action")
  rval$levels <- .getXlevels(mt, mf)
  rval$contrasts <- list(regressors = attr(X, "contrasts"), instruments = attr(Z, "contrasts"))
  if(model) rval$model <- mf
  if(y) rval$y <- Y
  if(x) rval$x <- list(regressors = X, instruments = Z, projected = rval$x)
    else rval$x <- NULL
      
  class(rval) <- "ivreg"
  return(rval)
}

ivreg.fit <- function(x, y, z, weights, offset, ...)
{
  ## model dimensions
  n <- NROW(y)
  p <- ncol(x)
  
  ## defaults
  if(missing(z)) z <- NULL
  if(missing(weights)) weights <- NULL
  if(missing(offset)) offset <- rep(0, n)
  
  ## sanity checks
  stopifnot(n == nrow(x))
  if(!is.null(z)) stopifnot(n == nrow(z))
  if(!is.null(weights)) stopifnot(n == NROW(weights))
  stopifnot(n == NROW(offset))
  
  ## project regressors x on image of instruments z
  if(!is.null(z)) {
    if(ncol(z) < ncol(x)) warning("more regressors than instruments")
    auxreg <- if(is.null(weights)) lm.fit(z, x, ...) else lm.wfit(z, x, weights, ...)
    xz <- as.matrix(auxreg$fitted.values)
    # pz <- z %*% chol2inv(auxreg$qr$qr) %*% t(z)
    colnames(xz) <- colnames(x)
  } else {
    xz <- x
    # pz <- diag(NROW(x))
    # colnames(pz) <- rownames(pz) <- rownames(x)
  }
  
  ## main regression
  fit <- if(is.null(weights)) lm.fit(xz, y, offset = offset, ...)
    else lm.wfit(xz, y, weights, offset = offset, ...)
 
  ## model fit information
  yhat <- drop(x %*% fit$coefficients)
  names(yhat) <- names(y)
  res <- y - yhat
  ucov <- chol2inv(fit$qr$qr[1:p, 1:p, drop = FALSE])
  colnames(ucov) <- rownames(ucov) <- names(fit$coefficients)
  rss <- if(is.null(weights)) sum(res^2) else sum(weights * res^2)
  ## hat <- diag(x %*% ucov %*% t(x) %*% pz)
  ## names(hat) <- rownames(x)

  rval <- list(
    coefficients = fit$coefficients,
    residuals = res,    
    fitted.values = yhat,
    weights = weights,
    offset = if(identical(offset, rep(0, n))) NULL else offset,
    n = n,
    nobs = if(is.null(weights)) n else sum(weights > 0),
    rank = fit$rank,
    df.residual = fit$df.residual,
    cov.unscaled = ucov,
    sigma = sqrt(rss/fit$df.residual), ## NOTE: Stata divides by n here and uses z tests rather than t tests...
    # hatvalues = hat,
    x = xz
  )
  
  return(rval)
}
   
vcov.ivreg <- function(object, ...)
  object$sigma^2 * object$cov.unscaled
    
bread.ivreg <- function (x, ...) 
    x$cov.unscaled * x$nobs

estfun.ivreg <- function (x, ...) 
{
    xmat <- model.matrix(x)
    wts <- weights(x)
    if(is.null(wts)) wts <- 1
    res <- residuals(x)
    rval <- as.vector(res) * wts * xmat
    attr(rval, "assign") <- NULL
    attr(rval, "contrasts") <- NULL
    return(rval)
}

hatvalues.ivreg <- function(model, ...) {
  xz <- model.matrix(model, component = "projected")
  x  <- model.matrix(model, component = "regressors")
  z  <- model.matrix(model, component = "instruments")
  solve_qr <- function(x) chol2inv(qr.R(qr(x)))
  diag(x %*% solve_qr(xz) %*% t(x) %*% z %*% solve_qr(z) %*% t(z))
}

terms.ivreg <- function(x, component = c("regressors", "instruments"), ...)
  x$terms[[match.arg(component)]]

model.matrix.ivreg <- function(object, component = c("projected", "regressors", "instruments"), ...) {
  component <- match.arg(component)
  if(!is.null(object$x)) rval <- object$x[[component]]
    else if(!is.null(object$model)) {
      X <- model.matrix(object$terms$regressors, object$model, contrasts = object$contrasts$regressors)
      Z <- if(is.null(object$terms$instruments)) NULL
        else model.matrix(object$terms$instruments, object$model, contrasts = object$contrasts$instruments)
      w <- weights(object)
      XZ <- if(is.null(Z)) X
        else if(is.null(w)) lm.fit(Z, X)$fitted.values else lm.wfit(Z, X, w)$fitted.values
      rval <- switch(component,
        "regressors" = X,
	"instruments" = Z,
	"projected" = XZ)
    } else stop("not enough information in fitted model to return model.matrix")
  return(rval)
}

predict.ivreg <- function(object, newdata, na.action = na.pass, ...)
{
  if(missing(newdata)) fitted(object)
  else {
    mf <- model.frame(delete.response(object$terms$full), newdata,
      na.action = na.action, xlev = object$levels)
    X <- model.matrix(delete.response(object$terms$regressors), mf,
      contrasts = object$contrasts$regressors)
    drop(X %*% object$coefficients)
  }
}

print.ivreg <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  cat("Coefficients:\n")
  print.default(format(coef(x), digits = digits), print.gap = 2, quote = FALSE)
  cat("\n")
  invisible(x)
}

summary.ivreg <- function(object, vcov. = NULL, df = NULL, diagnostics = FALSE, ...)
{
  ## weighted residuals
  res <- object$residuals
  y <- object$fitted.values + res
  n <- NROW(res)
  w <- object$weights
  if(is.null(w)) w <- rep(1, n)
  res <- res * sqrt(w)

  ## R-squared
  rss <- sum(res^2)
  if(attr(object$terms$regressors, "intercept")) {
    tss <- sum(w * (y - weighted.mean(y, w))^2)
    dfi <- 1    
  } else {
    tss <- sum(w * y^2)
    dfi <- 0
  }
  r.squared <- 1 - rss/tss
  adj.r.squared <- 1 - (1 - r.squared) * ((n - dfi)/object$df.residual)
  
  ## degrees of freedom (for z vs. t test)
  if(is.null(df)) df <- object$df.residual
  if(!is.finite(df)) df <- 0
  if(df > 0 & (df != object$df.residual)) {
    df <- object$df.residual
  }

  ## covariance matrix
  if(is.null(vcov.)) 
      vc <- vcov(object)
  else {
      if(is.function(vcov.)) vc <- vcov.(object)
        else vc <- vcov.
  }
  
  ## Wald test of each coefficient
  cf <- coeftest(object, vcov. = vc, df = df, ...)
  attr(cf, "method") <- NULL
  class(cf) <- "matrix"
  
  ## Wald test of all coefficients
  Rmat <- if(attr(object$terms$regressors, "intercept"))
    cbind(0, diag(length(coef(object))-1)) else diag(length(coef(object)))
  waldtest <- linearHypothesis(object, Rmat, vcov. = vcov., test = ifelse(df > 0, "F", "Chisq"))
  waldtest <- c(waldtest[2,3], waldtest[2,4], abs(waldtest[2,2]), if(df > 0) waldtest[2,1] else NULL)
  
  ## diagnostic tests
  diag <- if(diagnostics) ivdiag(object, vcov. = vcov.) else NULL
  
  rval <- list(
    call = object$call,
    terms = object$terms,
    residuals = res,
    weights <- object$weights,
    coefficients = cf,
    sigma = object$sigma,
    df = c(object$rank, if(df > 0) df else Inf, object$rank), ## aliasing not handled yet
    r.squared = r.squared,
    adj.r.squared = adj.r.squared,
    waldtest = waldtest,
    vcov = vc,
    diagnostics = diag)
    
  class(rval) <- "summary.ivreg"
  return(rval)
}
 
print.summary.ivreg <- function(x, digits = max(3, getOption("digits") - 3), 
    signif.stars = getOption("show.signif.stars"), ...)
{
  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat(if(!is.null(x$weights) && diff(range(x$weights))) "Weighted ", "Residuals:\n", sep = "")      
  if(NROW(x$residuals) > 5L) {
      nam <- c("Min", "1Q", "Median", "3Q", "Max")
      rq <- if(length(dim(x$residuals)) == 2) 
	  structure(apply(t(x$residuals), 1, quantile), dimnames = list(nam, dimnames(x$residuals)[[2]]))
      else structure(quantile(x$residuals), names = nam)
      print(rq, digits = digits, ...)
  } else {
      print(x$residuals, digits = digits, ...)
  }

  cat("\nCoefficients:\n")
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars,
    signif.legend = signif.stars & is.null(x$diagnostics), na.print = "NA", ...)

  if(!is.null(x$diagnostics)) {
    cat("\nDiagnostic tests:\n")
    printCoefmat(x$diagnostics, cs.ind = 1L:2L, tst.ind = 3L,
      has.Pvalue = TRUE, P.values = TRUE, digits = digits,
      signif.stars = signif.stars, na.print = "NA", ...)  
  }

  cat("\nResidual standard error:", format(signif(x$sigma, digits)),
    "on", x$df[2L], "degrees of freedom\n")

  cat("Multiple R-Squared:", formatC(x$r.squared, digits = digits))
  cat(",\tAdjusted R-squared:", formatC(x$adj.r.squared, digits = digits),
    "\nWald test:", formatC(x$waldtest[1L], digits = digits),
    "on", x$waldtest[3L], if(length(x$waldtest) > 3L) c("and", x$waldtest[4L]) else NULL,       
    "DF,  p-value:", format.pval(x$waldtest[2L], digits = digits), "\n\n")

  invisible(x)
}
    
anova.ivreg <- function(object, object2, test = "F", vcov = NULL, ...)
{
  rval <- waldtest(object, object2, test = test, vcov = vcov)
  if(is.null(vcov)) {
    head <- attr(rval, "heading")
    head[1] <- "Analysis of Variance Table\n"
    rss <- sapply(list(object, object2), function(x) sum(residuals(x)^2))
    dss <- c(NA, -diff(rss))
    rval <- cbind(rval, cbind("RSS" = rss, "Sum of Sq" = dss))[,c(1L, 5L, 2L, 6L, 3L:4L)]
    attr(rval, "heading") <- head
    class(rval) <- c("anova", "data.frame")
  }
  return(rval)
}

ivdiag <- function(obj, vcov. = NULL) {
  ## extract data
  y <- model.response(model.frame(obj))
  x <- model.matrix(obj, component = "regressors")
  z <- model.matrix(obj, component = "instruments")
  
  ## endogenous/instrument variables
  endo <- which(!(colnames(x) %in% colnames(z)))
  inst <- which(!(colnames(z) %in% colnames(x)))
  if((length(endo) <= 0L) | (length(inst) <= 0L))
    stop("no endogenous/instrument variables")

  ## return value
  rval <- matrix(NA, nrow = length(endo) + 2L, ncol = 4L)
  colnames(rval) <- c("df1", "df2", "statistic", "p-value")
  rownames(rval) <- c(if(length(endo) > 1L) paste0("Weak instruments (", colnames(x)[endo], ")") else "Weak instruments",
    "Wu-Hausman", "Sargan")
  
  ## convenience functions
  lmfit <- function(x, y) {
    rval <- lm.fit(x, y)
    rval$x <- x
    rval$y <- y
    return(rval)
  }
  rss <- function(obj) sum(obj$residuals^2)
  wald <- function(obj0, obj1, vcov. = NULL) {
    df <- c(obj1$rank - obj0$rank, obj1$df.residual)
    if(!is.function(vcov.)) {
      w <- ((rss(obj0) - rss(obj1)) / df[1L]) / (rss(obj1)/df[2L])
    } else {
      if(NCOL(obj0$coefficients) > 1L) {
        cf0 <- structure(as.vector(obj0$coefficients),
	  .Names = c(outer(rownames(obj0$coefficients), colnames(obj0$coefficients), paste, sep = ":")))
        cf1 <- structure(as.vector(obj1$coefficients),
	  .Names = c(outer(rownames(obj1$coefficients), colnames(obj1$coefficients), paste, sep = ":")))
      } else {
        cf0 <- obj0$coefficients
        cf1 <- obj1$coefficients
      }
      ovar <- which(!(names(cf1) %in% names(cf0)))
      vc <- vcov.(lm(obj1$y ~ 0 + obj1$x))
      w <- t(cf1[ovar]) %*% solve(vc[ovar,ovar]) %*% cf1[ovar]
      w <- w / df[1L]
    }
    pval <- pf(w, df[1L], df[2L], lower.tail = FALSE)
    c(df, w, pval)
  }
    
  # Test for weak instruments
  for(i in seq_along(endo)) {
  aux0 <- lmfit(z[, -inst, drop = FALSE], x[, endo[i]])
  aux1 <- lmfit(z,                        x[, endo[i]])
  rval[i, ] <- wald(aux0, aux1, vcov. = vcov.)
  }

  ## Wu-Hausman test for endogeneity
  if(length(endo) > 1L) aux1 <- lmfit(z, x[, endo])
  xfit <- as.matrix(aux1$fitted.values)
  colnames(xfit) <- paste("fit", colnames(xfit), sep = "_")
  auxo <- lmfit(      x,        y)
  auxe <- lmfit(cbind(x, xfit), y)
  rval[nrow(rval) - 1L, ] <- wald(auxo, auxe, vcov. = vcov.)

  ## Sargan test of overidentifying restrictions 
  r <- residuals(obj)  
  auxs <- lmfit(z, r)
  rval[nrow(rval), 1L] <- length(inst) - length(endo)
  if(rval[nrow(rval), 1L] > 0L) {
    rval[nrow(rval), 3L] <- length(r) * (1 - rss(auxs)/sum((r - mean(r))^2))
    rval[nrow(rval), 4L] <- pchisq(rval[nrow(rval), 3L], rval[nrow(rval), 1L], lower.tail = FALSE)
  }

  return(rval)
}

## If #Instruments = #Regressors then
##   b = (Z'X)^{-1} Z'y
## and solves the estimating equations
##   Z' (y - X beta) = 0
## For
##   cov(y) = Omega
## the following holds
##   cov(b) = (Z'X)^{-1} Z' Omega Z (X'Z)^{-1}
##   
## Generally:  
##   b = (X' P_Z X)^{-1} X' P_Z y
## with estimating equations
##   X' P_Z (y - X beta) = 0
## where P_Z is the usual projector (hat matrix wrt Z) and
##   cov(b) = (X' P_Z X)^{-1} X' P_Z Omega P_Z X (X' P_Z X)^{-1}
## Thus meat is X' P_Z Omega P_Z X and bread i (X' P_Z X)^{-1}
## 
## See
##   http://www.stata.com/support/faqs/stat/2sls.html
