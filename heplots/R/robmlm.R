### robust multivariate linear regression

## John Fox 2012-06-02
## revised: 2013-08-20 to avoid calling summary.mlm() directly in vcov.mlm()

robmlm <- function(X, ...){
  UseMethod("robmlm")
}

robmlm.default <- function(X, Y, w, P=2*pnorm(4.685, lower.tail=FALSE), 
                           tune, max.iter=100, psi=psi.bisquare, tol=1e-6, 
                           initialize, verbose=FALSE, ...){
  # args:
  #   X: model matrix, including constant (if present)
  #   Y: response matrix
  #   w: prior weights
  #   P: two-tail probability, to find cutoff quantile for chisq (tuning constant);
  #       default is set for bisquare weight function
  #   tune: tuning constant (if given directly)
  #   max.iter: iteration limit
  #   psi: robustness weight function; psi.bisquare (from MASS) is default
  #   tol: convergence tolerance, maximum relative change in coefficients
  #   initialize: modeling function to find start values for coefficients,
  #       equation-by-equation; if absent WLS is used
  #   verbose: show iteration history?
  #if (!require(MASS)) stop("MASS package missing")
  p <- ncol(Y)
  if (missing(w) || is.null(w)) w <- rep(1, nrow(Y))
  if (missing(tune)) tune <- qchisq(1 - P, df=p)
  if (missing(initialize)){
    fit.last <- lm.wfit(X, Y, w)
  }
  else{
    k <- ncol(X)
    n <- nrow(X)
    coef <- matrix(0, k, p)
    res <- matrix(0, n, p)
    for (j in 1:p){
      modj <- initialize(Y[, j] ~ X - 1, weights=w)
      coef[, j] <- coef(modj)
      res[, j] <- residuals(modj)
    }
    fit.last <- list(coefficients=coef, residuals=res)
  }
  B.last <- B.new <- fit.last$coefficients
  iter <- 0
  if (verbose){
    coefnames <- abbreviate(outer(rownames(B.new), colnames(B.new), function(x, y) paste(x, ":", y, sep="")), 10)
    b <- as.vector(B.new)
    names(b) <- coefnames
    cat("\n", iter, ":\n", sep="")
    print(b)
  }
  repeat {
    iter <- iter + 1
    if (iter > max.iter) break
    E <- fit.last$residuals
    S <- cov.trob(E, center=FALSE)$cov
    mahal <- mahalanobis(E, 0, S)
    wts <- psi(mahal, tune)
    fit.new <- lm.wfit(X, Y, w*wts)
    B.last <- B.new
    B.new <- fit.new$coefficients
    if (verbose){
      b <- as.vector(B.new)
      names(b) <- coefnames
      cat("\n", iter, ":\n", sep="")
      print(b)
    }
    if (max(abs((B.last - B.new)/(B.last + tol))) < tol) break
    fit.last <- fit.new
  }
  if (iter > max.iter) warning("maximum iterations exceeded")
  fit.new$iterations <- iter
  fit.new$weights <- wts
  fit.new$converged <- iter <= max.iter
  fit.new
}

robmlm.formula <- function(formula, data, subset, weights, na.action, model = TRUE,
                         contrasts = NULL, ...) {
  # ... passed to robmlm.default
  call <- match.call()  
  call[[1]] <- as.name("robmlm")
  mf <- match.call(expand.dots = FALSE)  
  args <- match(c("formula", "data", "subset", "weights", "na.action"),
                names(mf), 0)  
  mf <- mf[c(1, args)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval.parent(mf)  
  terms <- attr(mf, "terms")  
  Y <- model.response(mf, type="numeric")  
  X <- model.matrix(terms, mf, contrasts) 
  w <- as.vector(model.weights(mf))
  mod <- robmlm.default(X, Y, w, ...)
  mod$na.action <- attr(mf, "na.action")
  mod$contrasts <- attr(X, "contrasts")
  mod$xlevels <- .getXlevels(terms, mf)
  mod$call <- call
  mod$terms <- terms
  if (model)  mod$model <- mf
  class(mod) <- c("robmlm", "mlm", "lm")
  mod
}

print.robmlm <- function(x, ...){
  if (!x$converged) warning("failed to converge")
  NextMethod()
  cat("iterations = ", x$iterations)
  invisible(x)
}

summary.robmlm <- function(object, ...){
  res <- list()
  res[[1]] <- NextMethod()
  res$iterations <- object$iterations
  res$converged <- object$converged
  class(res) <- c("summary.robmlm", class(res[[1]]))
  res
}

print.summary.robmlm <- function(x, ...){
  if (!x$converged) warning("failed to converge")
  print(x[[1]])
  cat("iterations = ", x$iterations)
  invisible(x)
}

vcov.mlm <- function (object, ...) {
# override stats::vcov.mlm to allow weights
#   adapted from code provided by Michael Friendly
# For R 3.1.0, to avoid calling summary.mlm directly, change the object class
#   temporarily to c("mlm", "lm")
	SSD.mlm <- function (object, ...) {
		if (!is.null(object$weights)) { 
			SSD <- wcrossprod(residuals(object), w=object$weights)
			df <- sum(object$weights)
		}
		else {
			SSD <- crossprod(residuals(object))
			df <- object$df.residual
		}
		structure(list(SSD=SSD, call = object$call, 
						df = df), class = "SSD")
	}
	estVar.mlm <- function (object, ...) estVar(SSD(object))
	obj <- object
	class(obj) <- c("mlm", "lm")  # remove robmlm class temporarily
	so <- summary(obj)[[1L]]
	kronecker(estVar(object), so$cov.unscaled, make.dimnames = TRUE)
}
